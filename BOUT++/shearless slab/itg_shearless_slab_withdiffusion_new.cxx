/*************************************************************************  
 * 2-fluid equations for the ITG-mode in a shearless slab                      
          
 * This file evolves the normalised equations                                  
 
 * A diffusion term is added to surpess the numerical instability              
 *************************************************************************/

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <cmath>

class ITGslab : public PhysicsModel {
private:
  // equilibrium quantities                                                     
  Field3D Ni0, Ti0; // equilibrium density, temperature      

  // 3D involving fields  
  
  Field3D Ni, Ti, vi, phi, Te, vEx; // ion density, ion temperature, ion parallel veloci \
ty                                                                              
                                                                               \

  // 3D vectors                                                                
                                                                               \

  Vector3D vE, B, checkphi, checkTi; // drift velocity, magnetic field         
                                                                               \


  // constants                                                                  
  BoutReal mi, e, B_norm; // ion mass, elementary charge, B-field magnitude

  // reference scales                                                           
  BoutReal cs, rho_i, wci; // ion sound speed, ion Larmor radius, ion gyro-freq\
uency                                                                           

  
  BoutReal LrefLT, hypdif;

  // use this to check the size of the box
  //Field3D x=0,y=0,z=0;

  int init(bool UNUSED(restarting)) {

    Coordinates* coords = mesh->getCoordinates();

    // generate coordinate system                                               

    coords->g11 = 1.0;
    coords->g22 = 1.0;
    coords->g33 = 1.0;
    coords->g12 = 0.0;
    coords->g13 = 0.0;
    coords->g23 = 0.0;

    coords->g_11 = 1.0;
    coords->g_22 = 1.0;
    coords->g_33 = 1.0;
    coords->g_12 = 0.0;
    coords->g_13 = 0.0;
    coords->g_23 = 0.0;

    coords->geometry();


    // get options from the input file                                         

    auto& options = Options::root()["ITG-slab"];

    // Read parameters from input file                                          
                                                                              

    // The doc() provides some documentaion in BOUT.settings                    

    mi = options["mi"].doc("ion mass").withDefault(1.0);
    e = options["e"].doc("elementary charge").withDefault(1.0);

    Ni0 = options["Ni0"].doc("Initial Density").withDefault(1.0);
    //Te = options["Te"].doc("Electron temperature").withDefault(10.0);
    Ti0 = options["Ti0"].doc("Initial ion temperature").withDefault(1.0);

    LrefLT = options["LrefLT"].doc("reference length scale").withDefault(1.0);
    hypdif = options["hypdif"].doc("hyperdiffusion").withDefault(10.0);


    // write ion density and temperature as initial + perturbation
    //  Ni += Ni0;
    // Ti += Ti0;
    vi = 0.0;
    

    // tell BOUT which variable to evolve (normalised variables)                

    SOLVE_FOR(Ni); // normalised ion density                                
    SOLVE_FOR(Ti); // normalised ion temperature                            
    SOLVE_FOR(vi); // normalised ion velocity (parallel to B-field)         

    // Set background B-field                                                   
    B.x = 0.0;
    B.y = 1.0; // set B-field parallel to y-direction                           
    B.z = 0.0;                                                    

    // scales                                                                   
    // cs = sqrt(Te/mi);                                                        
    //B_norm = sqrt(B*B) // does not work for some reason                       
    // B_norm = 1;                                                              
    // rho_i = (mi*cs)/(e*B_norm);                                              
    // wci = cs/rho_i;                                                          

    // phi0 = Te/e; // assuming Ni=Ni0 so that they cancel out                  
    
    // set initial values for phi and vE
    phi = 0.0;
    vE = 0.0;

    Te = Ti0;
    

    // define the pow4 function

    #define pow4(f) (f*f*f*f)
    
    // use this to check the size of the box
    /*
    FieldFactory f(x->getMesh());
    x = f.create3D("x");
    y = f.create3D("y");
    z = f.create3D("z");
    */
    
    // save variables in the output file                                        
    SAVE_REPEAT(phi, vE);
    //SAVE_ONCE(x, y, z);
    // SAVE_ONCE(wci, rho_i, cs);                                               

   return 0;
  }

  int rhs(BoutReal UNUSED(time)) override{

    // communicate Ni, Ti, vi                                      
    mesh->communicate(Ni, Ti, vi);

    phi= Ni;//*Te calculate phi (Te=1 in this case)                               

    // communicate phi
    mesh->communicate(phi);

    vE = cross(Grad(phi), B); // calculate vE
                            

    // communicate vE                                                           
    mesh->communicate(vE);

    // write equations                                                          

    {
      TRACE("ddt(Ti)");
      ddt(Ti) = vE.x*LrefLT;
    }

    {
      TRACE("ddt(Ni)");
      ddt(Ni) = -Grad_par(vi);
    }

    {
      TRACE("ddt(vi)");
      ddt(vi) = -Grad_par(phi) - Grad_par(Ti);
    }
    
    // Add hyperdiffusion to surpress the numerical instability
    ddt(Ti) -= hypdif*pow4(mesh->getCoordinates()->dy)*D4DY4(Ti);
    ddt(Ni) -= hypdif*pow4(mesh->getCoordinates()->dy)*D4DY4(Ni);
    ddt(vi) -= hypdif*pow4(mesh->getCoordinates()->dy)*D4DY4(vi);
    

    return 0;
  }
};

BOUTMAIN(ITGslab)
