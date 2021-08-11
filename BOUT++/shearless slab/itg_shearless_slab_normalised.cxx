/*************************************************************************     \
                                                                                
 * 2-fluid equations for the ITG-mode in a shearless slab                      \
                                                                                
 * This file contains the normalised equations				       \
                                                                     
 *                                                                             \
                                                                                
 *                                                                             \
                                                                                
 *                                                                             \
                                                                                
 *************************************************************************/

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>


class ITGslab : public PhysicsModel {
private:
  // equilibrium quantities
  Field3D Ni0, Ti0, phi0; // equilibrium density, temperature, potential
  
  // 3D involving fields                                                       
                                                                                
  Field3D Ni, Ti, vi, phi; // ion density, ion temperature, ion parallel velocity                                        
                                                                                                                                                             
  // 3D vectors                                                                
                                                                                
  Vector3D vE, B, checkphi, checkTi; // drift velocity, magnetic field         
                                                                                

  // constants
  BoutReal mi, e, B_norm, Te; // ion mass, elementary charge, B-field magnitude

  // reference scales
  BoutReal cs, rho_i, wci; // ion sound speed, ion Larmor radius, ion gyro-frequency

  // normalised variables
  Field3D Ni_hat, Ti_hat, vi_hat, phi_hat, nabla_par_hat, Te_hat;
  Vector3D vE_hat, nabla_hat;

  BoutReal LrefLT; 
  
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

    Ni0 = options["Ni0"].doc("Initial Density").withDefault(1.0e20);
    Te = options["Te"].doc("Electron temperature").withDefault(100.0);
    Ti0 = options["Ti0"].doc("Initial ion temperature").withDefault(10.0);

    LrefLT = options["LrefLT"].doc("reference length scale").withDefault(1.0);

    
    // tell BOUT which variable to evolve (normalised variables)

    SOLVE_FOR(Ni_hat); // normalised ion density
    SOLVE_FOR(Ti_hat); // normalised ion temperature
    SOLVE_FOR(vi_hat); // normalised ion velocity (parallel to B-field)

    // Set background B-field
    B.x = 0.0;
    B.y = 1.0; // set B-field parallel to y-direction
    B.z = 0.0;


    // normalisations                                                                              
    Ni_hat = Ni_hat/Ni0;
    Ti_hat = Ti_hat/Ti0;
    vi_hat = 1.0;

    
    // write ion density and ion temperature as initial + perturbation
    //Ni_hat += Ni0;
    //Ti_hat += Ti0;

    // scales
    // cs = sqrt(Te/mi);
    //B_norm = sqrt(B*B) // does not work for some reason
    // B_norm = 1;
    // rho_i = (mi*cs)/(e*B_norm);
    // wci = cs/rho_i;

    // phi0 = Te/e; // assuming Ni=Ni0 so that they cancel out
    Te_hat = Te/Ti0;

    phi_hat = 0.0;
    vE = 0.0;
    
    // save variables in the output file
    SAVE_REPEAT(phi_hat, vE);
    // SAVE_ONCE(wci, rho_i, cs);

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override{

    // communicate Ni_hat, Ti_hat, vi_hat
    mesh->communicate(Ni_hat, Ti_hat, vi_hat, Te_hat);
    
    phi_hat = Ni_hat*Te_hat; // calculate phi_hat

    // communicate phi_hat
    mesh->communicate(phi_hat);

    vE = cross(-Grad(phi_hat), B)/(B*B); // calculate vE
   
    // communicate vE
    mesh->communicate(vE);

    // write equations

    {
      TRACE("ddt(Ti_hat)");
      ddt(Ti_hat) = -vE.x*LrefLT;
    }

    {
      TRACE("ddt(Ni_hat");
      ddt(Ni_hat) = -Grad_par(vi_hat);
    }

    {
      TRACE("ddt(vi_hat");
      ddt(vi_hat) = -Grad_par(phi_hat) - Grad_par(Ti_hat);
    }

    return 0;
  }
};

BOUTMAIN(ITGslab);


    

