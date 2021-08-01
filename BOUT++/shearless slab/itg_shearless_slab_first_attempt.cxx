/*************************************************************************
 * 2-fluid equations for the ITG-mode in a shearless slab
 *
 *
 *
 *
 *************************************************************************/

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>


class ITGslab : public PhysicsModel {
private:
  // 3D involving fields
  Field3D Ni, Ti, vi, phi, Te; // ion density, ion temperature, ion parallel velocity, potential, electron temperature

  // 2D fields
  Field2D Ni0, Ti0; // equilibrium ion density and ion temperature 

  // 3D vectors
  Vector3D vE, B, checkphi, checkTi; // drift velocity, magnetic field

  // parameters
  BoutReal mi; // ion mass
  BoutReal e; // elementary charge
  
  
  int init(bool UNUSED(restarting)) {
    auto *coords = mesh->getCoordinates();
   
    coords->g11 = 1.0;
    coords->g22 = 1.0;
    coords->g33 = 1.0;
    coords->g12 = 0.0;
    coords->g13 = 0.0;
    coords->g23 = 0.0;
   
    // coords->g_11 = 1.0;
    // coords->g_22 = 1.0;
    // coords->g_33 = 1.0;
    // coords->g_12 = 0.0;
    // coords->g_13 = 0.0;
    //coords->g_23 = 0.0;

    coords->geometry();
   

    // Get option from the input file
    auto& options = Options::root()["ITG-slab"];

    // Read from BOUT.inp, setting default to 1.0
    // The doc() provides some documentaion in BOUT.settings
    
    mi = options["mi"].doc("Ion mass").withDefault(1.0);
    e = options["e"].doc("Elementary charge").withDefault(1.0);

    Ni0 = options["Ni0"].doc("Initial Density").withDefault(1.0);
    Te = options["Te"].doc("Electron temperature").withDefault(1.0);
    Ti0 = options["Ti0"].doc("Ion temperature").withDefault(1.0);
    
    // tell BOUT which variables to evolve
    SOLVE_FOR(Ni);
    SOLVE_FOR(Ti);
    SOLVE_FOR(vi);

    // Background magnetic field
    B.x = 0.0;
    B.y = 1.0; // set the B-field parallel to the y-direction 
    B.z = 0.0;

    Ni += Ni0;
    Ti += Ti0;
    
    SAVE_REPEAT(phi,vE, B, checkphi.x, checkphi.y, checkphi.z, checkTi.x, checkTi.y, checkTi.z); // save variables in output file
 
  return 0;
  }

  int rhs(BoutReal UNUSED(time)) override{

    
    phi = Ni*Te/e; // Calculate phi

    vE = cross(-Grad(phi), B)/(B*B); // calculate vE

    checkphi  = Grad(phi);
    checkTi = Grad(Ti);
    
    mesh->communicate(Ni, Ti, vi, phi, vE);

    {
      TRACE("ddt(Ti)");
      ddt(Ti) = -V_dot_Grad(vE,Ti);
    }
    
    {
      TRACE("ddt(Ni)");
      ddt(Ni) = -Ni0*Grad_par(vi);
    }

    {
      TRACE("ddt(vi)");
      ddt(vi) = -(e/mi)*Grad_par(phi) - (1/mi)*Grad_par(Ti);
    }

    return 0;
  }
};

BOUTMAIN(ITGslab);

    


 
