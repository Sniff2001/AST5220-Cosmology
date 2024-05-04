#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  //compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array(n_k);
  k_array = Utils::linspace(std::log(k_min), std::log(k_max), n_k);
  std::for_each(k_array.begin(), k_array.end(), [](double &n) {n = std::exp(n);});
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  Vector3D y(Constants.n_scalars + 3, Vector2D(n_x, Vector(n_k, 0)));
  bool run_once = false;

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k, x_array);
    auto x_end_tight_idx = std::lower_bound(x_array.begin(), x_array.end(), x_end_tight) - x_array.begin();
    Vector sliced_x = Vector(x_array.begin(), x_array.begin() + x_end_tight_idx);

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    ODESolver solver_tight;
    solver_tight.solve(dydx_tight_coupling, sliced_x, y_tight_coupling_ini);
    Vector y_tight_coupling = solver_tight.get_final_data();
    for (int i = 0; i < y_tight_coupling.size(); i++) {
      if (run_once)
        break;
      std::cout << y_tight_coupling[i] << std::endl;
    }
    run_once = true;

    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector sliced_x_full = Vector(x_array.begin() + x_end_tight_idx, x_array.end());
    ODESolver solver_full;
    solver_full.solve(dydx_full, sliced_x_full, y_full_ini);

    //std::cout << y.size() << " " << y[0].size() << " " << y[0][0].size() << std::endl;

    for (int i = 0; i < (Constants.ind_start_theta + 2); i++) {
      Vector tc = solver_tight.get_data_by_component(i);
      Vector fl = solver_full.get_data_by_component(i);
      // Concatenates the two vectors
      tc.insert(tc.end(), fl.begin(), fl.end());
      for (int j = 0; j < n_x; j++) {
        y[i][j][ik] = tc[j];
      }
    }
    // Theta2 is needed, but is not included in the tight coupling
    Vector theta2 = solver_full.get_data_by_component(Constants.ind_start_theta + 2);
    for (int j = 0; j < x_end_tight_idx; j++)
      y[Constants.ind_start_theta + 2][j][ik] = -20 / (45. * rec->dtaudx_of_x(x_array[j]) * cosmo->Hp_of_x(x_array[j])) * Constants.c * k * y[Constants.ind_start_theta + 1][j][k];
    for (int j = x_end_tight_idx; j < n_x; j++)
      y[Constants.ind_start_theta + 2][j][ik] = theta2[j - x_end_tight_idx];

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
  delta_cdm_spline.create(x_array, k_array, y[0], "Function delta cdm");
  delta_b_spline.create(x_array, k_array, y[1], "Function delta b");
  v_cdm_spline.create(x_array, k_array, y[2], "Function vcdm");
  v_b_spline.create(x_array, k_array, y[3], "Function vb");
  Phi_spline.create(x_array, k_array, y[4], "Function Phi");
  Pi_spline.create(x_array, k_array, y[7], "Function Pi");
  Vector2D Psi(n_x, Vector(n_k, 0));
  for (int i = 0; i < n_x; i++)
    for (int j = 0; j < n_k; j++)
      Psi[i][j] = -y[4][i][j] - 12 * std::pow(cosmo->get_H0() / (Constants.c * k_array[j] * std::exp(x_array[i])), 2) * (cosmo->get_OmegaR(0) * y[7][i][j]);
  Psi_spline.create(x_array, k_array, Psi, "Function Psi");
  Theta0_spline.create(x_array, k_array, y[5], "Function Theta0");
  Theta1_spline.create(x_array, k_array, y[6], "Function Theta1");
  Theta2_spline.create(x_array, k_array, y[7], "Function Theta2");
  Theta_spline = {Theta0_spline, Theta1_spline, Theta2_spline};
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================
  // ...
  // ...

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  // ...
  // ...
  double ckHp = Constants.c * k / cosmo->Hp_of_x(x);
  double fv = 0.;
  double Psi = - 1 / (3 / 2. + 2 * fv / 5.);
  Phi = - (1 + 2 * fv / 5.) * Psi;
  delta_cdm = - 3 / 2. * Psi;
  delta_b = delta_cdm;
  v_cdm = - ckHp / 2. * Psi;
  v_b = v_cdm;

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  // ...
  Theta[0] = - 0.5 * Psi;
  Theta[1] = ckHp / 6. * Psi;

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  // ...
  // ...
  double ckHp = Constants.c * k / cosmo->Hp_of_x(x);
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  // ...
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -20 / (45. * rec->dtaudx_of_x(x)) * ckHp * Theta_tc[1];
  for (int l = 3; l < n_ell_theta; l++)
    Theta[l] = - l / (2. * l + 1.) * ckHp / rec->dtaudx_of_x(x) * Theta[l - 1];

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    // ...
    // ...
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k, Vector x_array) const{
  int i;
  for (i = 0; i < x_array.size(); i++) {
    if (std::abs(rec->dtaudx_of_x(x_array[i])) < 10. * Constants.c * k / cosmo->Hp_of_x(x_array[i]) ||
        std::abs(rec->dtaudx_of_x(x_array[i])) < 10. ||
        x_array[i] > -8.3) {
      break;
    }
  }
  double x_tight_coupling_end = x_array[i];

  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  double ckHp = Constants.c * k / cosmo->Hp_of_x(x);
  double a = std::exp(x);
  double Theta2 = - 20 / 45. * ckHp / rec->dtaudx_of_x(x) * Theta[1];
  double Psi = -Phi - 12 * std::pow(cosmo->get_H0() / (Constants.c * k * a), 2) * (cosmo->get_OmegaR(0) * Theta2);
  double R = 4 * cosmo->get_OmegaR(0) / (3 * cosmo->get_OmegaB(0) * a);
  dThetadx[0] = - ckHp * Theta[1] - dPhidx;
  double q = (-((1 - R) * rec->dtaudx_of_x(x) + (1 + R) * rec->ddtauddx_of_x(x)) * (3 * Theta[1] + v_b) + ckHp * (- Psi + (1 - cosmo->dHpdx_of_x(x)/cosmo->Hp_of_x(x)) * (-Theta[0] + 2 * Theta2) - (- ckHp * Theta[1] - dPhidx))) / ((1 + R) * rec->dtaudx_of_x(x) + cosmo->dHpdx_of_x(x) / cosmo->Hp_of_x(x) - 1);
  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - ckHp * ckHp / 3. * Phi + std::pow(cosmo->get_H0() / cosmo->Hp_of_x(x), 2) / 2. * (cosmo->get_OmegaCDM(0) / a * delta_cdm + cosmo->get_OmegaB(0) / a * delta_b + 4 * cosmo->get_OmegaR(0) / (a * a) * Theta[0]);
  ddelta_cdmdx = ckHp * v_cdm - 3 * dPhidx;
  dv_cdmdx = -v_cdm - ckHp * Psi;
  ddelta_bdx = ckHp * v_b - 3 * dPhidx;
  dv_bdx = 1 / (1 + R) * (-v_b - ckHp * Psi + R * (q + ckHp * (-Theta[0] + 2 * Theta2) - ckHp * Psi));

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = - ckHp * Theta[1] - dPhidx;
  dThetadx[1] = 1 / 3. * (q - dv_bdx);
  //dThetadx[1] = ckHp / 3. * (Theta[0] - 2 * Theta[2] + dPhidx) + rec->dtaudx_of_x(x) * (Theta[1] + 1/3. * v_b);

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  // double Hp = cosmo->Hp_of_x(x);
  // ...

  // Recombination variables
  // ...

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  double ckHp = Constants.c * k / cosmo->Hp_of_x(x);
  double a = std::exp(x);
  double Psi = -Phi - 12 * std::pow(cosmo->get_H0() / (Constants.c * k * a), 2) * (cosmo->get_OmegaR(0) * Theta[2]);
  double R = 4 * cosmo->get_OmegaR(0) / (3 * cosmo->get_OmegaB(0) * a);

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi - ckHp * ckHp / 3. * Phi + std::pow(cosmo->get_H0() / cosmo->Hp_of_x(x), 2) / 2. * (cosmo->get_OmegaCDM(0) / a * delta_cdm + cosmo->get_OmegaB(0) / a * delta_b + 4 * cosmo->get_OmegaR(0) / (a * a) * Theta[0]);
  ddelta_cdmdx = ckHp * v_cdm - 3 * dPhidx;
  dv_cdmdx = - v_cdm - ckHp * Psi;
  ddelta_bdx = ckHp * v_b - 3 * dPhidx;
  dv_bdx = - v_b - ckHp * Psi + rec->dtaudx_of_x(x) * R * (3 * Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = - ckHp * Theta[1] - dPhidx;
  dThetadx[1] = ckHp / 3. * (Theta[0] - 2 * Theta[2] + Psi) + rec->dtaudx_of_x(x) * (Theta[1] + 1 / 3. * v_b);
  // This only in non tight coupling
  for (int l = 2; l < (n_ell_theta - 1); l++) {
    dThetadx[l] = l / (2. * l + 1.) * ckHp * Theta[l-1] - (l - 1) / (2. * l + 1.) * ckHp * Theta[l+1] + rec->dtaudx_of_x(x) * (Theta[l] - 1 / 10. * Theta[2] * (l == 2));
  }
  int last_theta = n_ell_theta - 1; // Because we index from 0!
  dThetadx[last_theta] = ckHp * Theta[last_theta - 1] - Constants.c * (last_theta + 1) / (cosmo->Hp_of_x(x) * cosmo->eta_of_x(x)) * Theta[last_theta] + rec->dtaudx_of_x(x) * Theta[last_theta];

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // ...
    // ...
    // ...
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  std::cout << "Tight coupling time (k = " << k << "): " << this->get_tight_coupling_time(k, x_array) << std::endl;
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    /*fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";*/
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

