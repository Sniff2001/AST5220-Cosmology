#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}
//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector log_k_array = Utils::linspace(std::log(k_min), std::log(k_max), n_k);
  Vector k_array = exp(log_k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  const int z_max = 30000;
  Vector z_array(z_max);
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    Vector j_ell(z_max);
    for (int z = 0; z < z_max; z++) {
      z_array[z] = z;
      j_ell[z] = Utils::j_ell(ell, z);
    }
    
    j_ell_splines[i].create(z_array, j_ell);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");
  const int n_x        = 1000;
  const double x_start = Constants.x_start;
  const double x_end   = Constants.x_end;
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  for(size_t ik = 0; ik < k_array.size(); ik++){
    double k = k_array[ik];
    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    // ...
    // ...
    // ...
    for(size_t j = 0; j < ells.size(); j++){
      const int ell = ells[j];
      double theta_ini = 0;
      Vector theta_ic{theta_ini};

      ODEFunction dthetadx = [&](double x, const double *y, double *dydx){
        dydx[0] = j_ell_splines[j](k * (cosmo->eta_of_x(0) - cosmo->eta_of_x(x))) * source_function(x, k);
        return GSL_SUCCESS;
      };
      ODESolver solver_theta;
      solver_theta.solve(dthetadx, x_array, theta_ic);

    // Store the result for Source_ell(k) in results[ell][ik]
      result[j][ik] = solver_theta.get_final_data_by_component(0);
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  const int n_x        = 1000;
  const double x_start = Constants.x_start;
  const double x_end   = Constants.x_end;
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  for (int i = 0; i < nells; i++)
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i]);

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){

    // ...
    // ...
    // ...
    // ...

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  Vector result(nells);
  //Vector k_array = log_k_array;
  //std::for_each(k_array.begin(), k_array.end(), [](double &n) {n = std::exp(n);});
  Vector k_array = exp(log_k_array);
  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================


  for (int i = 0; i < nells; i++) {
    double theta_ini = 0;
    Vector theta_ic{theta_ini};

    ODEFunction cell = [&](double k, const double *y, double *dydx){
      dydx[0] = 4 * M_PI * primordial_power_spectrum(k) * f_ell_spline[i](k) * g_ell_spline[i](k) / k;
      return GSL_SUCCESS;
    };
    ODESolver solver_cell;
    solver_cell.solve(cell, k_array, theta_ic);

    result[i] = solver_cell.get_final_data_by_component(0);
  }


  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  double delta_M = std::pow(Constants.c * k_mpc / Constants.Mpc / cosmo->get_H0(), 2) * 2 * pert->get_Phi(x,k_mpc / Constants.Mpc) * std::exp(x) / (3. * (cosmo->get_OmegaB(0) + cosmo->get_OmegaCDM(0)));
  double pofk = delta_M * delta_M * primordial_power_spectrum(k_mpc / Constants.Mpc) * 2 * M_PI * M_PI / std::pow(k_mpc / Constants.Mpc, 3);
  return pofk * std::pow(cosmo->get_h() / Constants.Mpc, 3);
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                  << " ";
    fp << cell_TT_spline( ell ) * normfactor   << " ";
    fp << cell_TT_spline( ell )                << " ";
    /*
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }*/
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_ps(std::string filename) const{
  std::ofstream fp(filename.c_str());
  Vector log_k_array = Utils::linspace(std::log(k_min), std::log(k_max), 1e5);
  Vector k_array = exp(log_k_array);
  auto print_data = [&] (const double k) {
    fp << k                                 << " ";
    fp << thetaT_ell_of_k_spline[0](k)      << " ";
    fp << thetaT_ell_of_k_spline[3](k)      << " ";
    fp << thetaT_ell_of_k_spline[7](k)      << " ";
    fp << thetaT_ell_of_k_spline[19](k)      << " ";
    fp << thetaT_ell_of_k_spline[32](k)      << " ";
    fp << thetaT_ell_of_k_spline[42](k)      << " ";
    fp << thetaT_ell_of_k_spline[62](k)      << " ";
    fp << primordial_power_spectrum(k)      << " ";
    fp << get_matter_power_spectrum(Constants.x_end, k * Constants.Mpc) << " ";
    /*
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }*/
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

