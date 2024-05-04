#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts_rec_arrays);
  Vector Xe_arr;
  Vector ne_arr;
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };
  Vector peebles_array;
  int j = 0;
  std::cout << std::boolalpha;

  // Calculate recombination history
  bool saha_regime = true;
  bool calculated_spline = false;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    //std::cout << x_array[i] << " : " << Xe_current << std::endl;
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      Xe_arr.emplace_back(Xe_current);
      ne_arr.emplace_back(ne_current);

    } else {

      // Solve the ODE
      if (not calculated_spline) {
        double Xeini = Xe_current;
        Vector Xe_ic{Xeini};
        Vector sliced_x = Vector(x_array.begin() + i, x_array.end());
        peebles_Xe_ode.solve(dXedx, sliced_x, Xe_ic);
        peebles_array = peebles_Xe_ode.get_data_by_component(0);
        calculated_spline = true;
      }

      Xe_arr.emplace_back(peebles_array[j]);
      double scale_factor = std::exp(x_array[i]);
      double rho_crit = 3 * std::pow(cosmo->get_H0(), 2) / (8. * M_PI * Constants.G);
      double nb = cosmo->get_OmegaB() * rho_crit / (Constants.m_H * std::pow(scale_factor, 3));
      ne_arr.emplace_back(Xe_arr[i] * nb);
      j++;
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  Vector log_Xe;
  Vector log_ne;
  for (int i = 0; i < Xe_arr.size(); i++) {
    log_Xe.emplace_back(std::log(Xe_arr[i]));
    log_ne.emplace_back(std::log(ne_arr[i]));
  }

  log_Xe_of_x_spline.create(x_array, log_Xe, "Function Xe");
  log_ne_of_x_spline.create(x_array, log_ne, "Function ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  /*
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  //const double OmegaB      = cosmo->get_OmegaB();
  //...
  //...

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  */

  double scale_factor = std::exp(x);
  double rho_crit = 3 * std::pow(cosmo->get_H0(), 2) / (8. * M_PI * Constants.G);
  double nb = cosmo->get_OmegaB() * rho_crit / (Constants.m_H * std::pow(scale_factor, 3));
  double Tb = cosmo->get_TCMB() / scale_factor;
  double a = 1.;
  double b = 1/nb * std::pow(Constants.m_e * Tb / (2 * M_PI), 3/2.) * std::exp(-Constants.epsilon_0 / (Constants.k_b * Tb)) * std::pow(Constants.k_b / (Constants.hbar * Constants.hbar), 3/2.);
  double c = -b;
  //double Xe = (std::sqrt(b*b - 4*a*c) - b) / (2. * a);
  double Xe = -4*a*c/((std::sqrt(b*b - 4*a*c) + b) * (2. * a));
  double ne = Xe * nb;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double OmegaB = cosmo->get_OmegaB();
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...
  double scale_factor = std::exp(x);
  double rho_crit = 3 * std::pow(cosmo->get_H0(), 2) / (8. * M_PI * Constants.G);
  double Tb = cosmo->get_TCMB() / scale_factor;
  double fine_structure_constant = Constants.m_e * Constants.c / Constants.hbar * std::sqrt(3 * Constants.sigma_T / (8. * M_PI));
  double phi2 = 0.448 * std::log(Constants.epsilon_0 /(Constants.k_b * Tb));
  double alpha2 = 64 * M_PI / std::sqrt(27. * M_PI) * fine_structure_constant * fine_structure_constant / (Constants.m_e * Constants.m_e) * std::sqrt(Constants.epsilon_0 / (Constants.k_b * Tb)) * phi2 * Constants.hbar * Constants.hbar / Constants.c;
  double beta = alpha2 * std::pow(Constants.m_e * Tb / (2. * M_PI), 3/2.) * std::exp(-Constants.epsilon_0 /(Constants.k_b * Tb)) * std::pow(Constants.k_b / (Constants.hbar * Constants.hbar), 3/2.);
  double beta2 = alpha2 * std::pow(Constants.m_e * Tb / (2. * M_PI), 3/2.) * std::exp(-Constants.epsilon_0 /(4. * Constants.k_b * Tb)) * std::pow(Constants.k_b / (Constants.hbar * Constants.hbar), 3/2.);
  double nb = cosmo->get_OmegaB() * rho_crit / (Constants.m_H * std::pow(scale_factor, 3));
  double n1s = (1 - X_e) * nb;
  double Lambda_alpha = cosmo->H_of_x(x) * std::pow(3 * Constants.epsilon_0, 3) / (std::pow(8 * M_PI, 2) * n1s) / std::pow(Constants.hbar * Constants.c, 3);
  double Cr = (Constants.lambda_2s1s + Lambda_alpha) / (Constants.lambda_2s1s + Lambda_alpha + beta2);
  dXedx[0] = Cr / cosmo->H_of_x(x) * (beta * (1 - X_e) - nb * alpha2 * X_e * X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = - Constants.c * Constants.sigma_T * this->ne_of_x(x)/ cosmo->H_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //double tauini = 0.0;
  double tauini = - Constants.c * Constants.sigma_T * this->ne_of_x(x_end)/ cosmo->H_of_x(x_end);
  Vector tau_ic{tauini};
  ODESolver ode;
  ode.solve(dtaudx, x_array, tau_ic);
  auto tau = ode.get_data_by_component(0);
  std::reverse(tau.begin(), tau.end());
  std::reverse(x_array.begin(), x_array.end());
  tau_of_x_spline.create(x_array, tau, "Function tau");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Vector g_tilde;
  for (int i = 0; i < x_array.size(); i++)
    g_tilde.emplace_back(-this->dtaudx_of_x(x_array[i]) * std::exp(-this->tau_of_x(x_array[i])));

  g_tilde_of_x_spline.create(x_array, g_tilde, "Function visibility");

  double r = 4 * cosmo->get_OmegaR(0)/(3. * cosmo->get_OmegaB(0));
  double scale_factor;
  double shini = Constants.c * std::sqrt(r / std::exp(Constants.x_start) / (3 * (1 + r / std::exp(Constants.x_start)))) / cosmo->Hp_of_x(Constants.x_start);
  ODEFunction dshdx = [&](double x, const double *sh, double *dshdx){
    dshdx[0] = Constants.c * std::sqrt(r / std::exp(x) / (3 * (1 + r / std::exp(x)))) / cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };
  Vector sh_ic{shini};
  ODESolver ode2;
  ode2.solve(dshdx, x_array, tau_ic);
  auto sound_horizon = ode2.get_data_by_component(0);
  sound_horizon_of_x_spline.create(x_array, sound_horizon, "Function sound horizon");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return this->tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return this->tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return this->tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return this->g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return this->g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return this->g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  double Xe = std::exp(this->log_Xe_of_x_spline(x));
  return Xe;
}

double RecombinationHistory::ne_of_x(double x) const{
  double ne = std::exp(this->log_ne_of_x_spline(x));
  return ne;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

double RecombinationHistory::sound_horizon_of_x(double x) const{
  return this->sound_horizon_of_x_spline(x);
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 1e4;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                                                   << " ";
    fp << cosmo->t_of_x(x)                                    << " ";
    fp << this->Xe_of_x(x)                                    << " ";
    fp << this->ne_of_x(x)                                    << " ";
    fp << this->tau_of_x(x)                                   << " ";
    fp << this->dtaudx_of_x(x)                                << " ";
    fp << this->ddtauddx_of_x(x)                              << " ";
    fp << this->g_tilde_of_x(x)                               << " ";
    fp << this->dgdx_tilde_of_x(x)                            << " ";
    fp << this->ddgddx_tilde_of_x(x)                          << " ";
    fp << this->sound_horizon_of_x(x)                         << " ";
    fp << this->electron_fraction_from_saha_equation(x).first << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

