#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  H0 = Constants.H0_over_h * this->h;
  OmegaR = 2. * std::pow(M_PI, 2) / 30. * std::pow(Constants.k_b * this->TCMB, 4) / (std::pow(Constants.hbar, 3) * std::pow(Constants.c, 5)) * 8. * M_PI * Constants.G / (3. * std::pow(this->H0, 2));
  OmegaNu = this->Neff * 7. / 8. * std::pow(4. / 11., 4. / 3.) * this->OmegaR;
  OmegaLambda = 1 - (this->OmegaK + this->OmegaB + this->OmegaCDM + this->OmegaR + this->OmegaNu);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(int npts){
  Utils::StartTiming("Eta");
  
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    detadx[0] = Constants.c/this->Hp_of_x(x);

    return GSL_SUCCESS;
  };

  // Set the IC (with the given IC the solution is y(x) = c/H(x)
  double etaini = 0.0;
  Vector eta_ic{etaini};

  // Solve the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  // Get the solution (we only have one component so index 0 holds the solution)
  auto eta_array = ode.get_data_by_component(0);

  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    dtdx[0] = 1/this->H_of_x(x);

    return GSL_SUCCESS;
  };

  double tini = 1/(2*this->H_of_x(Constants.x_start));
  Vector t_ic{tini};

  ode.solve(dtdx, x_array, t_ic);
  auto t_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "Function eta");
  t_of_x_spline.create(x_array, t_array, "Function t");

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  if (x == 0.0) return this->H0;

  double a = std::exp(x);
  double H = this->H0 * std::sqrt((this->OmegaB + this->OmegaCDM) / std::pow(a, 3) + (this->OmegaR + this->OmegaNu) / std::pow(a, 4) + this->OmegaK / std::pow(a, 2) + this->OmegaLambda);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double a = std::exp(x);
  double Hp = a * H_of_x(x);
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double a = this->OmegaB + this->OmegaCDM;
  double b = this->OmegaR + this->OmegaNu;
  double c = this->OmegaK;
  double d = this->OmegaLambda;
  double dHpdx = (this->H0 * std::exp(-3*x) * (2 * d * std::exp(4 * x) - a * std::exp(x) - 2 * b)) 
                / (2 * std::sqrt(c * std::exp(-2 * x) + a * std::exp(-3 * x) + b * std::exp(-4 * x) + d));

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  double a = this->OmegaB + this->OmegaCDM;
  double b = this->OmegaR + this->OmegaNu;
  double c = this->OmegaK;
  double d = this->OmegaLambda;
  double ddHpddx = (this->H0 * std::exp(-7 * x) * (4 * std::pow(d, 2) * std::exp(8 * x) + 8 * c * d * std::exp(6 * x)
                  + 14 * a * d * std::exp(5 * x) + 24 * b * d * std::exp(4 * x) + 2 * a * c * std::exp(3 * x)
                  + (8 * b * c + std::pow(a, 2)) * std::exp(2 * x) + 6 * a * b * std::exp(x) + 4 * std::pow(b,2))) 
                  / (4 * std::pow(c * std::exp(-2 * x)+ a * std::exp(-3 * x) + b * std::exp(-4 * x) + d, (3./2.)));

  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  
  double a = std::exp(x);
  double OmegaB_of_a = this->OmegaB / (std::pow(a, 3) * std::pow(this->H_of_x(x) / this->H0, 2));
  return OmegaB_of_a;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  double a = std::exp(x);
  double OmegaR_of_a = this->OmegaR / (std::pow(a, 4) * std::pow(this->H_of_x(x) / this->H0, 2));
  return OmegaR_of_a;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  double a = std::exp(x);
  double OmegaNu_of_a = this->OmegaNu / (std::pow(a, 4) * std::pow(this->H_of_x(x) / this->H0, 2));
  return OmegaNu_of_a;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  double a = std::exp(x);
  double OmegaCDM_of_a = this->OmegaCDM / (std::pow(a, 3) * std::pow(this->H_of_x(x) / this->H0, 2));
  return OmegaCDM_of_a;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double a = std::exp(x);
  double OmegaLambda_of_a = this->OmegaLambda / (std::pow(this->H_of_x(x) / this->H0, 2));
  return OmegaLambda_of_a;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  double a = std::exp(x);
  double OmegaK_of_a = this->OmegaK / (std::pow(a, 2) * std::pow(this->H_of_x(x) / this->H0, 2));

  return OmegaK_of_a;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double r = this->get_comoving_distance_of_x(x);

  double a = std::exp(x);
  double luminosity_distance = r / a;

  return luminosity_distance;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  double chi = this->eta_of_x(Constants.x_end) - this->eta_of_x(x);
  double r;
  if (this->OmegaK == 0)
    r = chi;
  else {
    double k = std::sqrt(std::abs(this->OmegaK)) * this->H0 * chi / Constants.c;
    if (this->OmegaK < 0) {
      r = chi * std::sin(k)/k;
    }
    else {
      r = chi * std::sinh(k)/k;
    }
  }

  return r;
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  double r = this->get_comoving_distance_of_x(x);

  double a = std::exp(x);
  double angular_distance = a * r;

  return angular_distance;
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                                      << " ";
    fp << t_of_x(x)                              << " ";
    fp << eta_of_x(x)                            << " ";
    fp << Hp_of_x(x)                             << " ";
    fp << dHpdx_of_x(x)                          << " ";
    fp << ddHpddx_of_x(x)                        << " ";
    fp << get_OmegaB(x)                          << " ";
    fp << get_OmegaCDM(x)                        << " ";
    fp << get_OmegaLambda(x)                     << " ";
    fp << get_OmegaR(x)                          << " ";
    fp << get_OmegaNu(x)                         << " ";
    fp << get_OmegaK(x)                          << " ";
    fp << get_luminosity_distance_of_x(x)        << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

