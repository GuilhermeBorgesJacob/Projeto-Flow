#include <iostream>
#include <exception>
#include <iomanip>
#include "solverInviscidOgridV2.h"
#include "leitorVTK.h"
#include "Coeficientes.h"

int main() {
  try {
    CampoVTK c = lerVTK("malha3.vtk");
    std::cout << "Malha: " << c.columns << " x " << c.lines << "\n";

    double Pmin = 1e300, Pmax = -1e300;
    double umin = 1e300, umax = -1e300;
    double vmin = 1e300, vmax = -1e300;
    for (int i = 0; i < c.lines; ++i)
      for (int j = 0; j < c.columns; ++j) {
        Pmin = std::min(Pmin, c.P[i][j]); Pmax = std::max(Pmax, c.P[i][j]);
        umin = std::min(umin, c.u[i][j]); umax = std::max(umax, c.u[i][j]);
        vmin = std::min(vmin, c.v[i][j]); vmax = std::max(vmax, c.v[i][j]);
    }
    std::cout << "P  ∈ [" << Pmin << ", " << Pmax << "]\n";
    std::cout << "u  ∈ [" << umin << ", " << umax << "]\n";
    std::cout << "v  ∈ [" << vmin << ", " << vmax << "]\n";
    
    const double rho   = 1.2;
    const double Vinf  = 20.0;
    const double Pinf  = 1.0e5;
    const double corda = 1.0;

    Coeficientes coef = calcularCoeficientes(c, rho, Vinf, Pinf, corda);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Cl = " << coef.Cl << "\n";
    std::cout << "Cd = " << coef.Cd << "\n";
    std::cout << "Cm = " << coef.Cm << "\n";

  } catch (const std::exception& e) {
    std::cerr << "Erro: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
