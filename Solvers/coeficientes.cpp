#include "Coeficientes.h"
#include <iostream>
#include <cmath>
#include <vector>

Coeficientes calcularCoeficientes(const CampoVTK& c, 
                                  double rho,
                                  double Vinf,
                                  double Pinf, 
                                  double corda) {

  const double q_inf = 0.5 * rho * Vinf * Vinf; 

  // Fronteira sul
  const int i = c.lines - 1;

  // Intradorso: j = 0, 1, 2, ... , jLE
  // Extradorso: j = jLE, ..., columns-1
  
  const int jLE = (c.columns + 1) / 2;

  std::vector<double> xc, yc, Cp;

  // Intradorso: TE -> LE 
  for (int j = 0; j<=jLE; ++j) {
    double Cpv = (c.P[i][j] - Pinf) / q_inf;
    xc.push_back(c.x[i][j]);
    yc.push_back(c.y[i][j]);
    Cp.push_back(Cpv);
  }
  // Extradorso: LE -> TE 
  for (int j = jLE + 1; j < c.columns; ++j) {
    double Cpv = (c.P[i][j] - Pinf) / q_inf;
    xc.push_back(c.x[i][j]);
    yc.push_back(c.y[i][j]);
    Cp.push_back(Cpv);
  }
  // Fechar contorno
  xc.push_back(c.x[i][0]);
  yc.push_back(c.y[i][0]);
  Cp.push_back((c.P[i][0] - Pinf) / q_inf);

  const int n = static_cast<int>(xc.size());

  double Cl = 0.0, Cd = 0.0, Cm = 0.0;

  for (int k = 0; k < n - 1; ++k) {
    double dx = xc[k+1] - xc[k];
    double dy = yc[k+1] - yc[k];
    double Cp_med = 0.5 * (Cp[k] + Cp[k+1]);

    Cl += -Cp_med * dx / corda * (-1.0);
    Cd += -Cp_med * dy / corda * (-1.0);

    double dFy = Cp_med * dx;
    double dFx = Cp_med * dy;
    Cm += (xc[k] * dFy - yc[k] * dFx) / corda;
  }
  
  return {Cl, Cd, Cm};
}
