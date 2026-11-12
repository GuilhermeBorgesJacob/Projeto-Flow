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
    
    // Considerando profundidade unitária
    double ds = std::sqrt(dx*dx + dy*dy);
    if (ds == 0.0) continue; 

    double normal_x = -dy / ds;
    double normal_y =  dx / ds;

    // CP médio no painel
    double Cp_med = 0.5 * (Cp[k] + Cp[k+1]);

    double dFx = -Cp_med * normal_x * ds;
    double dFy = -Cp_med * normal_y * ds;

    // Assumindo ângulo de ataque = 0
    Cd += dFx / corda;
    Cl += dFy / corda;

    // Momento em relação à origem (usando o centroide do painel)
    double x_med = 0.5 * (xc[k] + xc[k+1]);
    double y_med = 0.5 * (yc[k] + yc[k+1]);
    
    Cm += (x_med * dFy - y_med * dFx) / pow(corda, 2.0); 
  } 
  return {Cl, Cd, Cm};
}
