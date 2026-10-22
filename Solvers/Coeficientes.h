#pragma once 
#include "leitorVTK.h"

struct Coeficientes { 
  double Cl;
  double Cd;
  double Cm; // momento em relação ao bordo de ataque (x=0)
};

Coeficientes calcularCoeficientes(const CampoVTK& c, 
                                  double rho,
                                  double Vinf,
                                  double Pinf, 
                                  double corda);
