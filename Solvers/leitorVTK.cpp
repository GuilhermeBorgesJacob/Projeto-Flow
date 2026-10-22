#include "leitorVTK.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;

CampoVTK lerVTK(const string& caminho) {
  ifstream arq(caminho);
  if (!arq.is_open())
    throw runtime_error("Não foi possível abrir");

  CampoVTK c;
  string linha;

  // Cabeçalho
  for (int i = 0; i < 4; ++i) getline(arq, linha);

  while (getline(arq, linha)) {
    if (linha.empty()) continue;
    if (linha.rfind("DIMENSIONS", 0) == 0) {
      istringstream ss(linha);
      string kw;
      int k;
      ss >> kw >> c.columns >> c.lines >> k;
      break;
    }
  }

  if (c.columns <= 0 || c.lines <= 0)  
    throw runtime_error("DIMENSIONS inválido no VTK");

  const int N = c.columns * c.lines;

  auto alocar = [&]() {
    return vector<vector<double>>(
        c.lines, vector<double>(c.columns, 0.0));
  };
  c.x = alocar(); c.y = alocar();
  c.P = alocar(); c.u = alocar(); c.v = alocar();

  // PONTOS
  while (getline(arq, linha)) {
    if (linha.empty()) continue;
    if (linha.rfind("POINTS", 0) == 0) break;
  }
  for (int idx = 0; idx < N; ++idx) {
    double xv, yv, zv;
    arq >> xv >> yv >> zv;
    int i = idx / c.columns;  // linha
    int j = idx % c.columns;  // coluna
    c.x[i][j] = xv;
    c.y[i][j] = yv;
  }

  // PRESSÃO ESCALAR
  while (getline(arq, linha)) {
    if (linha.rfind("SCALARS", 0) == 0 && 
      linha.rfind("Pressure") != std::string::npos) {
        getline(arq, linha);
        break;
    }
  }
  for (int idx = 0; idx < N; ++idx) {
    double pv;
    arq >> pv;
    c.P[idx / c.columns][idx % c.columns] = pv;
  }

  // VELOCIDADE VETORIAL
  while (getline(arq, linha)) {
    if (linha.rfind("VECTORS", 0) == 0) break;
  }
  for (int idx = 0; idx < N; ++idx) {
    double uv, vv, wv;
    arq >> uv >> vv >> wv;
    int i = idx / c.columns;
    int j = idx % c.columns;
    c.u[i][j] = uv;
    c.v[i][j] = vv;
  }

  return c;
}
