#pragma once
#include <vector>
#include <string>

struct CampoVTK {
  int lines   = 0;
  int columns = 0;
  std::vector<std::vector<double>> x, y;
  std::vector<std::vector<double>> P;
  std::vector<std::vector<double>> u, v;
};

CampoVTK lerVTK(const std::string& caminho);
