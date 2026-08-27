#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>

using namespace std;

// ============================================================
//  Estruturas
// ============================================================

struct NACAParameters
{
    int liftDigit;
    int camberPositionDigit;
    int camberTypeDigit;
    int thicknessPercent;
    double r = 0.0;
    double k1 = 0.0;
    double k2 = 0.0;
    bool isReflex = false;
    bool isValidCamber = false; // true se os 3 primeiros dígitos formam uma linha média tabelada
};

struct Point
{
    double x;
    double y;
};

// Parâmetros de geração da malha em H
struct MeshParameters
{
    double xInlet      = -8.0;  // posição x do plano de entrada (em cordas), corda = [0,1]
    double xOutlet      = 12.0;  // posição x do plano de saída
    double yFarfield    = 8.0;   // distância (em cordas) das fronteiras superior/inferior até y=0

    int nUpstream       = 40;    // pontos entre xInlet e o bordo de ataque (x=0)
    int nDownstream     = 60;    // pontos entre o bordo de fuga (x=1) e xOutlet
    int nAirfoil        = 100;   // pontos ao longo da corda (0 a 1) - mesma malha usada nas superfícies

    int nNormalHalf     = 40;    // pontos em j para CADA metade (inferior e superior) -> jmax = 2*nNormalHalf
    double betaAirfoil  = 2.5;   // fator de clustering perto da superfície do perfil (>1 concentra pontos)
    double betaInlet    = 2.0;   // fator de clustering perto do bordo de ataque, na direção i (upstream)
    double betaOutlet   = 2.0;   // fator de clustering perto do bordo de fuga, na direção i (downstream)
};

// Malha estruturada 2D (i = ao longo do perfil / esteira, j = direção normal)
struct StructuredMesh
{
    int imax = 0;
    int jmax = 0;
    vector<vector<double>> X; // X[i][j]
    vector<vector<double>> Y; // Y[i][j]
};

// ============================================================
//  Protótipos - geração do perfil (código original)
// ============================================================

bool isValidInputFormat(string nacaCode);
NACAParameters parseNACA(string nacaCode);
void getCamberParameters(int liftDigit, int camberPositionDigit, int camberTypeDigit,
                          double& r, double& k1, double& k2, bool& isReflex, bool& isValidCamber);
vector<double> generateXCoordinates(int numberOfPoints);
double calculateThickness(double x, double thickness);
double calculateStandardCamberLine(double x, double r, double k1);
double calculateStandardCamberSlope(double x, double r, double k1);
double calculateReflexCamberLine(double x, double r, double k1, double k2);
double calculateReflexCamberSlope(double x, double r, double k1, double k2);
Point calculateUpperSurfacePoint(double x, double yc, double yt, double theta);
Point calculateLowerSurfacePoint(double x, double yc, double yt, double theta);
vector<Point> generateSurface(vector<double>& xCoordinates, const NACAParameters& p, double thickness, bool upper);

// ============================================================
//  Protótipos - geração da malha em H (nova parte)
// ============================================================

vector<double> generateClusteredFraction(int n, double beta, bool clusterAtStart);
vector<double> buildSpineX(const MeshParameters& mp);
void surfaceYSimplified(double x, const NACAParameters& p, double thickness, double& yUpper, double& yLower);
StructuredMesh generateHMesh(const NACAParameters& p, double thickness, const MeshParameters& mp);
void exportMeshCSV(const StructuredMesh& mesh, const string& filename);
void exportMeshPlot3D(const StructuredMesh& mesh, const string& filename);
void exportMeshVTS(const StructuredMesh& mesh, const string& filename);
void exportMeshTecplotDAT(const StructuredMesh& mesh, const string& filename, const string& title);

// ============================================================
//  Main
// ============================================================

int main()
{
    string nacaCode;
    cout << "Digite o codigo NACA de 5 digitos: ";
    cin >> nacaCode;

    if (!isValidInputFormat(nacaCode))
    {
        cerr << "\nErro: informe exatamente 5 digitos numericos (ex: 23012).\n";
        return 1;
    }

    NACAParameters parameters = parseNACA(nacaCode);

    if (!parameters.isValidCamber)
    {
        cerr << "\nErro: o codigo NACA " << nacaCode << " nao corresponde a nenhuma linha media"
             << " da serie 5 digitos reconhecida.\n"
             << "Combinacoes validas: 210, 220, 230, 240, 250 (padrao) ou 221, 231, 241, 251 (reflexa).\n";
        return 1;
    }

    double thickness = parameters.thicknessPercent / 100.0;

    // ---------- (Mantido) geração e impressão das superfícies do perfil ----------
    vector<double> xCoordinates = generateXCoordinates(100);
    vector<Point> upperSurface = generateSurface(xCoordinates, parameters, thickness, true);
    vector<Point> lowerSurface = generateSurface(xCoordinates, parameters, thickness, false);

    cout << "\nNACA: " << nacaCode
         << "\nLift digit: " << parameters.liftDigit
         << "\nCamber position digit: " << parameters.camberPositionDigit
         << "\nCamber type digit: " << parameters.camberTypeDigit
         << "\nThickness: " << parameters.thicknessPercent << "%"
         << "\nReflex: " << (parameters.isReflex ? "Yes" : "No") << "\n";

    // ---------- (Novo) geração da malha em H ----------
    MeshParameters meshParams; // valores padrão; podem ser ajustados aqui se desejar

    StructuredMesh mesh = generateHMesh(parameters, thickness, meshParams);

    cout << "\nMalha em H gerada:"
         << "\n  imax = " << mesh.imax
         << "\n  jmax = " << mesh.jmax
         << "\n  total de pontos = " << (long long)mesh.imax * mesh.jmax << "\n";

    exportMeshCSV(mesh, "malha_h_" + nacaCode + ".csv");
    exportMeshPlot3D(mesh, "malha_h_" + nacaCode + ".p3d");
    exportMeshVTS(mesh, "malha_h_" + nacaCode + ".vts");
    exportMeshTecplotDAT(mesh, "malha_h_" + nacaCode + ".dat", "Malha H NACA " + nacaCode);

    cout << "\nArquivos gerados:"
         << "\n  malha_h_" << nacaCode << ".csv   (i,j,x,y  - somente para inspecao/plot manual)"
         << "\n  malha_h_" << nacaCode << ".p3d   (Plot3D ASCII completo - Tecplot, Pointwise, ICEM CFD)"
         << "\n  malha_h_" << nacaCode << ".vts   (VTK StructuredGrid - abre direto no ParaView)"
         << "\n  malha_h_" << nacaCode << ".dat   (Tecplot ASCII, formato POINT - visualizadores 2D)\n";

    return 0;
}

// ============================================================
//  Implementação - geração do perfil (código original, inalterado)
// ============================================================

bool isValidInputFormat(string nacaCode)
{
    return nacaCode.length() == 5 && nacaCode.find_first_not_of("0123456789") == string::npos;
}

NACAParameters parseNACA(string nacaCode)
{
    NACAParameters parameters;
    parameters.liftDigit = nacaCode[0] - '0';
    parameters.camberPositionDigit = nacaCode[1] - '0';
    parameters.camberTypeDigit = nacaCode[2] - '0';
    parameters.thicknessPercent = (nacaCode[3] - '0') * 10 + (nacaCode[4] - '0');

    getCamberParameters(parameters.liftDigit, parameters.camberPositionDigit, parameters.camberTypeDigit,
                         parameters.r, parameters.k1, parameters.k2,
                         parameters.isReflex, parameters.isValidCamber);
    return parameters;
}

void getCamberParameters(int liftDigit, int camberPositionDigit, int camberTypeDigit,
                          double& r, double& k1, double& k2, bool& isReflex, bool& isValidCamber)
{
    int camberCode = liftDigit * 100 + camberPositionDigit * 10 + camberTypeDigit;
    r = k1 = k2 = 0.0;
    isReflex = false;
    isValidCamber = true;

    switch (camberCode)
    {
        case 210: r = 0.0580; k1 = 361.400; break;
        case 220: r = 0.1260; k1 = 51.640;  break;
        case 230: r = 0.2025; k1 = 15.957;  break;
        case 240: r = 0.2900; k1 = 6.643;   break;
        case 250: r = 0.3910; k1 = 3.230;   break;

        case 221: r = 0.130; k1 = 51.990; k2 = k1 * 0.000764; isReflex = true; break;
        case 231: r = 0.217; k1 = 15.793; k2 = k1 * 0.00677;  isReflex = true; break;
        case 241: r = 0.318; k1 = 6.520;  k2 = k1 * 0.0303;   isReflex = true; break;
        case 251: r = 0.441; k1 = 3.191;  k2 = k1 * 0.1355;   isReflex = true; break;

        default: isValidCamber = false;
    }
}

vector<double> generateXCoordinates(int numberOfPoints)
{
    vector<double> coordinates;
    double pi = acos(-1.0);

    for (int i = 0; i < numberOfPoints; i++)
    {
        double beta = pi * i / (numberOfPoints - 1);
        coordinates.push_back(0.5 * (1 - cos(beta)));
    }
    return coordinates;
}

double calculateThickness(double x, double thickness)
{
    return 5 * thickness * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * pow(x, 2)
                             + 0.2843 * pow(x, 3) - 0.1015 * pow(x, 4));
}

double calculateStandardCamberLine(double x, double r, double k1)
{
    if (x <= r)
        return (k1 / 6.0) * (pow(x, 3) - 3 * r * pow(x, 2) + pow(r, 2) * (3 - r) * x);
    return (k1 * pow(r, 3) / 6.0) * (1 - x);
}

double calculateStandardCamberSlope(double x, double r, double k1)
{
    if (x <= r)
        return (k1 / 6.0) * (3 * pow(x, 2) - 6 * r * x + pow(r, 2) * (3 - r));
    return -(k1 * pow(r, 3)) / 6.0;
}

double calculateReflexCamberLine(double x, double r, double k1, double k2)
{
    double ratio = k2 / k1;
    if (x <= r)
        return (k1 / 6.0) * (pow(x - r, 3) - ratio * pow(1 - r, 3) * x - pow(r, 3) * x + pow(r, 3));
    return (k1 / 6.0) * (ratio * pow(x - r, 3) - ratio * pow(1 - r, 3) * x - pow(r, 3) * x + pow(r, 3));
}

double calculateReflexCamberSlope(double x, double r, double k1, double k2)
{
    double ratio = k2 / k1;
    if (x <= r)
        return (k1 / 6.0) * (3 * pow(x - r, 2) - ratio * pow(1 - r, 3) - pow(r, 3));
    return (k1 / 6.0) * (3 * ratio * pow(x - r, 2) - ratio * pow(1 - r, 3) - pow(r, 3));
}

Point calculateUpperSurfacePoint(double x, double yc, double yt, double theta)
{
    return { x - yt * sin(theta), yc + yt * cos(theta) };
}

Point calculateLowerSurfacePoint(double x, double yc, double yt, double theta)
{
    return { x + yt * sin(theta), yc - yt * cos(theta) };
}

vector<Point> generateSurface(vector<double>& xCoordinates, const NACAParameters& p, double thickness, bool upper)
{
    vector<Point> surface;

    for (double x : xCoordinates)
    {
        double yc, slope;
        if (p.isReflex)
        {
            yc = calculateReflexCamberLine(x, p.r, p.k1, p.k2);
            slope = calculateReflexCamberSlope(x, p.r, p.k1, p.k2);
        }
        else
        {
            yc = calculateStandardCamberLine(x, p.r, p.k1);
            slope = calculateStandardCamberSlope(x, p.r, p.k1);
        }

        double yt = calculateThickness(x, thickness);
        double theta = atan(slope);

        surface.push_back(upper ? calculateUpperSurfacePoint(x, yc, yt, theta)
                                 : calculateLowerSurfacePoint(x, yc, yt, theta));
    }
    return surface;
}

// ============================================================
//  Implementação - geração da malha em H (nova parte)
// ============================================================

// Gera n valores em [0,1] com clustering geometrico.
// beta > 1 concentra os pontos perto da extremidade escolhida (start ou end).
vector<double> generateClusteredFraction(int n, double beta, bool clusterAtStart)
{
    vector<double> t(n);
    for (int k = 0; k < n; k++)
    {
        double s = (n == 1) ? 0.0 : (double)k / (double)(n - 1); // s em [0,1] uniforme
        double v = pow(s, beta); // concentra perto de s=0 quando beta>1
        t[k] = clusterAtStart ? v : 1.0 - pow(1.0 - s, beta);
    }
    return t;
}

// Monta a "espinha" de coordenadas x ao longo de toda a direcao i:
//   [xInlet ... 0 (bordo de ataque) ... 1 (bordo de fuga) ... xOutlet]
// concentrando pontos perto do perfil (bordo de ataque e bordo de fuga).
vector<double> buildSpineX(const MeshParameters& mp)
{
    vector<double> spine;

    // Trecho upstream: xInlet -> 0, mais denso perto de x=0
    vector<double> tUp = generateClusteredFraction(mp.nUpstream, mp.betaInlet, false);
    for (double t : tUp)
        spine.push_back(mp.xInlet + t * (0.0 - mp.xInlet));

    // Trecho do perfil: 0 -> 1 (espacamento em cosseno, igual ao usado nas superficies)
    // Pula o primeiro ponto (x=0) pois ja foi incluido pelo trecho upstream.
    vector<double> xAirfoil = generateXCoordinates(mp.nAirfoil);
    for (size_t i = 1; i < xAirfoil.size(); i++)
        spine.push_back(xAirfoil[i]);

    // Trecho downstream: 1 -> xOutlet, mais denso perto de x=1
    // Pula o primeiro ponto (x=1) pois ja foi incluido pelo trecho do perfil.
    vector<double> tDown = generateClusteredFraction(mp.nDownstream, mp.betaOutlet, true);
    for (size_t i = 1; i < tDown.size(); i++)
        spine.push_back(1.0 + tDown[i] * (mp.xOutlet - 1.0));

    return spine;
}

// Espessura/camber "simplificados" (sem rotacao por theta) usados apenas
// como espinha da malha: yUpper(x) = yc(x) + yt(x), yLower(x) = yc(x) - yt(x).
// Fora do perfil (x<0 ou x>1) upper e lower colapsam em y=0 -> corte da esteira.
void surfaceYSimplified(double x, const NACAParameters& p, double thickness, double& yUpper, double& yLower)
{
    if (x < 0.0 || x > 1.0)
    {
        yUpper = 0.0;
        yLower = 0.0;
        return;
    }

    double yc;
    if (p.isReflex)
        yc = calculateReflexCamberLine(x, p.r, p.k1, p.k2);
    else
        yc = calculateStandardCamberLine(x, p.r, p.k1);

    double yt = calculateThickness(x, thickness);

    yUpper = yc + yt;
    yLower = yc - yt;
}

// Gera a malha estruturada em H completa.
StructuredMesh generateHMesh(const NACAParameters& p, double thickness, const MeshParameters& mp)
{
    StructuredMesh mesh;

    vector<double> spineX = buildSpineX(mp);
    int imax = (int)spineX.size();
    int nHalf = mp.nNormalHalf;
    int jmax = 2 * nHalf;

    mesh.imax = imax;
    mesh.jmax = jmax;
    mesh.X.assign(imax, vector<double>(jmax, 0.0));
    mesh.Y.assign(imax, vector<double>(jmax, 0.0));

    // Fracoes de clustering reutilizadas em cada coluna (concentradas perto do perfil)
    vector<double> tLower = generateClusteredFraction(nHalf, mp.betaAirfoil, false); // 0 -> 1, denso perto de 1 (perto do perfil)
    vector<double> tUpper = generateClusteredFraction(nHalf, mp.betaAirfoil, true);  // 0 -> 1, denso perto de 0 (perto do perfil)

    for (int i = 0; i < imax; i++)
    {
        double x = spineX[i];

        double yUpperAirfoil, yLowerAirfoil;
        surfaceYSimplified(x, p, thickness, yUpperAirfoil, yLowerAirfoil);

        double yBottom = -mp.yFarfield;
        double yTop = mp.yFarfield;

        // Sub-bloco inferior: de yBottom ate a superficie inferior (ou y=0 fora do perfil)
        for (int j = 0; j < nHalf; j++)
        {
            double t = tLower[j];
            double y = yBottom + t * (yLowerAirfoil - yBottom);
            mesh.X[i][j] = x;
            mesh.Y[i][j] = y;
        }

        // Sub-bloco superior: da superficie superior (ou y=0 fora do perfil) ate yTop
        for (int j = 0; j < nHalf; j++)
        {
            double t = tUpper[j];
            double y = yUpperAirfoil + t * (yTop - yUpperAirfoil);
            mesh.X[i][nHalf + j] = x;
            mesh.Y[i][nHalf + j] = y;
        }
    }

    return mesh;
}

// Exporta em CSV simples: i,j,x,y (facil de importar em Python/Excel/ParaView)
void exportMeshCSV(const StructuredMesh& mesh, const string& filename)
{
    ofstream file(filename);
    file << "i,j,x,y\n";
    for (int i = 0; i < mesh.imax; i++)
        for (int j = 0; j < mesh.jmax; j++)
            file << i << "," << j << "," << mesh.X[i][j] << "," << mesh.Y[i][j] << "\n";
    file.close();
}

// Exporta em formato Plot3D ASCII de bloco unico (3D "achatado": nk=1, Z=0),
// compativel com leitores estruturados como Tecplot, Pointwise e ICEM CFD.
void exportMeshPlot3D(const StructuredMesh& mesh, const string& filename)
{
    ofstream file(filename);
    int imax = mesh.imax, jmax = mesh.jmax, kmax = 1;

    file << 1 << "\n"; // 1 bloco
    file << imax << " " << jmax << " " << kmax << "\n";

    // Ordem: i varia mais rapido, depois j, depois k (padrao Plot3D)
    for (int k = 0; k < kmax; k++)
        for (int j = 0; j < jmax; j++)
            for (int i = 0; i < imax; i++)
                file << mesh.X[i][j] << " ";
    file << "\n";

    for (int k = 0; k < kmax; k++)
        for (int j = 0; j < jmax; j++)
            for (int i = 0; i < imax; i++)
                file << mesh.Y[i][j] << " ";
    file << "\n";

    for (int k = 0; k < kmax; k++)
        for (int j = 0; j < jmax; j++)
            for (int i = 0; i < imax; i++)
                file << 0.0 << " "; // Z = 0 (malha 2D "achatada")
    file << "\n";

    file.close();
}

// Exporta em VTK XML StructuredGrid (.vts) - abre direto no ParaView
// (Open Data File -> .vts), sem precisar configurar dimensoes manualmente.
void exportMeshVTS(const StructuredMesh& mesh, const string& filename)
{
    ofstream file(filename);
    int imax = mesh.imax, jmax = mesh.jmax;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <StructuredGrid WholeExtent=\"0 " << imax - 1 << " 0 " << jmax - 1 << " 0 0\">\n";
    file << "    <Piece Extent=\"0 " << imax - 1 << " 0 " << jmax - 1 << " 0 0\">\n";
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    // Ordem exigida pelo VTK: i varia mais rapido, depois j, depois k
    for (int j = 0; j < jmax; j++)
        for (int i = 0; i < imax; i++)
            file << mesh.X[i][j] << " " << mesh.Y[i][j] << " 0 ";

    file << "\n        </DataArray>\n";
    file << "      </Points>\n";
    file << "    </Piece>\n";
    file << "  </StructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}

// Exporta em formato Tecplot ASCII (.dat), zona estruturada "POINT" - o formato
// mais direto para visualizar a malha 2D (I x J) em Tecplot ou no leitor
// Tecplot do ParaView.
void exportMeshTecplotDAT(const StructuredMesh& mesh, const string& filename, const string& title)
{
    ofstream file(filename);
    int imax = mesh.imax, jmax = mesh.jmax;

    file << "TITLE = \"" << title << "\"\n";
    file << "VARIABLES = \"X\", \"Y\"\n";
    file << "ZONE T=\"malha_h\", I=" << imax << ", J=" << jmax << ", F=POINT\n";

    // Ordem POINT: i varia mais rapido, depois j
    for (int j = 0; j < jmax; j++)
        for (int i = 0; i < imax; i++)
            file << mesh.X[i][j] << " " << mesh.Y[i][j] << "\n";

    file.close();
}
