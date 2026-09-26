// Autor: Davi Guimarães Durval
// Data: 30 de Julho de 2026
// Código C++ que "lê" código NACA de 5 dígitos e fornece os valores de X e Y do aerofólio para plotar.

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

// ---------- Estruturas ----------

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

// ---------- Protótipos ----------

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

// ---------- Main ----------

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
    vector<double> xCoordinates = generateXCoordinates(100);

    vector<Point> upperSurface = generateSurface(xCoordinates, parameters, thickness, true);
    vector<Point> lowerSurface = generateSurface(xCoordinates, parameters, thickness, false);

    cout << "\nNACA: " << nacaCode
         << "\nLift digit: " << parameters.liftDigit
         << "\nCamber position digit: " << parameters.camberPositionDigit
         << "\nCamber type digit: " << parameters.camberTypeDigit
         << "\nThickness: " << parameters.thicknessPercent << "%"
         << "\nReflex: " << (parameters.isReflex ? "Yes" : "No") << "\n";

    cout << "\nUPPER SURFACE\n";
    for (const Point& p : upperSurface)
        cout << "x = " << p.x << " | y = " << p.y << "\n";

    cout << "\nLOWER SURFACE\n";
    for (const Point& p : lowerSurface)
        cout << "x = " << p.x << " | y = " << p.y << "\n";

    return 0;
}

// ---------- Validação da entrada ----------

bool isValidInputFormat(string nacaCode)
{
    return nacaCode.length() == 5 && nacaCode.find_first_not_of("0123456789") == string::npos;
}

// ---------- Parse do código NACA ----------

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

// ---------- Parâmetros da linha média (tabela NACA) ----------

void getCamberParameters(int liftDigit, int camberPositionDigit, int camberTypeDigit,
                          double& r, double& k1, double& k2, bool& isReflex, bool& isValidCamber)
{
    int camberCode = liftDigit * 100 + camberPositionDigit * 10 + camberTypeDigit;
    r = k1 = k2 = 0.0;
    isReflex = false;
    isValidCamber = true;

    switch (camberCode)
    {
        // Linhas médias padrão
        case 210: r = 0.0580; k1 = 361.400; break;
        case 220: r = 0.1260; k1 = 51.640;  break;
        case 230: r = 0.2025; k1 = 15.957;  break;
        case 240: r = 0.2900; k1 = 6.643;   break;
        case 250: r = 0.3910; k1 = 3.230;   break;

        // Linhas médias reflexas
        case 221: r = 0.130; k1 = 51.990; k2 = k1 * 0.000764; isReflex = true; break;
        case 231: r = 0.217; k1 = 15.793; k2 = k1 * 0.00677;  isReflex = true; break;
        case 241: r = 0.318; k1 = 6.520;  k2 = k1 * 0.0303;   isReflex = true; break;
        case 251: r = 0.441; k1 = 3.191;  k2 = k1 * 0.1355;   isReflex = true; break;

        default: isValidCamber = false; // combinação não tabelada pela NACA
    }
}

// ---------- Coordenadas x (espaçamento em cosseno) ----------

vector<double> generateXCoordinates(int numberOfPoints)
{
    vector<double> coordinates;
    double pi = acos(-1.0);

    for (int i = 0; i < numberOfPoints; i++)
    {
        double beta = pi * i / (numberOfPoints - 1);
        coordinates.push_back(0.5 * (1 - cos(beta))); // concentra pontos perto de LE/TE
    }
    return coordinates;
}

// ---------- Distribuição de espessura ----------

double calculateThickness(double x, double thickness)
{
    return 5 * thickness * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * pow(x, 2)
                             + 0.2843 * pow(x, 3) - 0.1015 * pow(x, 4));
}

// ---------- Linha média padrão ----------

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

// ---------- Linha média reflexa ----------

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

// ---------- Pontos das superfícies ----------

Point calculateUpperSurfacePoint(double x, double yc, double yt, double theta)
{
    return { x - yt * sin(theta), yc + yt * cos(theta) };
}

Point calculateLowerSurfacePoint(double x, double yc, double yt, double theta)
{
    return { x + yt * sin(theta), yc - yt * cos(theta) };
}

// ---------- Geração de superfície (superior ou inferior) ----------

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
