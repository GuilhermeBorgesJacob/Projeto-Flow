// ============================================================================
//  Solver 2D compressivel RANS k-omega SST-2003 (ou laminar) - VOLUME FINITO
//    * Fluxo convectivo: AUSM+ (Liou 1996) com reconstrucao MUSCL (van Albada)
//    * Viscoso: gradientes GG/LSQ com correcao na aresta; SST-2003 com parede resolvida
//    * Tempo: Runge-Kutta SSP de 3 estagios, passo global (transiente) ou local
//    * Malhas estruturadas de 1 bloco (tipos O, H, C, ...) - condicoes de contorno
//      por segmento de lado (imin/imax/jmin/jmax), ver parametros.txt
//    * Parametros em arquivo (nao precisa recompilar)
//    * Saida: um arquivo por instante (.vts XML ou .vtk legado) em 'saidas/' + solucao.pvd (indice da serie)
//
//  Compilar: g++ -O3 -march=native -fopenmp -std=c++17 rans_fv_sst.cpp -o rans_sst
//  Executar: OMP_NUM_THREADS=<n> ./rans_sst parametros.txt
//  NAO usar -ffast-math.
// ============================================================================
#define _USE_MATH_DEFINES
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <chrono>
#include <map>
#include <regex>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <omp.h>

constexpr int G = 2;     // camadas de celulas fantasma (MUSCL de 2a ordem)
constexpr int NV = 6;    // rho, rho u, rho v, rho E, rho k, rho omega
enum BCType : uint8_t { BC_NONE = 0, BC_WALL, BC_FARFIELD, BC_OUTFLOW, BC_SYMMETRY, BC_PERIODIC, BC_CUT };
enum Side { IMIN = 0, IMAX = 1, JMIN = 2, JMAX = 3 };
static const char* SIDE_NAME[4] = {"imin", "imax", "jmin", "jmax"};

// -------------------------- modelo Menter SST-2003 -------------------------
// Constantes internas (1: k-omega) e externas (2: k-epsilon), com blending F1.
struct SST {
    static constexpr double betaStar = 0.09, a1 = 0.31;
    static constexpr double beta1 = 0.075, beta2 = 0.0828;
    static constexpr double sigmaK1 = 0.85, sigmaK2 = 1.0;
    static constexpr double sigmaW1 = 0.5, sigmaW2 = 0.856;
    static constexpr double gamma1 = 5.0 / 9.0, gamma2 = 0.44;
    static constexpr double Prt = 0.9;
    static constexpr double kMin = 1e-10, omegaMin = 1e-10, muRatioMax = 1e5;
    static constexpr double productionLimit = 10.0; // SST-2003: Pk <= 10 beta* rho k omega
};

// ------------------------------ parametros ---------------------------------
struct Params {
    std::string grid_file = "malha.grid", grid_format = "auto", output_dir = "saidas", time_mode = "global";
    double mesh_scale = 1.0;
    double gamma = 1.4, R_gas = 287.0, mu_ref = 1.7894e-5, T_ref = 288.16, S_suth = 110.0, Pr = 0.71;
    double mach = 0.3, alpha_deg = 0.0, p_inf = 101325.0, T_inf = 288.2, T_wall = -1.0, p_back = -1.0;
    int turbulence = 1;
    double Tu = 1e-4, visc_ratio = 1e-2; // ponto de partida para aero externo; calibrar com Re e condicoes de tunel
    int order = 2;
    double cfl = 0.7, cfl_start = 0.2;
    long cfl_ramp = 1000;
    long max_iter = 100000;
    double t_end_star = -1.0;
    int report_every = 500;
    long write_every_iter = 0;
    double write_every_time = 0.1;     // [s] grava um arquivo a cada este tempo FISICO (time_mode=global)
    double write_every_tstar = 0.0;    // alternativa em tempo adimensional t* = t*U/l_ref
    std::string output_format = "vts"; // vts (XML, recomendado) | vtk (legado ASCII)
    double l_ref = -1.0;
    int test_freestream = 0;
    int turb_order = 1;        // ordem espacial de k e omega: 1 = upwind (recomendado: preserva positividade) | 2
    std::string gradient = "lsq";   // lsq = minimos quadrados (recomendado em malhas enviesadas) | gg = Green-Gauss
};

struct RawSeg { std::string side, range, type; };

static void validateParams(const Params& P) {
    auto bad = [](const std::string& m) { throw std::runtime_error("Parametro invalido: " + m); };
    if (!(P.gamma > 1.0)) bad("gamma deve ser > 1");
    if (!(P.R_gas > 0) || !(P.p_inf > 0) || !(P.T_inf > 0)) bad("R_gas, p_inf e T_inf devem ser > 0");
    if (!(P.mu_ref > 0) || !(P.T_ref > 0) || P.S_suth < 0 || !(P.Pr > 0)) bad("mu_ref, T_ref, Pr devem ser > 0 e S_suth >= 0");
    if (!(P.visc_ratio > 0) || !std::isfinite(P.visc_ratio)) bad("visc_ratio deve ser finito e > 0 (mu_t/mu no escoamento livre)");
    if (!(P.Tu >= 0.0) || !std::isfinite(P.Tu)) bad("Tu deve ser finito e >= 0");
    if (!(P.mach > 0)) bad("mach deve ser > 0");
    if (P.turbulence != 0 && P.turbulence != 1) bad("turbulence deve ser 0 ou 1");
    if (P.order != 1 && P.order != 2) bad("order deve ser 1 ou 2");
    if (P.turb_order != 1 && P.turb_order != 2) bad("turb_order deve ser 1 ou 2");
    if (P.gradient != "lsq" && P.gradient != "gg") bad("gradient deve ser lsq ou gg");
    if (!(P.cfl > 0) || !(P.cfl_start > 0) || P.cfl_ramp < 0) bad("cfl e cfl_start devem ser > 0");
    if (P.time_mode != "global" && P.time_mode != "local") bad("time_mode deve ser global ou local");
    if (P.max_iter < 1) bad("max_iter deve ser >= 1");
    if (P.report_every < 1) bad("report_every deve ser >= 1");
    if (P.write_every_time < 0 || P.write_every_tstar < 0 || P.write_every_iter < 0) bad("intervalos de gravacao devem ser >= 0");
}

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}
static BCType parseBC(const std::string& s) {
    if (s == "wall") return BC_WALL;
    if (s == "farfield") return BC_FARFIELD;
    if (s == "outflow") return BC_OUTFLOW;
    if (s == "symmetry") return BC_SYMMETRY;
    if (s == "periodic") return BC_PERIODIC;
    if (s == "cut") return BC_CUT;
    throw std::runtime_error("Tipo de contorno desconhecido: '" + s + "' (use wall, farfield, outflow, symmetry, periodic, cut)");
}
static int parseSide(const std::string& s) {
    for (int k = 0; k < 4; k++) if (s == SIDE_NAME[k]) return k;
    throw std::runtime_error("Lado desconhecido: '" + s + "' (use imin, imax, jmin, jmax)");
}

static Params readParams(const std::string& file, std::vector<RawSeg>& raw) {
    std::ifstream f(file);
    if (!f) throw std::runtime_error("Nao foi possivel abrir o arquivo de parametros: " + file);
    Params P;
    std::string line;
    int ln = 0;
    while (std::getline(f, line)) {
        ln++;
        size_t h = line.find('#'); if (h != std::string::npos) line = line.substr(0, h);
        line = trim(line); if (line.empty()) continue;
        if (line.rfind("boundary", 0) == 0) {
            std::istringstream ss(line); std::string w, a, b, c;
            ss >> w >> a >> b >> c;
            if (c.empty()) throw std::runtime_error("Linha " + std::to_string(ln) + ": use 'boundary <lado> <a:b|all> <tipo>'");
            raw.push_back({a, b, c});
            continue;
        }
        size_t eq = line.find('=');
        if (eq == std::string::npos) throw std::runtime_error("Linha " + std::to_string(ln) + " invalida: " + line);
        std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
        auto D = [&](double& x) { x = std::stod(v); };
        auto I = [&](int& x) { x = std::stoi(v); };
        auto L = [&](long& x) { x = std::stol(v); };
        if (k == "grid_file") P.grid_file = v;
        else if (k == "grid_format") P.grid_format = v;
        else if (k == "output_dir") P.output_dir = v;
        else if (k == "time_mode") P.time_mode = v;
        else if (k == "mesh_scale") D(P.mesh_scale);
        else if (k == "gamma") D(P.gamma);
        else if (k == "R_gas") D(P.R_gas);
        else if (k == "mu_ref") D(P.mu_ref);
        else if (k == "T_ref") D(P.T_ref);
        else if (k == "S_suth") D(P.S_suth);
        else if (k == "Pr") D(P.Pr);
        else if (k == "mach") D(P.mach);
        else if (k == "alpha_deg") D(P.alpha_deg);
        else if (k == "p_inf") D(P.p_inf);
        else if (k == "T_inf") D(P.T_inf);
        else if (k == "T_wall") D(P.T_wall);
        else if (k == "p_back") D(P.p_back);
        else if (k == "turbulence") I(P.turbulence);
        else if (k == "Tu") D(P.Tu);
        else if (k == "visc_ratio") D(P.visc_ratio);
        else if (k == "order") I(P.order);
        else if (k == "cfl") D(P.cfl);
        else if (k == "cfl_start") D(P.cfl_start);
        else if (k == "cfl_ramp") L(P.cfl_ramp);
        else if (k == "max_iter") L(P.max_iter);
        else if (k == "t_end_star") D(P.t_end_star);
        else if (k == "report_every") I(P.report_every);
        else if (k == "write_every_iter") L(P.write_every_iter);
        else if (k == "write_every_time") D(P.write_every_time);
        else if (k == "write_every_tstar") D(P.write_every_tstar);
        else if (k == "output_format") { P.output_format = v; if (v != "vts" && v != "vtk") throw std::runtime_error("output_format deve ser vts ou vtk"); }
        else if (k == "l_ref") D(P.l_ref);
        else if (k == "test_freestream") I(P.test_freestream);
        else if (k == "turb_order") I(P.turb_order);
        else if (k == "gradient") P.gradient = v;
        else std::cerr << "AVISO: parametro desconhecido ignorado: '" << k << "' (linha " << ln << ")\n";
    }
    return P;
}

// ------------------------------- leitura de malha --------------------------
// Formatos: Tecplot ASCII (POINT ou BLOCK; I varia mais rapido) e Plot3D ASCII 2D (1 bloco).
static void readGrid(const Params& P, int& NI, int& NJ, std::vector<double>& X, std::vector<double>& Y) {
    std::ifstream f(P.grid_file);
    if (!f) throw std::runtime_error("Nao foi possivel abrir a malha: " + P.grid_file);
    std::stringstream buf; buf << f.rdbuf();
    std::string text = buf.str(), up = text;
    for (auto& c : up) c = std::toupper(static_cast<unsigned char>(c));
    std::string fmt = P.grid_format;
    if (fmt == "auto") fmt = (up.find("ZONE") != std::string::npos) ? "tecplot" : "plot3d";
    size_t nn = 0;
    if (fmt == "tecplot") {
        size_t pz = up.find("ZONE");
        if (pz == std::string::npos) throw std::runtime_error("Malha Tecplot sem linha ZONE");
        size_t eol = up.find('\n', pz);
        std::string hdr = up.substr(pz, eol - pz);
        std::smatch m;
        std::regex ri(R"(I\s*=\s*(\d+))"), rj(R"(J\s*=\s*(\d+))");
        if (!std::regex_search(hdr, m, ri)) throw std::runtime_error("ZONE sem I=");
        NI = std::stoi(m[1]);
        if (!std::regex_search(hdr, m, rj)) throw std::runtime_error("ZONE sem J=");
        NJ = std::stoi(m[1]);
        bool block = hdr.find("BLOCK") != std::string::npos;
        std::istringstream ss(text.substr(eol + 1));
        nn = static_cast<size_t>(NI) * NJ;
        X.assign(nn, 0); Y.assign(nn, 0);
        if (block) { for (size_t k = 0; k < nn; k++) ss >> X[k]; for (size_t k = 0; k < nn; k++) ss >> Y[k]; }
        else { for (size_t k = 0; k < nn; k++) ss >> X[k] >> Y[k]; }
        if (!ss) throw std::runtime_error("Malha Tecplot: menos pontos que I*J");
    } else {
        std::istringstream ss(text);
        std::string first; std::getline(ss, first);
        std::istringstream l1(first); std::vector<long> tok; long v;
        while (l1 >> v) tok.push_back(v);
        if (tok.size() == 1) { ss >> NI >> NJ; }       // formato com numero de blocos na 1a linha
        else if (tok.size() >= 2) { NI = static_cast<int>(tok[0]); NJ = static_cast<int>(tok[1]); }
        else throw std::runtime_error("Plot3D: cabecalho invalido");
        nn = static_cast<size_t>(NI) * NJ;
        X.assign(nn, 0); Y.assign(nn, 0);
        for (size_t k = 0; k < nn; k++) ss >> X[k];
        for (size_t k = 0; k < nn; k++) ss >> Y[k];
        if (!ss) throw std::runtime_error("Malha Plot3D: menos pontos que NI*NJ");
    }
    if (NI < 4 || NJ < 4) throw std::runtime_error("Malha pequena demais");
    for (size_t k = 0; k < nn; k++) { X[k] *= P.mesh_scale; Y[k] *= P.mesh_scale; }
}

// --------------------------- AUSM+ (Liou 1996) -----------------------------
static inline double M1p(double M) { return 0.5 * (M + std::fabs(M)); }
static inline double M1m(double M) { return 0.5 * (M - std::fabs(M)); }
static inline double M2p(double M) { return 0.25 * (M + 1.0) * (M + 1.0); }
static inline double M2m(double M) { return -0.25 * (M - 1.0) * (M - 1.0); }
static inline double M4p(double M) { return std::fabs(M) >= 1.0 ? M1p(M) : M2p(M) * (1.0 - 2.0 * M2m(M)); }   // 16*beta = 2
static inline double M4m(double M) { return std::fabs(M) >= 1.0 ? M1m(M) : M2m(M) * (1.0 + 2.0 * M2p(M)); }
static inline double P5p(double M) { return std::fabs(M) >= 1.0 ? (M > 0 ? 1.0 : 0.0) : M2p(M) * ((2.0 - M) - 3.0 * M * M2m(M)); }   // 16*alpha = 3
static inline double P5m(double M) { return std::fabs(M) >= 1.0 ? (M > 0 ? 0.0 : 1.0) : M2m(M) * ((-2.0 - M) + 3.0 * M * M2p(M)); }
static inline double vanAlbada(double a, double b, double e2) { return (a * (b * b + e2) + b * (a * a + e2)) / (a * a + b * b + 2.0 * e2); }

// ================================ solver ====================================
struct Seg { int side, a, b; BCType type; int pside = -1, pa = 0, pb = 0; bool rev = false, paired = false; };
struct WallFace {
    int ic, jc, axis, fid;
    double sigma, nx, ny, A, y1, xf, yf;
    double xa = 0, ya = 0, xb = 0, yb = 0;
    double Fx = 0, Fy = 0, tau = 0, yplus = 0;
};
struct WallTreeNode {
    double xmin, xmax, ymin, ymax;
    int begin = 0, end = 0, left = -1, right = -1;
};

struct Solver {
    Params P;
    // --- malha ---
    int NI = 0, NJ = 0, nci = 0, ncj = 0, sx = 0, N = 0;
    std::vector<double> xn, yn, xc, yc, vol, Six, Siy, Sjx, Sjy, xfi, yfi, xfj, yfj;
    inline int c(int i, int j) const { return (j + G) * sx + (i + G); }
    inline int fi(int f, int j) const { return j * (nci + 1) + f; }
    inline int fj(int i, int f) const { return f * nci + i; }
    // --- contorno ---
    std::vector<Seg> segs;
    std::vector<uint8_t> btI, btJ;
    std::vector<int> wI, wJ;
    std::vector<WallFace> walls;
    // --- estado ---
    std::array<std::vector<double>, NV> U, U0, Ub, Lr, Fi, Fj;
    std::vector<double> rho, u, v, p, T, k, omega, mu, mut;
    std::array<std::vector<double>, 7> gx, gy;   // grad de u,v,T,k,omega,p,rho
    std::vector<double> F1, F2, sigmaK, sigmaOmega, gammaOmega, betaOmega, wallDist;
    std::vector<double> sinkK, sinkOmega, dtc;
    std::vector<int> wallTreeOrder;
    std::vector<WallTreeNode> wallTree;
    // --- constantes ---
    double Rg, gam, Cv, Cp, rhoI, uI, vI, pI, TI, kI, omegaI, Uinf, muI, rhoMin, eMin, lref;
    long clipsTurb = 0, backoffs = 0;
    double resRef[NV] = {0, 0, 0, 0, 0, 0}, resNow[NV] = {0, 0, 0, 0, 0, 0};
    double resScale[NV] = {1, 1, 1, 1, 1, 1}, tRef = 1.0;

    // ------------------------------------------------------------------ malha
    void buildMesh(int NI_, int NJ_, const std::vector<double>& X, const std::vector<double>& Y) {
        NI = NI_; NJ = NJ_; nci = NI - 1; ncj = NJ - 1; sx = nci + 2 * G; N = sx * (ncj + 2 * G);
        xn = X; yn = Y;
        auto nx_ = [&](int i, int j) { return xn[static_cast<size_t>(j) * NI + i]; };
        auto ny_ = [&](int i, int j) { return yn[static_cast<size_t>(j) * NI + i]; };
        double sumA = 0;
        std::vector<double> Asg(static_cast<size_t>(nci) * ncj);
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            double x0 = nx_(i, j), y0 = ny_(i, j), x1 = nx_(i + 1, j), y1 = ny_(i + 1, j);
            double x2 = nx_(i + 1, j + 1), y2 = ny_(i + 1, j + 1), x3 = nx_(i, j + 1), y3 = ny_(i, j + 1);
            double A = 0.5 * ((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
            Asg[static_cast<size_t>(j) * nci + i] = A; sumA += A;
        }
        double sgn = sumA >= 0 ? 1.0 : -1.0;
        xc.assign(N, 0); yc.assign(N, 0); vol.assign(N, 1.0);
        int folded = 0, nonconv = 0;
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            double A = Asg[static_cast<size_t>(j) * nci + i] * sgn;
            if (!(A > 0)) folded++;
            double px[4] = {nx_(i, j), nx_(i + 1, j), nx_(i + 1, j + 1), nx_(i, j + 1)};
            double py[4] = {ny_(i, j), ny_(i + 1, j), ny_(i + 1, j + 1), ny_(i, j + 1)};
            for (int m = 0; m < 4; m++) {            // celula nao-convexa: algum canto vira para o lado oposto
                int m0 = (m + 3) % 4, m2 = (m + 1) % 4;
                if (!(((px[m] - px[m0]) * (py[m2] - py[m]) - (py[m] - py[m0]) * (px[m2] - px[m])) * sgn > 0)) { nonconv++; break; }
            }
            int q = c(i, j);
            vol[q] = std::fabs(A);
            // centroide geometrico (ponderado pela area de 2 triangulos), nao a media dos vertices
            double a1 = 0.5 * ((px[1] - px[0]) * (py[2] - py[0]) - (py[1] - py[0]) * (px[2] - px[0]));
            double a2 = 0.5 * ((px[2] - px[0]) * (py[3] - py[0]) - (py[2] - py[0]) * (px[3] - px[0]));
            double at = a1 + a2;
            if (std::fabs(at) > 0) {
                xc[q] = (a1 * (px[0] + px[1] + px[2]) + a2 * (px[0] + px[2] + px[3])) / (3.0 * at);
                yc[q] = (a1 * (py[0] + py[1] + py[2]) + a2 * (py[0] + py[2] + py[3])) / (3.0 * at);
            } else { xc[q] = 0.25 * (px[0] + px[1] + px[2] + px[3]); yc[q] = 0.25 * (py[0] + py[1] + py[2] + py[3]); }
        }
        if (nonconv) std::cout << "AVISO: " << nonconv << " celula(s) nao-convexa(s) na malha (a precisao local diminui; confira a malha).\n";
        if (folded) throw std::runtime_error("Malha com " + std::to_string(folded) + " celula(s) dobrada(s)/degenerada(s) (area com sinal diferente). Corrija a malha.");
        size_t nfi = static_cast<size_t>(nci + 1) * ncj, nfj = static_cast<size_t>(nci) * (ncj + 1);
        Six.assign(nfi, 0); Siy.assign(nfi, 0); xfi.assign(nfi, 0); yfi.assign(nfi, 0);
        Sjx.assign(nfj, 0); Sjy.assign(nfj, 0); xfj.assign(nfj, 0); yfj.assign(nfj, 0);
        for (int j = 0; j < ncj; j++) for (int f = 0; f <= nci; f++) {
            double dx = nx_(f, j + 1) - nx_(f, j), dy = ny_(f, j + 1) - ny_(f, j);
            int id = fi(f, j);
            Six[id] = sgn * dy; Siy[id] = -sgn * dx;
            xfi[id] = 0.5 * (nx_(f, j + 1) + nx_(f, j)); yfi[id] = 0.5 * (ny_(f, j + 1) + ny_(f, j));
        }
        for (int f = 0; f <= ncj; f++) for (int i = 0; i < nci; i++) {
            double dx = nx_(i + 1, f) - nx_(i, f), dy = ny_(i + 1, f) - ny_(i, f);
            int id = fj(i, f);
            Sjx[id] = -sgn * dy; Sjy[id] = sgn * dx;
            xfj[id] = 0.5 * (nx_(i + 1, f) + nx_(i, f)); yfj[id] = 0.5 * (ny_(i + 1, f) + ny_(i, f));
        }
        // Centros fantasmas refletidos em relacao a LINHA da face (nao ao ponto medio).
        auto reflectLine = [](double x, double y, double ax, double ay, double bx, double by) {
            const double dx = bx - ax, dy = by - ay, l2 = dx * dx + dy * dy;
            const double t = l2 > 0.0 ? ((x - ax) * dx + (y - ay) * dy) / l2 : 0.0;
            const double px = ax + t * dx, py = ay + t * dy;
            return std::array<double, 2>{2.0 * px - x, 2.0 * py - y};
        };
        for (int t = 0; t < ncj; ++t) {
            auto L = reflectLine(xc[c(0, t)], yc[c(0, t)], nx_(0, t), ny_(0, t), nx_(0, t + 1), ny_(0, t + 1));
            auto R = reflectLine(xc[c(nci - 1, t)], yc[c(nci - 1, t)], nx_(nci, t), ny_(nci, t), nx_(nci, t + 1), ny_(nci, t + 1));
            xc[c(-1, t)] = L[0]; yc[c(-1, t)] = L[1];
            xc[c(nci, t)] = R[0]; yc[c(nci, t)] = R[1];
        }
        for (int t = 0; t < nci; ++t) {
            auto B = reflectLine(xc[c(t, 0)], yc[c(t, 0)], nx_(t, 0), ny_(t, 0), nx_(t + 1, 0), ny_(t + 1, 0));
            auto T = reflectLine(xc[c(t, ncj - 1)], yc[c(t, ncj - 1)], nx_(t, ncj), ny_(t, ncj), nx_(t + 1, ncj), ny_(t + 1, ncj));
            xc[c(t, -1)] = B[0]; yc[c(t, -1)] = B[1];
            xc[c(t, ncj)] = T[0]; yc[c(t, ncj)] = T[1];
        }
        // diagnostico GCL: soma das areas orientadas de cada celula deve ser ~0
        double gclMax = 0;
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            double sx_ = Six[fi(i + 1, j)] - Six[fi(i, j)] + Sjx[fj(i, j + 1)] - Sjx[fj(i, j)];
            double sy_ = Siy[fi(i + 1, j)] - Siy[fi(i, j)] + Sjy[fj(i, j + 1)] - Sjy[fj(i, j)];
            double ref = std::hypot(Six[fi(i, j)], Siy[fi(i, j)]) + 1e-300;
            gclMax = std::max(gclMax, std::hypot(sx_, sy_) / ref);
        }
        std::cout << "Malha: " << NI << " x " << NJ << " nos (" << nci << " x " << ncj << " celulas), orientacao "
                  << (sgn > 0 ? "direita" : "esquerda (tratada)") << ", erro GCL max = " << gclMax << "\n";
    }

    // -------------------------------------------------------------- contornos
    void interiorCell(int side, int t, int m, int& i, int& j) const {
        switch (side) {
            case IMIN: i = m; j = t; break;
            case IMAX: i = nci - 1 - m; j = t; break;
            case JMIN: i = t; j = m; break;
            default:   i = t; j = ncj - 1 - m; break;
        }
    }
    void ghostCell(int side, int t, int m, int& i, int& j) const {
        switch (side) {
            case IMIN: i = -1 - m; j = t; break;
            case IMAX: i = nci + m; j = t; break;
            case JMIN: i = t; j = -1 - m; break;
            default:   i = t; j = ncj + m; break;
        }
    }
    int sideLength(int side) const { return (side == IMIN || side == IMAX) ? ncj : nci; }
    void boundaryFace(int side, int t, double& Sx, double& Sy, int& fid, int& axis, double& sigma) const {
        switch (side) {
            case IMIN: fid = fi(0, t);   Sx = Six[fid]; Sy = Siy[fid]; axis = 0; sigma = -1; break;
            case IMAX: fid = fi(nci, t); Sx = Six[fid]; Sy = Siy[fid]; axis = 0; sigma = +1; break;
            case JMIN: fid = fj(t, 0);   Sx = Sjx[fid]; Sy = Sjy[fid]; axis = 1; sigma = -1; break;
            default:   fid = fj(t, ncj); Sx = Sjx[fid]; Sy = Sjy[fid]; axis = 1; sigma = +1; break;
        }
    }

    void setupBoundaries(const std::vector<RawSeg>& raw) {
        if (raw.empty()) throw std::runtime_error("Nenhuma linha 'boundary' no arquivo de parametros.");
        for (const auto& r : raw) {
            Seg s; s.side = parseSide(r.side); s.type = parseBC(r.type);
            int len = sideLength(s.side);
            if (r.range == "all") { s.a = 0; s.b = len - 1; }
            else {
                size_t col = r.range.find(':');
                if (col == std::string::npos) throw std::runtime_error("Faixa invalida '" + r.range + "' (use a:b, 'end' ou 'all')");
                auto val = [&](const std::string& x) { return x == "end" ? len - 1 : std::stoi(x); };
                s.a = val(r.range.substr(0, col)); s.b = val(r.range.substr(col + 1));
            }
            if (s.a < 0 || s.b >= len || s.a > s.b) throw std::runtime_error(std::string("Faixa fora do lado ") + SIDE_NAME[s.side] + " (celulas 0.." + std::to_string(len - 1) + ")");
            segs.push_back(s);
        }
        // cobertura completa de cada lado
        for (int sd = 0; sd < 4; sd++) {
            std::vector<int> cov(sideLength(sd), 0);
            for (auto& s : segs) if (s.side == sd) for (int t = s.a; t <= s.b; t++) cov[t]++;
            for (int t = 0; t < static_cast<int>(cov.size()); t++)
                if (cov[t] != 1) throw std::runtime_error(std::string("Lado ") + SIDE_NAME[sd] + ": celula " + std::to_string(t) + (cov[t] == 0 ? " sem contorno definido" : " definida mais de uma vez"));
        }
        // pareamento periodic (imin<->imax, jmin<->jmax) e cut (mesmo lado, reverso)
        for (auto& s : segs) {
            if (s.type == BC_PERIODIC && !s.paired) {
                int opp = s.side ^ 1;
                for (auto& q : segs) if (q.side == opp && q.type == BC_PERIODIC && !q.paired && (q.b - q.a) == (s.b - s.a)) {
                    s.pside = opp; s.pa = q.a; s.pb = q.b; s.rev = false; s.paired = true;
                    q.pside = s.side; q.pa = s.a; q.pb = s.b; q.rev = false; q.paired = true; break;
                }
                if (!s.paired) throw std::runtime_error(std::string("Segmento periodic em ") + SIDE_NAME[s.side] + " sem par no lado oposto (mesmo tamanho)");
            }
        }
        for (int sd = 0; sd < 4; sd++) {
            Seg* open = nullptr;
            for (auto& s : segs) if (s.side == sd && s.type == BC_CUT) {
                if (!open) open = &s;
                else {
                    if ((open->b - open->a) != (s.b - s.a)) throw std::runtime_error(std::string("Segmentos 'cut' pareados em ") + SIDE_NAME[sd] + " precisam ter o mesmo tamanho");
                    open->pside = sd; open->pa = s.a; open->pb = s.b; open->rev = true; open->paired = true;
                    s.pside = sd; s.pa = open->a; s.pb = open->b; s.rev = true; s.paired = true;
                    open = nullptr;
                }
            }
            if (open) throw std::runtime_error(std::string("Numero impar de segmentos 'cut' em ") + SIDE_NAME[sd]);
        }
        // Coordenadas fantasma coerentes nos pares periodicos/cut. Assume
        // pareamento por translacao (a mesma hipotese usada ao copiar os estados).
        // Sem isso, o LSQ usaria centros refletidos em vez dos centros do parceiro.
        for (const Seg& s : segs) if ((s.type == BC_PERIODIC || s.type == BC_CUT) && s.paired) {
            for (int t = s.a; t <= s.b; ++t) {
                const int tp = s.rev ? (s.pb - (t - s.a)) : (s.pa + (t - s.a));
                double Sx, Sy, sigma; int fid, axis;
                boundaryFace(s.side, t, Sx, Sy, fid, axis, sigma);
                const double xLocal = axis == 0 ? xfi[fid] : xfj[fid];
                const double yLocal = axis == 0 ? yfi[fid] : yfj[fid];
                double Spx, Spy, sigmap; int fidp, axisp;
                boundaryFace(s.pside, tp, Spx, Spy, fidp, axisp, sigmap);
                const double xPartner = axisp == 0 ? xfi[fidp] : xfj[fidp];
                const double yPartner = axisp == 0 ? yfi[fidp] : yfj[fidp];
                for (int m = 0; m < G; ++m) {
                    int pi_, pj_, gi, gj;
                    interiorCell(s.pside, tp, m, pi_, pj_);
                    ghostCell(s.side, t, m, gi, gj);
                    xc[c(gi, gj)] = xc[c(pi_, pj_)] + (xLocal - xPartner);
                    yc[c(gi, gj)] = yc[c(pi_, pj_)] + (yLocal - yPartner);
                }
            }
        }
        // marcacao de faces e lista de paredes
        btI.assign(static_cast<size_t>(nci + 1) * ncj, BC_NONE); btJ.assign(static_cast<size_t>(nci) * (ncj + 1), BC_NONE);
        wI.assign(btI.size(), -1); wJ.assign(btJ.size(), -1);
        for (auto& s : segs) for (int t = s.a; t <= s.b; t++) {
            double Sx, Sy, sigma; int fid, axis;
            boundaryFace(s.side, t, Sx, Sy, fid, axis, sigma);
            (axis == 0 ? btI : btJ)[fid] = s.type;
            if (s.type == BC_WALL) {
                WallFace w; w.axis = axis; w.fid = fid; w.sigma = sigma;
                interiorCell(s.side, t, 0, w.ic, w.jc);
                w.A = std::hypot(Sx, Sy); w.nx = sigma * Sx / w.A; w.ny = sigma * Sy / w.A;
                w.xf = axis == 0 ? xfi[fid] : xfj[fid]; w.yf = axis == 0 ? yfi[fid] : yfj[fid];
                if (axis == 0) {
                    int f = fid % (nci + 1), j = fid / (nci + 1);
                    w.xa = xn[static_cast<size_t>(j) * NI + f];     w.ya = yn[static_cast<size_t>(j) * NI + f];
                    w.xb = xn[static_cast<size_t>(j + 1) * NI + f]; w.yb = yn[static_cast<size_t>(j + 1) * NI + f];
                } else {
                    int i = fid % nci, f = fid / nci;
                    w.xa = xn[static_cast<size_t>(f) * NI + i];     w.ya = yn[static_cast<size_t>(f) * NI + i];
                    w.xb = xn[static_cast<size_t>(f) * NI + i + 1]; w.yb = yn[static_cast<size_t>(f) * NI + i + 1];
                }
                int q = c(w.ic, w.jc);
                w.y1 = std::max((w.xf - xc[q]) * w.nx + (w.yf - yc[q]) * w.ny, 1e-14);
                (axis == 0 ? wI : wJ)[fid] = static_cast<int>(walls.size());
                walls.push_back(w);
            }
        }
        buildWallDistanceTree();
        std::cout << "Contornos: " << segs.size() << " segmento(s), " << walls.size() << " face(s) de parede\n";
    }

    // ------------------------- distancia ate a parede ------------------------
    // BVH sobre os segmentos de parede: distancia euclidiana exata sem o custo
    // O(numero_de_celulas * numero_de_faces_de_parede).
    static double pointSegmentDistance2(double x, double y, const WallFace& w) {
        const double dx = w.xb - w.xa, dy = w.yb - w.ya;
        const double l2 = dx * dx + dy * dy;
        double t = l2 > 0.0 ? ((x - w.xa) * dx + (y - w.ya) * dy) / l2 : 0.0;
        t = std::max(0.0, std::min(1.0, t));
        const double ex = x - (w.xa + t * dx), ey = y - (w.ya + t * dy);
        return ex * ex + ey * ey;
    }
    int buildWallTreeNode(int begin, int end) {
        WallTreeNode node;
        node.xmin = node.ymin = std::numeric_limits<double>::max();
        node.xmax = node.ymax = -std::numeric_limits<double>::max();
        node.begin = begin; node.end = end;
        for (int n = begin; n < end; ++n) {
            const WallFace& w = walls[wallTreeOrder[n]];
            node.xmin = std::min({node.xmin, w.xa, w.xb}); node.xmax = std::max({node.xmax, w.xa, w.xb});
            node.ymin = std::min({node.ymin, w.ya, w.yb}); node.ymax = std::max({node.ymax, w.ya, w.yb});
        }
        const int id = static_cast<int>(wallTree.size());
        wallTree.push_back(node);
        if (end - begin > 8) {
            const bool splitX = (node.xmax - node.xmin) >= (node.ymax - node.ymin);
            const int mid = begin + (end - begin) / 2;
            std::nth_element(wallTreeOrder.begin() + begin, wallTreeOrder.begin() + mid, wallTreeOrder.begin() + end,
                [&](int a, int b) {
                    const WallFace& A = walls[a]; const WallFace& B = walls[b];
                    const double ca = splitX ? 0.5 * (A.xa + A.xb) : 0.5 * (A.ya + A.yb);
                    const double cb = splitX ? 0.5 * (B.xa + B.xb) : 0.5 * (B.ya + B.yb);
                    return ca < cb;
                });
            const int left = buildWallTreeNode(begin, mid);
            const int right = buildWallTreeNode(mid, end);
            wallTree[id].left = left; wallTree[id].right = right;
        }
        return id;
    }
    void buildWallDistanceTree() {
        wallTree.clear(); wallTreeOrder.resize(walls.size());
        for (size_t n = 0; n < walls.size(); ++n) wallTreeOrder[n] = static_cast<int>(n);
        if (!walls.empty()) buildWallTreeNode(0, static_cast<int>(walls.size()));
    }
    static double pointAabbDistance2(double x, double y, const WallTreeNode& n) {
        const double dx = (x < n.xmin) ? n.xmin - x : ((x > n.xmax) ? x - n.xmax : 0.0);
        const double dy = (y < n.ymin) ? n.ymin - y : ((y > n.ymax) ? y - n.ymax : 0.0);
        return dx * dx + dy * dy;
    }
    void queryWallTree(int id, double x, double y, double& best2) const {
        const WallTreeNode& n = wallTree[id];
        if (pointAabbDistance2(x, y, n) >= best2) return;
        if (n.left < 0) {
            for (int m = n.begin; m < n.end; ++m)
                best2 = std::min(best2, pointSegmentDistance2(x, y, walls[wallTreeOrder[m]]));
            return;
        }
        const double dl = pointAabbDistance2(x, y, wallTree[n.left]);
        const double dr = pointAabbDistance2(x, y, wallTree[n.right]);
        if (dl < dr) { queryWallTree(n.left, x, y, best2); queryWallTree(n.right, x, y, best2); }
        else         { queryWallTree(n.right, x, y, best2); queryWallTree(n.left, x, y, best2); }
    }
    void computeWallDistances() {
        if (wallDist.size() != static_cast<size_t>(N)) wallDist.assign(N, 0.0);
        const double noWallDistance = 1.0e6 * std::max(lref, 1e-30);
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < ncj; ++j) for (int i = 0; i < nci; ++i) {
            const int q = c(i, j);
            if (wallTree.empty()) wallDist[q] = noWallDistance;
            else {
                double d2 = std::numeric_limits<double>::max();
                queryWallTree(0, xc[q], yc[q], d2);
                wallDist[q] = std::max(std::sqrt(std::max(d2, 0.0)), 1e-14 * std::max(lref, 1e-30));
            }
        }
    }
    double wallOmegaValue(const WallFace& W, int q) const {
        const double Tw = (P.T_wall > 0.0) ? P.T_wall : T[q];
        const double rhoW = p[q] / (Rg * Tw);
        const double nuW = viscosity(Tw) / std::max(rhoW, rhoMin);
        // Condicao SST-2003 de parede: omega_w = 60 nu_w / (beta_1 d_1^2).
        return std::max(60.0 * nuW / (SST::beta1 * W.y1 * W.y1), SST::omegaMin);
    }

    // ----------------------------------------------------------- estado base
    void allocate() {
        for (int v = 0; v < NV; v++) { U[v].assign(N, 0); U0[v].assign(N, 0); Ub[v].assign(N, 0); Lr[v].assign(N, 0);
            Fi[v].assign(static_cast<size_t>(nci + 1) * ncj, 0); Fj[v].assign(static_cast<size_t>(nci) * (ncj + 1), 0); }
        for (auto* a : {&rho, &u, &v, &p, &T, &k, &omega, &mu, &mut, &F1, &F2, &sigmaK, &sigmaOmega, &gammaOmega, &betaOmega, &wallDist, &sinkK, &sinkOmega, &dtc}) a->assign(N, 0);
        for (int q = 0; q < 7; q++) { gx[q].assign(N, 0); gy[q].assign(N, 0); }
    }
    void setConstants() {
        gam = P.gamma; Rg = P.R_gas; Cv = Rg / (gam - 1.0); Cp = gam * Cv;
        TI = P.T_inf; pI = P.p_inf; rhoI = pI / (Rg * TI);
        double cI = std::sqrt(gam * Rg * TI);
        Uinf = P.mach * cI;
        double al = P.alpha_deg * M_PI / 180.0;
        uI = Uinf * std::cos(al); vI = Uinf * std::sin(al);
        muI = viscosity(TI);
        kI = std::max(1.5 * std::pow(P.Tu * Uinf, 2), SST::kMin);
        // Para mu_t/mu especificado no escoamento livre, mu_t ~= rho*k/omega.
        omegaI = std::max(rhoI * kI / (P.visc_ratio * muI), SST::omegaMin);
        rhoMin = 1e-4 * rhoI; eMin = Cv * 0.05 * TI;
        // comprimento de referencia
        lref = P.l_ref;
        if (lref <= 0) {
            double xmin = 1e300, xmax = -1e300;
            for (auto& w : walls) { xmin = std::min(xmin, w.xf); xmax = std::max(xmax, w.xf); }
            lref = walls.empty() ? 1.0 : std::max(xmax - xmin, 1e-30);
        }
        // escalas FIXAS dos residuos adimensionais: res = rms(dU/dt) * tRef / escala,  tRef = h_min / (U + c)
        double hmin = 1e300;
        for (int jj = 0; jj < ncj; ++jj) for (int ii = 0; ii < nci; ++ii) {
            // Altura minima aproximada da celula: A / maior comprimento de aresta.
            // sqrt(A) superestima a escala em celulas muito alongadas.
            const double x0 = xn[static_cast<size_t>(jj) * NI + ii],     y0 = yn[static_cast<size_t>(jj) * NI + ii];
            const double x1 = xn[static_cast<size_t>(jj) * NI + ii + 1], y1 = yn[static_cast<size_t>(jj) * NI + ii + 1];
            const double x2 = xn[static_cast<size_t>(jj + 1) * NI + ii + 1], y2 = yn[static_cast<size_t>(jj + 1) * NI + ii + 1];
            const double x3 = xn[static_cast<size_t>(jj + 1) * NI + ii],   y3 = yn[static_cast<size_t>(jj + 1) * NI + ii];
            const double e0 = std::hypot(x1 - x0, y1 - y0), e1 = std::hypot(x2 - x1, y2 - y1);
            const double e2 = std::hypot(x3 - x2, y3 - y2), e3 = std::hypot(x0 - x3, y0 - y3);
            const double emax = std::max({e0, e1, e2, e3, 1e-300});
            hmin = std::min(hmin, vol[c(ii, jj)] / emax);
        }
        hmin = std::max(hmin, 1e-30);
        tRef = hmin / (Uinf + cI);
        const double ks = std::max(kI, 1e-4 * cI * cI);
        const double omegaScale = std::max({omegaI, Uinf / std::max(lref, 1e-12), 1e-12});
        resScale[0] = rhoI; resScale[1] = resScale[2] = rhoI * std::max(Uinf, 0.1 * cI); resScale[3] = rhoI * cI * cI;
        resScale[4] = P.turbulence ? rhoI * ks : 1.0; resScale[5] = P.turbulence ? rhoI * omegaScale : 1.0;
        computeWallDistances();
    }
    double viscosity(double Tt) const { return P.mu_ref * std::pow(Tt / P.T_ref, 1.5) * (P.T_ref + P.S_suth) / (Tt + P.S_suth); }

    void initUniform() {
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            int q = c(i, j);
            U[0][q] = rhoI; U[1][q] = rhoI * uI; U[2][q] = rhoI * vI;
            U[3][q] = rhoI * (Cv * TI + 0.5 * (uI * uI + vI * vI));
            U[4][q] = rhoI * kI; U[5][q] = rhoI * omegaI;
        }
    }

    // ------------------------------------------------------------ primitivas
    // Converte o estado conservativo em primitivo. mut e calculado depois dos
    // gradientes, pois o SST necessita de F2 e do invariante de deformacao S.
    int primitives() {
        int bad = 0; long ck = 0;
        #pragma omp parallel for reduction(+ : bad, ck) schedule(static)
        for (int j = 0; j < ncj; ++j) for (int i = 0; i < nci; ++i) {
            const int q = c(i, j);
            double r = U[0][q];
            if (!(r > rhoMin) || !std::isfinite(r)) { ++bad; r = rhoI; }
            double uu = U[1][q] / r, vv = U[2][q] / r;
            double e = U[3][q] / r - 0.5 * (uu * uu + vv * vv);
            if (!(e > eMin) || !std::isfinite(uu) || !std::isfinite(vv) || !std::isfinite(e)) {
                ++bad; e = Cv * TI; uu = uI; vv = vI;
            }
            double kk = U[4][q] / r, ww = U[5][q] / r;
            if (P.turbulence) {
                if (!(kk > 0.0) || !(ww > 0.0) || !std::isfinite(kk) || !std::isfinite(ww)) {
                    ++bad; kk = kI; ww = omegaI;
                }
            }
            if (kk < SST::kMin) { kk = SST::kMin; U[4][q] = r * kk; ++ck; }
            if (ww < SST::omegaMin) { ww = SST::omegaMin; U[5][q] = r * ww; ++ck; }
            rho[q] = r; u[q] = uu; v[q] = vv;
            T[q] = e / Cv; p[q] = r * Rg * T[q];
            k[q] = kk; omega[q] = ww;
            mu[q] = viscosity(T[q]);
            mut[q] = 0.0;
        }
        clipsTurb += ck;
        return bad;
    }

    // ----------------------------------------------------------- contorno BC
    void farfieldState(int q0, double nx, double ny, double* Wb) const {
        double ri = rho[q0], ui = u[q0], vi = v[q0], pi_ = p[q0];
        double ci = std::sqrt(gam * pi_ / ri), uni = ui * nx + vi * ny;
        double cinf = std::sqrt(gam * pI / rhoI), unf = uI * nx + vI * ny;
        double rb, ub, vb, pb, kb, eb;
        if (uni / ci >= 1.0) { rb = ri; ub = ui; vb = vi; pb = pi_; kb = k[q0]; eb = omega[q0]; }      // saida supersonica
        else if (uni / ci <= -1.0) { rb = rhoI; ub = uI; vb = vI; pb = pI; kb = kI; eb = omegaI; }       // entrada supersonica
        else {                                                                                           // subsonico: invariantes de Riemann
            double Rp = uni + 2.0 * ci / (gam - 1.0), Rm = unf - 2.0 * cinf / (gam - 1.0);
            double unb = 0.5 * (Rp + Rm), cb = 0.25 * (gam - 1.0) * (Rp - Rm);
            double s, utx, uty;
            if (unb > 0) { s = pi_ / std::pow(ri, gam); utx = ui - uni * nx; uty = vi - uni * ny; kb = k[q0]; eb = omega[q0]; }
            else { s = pI / std::pow(rhoI, gam); utx = uI - unf * nx; uty = vI - unf * ny; kb = kI; eb = omegaI; }
            rb = std::pow(cb * cb / (gam * s), 1.0 / (gam - 1.0));
            pb = rb * cb * cb / gam;
            ub = utx + unb * nx; vb = uty + unb * ny;
        }
        Wb[0] = rb; Wb[1] = ub; Wb[2] = vb; Wb[3] = pb; Wb[4] = kb; Wb[5] = eb;
    }

    void applyBC() {
        for (const Seg& s : segs) {
            for (int t = s.a; t <= s.b; ++t) {
                double Sx, Sy, sigma; int fid, axis;
                boundaryFace(s.side, t, Sx, Sy, fid, axis, sigma);
                const double A = std::hypot(Sx, Sy), nx = sigma * Sx / A, ny = sigma * Sy / A;
                int i0, j0; interiorCell(s.side, t, 0, i0, j0);
                const int q0 = c(i0, j0);
                double Wb[6] = {0, 0, 0, 0, 0, 0};
                if (s.type == BC_FARFIELD) farfieldState(q0, nx, ny, Wb);
                double omegaWall = SST::omegaMin;
                if (s.type == BC_WALL) {
                    const int wid = (axis == 0 ? wI : wJ)[fid];
                    if (wid >= 0) omegaWall = wallOmegaValue(walls[wid], q0);
                }
                for (int m = 0; m < G; ++m) {
                    int ii, jj, gi, gj;
                    interiorCell(s.side, t, m, ii, jj); ghostCell(s.side, t, m, gi, gj);
                    const int qi = c(ii, jj), qg = c(gi, gj);
                    switch (s.type) {
                        case BC_WALL: {
                            // SST wall-resolved: no-slip, k_wall=0 e omega_wall prescrita.
                            u[qg] = -u[qi]; v[qg] = -v[qi]; p[qg] = p[qi];
                            T[qg] = (P.T_wall > 0.0) ? std::max(2.0 * P.T_wall - T[qi], 0.2 * P.T_wall) : T[qi];
                            rho[qg] = p[qg] / (Rg * T[qg]);
                            k[qg] = -k[qi];
                            omega[qg] = 2.0 * omegaWall - omega[qi];
                        } break;
                        case BC_SYMMETRY: {
                            const double un = u[qi] * nx + v[qi] * ny;
                            u[qg] = u[qi] - 2.0 * un * nx; v[qg] = v[qi] - 2.0 * un * ny;
                            p[qg] = p[qi]; T[qg] = T[qi]; rho[qg] = rho[qi];
                            k[qg] = k[qi]; omega[qg] = omega[qi];
                        } break;
                        case BC_FARFIELD: {
                            rho[qg] = Wb[0]; u[qg] = Wb[1]; v[qg] = Wb[2]; p[qg] = Wb[3];
                            k[qg] = Wb[4]; omega[qg] = Wb[5]; T[qg] = p[qg] / (Rg * rho[qg]);
                        } break;
                        case BC_OUTFLOW: {
                            u[qg] = u[qi]; v[qg] = v[qi]; k[qg] = k[qi]; omega[qg] = omega[qi];
                            if (P.p_back > 0.0) { p[qg] = P.p_back; T[qg] = T[qi]; rho[qg] = p[qg] / (Rg * T[qg]); }
                            else { p[qg] = p[qi]; T[qg] = T[qi]; rho[qg] = rho[qi]; }
                        } break;
                        case BC_PERIODIC: case BC_CUT: {
                            const int tp = s.rev ? (s.pb - (t - s.a)) : (s.pa + (t - s.a));
                            int pi_, pj_; interiorCell(s.pside, tp, m, pi_, pj_);
                            const int qp = c(pi_, pj_);
                            rho[qg] = rho[qp]; u[qg] = u[qp]; v[qg] = v[qp]; p[qg] = p[qp]; T[qg] = T[qp];
                            k[qg] = k[qp]; omega[qg] = omega[qp];
                        } break;
                        default: break;
                    }
                    mu[qg] = viscosity(T[qg]);
                }
            }
        }
    }

    // -------------------------------------------------------------- gradientes
    void computeGradients() {
        const std::vector<double>* Q[7] = {&u, &v, &T, &k, &omega, &p, &rho};
        const bool lsq = (P.gradient == "lsq");
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < ncj; ++j) for (int i = 0; i < nci; ++i) {
            const int q = c(i, j), qe = c(i + 1, j), qw = c(i - 1, j), qn = c(i, j + 1), qs = c(i, j - 1);
            bool useLsq = false;
            double dxn[4], dyn[4], wn[4], a11 = 0.0, a12 = 0.0, a22 = 0.0, det = 0.0;
            if (lsq) {
                const int nb[4] = {qe, qw, qn, qs};
                for (int m = 0; m < 4; ++m) {
                    dxn[m] = xc[nb[m]] - xc[q]; dyn[m] = yc[nb[m]] - yc[q];
                    const double d2 = dxn[m] * dxn[m] + dyn[m] * dyn[m];
                    wn[m] = 1.0 / std::max(d2, 1e-300);
                    a11 += wn[m] * dxn[m] * dxn[m]; a12 += wn[m] * dxn[m] * dyn[m]; a22 += wn[m] * dyn[m] * dyn[m];
                }
                det = a11 * a22 - a12 * a12;
                const double scale = (a11 + a22) * (a11 + a22);
                useLsq = std::isfinite(det) && det > 1e-12 * std::max(scale, 1e-300);
                if (useLsq) {
                    for (int a = 0; a < 7; ++a) {
                        const std::vector<double>& phi = *Q[a];
                        double bx = 0.0, by = 0.0;
                        for (int m = 0; m < 4; ++m) {
                            const double dphi = phi[nb[m]] - phi[q];
                            bx += wn[m] * dxn[m] * dphi; by += wn[m] * dyn[m] * dphi;
                        }
                        gx[a][q] = (a22 * bx - a12 * by) / det;
                        gy[a][q] = (a11 * by - a12 * bx) / det;
                    }
                }
            }
            if (!useLsq) { // Green-Gauss como alternativa para LSQ singular/mal condicionado.
                const int fe = fi(i + 1, j), fw = fi(i, j), fn = fj(i, j + 1), fs = fj(i, j);
                const double V = vol[q];
                for (int a = 0; a < 7; ++a) {
                    const std::vector<double>& phi = *Q[a];
                    const double pe = 0.5 * (phi[q] + phi[qe]), pw = 0.5 * (phi[q] + phi[qw]);
                    const double pn = 0.5 * (phi[q] + phi[qn]), ps = 0.5 * (phi[q] + phi[qs]);
                    gx[a][q] = (pe * Six[fe] - pw * Six[fw] + pn * Sjx[fn] - ps * Sjx[fs]) / V;
                    gy[a][q] = (pe * Siy[fe] - pw * Siy[fw] + pn * Sjy[fn] - ps * Sjy[fs]) / V;
                }
            }
        }
        // Gradientes nas celulas fantasma da primeira camada.
        for (const Seg& s : segs) for (int t = s.a; t <= s.b; ++t) {
            int i0, j0, gi, gj; interiorCell(s.side, t, 0, i0, j0); ghostCell(s.side, t, 0, gi, gj);
            int src = c(i0, j0);
            if (s.type == BC_PERIODIC || s.type == BC_CUT) {
                const int tp = s.rev ? (s.pb - (t - s.a)) : (s.pa + (t - s.a));
                int pi_, pj_; interiorCell(s.pside, tp, 0, pi_, pj_); src = c(pi_, pj_);
            }
            const int qg = c(gi, gj);
            for (int a = 0; a < 7; ++a) { gx[a][qg] = gx[a][src]; gy[a][qg] = gy[a][src]; }
        }
    }

    // Atualiza F1/F2, coeficientes blended e viscosidade turbulenta SST.
    void updateSSTProperties() {
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < ncj; ++j) for (int i = 0; i < nci; ++i) {
            const int q = c(i, j);
            if (!P.turbulence) {
                F1[q] = F2[q] = 0.0; sigmaK[q] = SST::sigmaK2; sigmaOmega[q] = SST::sigmaW2;
                gammaOmega[q] = SST::gamma2; betaOmega[q] = SST::beta2; mut[q] = 0.0;
                continue;
            }
            const double kk = std::max(k[q], SST::kMin), ww = std::max(omega[q], SST::omegaMin);
            const double d = std::max(wallDist[q], 1e-14 * std::max(lref, 1e-30));
            const double nu = mu[q] / rho[q];
            const double ux = gx[0][q], uy = gy[0][q], vx = gx[1][q], vy = gy[1][q];
            const double strain = std::sqrt(std::max(2.0 * ux * ux + 2.0 * vy * vy + (uy + vx) * (uy + vx), 0.0));
            const double gradDot = gx[3][q] * gx[4][q] + gy[3][q] * gy[4][q];
            const double cdkw = std::max(2.0 * rho[q] * SST::sigmaW2 * gradDot / ww, 1e-10);
            const double arg1a = std::sqrt(kk) / (SST::betaStar * ww * d);
            const double arg1b = 500.0 * nu / (d * d * ww);
            const double arg1c = 4.0 * rho[q] * SST::sigmaW2 * kk / (cdkw * d * d);
            const double arg1 = std::min(std::max(arg1a, arg1b), arg1c);
            const double arg2 = std::max(2.0 * std::sqrt(kk) / (SST::betaStar * ww * d), arg1b);
            F1[q] = std::tanh(std::pow(std::max(arg1, 0.0), 4.0));
            F2[q] = std::tanh(std::max(arg2, 0.0) * std::max(arg2, 0.0));
            sigmaK[q] = F1[q] * SST::sigmaK1 + (1.0 - F1[q]) * SST::sigmaK2;
            sigmaOmega[q] = F1[q] * SST::sigmaW1 + (1.0 - F1[q]) * SST::sigmaW2;
            gammaOmega[q] = F1[q] * SST::gamma1 + (1.0 - F1[q]) * SST::gamma2;
            betaOmega[q] = F1[q] * SST::beta1 + (1.0 - F1[q]) * SST::beta2;
            const double denom = std::max(SST::a1 * ww, strain * F2[q]);
            mut[q] = std::min(rho[q] * SST::a1 * kk / std::max(denom, 1e-30), SST::muRatioMax * mu[q]);
        }
        // Ghost viscosities and blends: wall ghosts are set to mu_t=0; other
        // boundary ghosts use the adjacent (or periodic partner) blending state.
        for (const Seg& s : segs) for (int t = s.a; t <= s.b; ++t) {
            for (int m = 0; m < G; ++m) {
                int ii, jj, gi, gj; interiorCell(s.side, t, m, ii, jj); ghostCell(s.side, t, m, gi, gj);
                const int qi = c(ii, jj), qg = c(gi, gj);
                int src = qi;
                if (s.type == BC_PERIODIC || s.type == BC_CUT) {
                    const int tp = s.rev ? (s.pb - (t - s.a)) : (s.pa + (t - s.a));
                    int pi_, pj_; interiorCell(s.pside, tp, m, pi_, pj_); src = c(pi_, pj_);
                }
                F1[qg] = F1[src]; F2[qg] = F2[src]; sigmaK[qg] = sigmaK[src]; sigmaOmega[qg] = sigmaOmega[src];
                gammaOmega[qg] = gammaOmega[src]; betaOmega[qg] = betaOmega[src]; wallDist[qg] = wallDist[src];
                if (!P.turbulence || s.type == BC_WALL) { mut[qg] = 0.0; continue; }
                const double ux = gx[0][qg], uy = gy[0][qg], vx = gx[1][qg], vy = gy[1][qg];
                const double strain = std::sqrt(std::max(2.0 * ux * ux + 2.0 * vy * vy + (uy + vx) * (uy + vx), 0.0));
                const double denom = std::max(SST::a1 * std::max(omega[qg], SST::omegaMin), strain * F2[qg]);
                mut[qg] = std::min(rho[qg] * SST::a1 * std::max(k[qg], 0.0) / std::max(denom, 1e-30), SST::muRatioMax * mu[qg]);
            }
        }
    }

    // --------------------------------------------------------------- reconstrucao
    inline void reconstruct(int cLL, int cL, int cR, int cRR, double* WL, double* WR) const {
        const std::vector<double>* Q[6] = {&rho, &u, &v, &p, &k, &omega};
        for (int a = 0; a < 6; a++) {
            const std::vector<double>& q = *Q[a];
            double qLL = q[cLL], qL = q[cL], qR = q[cR], qRR = q[cRR];
            if (P.order == 1 || (a >= 4 && P.turb_order == 1)) { WL[a] = qL; WR[a] = qR; continue; }   // k, omega: upwind de 1a ordem
            double sc = std::fabs(qL) + std::fabs(qR), e2 = 1e-6 * sc * sc + 1e-30;   // piso absoluto evita 0/0 (ex.: v=0 uniforme)
            WL[a] = qL + 0.5 * vanAlbada(qL - qLL, qR - qL, e2);
            WR[a] = qR - 0.5 * vanAlbada(qR - qL, qRR - qR, e2);
        }
        if (!(WL[0] > 0 && WL[3] > 0 && WR[0] > 0 && WR[3] > 0 && WL[4] >= 0 && WR[4] >= 0 && WL[5] >= 0 && WR[5] >= 0))
            for (int a = 0; a < 6; a++) { WL[a] = (*Q[a])[cL]; WR[a] = (*Q[a])[cR]; }   // fallback de 1a ordem (preserva positividade)
    }

    // --------------------------------------------------------------- fluxo AUSM+
    inline void ausmPlus(const double* WL, const double* WR, double nx, double ny, double A, double* F) const {
        double rL = WL[0], uL = WL[1], vL = WL[2], pL = WL[3];
        double rR = WR[0], uR = WR[1], vR = WR[2], pR = WR[3];
        double qL2 = uL * uL + vL * vL, qR2 = uR * uR + vR * vR;
        double HL = gam / (gam - 1.0) * pL / rL + 0.5 * qL2, HR = gam / (gam - 1.0) * pR / rR + 0.5 * qR2;
        double unL = uL * nx + vL * ny, unR = uR * nx + vR * ny;
        double as2L = 2.0 * (gam - 1.0) / (gam + 1.0) * HL, as2R = 2.0 * (gam - 1.0) / (gam + 1.0) * HR;
        double asL = std::sqrt(as2L), asR = std::sqrt(as2R);
        double atL = as2L / std::max(asL, unL), atR = as2R / std::max(asR, -unR);
        double a12 = std::min(atL, atR);
        double ML = unL / a12, MR = unR / a12;
        double M12 = M4p(ML) + M4m(MR);
        double p12 = P5p(ML) * pL + P5m(MR) * pR;
        double mdot = a12 * M12 * (M12 > 0 ? rL : rR);
        const double* W = (mdot >= 0) ? WL : WR;
        double H = (mdot >= 0) ? HL : HR;
        F[0] = A * mdot;
        F[1] = A * (mdot * W[1] + p12 * nx);
        F[2] = A * (mdot * W[2] + p12 * ny);
        F[3] = A * mdot * H;
        F[4] = A * mdot * W[4];
        F[5] = A * mdot * W[5];
    }

    // ------------------------------------------------------- fluxo de uma face
    inline void faceFlux(int cLL, int cL, int cR, int cRR, double Sx, double Sy, double* F) const {
        double A = std::hypot(Sx, Sy), nx = Sx / A, ny = Sy / A;
        double WL[6], WR[6];
        reconstruct(cLL, cL, cR, cRR, WL, WR);
        ausmPlus(WL, WR, nx, ny, A, F);
        // ---- parte viscosa (face: media + gradiente corrigido na aresta) ----
        double dx = xc[cR] - xc[cL], dy = yc[cR] - yc[cL], dd = std::hypot(dx, dy);
        if (dd <= 0) return;
        double tx = dx / dd, ty = dy / dd;
        auto fg = [&](int a, const std::vector<double>& phi, double& gxf, double& gyf) {
            double gxm = 0.5 * (gx[a][cL] + gx[a][cR]), gym = 0.5 * (gy[a][cL] + gy[a][cR]);
            double corr = (phi[cR] - phi[cL]) / dd - (gxm * tx + gym * ty);
            gxf = gxm + corr * tx; gyf = gym + corr * ty;
        };
        double uxx, uxy, vxx, vxy, Txx, Txy;
        fg(0, u, uxx, uxy); fg(1, v, vxx, vxy); fg(2, T, Txx, Txy);
        double uf = 0.5 * (u[cL] + u[cR]), vf = 0.5 * (v[cL] + v[cR]);
        double muf = 0.5 * (mu[cL] + mu[cR]), mutf = 0.5 * (mut[cL] + mut[cR]);
        double muE = muf + mutf, kE = Cp * (muf / P.Pr + mutf / SST::Prt);
        double div = uxx + vxy;
        double rk = 0.0;
        if (P.turbulence) rk = (2.0 / 3.0) * 0.25 * (rho[cL] + rho[cR]) * (k[cL] + k[cR]);   // 2/3 rho k (media de rho e k)
        double txx = muE * (2.0 * uxx - (2.0 / 3.0) * div) - rk;
        double tyy = muE * (2.0 * vxy - (2.0 / 3.0) * div) - rk;
        double txy = muE * (uxy + vxx);
        double tSx = txx * Sx + txy * Sy, tSy = txy * Sx + tyy * Sy;
        F[1] -= tSx; F[2] -= tSy;
        F[3] += -(uf * tSx + vf * tSy) - kE * (Txx * Sx + Txy * Sy);
        if (P.turbulence) {
            double kx, ky, wx, wy; fg(3, k, kx, ky); fg(4, omega, wx, wy);
            const double skf = 0.5 * (sigmaK[cL] + sigmaK[cR]);
            const double swf = 0.5 * (sigmaOmega[cL] + sigmaOmega[cR]);
            F[4] -= (muf + mutf * skf) * (kx * Sx + ky * Sy);
            F[5] -= (muf + mutf * swf) * (wx * Sx + wy * Sy);
        }
    }

    // ------------------------------ fluxo de parede SST ----------------------
    // Parede resolvida: k_wall=0, omega_wall=60 nu/(beta1 d1^2), sem lei log.
    inline void wallFlux(int w, double Sx, double Sy, double* F) {
        WallFace& W = walls[w];
        const int q = c(W.ic, W.jc);
        double pw = p[q] + gx[5][q] * (W.xf - xc[q]) + gy[5][q] * (W.yf - yc[q]);
        pw = std::min(std::max(pw, 0.5 * p[q]), 2.0 * p[q]);
        const double un = u[q] * W.nx + v[q] * W.ny;
        const double utx = u[q] - un * W.nx, uty = v[q] - un * W.ny;
        const double Ut = std::hypot(utx, uty);
        const double tx = (Ut > 1e-14) ? utx / Ut : 0.0;
        const double ty = (Ut > 1e-14) ? uty / Ut : 0.0;
        const double Tw = (P.T_wall > 0.0) ? P.T_wall : T[q];
        const double muWall = viscosity(Tw);
        const double muFace = 0.5 * (mu[q] + muWall);
        const double tw = muFace * Ut / W.y1; // tau_w resolvida pela primeira celula
        const double uTau = std::sqrt(std::max(tw / rho[q], 0.0));
        const double yp = rho[q] * uTau * W.y1 / mu[q];
        const double omegaW = wallOmegaValue(W, q);
        const double gammaThermal = Cp * (muFace / P.Pr + 0.5 * mut[q] / SST::Prt);
        const double gammaK = muFace + 0.5 * sigmaK[q] * mut[q];
        const double gammaW = muFace + 0.5 * sigmaOmega[q] * mut[q];
        F[0] = 0.0;
        F[1] = pw * Sx + W.sigma * tw * tx * W.A;
        F[2] = pw * Sy + W.sigma * tw * ty * W.A;
        F[3] = (P.T_wall > 0.0) ? -W.sigma * gammaThermal * (P.T_wall - T[q]) / W.y1 * W.A : 0.0;
        // F = -Gamma grad(phi).S, com gradiente normal entre centro e parede.
        F[4] = P.turbulence ? W.sigma * gammaK * k[q] / W.y1 * W.A : 0.0;
        F[5] = P.turbulence ? W.sigma * gammaW * (omega[q] - omegaW) / W.y1 * W.A : 0.0;
        W.Fx = pw * W.nx * W.A + tw * tx * W.A;
        W.Fy = pw * W.ny * W.A + tw * ty * W.A;
        W.tau = tw; W.yplus = yp;
    }

    // ------------------------------------------------------------------ residuo
    // SST-2003: fontes de k/omega; destruicoes tratadas implicitamente no RK.
    int computeResidual() {
        int bad = primitives();
        applyBC();
        computeGradients();
        updateSSTProperties();
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < ncj; ++j) for (int f = 0; f <= nci; ++f) {
            const int id = fi(f, j); double F[NV];
            if (btI[id] == BC_WALL) wallFlux(wI[id], Six[id], Siy[id], F);
            else faceFlux(c(f - 2, j), c(f - 1, j), c(f, j), c(f + 1, j), Six[id], Siy[id], F);
            for (int a = 0; a < NV; ++a) Fi[a][id] = F[a];
        }
        #pragma omp parallel for schedule(static)
        for (int f = 0; f <= ncj; ++f) for (int i = 0; i < nci; ++i) {
            const int id = fj(i, f); double F[NV];
            if (btJ[id] == BC_WALL) wallFlux(wJ[id], Sjx[id], Sjy[id], F);
            else faceFlux(c(i, f - 2), c(i, f - 1), c(i, f), c(i, f + 1), Sjx[id], Sjy[id], F);
            for (int a = 0; a < NV; ++a) Fj[a][id] = F[a];
        }
        double r2[NV] = {0, 0, 0, 0, 0, 0};
        #pragma omp parallel for schedule(static) reduction(+ : r2[:NV])
        for (int j = 0; j < ncj; ++j) for (int i = 0; i < nci; ++i) {
            const int q = c(i, j), f1 = fi(i, j), f2 = fi(i + 1, j), g1 = fj(i, j), g2 = fj(i, j + 1);
            const double V = vol[q];
            for (int a = 0; a < NV; ++a)
                Lr[a][q] = (Fi[a][f1] - Fi[a][f2] + Fj[a][g1] - Fj[a][g2]) / V;
            sinkK[q] = 0.0; sinkOmega[q] = 0.0;
            if (P.turbulence) {
                const double ux = gx[0][q], uy = gy[0][q], vx = gx[1][q], vy = gy[1][q];
                const double div = ux + vy;
                const double strain2 = 2.0 * (ux * ux + vy * vy) + (uy + vx) * (uy + vx);
                // Producao exata tau_ij * dU_i/dx_j para a aproximacao 2D compressivel.
                double Pk = mut[q] * (strain2 - (2.0 / 3.0) * div * div)
                          - (2.0 / 3.0) * rho[q] * k[q] * div;
                // SST-2003: limita superiormente, mas preserva a producao exata
                // (inclusive valores negativos em regioes de dilatacao forte).
                Pk = std::min(Pk,
                    SST::productionLimit * SST::betaStar * rho[q] * omega[q] * k[q]);
                const double gradKw = gx[3][q] * gx[4][q] + gy[3][q] * gy[4][q];
                const double crossDiff = 2.0 * (1.0 - F1[q]) * rho[q] * SST::sigmaW2
                                       * gradKw / std::max(omega[q], SST::omegaMin);
                const double mutSafe = std::max(mut[q], 1e-20 * mu[q]);
                Lr[4][q] += Pk;
                Lr[5][q] += gammaOmega[q] * rho[q] * Pk / mutSafe + crossDiff;
                sinkK[q] = SST::betaStar * omega[q];
                sinkOmega[q] = betaOmega[q] * omega[q];
            }
            for (int a = 0; a < NV; ++a) {
                double Rfull = Lr[a][q];
                if (a == 4) Rfull -= sinkK[q] * U[4][q];
                else if (a == 5) Rfull -= sinkOmega[q] * U[5][q];
                r2[a] += Rfull * Rfull;
            }
        }
        const double n = static_cast<double>(nci) * ncj;
        for (int a = 0; a < NV; ++a) resNow[a] = std::sqrt(r2[a] / n);
        return bad;
    }

    // ------------------------------------------------------------- passo de tempo
    double computeDt(double cfl) {
        double dmin = std::numeric_limits<double>::max();
        #pragma omp parallel for schedule(static) reduction(min : dmin)
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            int q = c(i, j);
            double S1x = 0.5 * (Six[fi(i, j)] + Six[fi(i + 1, j)]), S1y = 0.5 * (Siy[fi(i, j)] + Siy[fi(i + 1, j)]);
            double S2x = 0.5 * (Sjx[fj(i, j)] + Sjx[fj(i, j + 1)]), S2y = 0.5 * (Sjy[fj(i, j)] + Sjy[fj(i, j + 1)]);
            double A1 = std::hypot(S1x, S1y), A2 = std::hypot(S2x, S2y);
            double cs = std::sqrt(gam * p[q] / rho[q]);
            double l1 = std::fabs(u[q] * S1x + v[q] * S1y) + cs * A1, l2 = std::fabs(u[q] * S2x + v[q] * S2y) + cs * A2;
            double nu = std::max((mu[q] + mut[q]) / rho[q], gam * (mu[q] / P.Pr + mut[q] / SST::Prt) / rho[q]);   // inclui a difusao termica turbulenta
            double visc = 2.0 * nu * (A1 * A1 + A2 * A2) / vol[q];
            double dt = cfl * vol[q] / (l1 + l2 + visc);
            dtc[q] = dt; dmin = std::min(dmin, dt);
        }
        return dmin;
    }

    // ------------------------------------------------ estagio RK3 com sumidouro implicito
    void stageUpdate(double a, double b, double dtg, bool local) {
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            int q = c(i, j);
            double dt = local ? dtc[q] : dtg;
            for (int s = 0; s < NV; s++) {
                double val = a * U0[s][q] + b * (U[s][q] + dt * Lr[s][q]);
                if (s == 4) val /= (1.0 + b * dt * sinkK[q]);
                else if (s == 5) val /= (1.0 + b * dt * sinkOmega[q]);
                U[s][q] = val;
            }
        }
    }

    // ----------------------------------------------------------------- um passo
    // devolve 0 se ok; >0 se o estado ficou nao-fisico (passo deve ser rejeitado)
    int rk3Step(double dtg, bool local) {
        int bad = 0;
        for (int s = 0; s < NV; s++) U0[s] = U[s];
        bad += computeResidual(); stageUpdate(0.0, 1.0, dtg, local);
        bad += computeResidual(); stageUpdate(0.75, 0.25, dtg, local);
        bad += computeResidual(); stageUpdate(1.0 / 3.0, 2.0 / 3.0, dtg, local);
        // validacao do estado final
        bad += primitives();
        for (int j = 0; j < ncj && !bad; j++) for (int i = 0; i < nci; i++) {
            int q = c(i, j);
            for (int s = 0; s < NV; s++) if (!std::isfinite(U[s][q])) { bad++; break; }
        }
        return bad;
    }

    // ------------------------------------------------------------------ saida
    void refreshDerived() { primitives(); applyBC(); computeGradients(); updateSSTProperties(); }

    struct Arr { const char* name; int nc; std::vector<double> d; };

    // Monta os campos por celula (atualiza contornos e gradientes antes). Nunca devolve nan/inf.
    std::vector<Arr> buildArrays() {
        refreshDerived();
        double q0 = 0.5 * rhoI * Uinf * Uinf + 1e-300;
        std::vector<Arr> arr; arr.reserve(16);          // sem realocacao: os ponteiros abaixo permanecem validos
        size_t nc = static_cast<size_t>(nci) * ncj;
        auto mk = [&](const char* n, int ncomp) { arr.push_back({n, ncomp, std::vector<double>(nc * ncomp, 0.0)}); return &arr.back().d; };
        auto* aP = mk("Pressure", 1); auto* aR = mk("Density", 1); auto* aT = mk("Temperature", 1); auto* aV = mk("Velocity", 3);
        auto* aM = mk("Mach", 1); auto* aK = mk("TKE", 1); auto* aE = mk("Omega", 1); auto* aMt = mk("MutOverMu", 1);
        auto* aF1 = mk("SST_F1", 1); auto* aF2 = mk("SST_F2", 1); auto* aD = mk("WallDistance", 1);
        auto* aW = mk("Vorticity", 1); auto* aS = mk("Schlieren", 1); auto* aC = mk("Cp", 1);
        for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) {
            int q = c(i, j); size_t id = static_cast<size_t>(j) * nci + i;
            double cs = std::sqrt(gam * p[q] / rho[q]);
            (*aP)[id] = p[q]; (*aR)[id] = rho[q]; (*aT)[id] = T[q];
            (*aV)[3 * id] = u[q]; (*aV)[3 * id + 1] = v[q]; (*aV)[3 * id + 2] = 0;
            (*aM)[id] = std::hypot(u[q], v[q]) / cs; (*aK)[id] = k[q]; (*aE)[id] = omega[q];
            (*aMt)[id] = mut[q] / mu[q]; (*aF1)[id] = F1[q]; (*aF2)[id] = F2[q]; (*aD)[id] = wallDist[q];
            (*aW)[id] = gx[1][q] - gy[0][q];
            (*aS)[id] = std::hypot(gx[6][q], gy[6][q]) / rho[q];
            (*aC)[id] = (p[q] - pI) / q0;
        }
        for (auto& a : arr) for (auto& x : a.d) if (!std::isfinite(x)) x = 0.0;
        return arr;
    }
    // valor no no (i,j): media das celulas vizinhas
    double pointValue(const Arr& a, int i, int j, int cc) const {
        double s = 0; int n = 0;
        for (int dj = -1; dj <= 0; dj++) for (int di = -1; di <= 0; di++) {
            int ci_ = i + di, cj_ = j + dj;
            if (ci_ >= 0 && ci_ < nci && cj_ >= 0 && cj_ < ncj) { s += a.d[(static_cast<size_t>(cj_) * nci + ci_) * a.nc + cc]; n++; }
        }
        return s / n;
    }

    // ---- formato XML (.vts): recomendado ----
    void writeVTS(const std::string& path) {
        std::vector<Arr> arr = buildArrays();
        std::ofstream o(path);
        if (!o) { std::cerr << "Nao foi possivel escrever " << path << "\n"; return; }
        o << std::setprecision(7);
        size_t nc = static_cast<size_t>(nci) * ncj;
        o << "<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        o << "<StructuredGrid WholeExtent=\"0 " << NI - 1 << " 0 " << NJ - 1 << " 0 0\">\n";
        o << "<Piece Extent=\"0 " << NI - 1 << " 0 " << NJ - 1 << " 0 0\">\n";
        o << "<PointData Scalars=\"Mach\">\n";
        for (auto& a : arr) {
            o << "<DataArray type=\"Float32\" Name=\"" << a.name << "\" NumberOfComponents=\"" << a.nc << "\" format=\"ascii\">\n";
            for (int j = 0; j < NJ; j++) for (int i = 0; i < NI; i++)
                for (int cc = 0; cc < a.nc; cc++) o << pointValue(a, i, j, cc) << (cc + 1 < a.nc ? " " : "\n");
            o << "</DataArray>\n";
        }
        o << "</PointData>\n<CellData Scalars=\"Mach\">\n";
        for (auto& a : arr) {
            o << "<DataArray type=\"Float32\" Name=\"" << a.name << "\" NumberOfComponents=\"" << a.nc << "\" format=\"ascii\">\n";
            for (size_t id = 0; id < nc; id++) for (int cc = 0; cc < a.nc; cc++) o << a.d[id * a.nc + cc] << (cc + 1 < a.nc ? " " : "\n");
            o << "</DataArray>\n";
        }
        o << "</CellData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < NJ; j++) for (int i = 0; i < NI; i++)
            o << xn[static_cast<size_t>(j) * NI + i] << " " << yn[static_cast<size_t>(j) * NI + i] << " 0\n";
        o << "</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
    }

    // ---- formato legado (.vtk, ASCII): opcional (output_format = vtk) ----
    void writeLegacyVTK(const std::string& path, double tval) {
        std::vector<Arr> arr = buildArrays();
        std::ofstream o(path);
        if (!o) { std::cerr << "Nao foi possivel escrever " << path << "\n"; return; }
        o << std::setprecision(7);
        size_t nc = static_cast<size_t>(nci) * ncj, np = static_cast<size_t>(NI) * NJ;
        o << "# vtk DataFile Version 3.0\nCampo 2D RANS k-omega SST-2003 t=" << tval << "\nASCII\nDATASET STRUCTURED_GRID\n";
        o << "DIMENSIONS " << NI << " " << NJ << " 1\nPOINTS " << np << " double\n";
        for (int j = 0; j < NJ; j++) for (int i = 0; i < NI; i++)
            o << xn[static_cast<size_t>(j) * NI + i] << " " << yn[static_cast<size_t>(j) * NI + i] << " 0\n";
        o << "POINT_DATA " << np << "\n";
        for (auto& a : arr) {
            if (a.nc == 1) o << "SCALARS " << a.name << " float 1\nLOOKUP_TABLE default\n"; else o << "VECTORS " << a.name << " float\n";
            for (int j = 0; j < NJ; j++) for (int i = 0; i < NI; i++)
                for (int cc = 0; cc < a.nc; cc++) o << pointValue(a, i, j, cc) << (cc + 1 < a.nc ? " " : "\n");
        }
        o << "CELL_DATA " << nc << "\n";
        for (auto& a : arr) {
            if (a.nc == 1) o << "SCALARS " << a.name << " float 1\nLOOKUP_TABLE default\n"; else o << "VECTORS " << a.name << " float\n";
            for (size_t id = 0; id < nc; id++) for (int cc = 0; cc < a.nc; cc++) o << a.d[id * a.nc + cc] << (cc + 1 < a.nc ? " " : "\n");
        }
    }

    std::string ext() const { return P.output_format == "vtk" ? ".vtk" : ".vts"; }
    void writeField(const std::string& base, double tval) {          // base sem extensao
        if (P.output_format == "vtk") writeLegacyVTK(base + ".vtk", tval); else writeVTS(base + ".vts");
    }

    // ---- indice da serie temporal (solucao.pvd): so uma lista "tempo -> arquivo", sem dados ----
    std::vector<std::pair<double, std::string>> pvd;
    void writePVD() {
        std::ofstream o(P.output_dir + "/solucao.pvd");
        o << std::setprecision(10);
        o << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n";
        for (auto& e : pvd) o << "<DataSet timestep=\"" << e.first << "\" group=\"\" part=\"0\" file=\"" << e.second << "\"/>\n";
        o << "</Collection>\n</VTKFile>\n";
    }
    int frame = 0;
    // tval = tempo fisico [s] (time_mode=global) ou numero da iteracao (time_mode=local)
    void writeFrame(double tval) {
        char nm[64]; std::snprintf(nm, sizeof(nm), "campo_%05d", frame++);
        writeField(P.output_dir + "/" + nm, tval);
        pvd.push_back({tval, std::string(nm) + ext()});
        writePVD();
    }

    // reavalia os dados de parede no estado atual (para Cd, Cl e y+ do estado final do passo)
    void evalWalls() {
        primitives(); applyBC(); computeGradients(); updateSSTProperties();
        for (size_t w = 0; w < walls.size(); w++) {
            const WallFace& W = walls[w];
            double Sx = W.axis == 0 ? Six[W.fid] : Sjx[W.fid], Sy = W.axis == 0 ? Siy[W.fid] : Sjy[W.fid];
            double F[NV]; wallFlux(static_cast<int>(w), Sx, Sy, F);
        }
    }

    // forcas e diagnosticos de parede
    void forces(double& Cd, double& Cl, double& ypMin, double& ypMean, double& ypMax) const {
        double Fx = 0, Fy = 0; ypMin = 1e300; ypMax = 0; ypMean = 0;
        for (auto& w : walls) { Fx += w.Fx; Fy += w.Fy; ypMin = std::min(ypMin, w.yplus); ypMax = std::max(ypMax, w.yplus); ypMean += w.yplus; }
        if (!walls.empty()) ypMean /= walls.size(); else ypMin = 0;
        double al = P.alpha_deg * M_PI / 180.0, q = 0.5 * rhoI * Uinf * Uinf * lref;
        Cd = (Fx * std::cos(al) + Fy * std::sin(al)) / q;
        Cl = (-Fx * std::sin(al) + Fy * std::cos(al)) / q;
    }

    // ------------------------------------------------------------------- execucao
    int run() {
        std::filesystem::create_directories(P.output_dir);
        const bool local = (P.time_mode == "local");
        initUniform();
        std::cout << std::scientific << std::setprecision(4)
                  << "Mach=" << P.mach << "  U=" << Uinf << " m/s  Re/m=" << rhoI * Uinf / muI << "  Re_ref=" << rhoI * Uinf * lref / muI
                  << "  l_ref=" << lref << " m  tempo_convectivo(l_ref/U)=" << lref / Uinf << " s  turbulencia=" << P.turbulence << "  ordem=" << P.order
                  << "  modo_tempo=" << P.time_mode << "  threads=" << omp_get_max_threads() << "\n";
        if (P.turbulence) {
            const double ReL = rhoI * Uinf * lref / muI;
            const double kStar = kI * ReL / (Uinf * Uinf);
            const double omegaStar = omegaI * lref / Uinf;
            std::cout << "Modelo=Menter SST-2003 (parede resolvida, requer y+~1)  k_inf=" << kI
                      << "  omega_inf=" << omegaI << "  mut/mu_inf~=" << P.visc_ratio << "\n";
            if (P.visc_ratio < 1e-5 || P.visc_ratio > 1e-2)
                std::cout << "AVISO: mut/mu_inf fora da faixa de referencia 1e-5..1e-2 da NASA; calibre Tu e visc_ratio.\n";
            if (kStar < 1e-5 || kStar > 0.1)
                std::cout << "AVISO: k_inf*Re_L/U_inf^2=" << kStar << " fora da faixa de referencia 1e-5..0.1 da NASA.\n";
            if (omegaStar < 1.0 || omegaStar > 10.0)
                std::cout << "AVISO: omega_inf*L/U_inf=" << omegaStar << " fora da faixa de referencia 1..10 da NASA.\n";
        }
        std::ofstream hist(P.output_dir + "/historico.csv");
        hist << "it,t,tstar,dt,cfl_scale,Cd,Cl,res_rho,res_rhou,res_rhov,res_rhoE,res_rhok,res_rhoomega,y+min,y+med,y+max,Mach_max\n";
        double t = 0, cflScale = 1.0;
        double nextWriteT = P.write_every_time, nextWriteS = P.write_every_tstar;   // proximos instantes de gravacao
        long it = 0;
        auto T0 = std::chrono::high_resolution_clock::now();
        if (local && (P.write_every_time > 0 || P.write_every_tstar > 0))
            std::cout << "AVISO: em time_mode=local nao ha tempo fisico; use write_every_iter (ou time_mode=global).\n";
        writeFrame(0.0);
        while (true) {
            double tstar = t * Uinf / lref;
            if (it >= P.max_iter) break;
            if (P.t_end_star > 0 && !local && tstar >= P.t_end_star) break;
            double cflT = (it < P.cfl_ramp) ? P.cfl_start + (P.cfl - P.cfl_start) * static_cast<double>(it) / P.cfl_ramp : P.cfl;
            int tries = 0; double dt = 0;
            for (;;) {
                for (int s = 0; s < NV; s++) Ub[s] = U[s];
                primitives();
                dt = computeDt(cflT * cflScale);
                int bad = rk3Step(dt, local);
                if (!bad) break;
                for (int s = 0; s < NV; s++) U[s] = Ub[s];       // rejeita o passo: volta ao estado anterior
                cflScale *= 0.5; backoffs++; tries++;
                if (tries > 12 || cflScale < 1e-4) {
                    std::cout << "\nFALHA: estado nao-fisico persistente na iteracao " << it << " (t*=" << tstar << "). Salvando falha"+ext()+"\n"
                              << "Sugestoes: reduza cfl; use order=1 para testar; confira a malha (celulas muito esticadas) e os contornos.\n";
                    writeField(P.output_dir + "/falha", t);
                    return 2;
                }
            }
            cflScale = std::min(1.0, cflScale * 1.05);
            t += local ? 0.0 : dt; it++;
            tstar = t * Uinf / lref;
            if (!local && P.write_every_time > 0) {
                if (t >= nextWriteT) { writeFrame(t); while (nextWriteT <= t) nextWriteT += P.write_every_time; }
            } else if (!local && P.write_every_tstar > 0) {
                if (tstar >= nextWriteS) { writeFrame(t); while (nextWriteS <= tstar) nextWriteS += P.write_every_tstar; }
            } else if (P.write_every_iter > 0 && it % P.write_every_iter == 0) writeFrame(local ? static_cast<double>(it) : t);
            if (it % P.report_every == 0 || it == 1) {
                evalWalls(); double Cd, Cl, a, b, cc; forces(Cd, Cl, a, b, cc);
                double mmax = 0;
                for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) { int q = c(i, j); mmax = std::max(mmax, std::hypot(u[q], v[q]) / std::sqrt(gam * p[q] / rho[q])); }
                double secs = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - T0).count();
                std::cout << "it=" << it << " t*=" << std::fixed << std::setprecision(3) << tstar << std::scientific << std::setprecision(3)
                          << " dt=" << dt << " cflx=" << cflScale << " res=";
                for (int s = 0; s < NV; s++) std::cout << resNow[s] * tRef / resScale[s] << (s < NV - 1 ? "," : "");
                std::cout << std::fixed << std::setprecision(4) << "  Cd=" << Cd << " Cl=" << Cl << std::setprecision(2)
                          << "  y+=" << a << "/" << b << "/" << cc << "  Mmax=" << mmax
                          << "  clipsTurb=" << clipsTurb << " backoffs=" << backoffs << "  " << std::setprecision(0) << secs << "s\n" << std::scientific;
                hist << it << "," << t << "," << tstar << "," << dt << "," << cflScale << "," << Cd << "," << Cl;
                for (int s = 0; s < NV; s++) hist << "," << resNow[s] * tRef / resScale[s];
                hist << "," << a << "," << b << "," << cc << "," << mmax << "\n"; hist.flush();
            }
        }
        double tstar = t * Uinf / lref;
        { double tv = local ? static_cast<double>(it) : t; if (pvd.empty() || pvd.back().first != tv) writeFrame(tv); }   // evita quadro duplicado
        if (local) writeField(P.output_dir + "/campo_final", t);   // em modo global o ultimo quadro numerado ja e o final
        if (P.test_freestream) {
            double dmax = 0;
            for (int j = 0; j < ncj; j++) for (int i = 0; i < nci; i++) { int q = c(i, j);
                dmax = std::max({dmax, std::fabs(u[q] - uI) / Uinf, std::fabs(v[q] - vI) / Uinf, std::fabs(p[q] - pI) / pI, std::fabs(rho[q] - rhoI) / rhoI}); }
            std::cout << "\nTESTE DE ESCOAMENTO UNIFORME: desvio maximo relativo = " << dmax << "  (esperado ~1e-12 ou menor)\n";
        }
        std::cout << "Fim: " << it << " iteracoes, t*=" << std::fixed << std::setprecision(3) << tstar << ". Arquivos em '" << P.output_dir << "/' (abra solucao.pvd no ParaView)\n";
        return 0;
    }
};

int main(int argc, char** argv) {
    std::cout << std::unitbuf;   // saida sem buffer (funciona com 'tee' e nohup)
    try {
        std::string pf = (argc > 1) ? argv[1] : "parametros.txt";
        std::vector<RawSeg> raw;
        Solver S;
        S.P = readParams(pf, raw);
        validateParams(S.P);
        int NI, NJ; std::vector<double> X, Y;
        readGrid(S.P, NI, NJ, X, Y);
        S.buildMesh(NI, NJ, X, Y);
        S.setupBoundaries(raw);
        S.allocate();
        S.setConstants();
        return S.run();
    } catch (const std::exception& e) {
        std::cerr << "Erro: " << e.what() << "\n";
        return 1;
    }
}
