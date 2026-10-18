// ============================================================================
//  Solver RANS k-epsilon compressivel 2D - MacCormack Navier-Stokes
//  COM "Soft Start" de Inicializacao, Protecao NaN e Saida de Video (ParaView)
// ============================================================================
#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <limits>

struct KE {
    static constexpr double Cmu = 0.09, C1 = 1.44, C2 = 1.92;
    static constexpr double sigk = 1.0, sige = 1.3, Prt = 0.9;
    static constexpr double kappa = 0.41, E = 9.793, yplusLam = 11.225;
    static constexpr double tkeMin = 1.0e-10, epsMin = 1.0e-10;
    static constexpr double muRatioMax = 1.0e5;   
    static constexpr double prodLimit = 10.0;     
};

using Grid2D = std::vector<std::vector<double>>;

class scalarField0 {
public:
    int lines = 0, columns = 0;
    Grid2D F, Ddxi, Ddeta;

    void constructor(int lin, int col) {
        lines = lin; columns = col;
        F.assign(lin, std::vector<double>(col, 0.0));
        Ddxi.assign(lin, std::vector<double>(col, 0.0));
        Ddeta.assign(lin, std::vector<double>(col, 0.0));
    }
    void calculateMetrics() {
        for (int i = 0; i < lines; i++) {
            for (int j = 1; j < columns - 1; j++) Ddxi[i][j] = (F[i][j + 1] - F[i][j - 1]) / 2.0;
            Ddxi[i][0] = (F[i][1] - F[i][columns - 2]) / 2.0;
            Ddxi[i][columns - 1] = Ddxi[i][0];
        }
        for (int j = 0; j < columns; j++) {
            for (int i = 1; i < lines - 1; i++) Ddeta[i][j] = (F[i - 1][j] - F[i + 1][j]) / 2.0;
            Ddeta[0][j] = (F[2][j] - 4.0 * F[1][j] + 3.0 * F[0][j]) / 2.0;
            Ddeta[lines - 1][j] = (-F[lines - 3][j] + 4.0 * F[lines - 2][j] - 3.0 * F[lines - 1][j]) / 2.0;
        }
    }
};

struct gridMetrics {
    int lines = 0, columns = 0;
    scalarField0 x, y;
    Grid2D J;

    void build(const scalarField0& x0, const scalarField0& y0) {
        x = x0; y = y0; lines = x.lines; columns = x.columns;
        J.assign(lines, std::vector<double>(columns, 0.0));
        for (int i = 0; i < lines; i++)
            for (int j = 0; j < columns; j++)
                J[i][j] = x.Ddxi[i][j] * y.Ddeta[i][j] - y.Ddxi[i][j] * x.Ddeta[i][j];
    }
};

class scalarField {
public:
    const gridMetrics* g = nullptr;
    int lines = 0, columns = 0;
    Grid2D F;

    void constructor(const gridMetrics& grid) {
        g = &grid; lines = grid.lines; columns = grid.columns;
        F.assign(lines, std::vector<double>(columns, 0.0));
    }

    double Ddxi(int i, int j) const {
        if (j == 0 || j == columns - 1) return (F[i][1] - F[i][columns - 2]) / 2.0;
        return (F[i][j + 1] - F[i][j - 1]) / 2.0;
    }
    double Ddxi_forward(int i, int j) const {
        if (j == columns - 1) return F[i][1] - F[i][j];
        return F[i][j + 1] - F[i][j];
    }
    double Ddxi_rearward(int i, int j) const {
        if (j == 0) return F[i][j] - F[i][columns - 2];
        return F[i][j] - F[i][j - 1];
    }
    double Ddeta(int i, int j) const {
        if (i == 0) return (F[i + 2][j] - 4.0 * F[i + 1][j] + 3.0 * F[i][j]) / 2.0;
        if (i == lines - 1) return (-F[i - 2][j] + 4.0 * F[i - 1][j] - 3.0 * F[i][j]) / 2.0;
        return (F[i - 1][j] - F[i + 1][j]) / 2.0;
    }
    double Ddeta_forward(int i, int j) const {
        if (i == 0) return F[i][j] - F[i + 1][j];
        return F[i - 1][j] - F[i][j];
    }
    double Ddeta_rearward(int i, int j) const {
        if (i == lines - 1) return F[i - 1][j] - F[i][j];
        return F[i][j] - F[i + 1][j];
    }

    double Ddx(int i, int j) const { return (Ddxi(i, j) * g->y.Ddeta[i][j] - Ddeta(i, j) * g->y.Ddxi[i][j]) / g->J[i][j]; }
    double Ddy(int i, int j) const { return (Ddeta(i, j) * g->x.Ddxi[i][j] - Ddxi(i, j) * g->x.Ddeta[i][j]) / g->J[i][j]; }
    double Ddx_forward(int i, int j) const { return (Ddxi_forward(i, j) * g->y.Ddeta[i][j] - Ddeta_forward(i, j) * g->y.Ddxi[i][j]) / g->J[i][j]; }
    double Ddy_forward(int i, int j) const { return (Ddeta_forward(i, j) * g->x.Ddxi[i][j] - Ddxi_forward(i, j) * g->x.Ddeta[i][j]) / g->J[i][j]; }
    double Ddx_rearward(int i, int j) const { return (Ddxi_rearward(i, j) * g->y.Ddeta[i][j] - Ddeta_rearward(i, j) * g->y.Ddxi[i][j]) / g->J[i][j]; }
    double Ddy_rearward(int i, int j) const { return (Ddeta_rearward(i, j) * g->x.Ddxi[i][j] - Ddxi_rearward(i, j) * g->x.Ddeta[i][j]) / g->J[i][j]; }
};

class flowSim {
public:
    scalarField rho, u, v, p, e, T;
    scalarField mu, mut, muE, kE;          
    scalarField tke, eps;                  
    scalarField TAUxx, TAUxy, TAUyy, qx, qy;
    scalarField U1, U2, U3, U4, U5, U6;
    scalarField E1, E2, E3, E4, E5, E6;
    scalarField F1, F2, F3, F4, F5, F6;

    double T0 = 0, mu0 = 0, Pr = 0, Cp = 0, Cv = 0, R = 0;
    bool useTurbulence = true;
    std::vector<double> wTau, wTx, wTy, wNx, wNy, wYplus;

    void constructor(double T1, double mu1, double Pr0, double Cp0, double Cv0, const gridMetrics& grid, bool turb) {
        T0 = T1; mu0 = mu1; Pr = Pr0; Cp = Cp0; Cv = Cv0; R = Cp - Cv; useTurbulence = turb;
        for (scalarField* f : {&rho, &u, &v, &p, &e, &T, &mu, &mut, &muE, &kE, &tke, &eps, &TAUxx, &TAUxy, &TAUyy, &qx, &qy,
                               &U1, &U2, &U3, &U4, &U5, &U6, &E1, &E2, &E3, &E4, &E5, &E6, &F1, &F2, &F3, &F4, &F5, &F6})
            f->constructor(grid);
        for (auto* w : {&wTau, &wTx, &wTy, &wNx, &wNy, &wYplus}) w->assign(grid.columns, 0.0);
    }

    void dynamicViscosity() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++)
                mu.F[i][j] = mu0 * std::pow(std::max(T.F[i][j], 10.0) / T0, 1.5) * (T0 + 110.0) / (std::max(T.F[i][j], 10.0) + 110.0);
    }

    void turbulenceProperties() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double mt = 0.0;
                if (useTurbulence) {
                    mt = rho.F[i][j] * KE::Cmu * tke.F[i][j] * tke.F[i][j] / std::max(eps.F[i][j], KE::epsMin);
                    mt = std::min(mt, KE::muRatioMax * mu.F[i][j]);
                }
                mut.F[i][j] = mt; muE.F[i][j] = mu.F[i][j] + mt; kE.F[i][j] = Cp * (mu.F[i][j] / Pr + mt / KE::Prt);
            }
    }

    void calculateU() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double r = rho.F[i][j], uu = u.F[i][j], vv = v.F[i][j];
                U1.F[i][j] = r; U2.F[i][j] = r * uu; U3.F[i][j] = r * vv;
                U4.F[i][j] = r * (e.F[i][j] + 0.5 * (uu * uu + vv * vv));
                U5.F[i][j] = r * tke.F[i][j]; U6.F[i][j] = r * eps.F[i][j];
            }
    }

    void decodeU() {
        const double MIN_RHO = 1e-4, MIN_T = 10.0;
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double r = U1.F[i][j];
                if (!std::isfinite(r) || r < MIN_RHO) { r = MIN_RHO; U1.F[i][j] = r; }
                rho.F[i][j] = r;

                u.F[i][j] = U2.F[i][j] / r; v.F[i][j] = U3.F[i][j] / r;
                if (!std::isfinite(u.F[i][j])) u.F[i][j] = 0.0;
                if (!std::isfinite(v.F[i][j])) v.F[i][j] = 0.0;

                double e_int = U4.F[i][j] / r - 0.5 * (u.F[i][j] * u.F[i][j] + v.F[i][j] * v.F[i][j]);
                if (!std::isfinite(e_int) || e_int < MIN_T * Cv) {
                    e_int = MIN_T * Cv;
                    U4.F[i][j] = r * (e_int + 0.5 * (u.F[i][j] * u.F[i][j] + v.F[i][j] * v.F[i][j]));
                }
                e.F[i][j] = e_int; T.F[i][j] = e_int / Cv; p.F[i][j] = r * R * T.F[i][j];

                double k_val = U5.F[i][j] / r;
                if (!std::isfinite(k_val) || k_val < KE::tkeMin) k_val = KE::tkeMin;
                tke.F[i][j] = k_val; U5.F[i][j] = r * k_val;

                double eps_val = U6.F[i][j] / r;
                if (!std::isfinite(eps_val) || eps_val < KE::epsMin) eps_val = KE::epsMin;
                eps.F[i][j] = eps_val; U6.F[i][j] = r * eps_val;
            }
    }

    void computeWallShear() {
        const int i = lines() - 1, i1 = lines() - 2;
        const Grid2D& gx = u.g->x.F; const Grid2D& gy = u.g->y.F;
        #pragma omp parallel for
        for (int j = 0; j < columns(); j++) {
            double nx = gx[i1][j] - gx[i][j], ny = gy[i1][j] - gy[i][j], y1 = std::hypot(nx, ny);
            nx /= y1; ny /= y1;
            double u1 = u.F[i1][j], v1 = v.F[i1][j], un = u1 * nx + v1 * ny, utx = u1 - un * nx, uty = v1 - un * ny, Ut = std::hypot(utx, uty);
            double r1 = rho.F[i1][j], mu1 = mu.F[i1][j];
            double ustar = useTurbulence ? std::pow(KE::Cmu, 0.25) * std::sqrt(std::max(tke.F[i1][j], KE::tkeMin)) : 0.0;
            double ystar = r1 * ustar * y1 / mu1;
            double tw = 0.0, tx = 0.0, ty = 0.0;
            if (Ut > 1.0e-12) {
                tx = utx / Ut; ty = uty / Ut;
                if (ystar > KE::yplusLam) tw = r1 * ustar * Ut * KE::kappa / std::log(KE::E * ystar);
                else                      tw = mu1 * Ut / y1;
            }
            wTau[j] = tw; wTx[j] = tx; wTy[j] = ty; wNx[j] = nx; wNy[j] = ny; wYplus[j] = ystar;
        }
    }

    void calculateE(const std::string& step) {
        const bool pred = (step == "predictor");
        const double turb = useTurbulence ? 1.0 : 0.0;
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double ux = pred ? u.Ddx_rearward(i, j) : u.Ddx_forward(i, j), vx = pred ? v.Ddx_rearward(i, j) : v.Ddx_forward(i, j);
                double Tx = pred ? T.Ddx_rearward(i, j) : T.Ddx_forward(i, j);
                double kx = pred ? tke.Ddx_rearward(i, j) : tke.Ddx_forward(i, j), ex = pred ? eps.Ddx_rearward(i, j) : eps.Ddx_forward(i, j);
                double uy = u.Ddy(i, j), vy = v.Ddy(i, j);
                double me = muE.F[i][j], rk = (2.0 / 3.0) * rho.F[i][j] * tke.F[i][j] * turb;
                TAUxx.F[i][j] = (2.0 / 3.0) * me * (2.0 * ux - vy) - rk;
                TAUxy.F[i][j] = me * (uy + vx);
                TAUyy.F[i][j] = (2.0 / 3.0) * me * (2.0 * vy - ux) - rk;
                qx.F[i][j] = -kE.F[i][j] * Tx;
                double muk = mu.F[i][j] + mut.F[i][j] / KE::sigk, mue = mu.F[i][j] + mut.F[i][j] / KE::sige;
                E5.F[i][j] = rho.F[i][j] * u.F[i][j] * tke.F[i][j] - muk * kx;
                E6.F[i][j] = rho.F[i][j] * u.F[i][j] * eps.F[i][j] - mue * ex;
            }
        applyWallStress();
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double r = rho.F[i][j], uu = u.F[i][j], vv = v.F[i][j];
                E1.F[i][j] = r * uu;
                E2.F[i][j] = r * uu * uu + p.F[i][j] - TAUxx.F[i][j];
                E3.F[i][j] = r * uu * vv - TAUxy.F[i][j];
                E4.F[i][j] = (U4.F[i][j] + p.F[i][j] - TAUxx.F[i][j]) * uu + qx.F[i][j] - vv * TAUxy.F[i][j];
            }
    }

    void calculateF(const std::string& step) {
        const bool pred = (step == "predictor");
        const double turb = useTurbulence ? 1.0 : 0.0;
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double uy = pred ? u.Ddy_rearward(i, j) : u.Ddy_forward(i, j), vy = pred ? v.Ddy_rearward(i, j) : v.Ddy_forward(i, j);
                double Ty = pred ? T.Ddy_rearward(i, j) : T.Ddy_forward(i, j);
                double ky = pred ? tke.Ddy_rearward(i, j) : tke.Ddy_forward(i, j), ey = pred ? eps.Ddy_rearward(i, j) : eps.Ddy_forward(i, j);
                double ux = u.Ddx(i, j), vx = v.Ddx(i, j);
                double me = muE.F[i][j], rk = (2.0 / 3.0) * rho.F[i][j] * tke.F[i][j] * turb;
                TAUxx.F[i][j] = (2.0 / 3.0) * me * (2.0 * ux - vy) - rk;
                TAUxy.F[i][j] = me * (uy + vx);
                TAUyy.F[i][j] = (2.0 / 3.0) * me * (2.0 * vy - ux) - rk;
                qy.F[i][j] = -kE.F[i][j] * Ty;
                double muk = mu.F[i][j] + mut.F[i][j] / KE::sigk, mue = mu.F[i][j] + mut.F[i][j] / KE::sige;
                F5.F[i][j] = rho.F[i][j] * v.F[i][j] * tke.F[i][j] - muk * ky;
                F6.F[i][j] = rho.F[i][j] * v.F[i][j] * eps.F[i][j] - mue * ey;
            }
        applyWallStress();
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++)
            for (int j = 0; j < u.columns; j++) {
                double r = rho.F[i][j], uu = u.F[i][j], vv = v.F[i][j];
                F1.F[i][j] = r * vv;
                F2.F[i][j] = r * uu * vv - TAUxy.F[i][j];
                F3.F[i][j] = r * vv * vv + p.F[i][j] - TAUyy.F[i][j];
                F4.F[i][j] = (U4.F[i][j] + p.F[i][j] - TAUyy.F[i][j]) * vv + qy.F[i][j] - uu * TAUxy.F[i][j];
            }
    }

    int lines() const { return u.lines; }
    int columns() const { return u.columns; }

private:
    void applyWallStress() {
        if (!useTurbulence) return;
        const int i = lines() - 1;
        for (int j = 0; j < columns(); j++) {
            TAUxx.F[i][j] = 2.0 * wTau[j] * wNx[j] * wTx[j];
            TAUyy.F[i][j] = 2.0 * wTau[j] * wNy[j] * wTy[j];
            TAUxy.F[i][j] = wTau[j] * (wNx[j] * wTy[j] + wNy[j] * wTx[j]);
        }
    }
};

class boundaryConditionsClass {
public:
    double rho_inf, u_inf, v_inf, p_inf, T_wall, tke_inf, eps_inf;
    boundaryConditionsClass(double rho, double u, double v, double p, double Tw, double k, double e)
        : rho_inf(rho), u_inf(u), v_inf(v), p_inf(p), T_wall(Tw), tke_inf(k), eps_inf(e) {}

    void apply(flowSim& f) const {
        const int L = f.lines(), C = f.columns();
        const Grid2D& gx = f.u.g->x.F; const Grid2D& gy = f.u.g->y.F;
        const double Tinf = p_inf / (f.R * rho_inf);

        for (int j = 0; j < C; j++) {
            f.rho.F[0][j] = rho_inf; f.u.F[0][j] = u_inf; f.v.F[0][j] = v_inf;
            f.p.F[0][j] = p_inf;     f.T.F[0][j] = Tinf;  f.e.F[0][j] = f.Cv * Tinf;
            double nx = gx[0][j] - gx[1][j], ny = gy[0][j] - gy[1][j];
            if ((u_inf * nx + v_inf * ny) >= 0.0) { f.tke.F[0][j] = f.tke.F[1][j]; f.eps.F[0][j] = f.eps.F[1][j]; }
            else { f.tke.F[0][j] = tke_inf; f.eps.F[0][j] = eps_inf; }
        }

        const int iw = L - 1, i1 = L - 2;
        for (int j = 0; j < C; j++) {
            f.u.F[iw][j] = 0.0; f.v.F[iw][j] = 0.0;
            f.T.F[iw][j] = T_wall; f.e.F[iw][j] = f.Cv * T_wall;
            
            // Extrapolação de Pressão estável de 1a Ordem
            f.p.F[iw][j] = std::max(f.p.F[i1][j], 10.0);
            f.rho.F[iw][j] = f.p.F[iw][j] / (f.R * T_wall);
            
            if (f.useTurbulence) {
                double y1 = std::hypot(gx[i1][j] - gx[iw][j], gy[i1][j] - gy[iw][j]);
                double k1 = std::max(f.tke.F[i1][j], KE::tkeMin);
                f.eps.F[i1][j] = std::pow(KE::Cmu, 0.75) * std::pow(k1, 1.5) / (KE::kappa * y1);
                f.tke.F[iw][j] = k1; f.eps.F[iw][j] = f.eps.F[i1][j];
            }
        }
    }
};

struct SimParams {
    double cfl = 0.4;
    double dissipationCoeff = 0.002;
    bool   useTurbulence = true;
    int    maxIterations = 500000;
    int    reportEvery = 100;
    int    writeEvery = 200;
};

struct StepInfo { double dt = 0; double res[6] = {0, 0, 0, 0, 0, 0}; };

double computeTimeStep(flowSim& f, double cfl, Grid2D& dtField) {
    const int L = f.lines(), C = f.columns();
    const double gamma = f.Cp / f.Cv;
    const Grid2D& gx = f.u.g->x.F; const Grid2D& gy = f.u.g->y.F;
    double dtMin = std::numeric_limits<double>::max();
    #pragma omp parallel for reduction(min : dtMin)
    for (int i = 0; i < L; i++)
        for (int j = 0; j < C; j++) {
            int jp = (j == C - 1) ? 1 : j + 1, jm = (j == 0) ? C - 2 : j - 1;
            int ip = std::min(i + 1, L - 1), im = std::max(i - 1, 0);
            double ds1 = std::max(0.5 * std::hypot(gx[i][jp] - gx[i][jm], gy[i][jp] - gy[i][jm]), 1e-12);
            double ds2 = std::max(0.5 * std::hypot(gx[ip][j] - gx[im][j], gy[ip][j] - gy[im][j]), 1e-12);
            double c = std::sqrt(gamma * f.R * std::max(f.T.F[i][j], 1.0));
            double nu = f.muE.F[i][j] / std::max(f.rho.F[i][j], 1e-8);
            double den = (std::fabs(f.u.F[i][j]) + c) / ds1 + (std::fabs(f.v.F[i][j]) + c) / ds2 + 2.0 * nu * (1.0 / (ds1 * ds1) + 1.0 / (ds2 * ds2));
            double dtl = cfl / (den + 1e-30);
            dtField[i][j] = dtl;
            dtMin = std::min(dtMin, dtl);
        }
    return dtMin;
}

void addArtificialDissipation(flowSim& f, double eps) {
    const int L = f.lines(), C = f.columns();
    for (scalarField* Uk : {&f.U1, &f.U2, &f.U3, &f.U4, &f.U5, &f.U6}) {
        const Grid2D ref = Uk->F;
        #pragma omp parallel for
        for (int i = 1; i < L - 1; i++)
            for (int j = 0; j < C; j++) {
                int jp = (j == C - 1) ? 1 : j + 1, jm = (j == 0) ? C - 2 : j - 1;
                double lap = ref[i + 1][j] + ref[i - 1][j] + ref[i][jp] + ref[i][jm] - 4.0 * ref[i][j];
                Uk->F[i][j] += eps * lap;
            }
    }
}

void applyTurbulenceSources(flowSim& old, flowSim& target, Grid2D& dtField) {
    const int L = old.lines(), C = old.columns();
    #pragma omp parallel for
    for (int i = 0; i < L; i++)
        for (int j = 0; j < C; j++) {
            double dtc = dtField[i][j];
            double r = old.rho.F[i][j], k = old.tke.F[i][j], ep = old.eps.F[i][j];
            double ux = old.u.Ddx(i, j), uy = old.u.Ddy(i, j), vx = old.v.Ddx(i, j), vy = old.v.Ddy(i, j), div = ux + vy;
            double S2 = 2.0 * (ux * ux + vy * vy) + (uy + vx) * (uy + vx) - (2.0 / 3.0) * div * div;
            double Pk = old.mut.F[i][j] * S2 - (2.0 / 3.0) * r * k * div;
            Pk = std::min(std::max(Pk, 0.0), KE::prodLimit * r * ep);
            double invT = ep / k;
            target.U5.F[i][j] = (target.U5.F[i][j] + dtc * Pk) / (1.0 + dtc * invT);
            target.U6.F[i][j] = (target.U6.F[i][j] + dtc * KE::C1 * invT * Pk) / (1.0 + dtc * KE::C2 * invT);
        }
}

StepInfo macCormack(flowSim& flow, flowSim& pred, flowSim& corr, const boundaryConditionsClass& bc, const SimParams& sp, Grid2D& dtField) {
    StepInfo info;
    const int L = flow.lines(), C = flow.columns();

    flow.dynamicViscosity(); flow.turbulenceProperties(); flow.calculateU(); flow.computeWallShear();
    flow.calculateE("predictor"); flow.calculateF("predictor");
    info.dt = computeTimeStep(flow, sp.cfl, dtField);

    #pragma omp parallel for
    for (int i = 0; i < L; i++)
        for (int j = 0; j < C; j++) {
            double d = dtField[i][j];
            pred.U1.F[i][j] = flow.U1.F[i][j] - (flow.E1.Ddx_forward(i, j) + flow.F1.Ddy_forward(i, j)) * d;
            pred.U2.F[i][j] = flow.U2.F[i][j] - (flow.E2.Ddx_forward(i, j) + flow.F2.Ddy_forward(i, j)) * d;
            pred.U3.F[i][j] = flow.U3.F[i][j] - (flow.E3.Ddx_forward(i, j) + flow.F3.Ddy_forward(i, j)) * d;
            pred.U4.F[i][j] = flow.U4.F[i][j] - (flow.E4.Ddx_forward(i, j) + flow.F4.Ddy_forward(i, j)) * d;
            pred.U5.F[i][j] = flow.U5.F[i][j] - (flow.E5.Ddx_forward(i, j) + flow.F5.Ddy_forward(i, j)) * d;
            pred.U6.F[i][j] = flow.U6.F[i][j] - (flow.E6.Ddx_forward(i, j) + flow.F6.Ddy_forward(i, j)) * d;
        }
    pred.decodeU(); bc.apply(pred);

    pred.dynamicViscosity(); pred.turbulenceProperties(); pred.calculateU(); pred.computeWallShear();
    pred.calculateE("corrector"); pred.calculateF("corrector");

    #pragma omp parallel for
    for (int i = 0; i < L; i++)
        for (int j = 0; j < C; j++) {
            double d = dtField[i][j];
            corr.U1.F[i][j] = 0.5 * (pred.U1.F[i][j] + flow.U1.F[i][j] - (pred.E1.Ddx_rearward(i, j) + pred.F1.Ddy_rearward(i, j)) * d);
            corr.U2.F[i][j] = 0.5 * (pred.U2.F[i][j] + flow.U2.F[i][j] - (pred.E2.Ddx_rearward(i, j) + pred.F2.Ddy_rearward(i, j)) * d);
            corr.U3.F[i][j] = 0.5 * (pred.U3.F[i][j] + flow.U3.F[i][j] - (pred.E3.Ddx_rearward(i, j) + pred.F3.Ddy_rearward(i, j)) * d);
            corr.U4.F[i][j] = 0.5 * (pred.U4.F[i][j] + flow.U4.F[i][j] - (pred.E4.Ddx_rearward(i, j) + pred.F4.Ddy_rearward(i, j)) * d);
            corr.U5.F[i][j] = 0.5 * (pred.U5.F[i][j] + flow.U5.F[i][j] - (pred.E5.Ddx_rearward(i, j) + pred.F5.Ddy_rearward(i, j)) * d);
            corr.U6.F[i][j] = 0.5 * (pred.U6.F[i][j] + flow.U6.F[i][j] - (pred.E6.Ddx_rearward(i, j) + pred.F6.Ddy_rearward(i, j)) * d);
        }

    if (sp.dissipationCoeff > 0.0) addArtificialDissipation(corr, sp.dissipationCoeff);
    if (flow.useTurbulence) applyTurbulenceSources(flow, corr, dtField);

    double r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0, r6 = 0;
    #pragma omp parallel for reduction(+ : r1, r2, r3, r4, r5, r6)
    for (int i = 1; i < L - 1; i++)
        for (int j = 0; j < C; j++) {
            double d = dtField[i][j];
            double a1 = (corr.U1.F[i][j] - flow.U1.F[i][j]) / d, a2 = (corr.U2.F[i][j] - flow.U2.F[i][j]) / d;
            double a3 = (corr.U3.F[i][j] - flow.U3.F[i][j]) / d, a4 = (corr.U4.F[i][j] - flow.U4.F[i][j]) / d;
            double a5 = (corr.U5.F[i][j] - flow.U5.F[i][j]) / d, a6 = (corr.U6.F[i][j] - flow.U6.F[i][j]) / d;
            r1 += a1*a1; r2 += a2*a2; r3 += a3*a3; r4 += a4*a4; r5 += a5*a5; r6 += a6*a6;
        }
    const double N = static_cast<double>(L - 2) * C;
    info.res[0] = std::sqrt(r1 / N); info.res[1] = std::sqrt(r2 / N); info.res[2] = std::sqrt(r3 / N);
    info.res[3] = std::sqrt(r4 / N); info.res[4] = std::sqrt(r5 / N); info.res[5] = std::sqrt(r6 / N);

    corr.decodeU(); bc.apply(corr);
    return info;
}

// Leitor de Malha via Regex
void readGridFile(const std::string& fileName, scalarField0& x, scalarField0& y) {
    std::ifstream malha(fileName);
    if (!malha.is_open()) throw std::runtime_error("Nao foi possivel abrir a malha: " + fileName);
    
    int lin = 0, col = 0;
    std::string linha;
    while (std::getline(malha, linha)) {
        std::regex pattern(R"([Ii]\s*=\s*(\d+).*?[Jj]\s*=\s*(\d+))");
        std::smatch match;
        if (std::regex_search(linha, match, pattern)) { 
            col = std::stoi(match[1]); lin = std::stoi(match[2]); break; 
        }
    }
    if (lin < 4 || col < 4) throw std::runtime_error("Impossivel determinar I e J na malha.");
    std::cout << "Malha Lida Com Sucesso: " << lin << " (I/Radial) x " << col << " (J/Circunf.)\n";
    
    x.constructor(lin, col); y.constructor(lin, col);
    int i1 = 0, i2 = 0;
    while (std::getline(malha, linha)) {
        std::stringstream ss(linha); double tx, ty;
        while (ss >> tx >> ty) {
            x.F[i1][i2] = tx; y.F[i1][i2] = ty;
            if (++i2 == col) { i2 = 0; if (++i1 == lin) break; }
        }
        if (i1 == lin) break;
    }
    x.calculateMetrics(); y.calculateMetrics();
}

// Escrita VTK Blindada em Vídeo
void writeVTKfile(flowSim& f, const std::string& fileName) {
    f.dynamicViscosity(); f.turbulenceProperties();
    std::ofstream out(fileName, std::ios::out | std::ios::trunc);
    if (!out.is_open()) return;
    
    const int L = f.lines(), C = f.columns();
    out << "# vtk DataFile Version 3.0\nVideo 2D RANS MacCormack\nASCII\nDATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << C << " " << L << " 1\nPOINTS " << C * L << " float\n";
    
    for (int i = 0; i < L; i++)
        for (int j = 0; j < C; j++) out << f.u.g->x.F[i][j] << " " << f.u.g->y.F[i][j] << " 0.0\n";
            
    out << "\nPOINT_DATA " << C * L << "\n";
    auto sanitize = [](double v) { return (std::isfinite(v)) ? static_cast<float>(v) : 0.0f; };
    
    auto scalar = [&](const char* name, auto getter) {
        out << "SCALARS " << name << " float 1\nLOOKUP_TABLE default\n";
        for (int i = 0; i < L; i++)
            for (int j = 0; j < C; j++) out << sanitize(getter(i, j)) << "\n";
    };
    
    scalar("Pressure", [&](int i, int j) { return f.p.F[i][j]; });
    scalar("Density", [&](int i, int j) { return f.rho.F[i][j]; });
    scalar("Temperature", [&](int i, int j) { return f.T.F[i][j]; });
    scalar("Mach", [&](int i, int j) { return std::hypot(f.u.F[i][j], f.v.F[i][j]) / std::sqrt(1.4 * f.R * std::max(f.T.F[i][j], 1.0)); });
    scalar("TKE", [&](int i, int j) { return f.tke.F[i][j]; });
    scalar("Epsilon", [&](int i, int j) { return f.eps.F[i][j]; });
    
    out << "VECTORS Velocity float\n";
    for (int i = 0; i < L; i++)
        for (int j = 0; j < C; j++) out << sanitize(f.u.F[i][j]) << " " << sanitize(f.v.F[i][j]) << " 0.0\n";
}

int main(int argc, char** argv) {
    std::string gridFileName = (argc > 1) ? argv[1] : "cilindro.grid";
    SimParams sp;
    if (argc > 2) sp.maxIterations = std::stoi(argv[2]);
    if (argc > 3) sp.dissipationCoeff = std::stod(argv[3]);
    if (argc > 4) sp.cfl = std::stod(argv[4]);
    if (argc > 5) sp.useTurbulence = (std::stoi(argv[5]) != 0);

    const double target_cfl = sp.cfl;
    const double target_dissip = sp.dissipationCoeff;

    const double T0 = 288.16, mu0 = 1.7894e-5, Pr = 0.71, Cp = 1004.0, Cv = 717.0, R = Cp - Cv;
    const double rho_inf = 1.225, u_inf = 102.0, v_inf = 0.0, p_inf = 1.01325e5;
    const double T_inf = p_inf / (R * rho_inf), T_wall = T_inf;
    const double mu_inf = mu0 * std::pow(T_inf / T0, 1.5) * (T0 + 110.0) / (T_inf + 110.0);
    const double U_inf = std::hypot(u_inf, v_inf);
    const double tke_inf = 1.5 * std::pow(0.01 * U_inf, 2);                        
    const double eps_inf = rho_inf * KE::Cmu * tke_inf * tke_inf / (10.0 * mu_inf);

    try {
        scalarField0 xr, yr;
        readGridFile(gridFileName, xr, yr);
        gridMetrics grid; grid.build(xr, yr);

        flowSim initialFlow, predictedFlow, correctedFlow;
        for (flowSim* f : {&initialFlow, &predictedFlow, &correctedFlow}) f->constructor(T0, mu0, Pr, Cp, Cv, grid, sp.useTurbulence);

        boundaryConditionsClass bc(rho_inf, u_inf, v_inf, p_inf, T_wall, tke_inf, eps_inf);

        for (int i = 0; i < grid.lines; i++)
            for (int j = 0; j < grid.columns; j++) {
                correctedFlow.rho.F[i][j] = rho_inf; correctedFlow.u.F[i][j] = u_inf;
                correctedFlow.v.F[i][j] = v_inf;     correctedFlow.p.F[i][j] = p_inf;
                correctedFlow.T.F[i][j] = T_inf;     correctedFlow.e.F[i][j] = Cv * T_inf;
                correctedFlow.tke.F[i][j] = tke_inf; correctedFlow.eps.F[i][j] = eps_inf;
            }
        bc.apply(correctedFlow);

        Grid2D dtField(grid.lines, std::vector<double>(grid.columns, 0.0));
        std::cout << "MacCormack Viscoso Restaurado | Mach = " << U_inf / std::sqrt(1.4 * R * T_inf) << "\n";

        double resRef[6] = {0};
        auto tStart = std::chrono::high_resolution_clock::now();
        
        for (int it = 1; it <= sp.maxIterations; it++) {
            std::swap(initialFlow, correctedFlow);     
            
            // -------------------------------------------------------------
            // SOFT START DUPLO (Rampa de CFL + Rampa de Dissipação)
            // Impede a explosão térmica transônica nas primeiras iterações
            // -------------------------------------------------------------
            if (it <= 3000) {
                sp.cfl = 0.01 + (target_cfl - 0.01) * ((double)it / 3000.0);
                sp.dissipationCoeff = 0.05 - (0.05 - target_dissip) * ((double)it / 3000.0); // Começa em 0.05 e decai
            } else {
                sp.cfl = target_cfl;
                sp.dissipationCoeff = target_dissip;
            }

            StepInfo info = macCormack(initialFlow, predictedFlow, correctedFlow, bc, sp, dtField);

            for (int k = 0; k < 6; k++) resRef[k] = std::max(resRef[k], std::max(info.res[k], 1e-300));

            if (it % sp.reportEvery == 0 || it == 1) {
                auto t = std::chrono::high_resolution_clock::now();
                double secs = std::chrono::duration<double>(t - tStart).count();
                std::cout << std::scientific << std::setprecision(3) << "it=" << it 
                          << " cfl=" << sp.cfl << " diss=" << sp.dissipationCoeff << " res=";
                for (int k = 0; k < 6; k++) std::cout << info.res[k] / resRef[k] << (k < 5 ? "," : "");
                std::cout << "  t=" << std::fixed << std::setprecision(2) << secs << "s\n" << std::defaultfloat;
            }
            if (it % sp.writeEvery == 0 || it == 1) {
                char buffer[256];
                snprintf(buffer, sizeof(buffer), "escoamento_%08d.vtk", it);
                writeVTKfile(correctedFlow, std::string(buffer));
            }
        }
    } catch (const std::exception& ex) {
        std::cerr << "Erro fatal: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
