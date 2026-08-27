// A equacao governante para o potencial de velocidade phi (escoamento
// incompressivel e irrotacional) e a equacao de Laplace:
//
//        div( grad(phi) ) = 0
//
// Em coordenadas polares (r, theta), com o cilindro de raio R centrado
// na origem e o dominio se estendendo ate um raio externo Rmax:
//
//   d2phi/dr2 + (1/r) dphi/dr + (1/r^2) d2phi/dtheta2 = 0
//
// CONDICOES DE CONTORNO
// ----------------------
//  - Contorno externo (r = Rmax): Dirichlet, phi = solucao analitica
//    (escoamento uniforme + dipolo). Isso e consistente pois a solucao
//    exata satisfaz Laplace em todo o dominio r > R; assim o unico erro
//    introduzido no interior do dominio vem da discretizacao, nao do
//    truncamento do dominio. Isso permite avaliar a ordem de convergencia
//    do esquema numerico de forma limpa.
//  - Superficie do cilindro (r = R): Neumann, dphi/dr = 0 (impenetrabilidade,
//    velocidade normal nula), imposta via ponto fantasma phi(-1,j)=phi(1,j).
//  - Direcao theta: periodica (o escoamento da volta completa no cilindro).
//
// SOLUCAO ANALITICA (escoamento uniforme Uinf + dipolo de intensidade
// mu = 2*pi*Uinf*R^2, resultando no cilindro de raio R como linha de
// corrente):
//
//   phi(r,theta)  =  Uinf*(r + R^2/r)*cos(theta)
//   u_r(r,theta)  =  Uinf*(1 - R^2/r^2)*cos(theta)
//   u_t(r,theta)  = -Uinf*(1 + R^2/r^2)*sin(theta)
//
// Na superficie (r=R): u_r = 0 e u_theta = -2*Uinf*sin(theta), logo
//   Cp(theta) = 1 - 4*sin^2(theta)
//
// SAIDAS
// -------
//  - Campo completo (r,theta,x,y,phi,u,v,Cp) numerico x analitico, para a
//    malha escolhida como referencia.
//  - Cp ao longo da superficie do cilindro, numerico x analitico.
//  - Estudo de convergencia de malha: para uma lista de malhas (Nr x Nt),
//    calcula os erros L2 e Linf de phi e de Cp na superficie.
//
// COMPILACAO
// -----------
//   g++ -O2 -std=c++17 potential_flow_cylinder.cpp -o potential_flow_cylinder
//
// =============================================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>

static const double PI = 3.14159265358979323846;

// -----------------------------------------------------------------------
// Estruturas de dados
// -----------------------------------------------------------------------
struct Grid {
    int Nr, Nt;                 // numero de pontos em r e theta
    double R, Rmax, Uinf;        // raio do cilindro, raio externo, veloc. de escoamento livre
    std::vector<double> r, theta;
    double dr, dtheta;
};

inline int IDX(int i, int j, int Nt) { return i * Nt + j; }

// -----------------------------------------------------------------------
// Solucao analitica
// -----------------------------------------------------------------------
double phi_exact(double r, double theta, double R, double Uinf) {
    return Uinf * (r + R * R / r) * std::cos(theta);
}

void vel_exact(double r, double theta, double R, double Uinf, double &u, double &v) {
    double ur = Uinf * (1.0 - R * R / (r * r)) * std::cos(theta);
    double ut = -Uinf * (1.0 + R * R / (r * r)) * std::sin(theta);
    u = ur * std::cos(theta) - ut * std::sin(theta);
    v = ur * std::sin(theta) + ut * std::cos(theta);
}

double cp_exact(double r, double theta, double R, double Uinf) {
    double u, v;
    vel_exact(r, theta, R, Uinf, u, v);
    return 1.0 - (u * u + v * v) / (Uinf * Uinf);
}

// -----------------------------------------------------------------------
// Solver: Laplace 2D em coordenadas polares via SOR
// -----------------------------------------------------------------------
struct Result {
    Grid g;
    std::vector<double> phi;
    int iters;
    double residual;
};

Result solve(int Nr, int Nt, double R, double Rmax, double Uinf,
             double omega = 1.85, double tol = 1e-11, int maxit = 400000) {
    Grid g;
    g.Nr = Nr; g.Nt = Nt; g.R = R; g.Rmax = Rmax; g.Uinf = Uinf;
    g.r.resize(Nr);
    g.theta.resize(Nt);
    g.dr = (Rmax - R) / (Nr - 1);
    g.dtheta = 2.0 * PI / Nt;
    for (int i = 0; i < Nr; ++i) g.r[i] = R + i * g.dr;
    for (int j = 0; j < Nt; ++j) g.theta[j] = j * g.dtheta;

    std::vector<double> phi(Nr * Nt, 0.0);

    // Chute inicial: potencial de escoamento uniforme
    for (int i = 0; i < Nr; ++i)
        for (int j = 0; j < Nt; ++j)
            phi[IDX(i, j, Nt)] = Uinf * g.r[i] * std::cos(g.theta[j]);

    // Contorno externo: Dirichlet = solucao analitica
    for (int j = 0; j < Nt; ++j)
        phi[IDX(Nr - 1, j, Nt)] = phi_exact(Rmax, g.theta[j], R, Uinf);

    double dr2 = g.dr * g.dr;
    double dt2 = g.dtheta * g.dtheta;

    int it = 0;
    double res = 1e30;
    for (it = 0; it < maxit; ++it) {
        res = 0.0;

        // --- i = 0 : superficie do cilindro, Neumann dphi/dr = 0 ---
        // ponto fantasma: phi(-1,j) = phi(1,j)
        {
            double ri = g.r[0];
            double diag = 2.0 / dr2 + 2.0 / (ri * ri * dt2);
            for (int j = 0; j < Nt; ++j) {
                int jp = (j + 1) % Nt, jm = (j - 1 + Nt) % Nt;
                double rhs = (2.0 * phi[IDX(1, j, Nt)]) / dr2
                           + (phi[IDX(0, jp, Nt)] + phi[IDX(0, jm, Nt)]) / (ri * ri * dt2);
                double newval = rhs / diag;
                double old = phi[IDX(0, j, Nt)];
                double updated = old + omega * (newval - old);
                phi[IDX(0, j, Nt)] = updated;
                res = std::max(res, std::fabs(updated - old));
            }
        }

        // --- pontos internos 1 .. Nr-2 ---
        for (int i = 1; i < Nr - 1; ++i) {
            double ri = g.r[i];
            double A = 1.0 / dr2 + 1.0 / (2.0 * ri * g.dr);   // coef. phi(i+1,j)
            double B = 1.0 / dr2 - 1.0 / (2.0 * ri * g.dr);   // coef. phi(i-1,j)
            double C = 1.0 / (ri * ri * dt2);                  // coef. phi(i,j+-1)
            double diag = 2.0 / dr2 + 2.0 * C;
            for (int j = 0; j < Nt; ++j) {
                int jp = (j + 1) % Nt, jm = (j - 1 + Nt) % Nt;
                double rhs = A * phi[IDX(i + 1, j, Nt)] + B * phi[IDX(i - 1, j, Nt)]
                           + C * (phi[IDX(i, jp, Nt)] + phi[IDX(i, jm, Nt)]);
                double newval = rhs / diag;
                double old = phi[IDX(i, j, Nt)];
                double updated = old + omega * (newval - old);
                phi[IDX(i, j, Nt)] = updated;
                res = std::max(res, std::fabs(updated - old));
            }
        }

        // i = Nr-1 permanece fixo (Dirichlet)

        if (res < tol) { ++it; break; }
    }

    Result out;
    out.g = g;
    out.phi = phi;
    out.iters = it;
    out.residual = res;
    return out;
}

// -----------------------------------------------------------------------
// Pos-processamento: velocidades (u,v) e Cp 
// -----------------------------------------------------------------------
void compute_velocities(const Result &res, std::vector<double> &U, std::vector<double> &V) {
    const Grid &g = res.g;
    int Nr = g.Nr, Nt = g.Nt;
    U.assign(Nr * Nt, 0.0);
    V.assign(Nr * Nt, 0.0);

    for (int i = 0; i < Nr; ++i) {
        double ri = g.r[i];
        for (int j = 0; j < Nt; ++j) {
            int jp = (j + 1) % Nt, jm = (j - 1 + Nt) % Nt;

            double dphidr;
            if (i == 0) {
                // dphi/dr = 0 por construcao (condicao de parede)
                dphidr = 0.0;
            } else if (i == Nr - 1) {
                // diferenca atrasada (contorno externo)
                dphidr = (res.phi[IDX(i, j, Nt)] - res.phi[IDX(i - 1, j, Nt)]) / g.dr;
            } else {
                dphidr = (res.phi[IDX(i + 1, j, Nt)] - res.phi[IDX(i - 1, j, Nt)]) / (2.0 * g.dr);
            }

            double dphidtheta = (res.phi[IDX(i, jp, Nt)] - res.phi[IDX(i, jm, Nt)]) / (2.0 * g.dtheta);

            double ur = dphidr;
            double ut = dphidtheta / ri;

            double theta = g.theta[j];
            U[IDX(i, j, Nt)] = ur * std::cos(theta) - ut * std::sin(theta);
            V[IDX(i, j, Nt)] = ur * std::sin(theta) + ut * std::cos(theta);
        }
    }
}

// -----------------------------------------------------------------------
// Escrita de arquivos CSV
// -----------------------------------------------------------------------
void write_field_csv(const std::string &path, const Result &res,
                      const std::vector<double> &U, const std::vector<double> &V) {
    const Grid &g = res.g;
    std::ofstream f(path);
    f << std::setprecision(10);
    f << "i,j,r,theta,x,y,phi_num,phi_exact,u_num,v_num,u_exact,v_exact,cp_num,cp_exact\n";
    for (int i = 0; i < g.Nr; ++i) {
        for (int j = 0; j < g.Nt; ++j) {
            double r = g.r[i], theta = g.theta[j];
            double x = r * std::cos(theta), y = r * std::sin(theta);
            double phin = res.phi[IDX(i, j, g.Nt)];
            double phie = phi_exact(r, theta, g.R, g.Uinf);
            double un = U[IDX(i, j, g.Nt)], vn = V[IDX(i, j, g.Nt)];
            double ue, ve;
            vel_exact(r, theta, g.R, g.Uinf, ue, ve);
            double cpn = 1.0 - (un * un + vn * vn) / (g.Uinf * g.Uinf);
            double cpe = cp_exact(r, theta, g.R, g.Uinf);
            f << i << "," << j << "," << r << "," << theta << "," << x << "," << y << ","
              << phin << "," << phie << "," << un << "," << vn << "," << ue << "," << ve << ","
              << cpn << "," << cpe << "\n";
        }
    }
}

void write_surface_csv(const std::string &path, const Result &res,
                        const std::vector<double> &U, const std::vector<double> &V) {
    const Grid &g = res.g;
    std::ofstream f(path);
    f << std::setprecision(10);
    f << "theta,cp_num,cp_exact,u_num,v_num,u_exact,v_exact\n";
    int i = 0; // superficie do cilindro
    for (int j = 0; j < g.Nt; ++j) {
        double theta = g.theta[j];
        double un = U[IDX(i, j, g.Nt)], vn = V[IDX(i, j, g.Nt)];
        double ue, ve;
        vel_exact(g.R, theta, g.R, g.Uinf, ue, ve);
        double cpn = 1.0 - (un * un + vn * vn) / (g.Uinf * g.Uinf);
        double cpe = cp_exact(g.R, theta, g.R, g.Uinf);
        f << theta << "," << cpn << "," << cpe << "," << un << "," << vn << "," << ue << "," << ve << "\n";
    }
}

// -----------------------------------------------------------------------
// Erros (norma L2 e Linf) entre solucao numerica e analitica
// -----------------------------------------------------------------------
struct Errors { double L2_phi, Linf_phi, L2_cp_surf, Linf_cp_surf; };

Errors compute_errors(const Result &res, const std::vector<double> &U, const std::vector<double> &V) {
    const Grid &g = res.g;
    double sum2phi = 0.0, maxphi = 0.0;
    int n = g.Nr * g.Nt;
    for (int i = 0; i < g.Nr; ++i)
        for (int j = 0; j < g.Nt; ++j) {
            double e = res.phi[IDX(i, j, g.Nt)] - phi_exact(g.r[i], g.theta[j], g.R, g.Uinf);
            sum2phi += e * e;
            maxphi = std::max(maxphi, std::fabs(e));
        }
    double L2_phi = std::sqrt(sum2phi / n);

    double sum2cp = 0.0, maxcp = 0.0;
    int i = 0;
    for (int j = 0; j < g.Nt; ++j) {
        double un = U[IDX(i, j, g.Nt)], vn = V[IDX(i, j, g.Nt)];
        double cpn = 1.0 - (un * un + vn * vn) / (g.Uinf * g.Uinf);
        double cpe = cp_exact(g.R, g.theta[j], g.R, g.Uinf);
        double e = cpn - cpe;
        sum2cp += e * e;
        maxcp = std::max(maxcp, std::fabs(e));
    }
    double L2_cp = std::sqrt(sum2cp / g.Nt);

    return {L2_phi, maxphi, L2_cp, maxcp};
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main() {
    const double R = 1.0;      // raio do cilindro
    const double Rmax = 20.0;  // raio externo do dominio (~20R, "infinito" numerico)
    const double Uinf = 1.0;   // velocidade do escoamento nao perturbado

    // ---- 1) Malha de referencia: campo completo (P,u,v) + comparacao ----
    int Nr_ref = 60, Nt_ref = 120;
    std::cout << "Resolvendo malha de referencia " << Nr_ref << " x " << Nt_ref << "...\n";
    Result ref = solve(Nr_ref, Nt_ref, R, Rmax, Uinf);
    std::vector<double> Uref, Vref;
    compute_velocities(ref, Uref, Vref);
    write_field_csv("field_reference.csv", ref, Uref, Vref);
    write_surface_csv("surface_reference.csv", ref, Uref, Vref);
    Errors eref = compute_errors(ref, Uref, Vref);
    std::cout << "  iteracoes SOR: " << ref.iters << ", residuo final: " << ref.residual << "\n";
    std::cout << "  L2(phi)=" << eref.L2_phi << "  Linf(phi)=" << eref.Linf_phi << "\n";
    std::cout << "  L2(Cp,superficie)=" << eref.L2_cp_surf << "  Linf(Cp,superficie)=" << eref.Linf_cp_surf << "\n\n";

    // ---- 2) Estudo de convergencia de malha ----
    std::vector<std::pair<int,int>> meshes = {
        {10, 20}, {15, 30}, {20, 40}, {30, 60}, {40, 80}, {60, 120}, {80, 160}, {120, 240}
    };

    std::ofstream conv("convergence.csv");
    conv << "Nr,Nt,dr,dtheta,h,iters,L2_phi,Linf_phi,L2_cp_surf,Linf_cp_surf\n";
    std::cout << "Estudo de convergencia de malha:\n";
    std::cout << std::left << std::setw(6) << "Nr" << std::setw(6) << "Nt"
              << std::setw(12) << "dr" << std::setw(14) << "L2(phi)"
              << std::setw(14) << "L2(Cp_surf)" << "iters\n";
    for (auto &m : meshes) {
        int Nr = m.first, Nt = m.second;
        Result r = solve(Nr, Nt, R, Rmax, Uinf);
        std::vector<double> U, V;
        compute_velocities(r, U, V);
        Errors e = compute_errors(r, U, V);
        double h = r.g.dr; // metrica de referencia da malha (espacamento radial)
        conv << Nr << "," << Nt << "," << r.g.dr << "," << r.g.dtheta << "," << h << ","
             << r.iters << "," << e.L2_phi << "," << e.Linf_phi << ","
             << e.L2_cp_surf << "," << e.Linf_cp_surf << "\n";
        std::cout << std::left << std::setw(6) << Nr << std::setw(6) << Nt
                  << std::setw(12) << r.g.dr << std::setw(14) << e.L2_phi
                  << std::setw(14) << e.L2_cp_surf << r.iters << "\n";
    }

    std::cout << "\nArquivos gerados: field_reference.csv, surface_reference.csv, convergence.csv\n";
    return 0;
}
