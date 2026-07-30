//CFD solver baseado no método de correção de pressão

#define _USE_MATH_DEFINES //utilizado para obter o valor de pi
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <omp.h>
#include <chrono>
#include "Eigen/Sparse"


class scalarField0{
    //estrutura que armazena os dados referentes aos campos escalares da malha e funções que representão suas derivardas
    public:
        int lines;
        int columns;
        std::vector<std::vector<double>> F;
        
        std::vector<std::vector<double>> Ddxi;
        std::vector<std::vector<double>> Ddeta;
        std::vector<std::vector<double>> D2dxi2;
        std::vector<std::vector<double>> D2deta2;
        std::vector<std::vector<double>> D2dxideta;
        
        void constructor (int lin, int col) {
            lines = lin;
            columns = col;
            std::vector<double> temp; 
            for (int j = 0; j<col; j++) {
                temp.push_back(0.0);
            } 
            for (int i = 0; i<lin; i++) {
                F.push_back(temp);
                Ddxi.push_back(temp);
                Ddeta.push_back(temp);
                D2dxi2.push_back(temp);
                D2deta2.push_back(temp);
                D2dxideta.push_back(temp);
            }
             
        }
        void calculateMetrics() {
            for (int i = 0; i <lines; i++) {
                for (int j = 1; j<columns-1; j++) {
                    Ddxi[i][j] = (F[i][j+1] - F[i][j-1])/2.0;
                    D2dxi2[i][j] = (F[i][j+1] -2.0*F[i][j] + F[i][j-1]);
                }
            }
            for (int i = 1; i <lines-1; i++) {
                for (int j = 0; j<columns; j++) {
                    Ddeta[i][j] = (F[i-1][j] - F[i+1][j])/2.0;
                    D2deta2[i][j] = (F[i-1][j] -2.0*F[i][j] + F[i+1][j]);
                }
            }
            for (int i = 1; i <lines-1; i++) {
                for (int j = 1; j<columns-1; j++) {
                    D2dxideta[i][j] = (F[i-1][j+1] + F[i+1][j-1] - F[i-1][j-1] - F[i+1][j+1])/4.0;
                }
            }

            //Norte
            for (int j = 0; j<columns; j++) {
                int i = 0;
                Ddeta[i][j] = (F[i+2][j] - 4.0*F[i+1][j] + 3.0*F[i][j])/2.0;
                D2deta2[i][j] = (F[i+3][j] - 4.0*F[i+2][j] + 5.0*F[i+1][j] - 2.0*F[i][j]);
            }
            for (int j = 1; j<columns-1;j++) {
                int i = 0;
                D2dxideta[i][j] = (Ddeta[i][j+1] - Ddeta[i][j-1])/2.0;
            }

            //Sul 
            for (int j = 0; j<columns;j++) {
                int i = lines-1;
                Ddeta[i][j] = (-F[i-2][j] + 4.0*F[i-1][j] - 3.0*F[i][j])/2.0;
                D2deta2[i][j] = (-F[i-3][j] + 4.0*F[i-2][j] - 5.0*F[i-1][j] + 2.0*F[i][j]);
            }
            for (int j = 1; j<columns-1;j++) {
                int i = lines-1;
                D2dxideta[i][j] = (Ddeta[i][j+1] - Ddeta[i][j-1])/2.0;
            }

            //Leste 
            for (int i = 0; i <lines; i++) {
                int j = columns-1;
                Ddxi[i][j] = (F[i][1] - F[i][j-1])/2.0;
                D2dxi2[i][j] = (F[i][1] -2.0*F[i][j] + F[i][j-1]);
            }    
            for (int i = 1; i <lines-1; i++) {
                int j = columns-1;
                D2dxideta[i][j] = (F[i-1][1] + F[i+1][j-1] - F[i-1][j-1] - F[i+1][1])/4.0;
            }

            //Oeste
            for (int i = 0; i <lines; i++) {
                int j = 0;
                Ddxi[i][j] = (F[i][j+1] - F[i][columns-2])/2.0;
                D2dxi2[i][j] = (F[i][j+1] -2.0*F[i][j] + F[i][columns-2]);
            }
            for (int i = 1; i <lines-1; i++) {
                int j = 0;
                D2dxideta[i][j] = (F[i-1][j+1] + F[i+1][columns-2] - F[i-1][columns-2] - F[i+1][j+1])/4.0;
            }
            //Noroeste
            D2dxideta[0][0] = (Ddeta[0][1] - Ddeta[0][columns-2])/2.0;
            //Nordeste
            D2dxideta[0][columns-1] = (Ddeta[0][1] - Ddeta[0][columns-2])/2.0;
            //Sudeste
            D2dxideta[lines-1][columns-1] = (Ddeta[lines-1][1] - Ddeta[lines-1][columns-2])/2.0;
            //Sudoeste
            D2dxideta[lines-1][0] = (Ddeta[lines-1][1] - Ddeta[lines-1][columns-2])/2.0;
        }
};
class scalarField{
    //estrutura que armazena os dados referentes aos campos escalares do fluído e funções que representão suas derivardas
    public:
        int lines;
        int columns;
        
        scalarField0 x;    
        scalarField0 y;

        std::vector<std::vector<double>> F;

        std::vector<std::vector<double>> J;
        std::vector<std::vector<double>> DJdxi;
        std::vector<std::vector<double>> DJdeta;

        scalarField (int lin, int col, scalarField0 x0,scalarField0 y0){
            x=x0;
            y=y0;
            lines = lin;
            columns = col;
            std::vector<double> temp; 
            for (int j = 0; j<col; j++) {
                temp.push_back(0.0);
            } 
            for (int i = 0; i<lin; i++) {
                F.push_back(temp);
                J.push_back(temp);
                DJdxi.push_back(temp);
                DJdeta.push_back(temp);
            }
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < columns; j++) {
                    J[i][j] = x.Ddxi[i][j]*y.Ddeta[i][j] - y.Ddxi[i][j]*x.Ddeta[i][j];
                    DJdxi[i][j] = (x.D2dxi2[i][j]*y.Ddeta[i][j] + x.Ddxi[i][j]*y.D2dxideta[i][j]
                                - y.D2dxi2[i][j]*x.Ddeta[i][j] - y.Ddxi[i][j]*x.D2dxideta[i][j]);
                    DJdeta[i][j] = (x.D2dxideta[i][j]*y.Ddeta[i][j] + x.Ddxi[i][j]*y.D2deta2[i][j]
                                - y.D2dxideta[i][j]*x.Ddeta[i][j] - y.Ddxi[i][j]*x.D2deta2[i][j]);
                }
            }
        }
        double Ddx(int i,int j) {
            //Derivada total em relação à x
            return (Ddxi(i,j)*y.Ddeta[i][j] - Ddeta(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy(int i,int j) {
            //Derivada total em relação à y
            return (Ddeta(i,j)*x.Ddxi[i][j] - Ddxi(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2(int i, int j) {
            //Derivada total de segunda ordem em relação à x
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            
            return (Ddxi(i,j)*A + Ddeta(i,j)*B + D2dxideta(i,j)*C + D2dxi2(i,j)*D + D2deta2(i,j)*E)/std::pow(J[i][j],3);   
        }
        double D2dy2(int i, int j) {
            //Derivada total de segunda ordem em relação à y
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi(i,j)*A + Ddeta(i,j)*B + D2dxideta(i,j)*C + D2dxi2(i,j)*D + D2deta2(i,j)*E)/std::pow(J[i][j],3);
        }
        
        double Ddxi(int i,int j) {
            if (j == 0 || j == columns-1) {
                return (F[i][1] - F[i][columns-2])/2.0;
            }
            return (F[i][j+1] - F[i][j-1])/2.0;
        }
        double Ddeta(int i,int j) {
            if (i == lines-1) {
                return (-F[i-2][j] + 4.0*F[i-1][j] - 3.0*F[i][j])/2.0;
            }
            if (i == 0) {
                return (+F[i+2][j] - 4.0*F[i+1][j] + 3.0*F[i][j])/2.0;
            }
            return (F[i-1][j] - F[i+1][j])/2.0;
        }
        double D2dxi2(int i,int j) {
            if (j == 0 || j == columns-1) {
                return (F[i][1] - 2.0*F[i][j] + F[i][columns-2]);
            }
            return (F[i][j+1] -2.0*F[i][j] + F[i][j-1]);
        }
        double D2deta2(int i,int j) {
            if (i == lines-1) {
                return (-F[i-3][j] + 4.0*F[i-2][j] - 5.0*F[i-1][j] + 2.0*F[i][j]);
            }
            if (i == 0) {
                return (F[i+3][j] - 4.0*F[i+2][j] + 5.0*F[i+1][j] - 2.0*F[i][j]);
            }
            return (F[i-1][j] -2.0*F[i][j] + F[i+1][j]);
        }
        double D2dxideta(int i,int j) {
            return (Ddxi(i-1,j) - Ddxi(i+1,j))/2.0;
        }
};

bool checkConvergence(std::vector<std::vector<double>> M, std::vector<std::vector<double>> M0, int lin, int col);
double NACAsup(int m, int p, int t, double c, double x);
double NACAinf(int m, int p, int t, double c, double x);
double Dydx(double x, std::string pos, int m, int p, int t, double c);
int matrixID(int i, int j, int lin, int col);

int main() {
    std::cout << "esta rodando\n";
    std::cout << "OpenMP Threads available: " << omp_get_max_threads() << "\n";
    int lines;
    int columns;
    std::ifstream malha;
    std::ofstream tecplot;
    std::string linha;
    // ler o valor de "I" e armazenar em "columns" e "J" em "lines" do arquivo
    //malha.open("OgridNACA20011.grid", std::ios::in);//Arquivo 1 .grid
    malha.open("OgridNACA20011.grid", std::ios::in);
    if (malha.is_open()) {
        std::getline(malha, linha); 

        if (std::getline(malha, linha)) {
            
            std::regex pattern(R"(I=\s*(\d+).*J=\s*(\d+))");
            std::smatch match;

            if (std::regex_search(linha, match, pattern)) {
                columns = std::stoi(match[1]);
                
                lines = std::stoi(match[2]);
                
                //i=columns;
                //j=lines;
            }
        }
    }
    std::cout << lines << '\n';
    std::cout << columns << '\n';
    //check
    scalarField0 x;
    scalarField0 y;
    x.constructor(lines,columns);
    y.constructor(lines,columns);
    
    // ler os dados dos pontos do arquivo .grid e armazenar em x.F[i][j] & y.F[i][j]
    int index1 = 0;
    int index2 = 0;
    while (std::getline(malha, linha)) {
        std::stringstream ss(linha);
        double tempx, tempy;
    
        while (ss >> tempx >> tempy) {
            x.F[index1][index2] = tempx;
            y.F[index1][index2] = tempy;
            index2++;
            if (index2==columns) {
                index2=0;
                index1++;
                break;
            }
        }
        if (index1 == lines) {
            break;
        } 
    }
    std::cout << "malha lida" << '\n';
    
    x.calculateMetrics();
    y.calculateMetrics();

    
    
    double density = 1.2;
    
    double Vinfinity = 20;
    double Pinfinity = 1.0e5;
    double Alpha = 10.0*M_PI/180.0;
    
    int m = 0;
    int p = 0;
    int tt = 11;
    double c = 1.0; //tamanho da corda

    scalarField u(lines,columns,x,y);
    scalarField v(lines,columns,x,y);
    scalarField phi(lines,columns,x,y);

    for (int i = 0; i<lines; i++) {
        for(int j = 0; j<columns; j++) {
            if (std::isnan(u.J[i][j]) || std::isnan(u.DJdeta[i][j]) || std::isnan(u.DJdxi[i][j])) {
                std::cout << "deu nan i="<<i<<", j="<<j;
            }
        }
    }

    for (int j = 0; j < columns; j++) {
        u.F[0][j] = Vinfinity*std::cos(Alpha);
    }
    for (int j = 0; j < columns; j++) {
        v.F[0][j] = Vinfinity*std::sin(Alpha);
    }
    
    for (int j = 0; j < columns; j++) {
        int i = 0;
        phi.F[i][j] = u.F[i][j]*x.F[i][j] + v.F[i][j]*y.F[i][j];
    }
    
    int iterations=100000;
    double dt=0.000001;

    std::vector<double> tempVector;
    std::vector<std::vector<double>> A,B,C,D,E;
    std::vector<std::vector<double>> Dphi2dxi2, Dphi2deta2, Dphi2dxideta, Dphidxi, Dphideta;
    std::vector<std::vector<double>> phi_N, phi_NE, phi_E, phi_SE, phi_S, phi_SW, phi_W, phi_NW, phi_C, phi_NN;
    std::vector<double> Pconstants;
    for (int i = 0; i < (lines-1)*(columns-1)+1; i++) {
        Pconstants.push_back(0.0);
    }
    for (int j = 0; j < columns; j++) {
        tempVector.push_back(0.0);
    }
    for (int i = 0; i < lines+1; i++) {
        A.push_back(tempVector);
        B.push_back(tempVector);
        C.push_back(tempVector);
        D.push_back(tempVector);
        E.push_back(tempVector);
        Dphi2dxi2.push_back(tempVector);
        Dphi2deta2.push_back(tempVector);
        Dphi2dxideta.push_back(tempVector);
        Dphidxi.push_back(tempVector);
        Dphideta.push_back(tempVector);
    }
    for (int i = 0; i < lines-1; i++) {
        for (int j = 0; j< columns; j++) {
            A[i][j] = std::pow(x.Ddeta[i][j],2.0) + std::pow(y.Ddeta[i][j],2.0);
            B[i][j] = std::pow(x.Ddxi[i][j],2.0) + std::pow(y.Ddxi[i][j],2.0);
            C[i][j] = x.Ddxi[i][j]*x.Ddeta[i][j] + y.Ddxi[i][j]*y.Ddeta[i][j];
            D[i][j] = x.D2dxideta[i][j]*x.Ddeta[i][j] - x.D2deta2[i][j]*x.Ddxi[i][j] - y.D2deta2[i][j]*y.Ddxi[i][j] + y.D2dxideta[i][j]*y.Ddeta[i][j];
            E[i][j] = x.D2dxideta[i][j]*x.Ddxi[i][j] - x.D2dxi2[i][j]*x.Ddeta[i][j] - y.D2dxi2[i][j]*y.Ddeta[i][j] + y.D2dxideta[i][j]*y.Ddxi[i][j];
        }
    }
    for (int i = 0; i < lines-1; i++) {
        for (int j = 0; j< columns; j++) {
            Dphi2dxi2[i][j] = phi.J[i][j]*A[i][j]/std::pow(phi.J[i][j], 3.0);
            Dphi2deta2[i][j] = phi.J[i][j]*B[i][j]/std::pow(phi.J[i][j], 3.0);
            Dphi2dxideta[i][j] = - 2.0*phi.J[i][j]*C[i][j]/std::pow(phi.J[i][j], 3.0);
            Dphidxi[i][j] = (phi.J[i][j]*D[i][j] - phi.DJdxi[i][j]*A[i][j] + phi.DJdeta[i][j]*C[i][j])/std::pow(phi.J[i][j], 3.0);
            Dphideta[i][j] = (phi.J[i][j]*E[i][j] + phi.DJdxi[i][j]*C[i][j] - phi.DJdeta[i][j]*B[i][j])/std::pow(phi.J[i][j], 3.0);
        }
    }
    for (int j = 1; j < columns-1; j++) {
        int i = lines -1;
        if (j < (columns+1)/2) {
            Dphidxi[i][j] = (y.Ddeta[i][j]*Dydx(x.F[i][j], "NACAinf", m, p, tt, c) + x.Ddeta[i][j])/phi.J[i][j];;
            Dphideta[i][j] = - (y.Ddxi[i][j]*Dydx(x.F[i][j], "NACAinf", m, p, tt, c) + x.Ddxi[i][j])/phi.J[i][j];;
        }
        else {
            Dphidxi[i][j] = (y.Ddeta[i][j]*Dydx(x.F[i][j], "NACAsup", m, p, tt, c) + x.Ddeta[i][j])/phi.J[i][j];;
            Dphideta[i][j] = - (y.Ddxi[i][j]*Dydx(x.F[i][j], "NACAsup", m, p, tt, c) + x.Ddxi[i][j])/phi.J[i][j];;
        }
    }
    for (int j = 0; j< 1; j++) {
        int i = lines-1;
        Dphidxi[i][j] = y.Ddeta[i][j]/phi.J[i][j];;
        Dphideta[i][j] =  - y.Ddxi[i][j]/phi.J[i][j];
    }

    int SacrificeID1 = lines;
    int SacrificeID2 = 0;
    for (int j = SacrificeID2; j< SacrificeID2 + 1; j++) {
        int i = SacrificeID1;
        Dphidxi[i][j] = - x.Ddeta[lines-1][0]/phi.J[lines-1][0];
        Dphideta[i][j] = x.Ddxi[lines-1][0]/phi.J[lines-1][0];
    }

    std::vector<std::vector<std::vector<std::vector<double>>>> matrixCoefficients(
        lines+1, std::vector<std::vector<std::vector<double>>>(
               columns, std::vector<std::vector<double>>(
                      3, std::vector<double>(3, 0.0))));
    
    for (int i = 0; i < lines-1; i++) {
        for (int j = 0; j < columns; j++) {
            if (i == SacrificeID1 && j == SacrificeID2) {
                continue;
            }
            double C = -2.0*Dphi2dxi2[i][j] -2.0*Dphi2deta2[i][j];
            double E = Dphi2dxi2[i][j] + Dphidxi[i][j]/2.0;
            double W = Dphi2dxi2[i][j] - Dphidxi[i][j]/2.0;
            double N = Dphi2deta2[i][j] + Dphideta[i][j]/2.0;
            double S = Dphi2deta2[i][j] - Dphideta[i][j]/2.0;
            double NE = Dphi2dxideta[i][j]/4.0;
            double NW = - Dphi2dxideta[i][j]/4.0;
            double SE = - Dphi2dxideta[i][j]/4.0;
            double SW = Dphi2dxideta[i][j]/4.0;

            matrixCoefficients[i][j] = {{NW, N, NE},
                                        {W, C, E},
                                        {SW, S, SE}};
        }
    }
    for (int j = 0; j< columns-1; j++) {
        int i = lines-1;
        if (i == SacrificeID1 && j == SacrificeID2) {
            continue;
        }
        double C = (-3.0/2.0)*Dphideta[i][j];
        double E = (1.0/2.0)*Dphidxi[i][j];
        double W = (-1.0/2.0)*Dphidxi[i][j];
        double N = (4.0/2.0)*Dphideta[i][j];
        double NN = (-1.0/2.0)*Dphideta[i][j];

        matrixCoefficients[i][j] = {{0.0, NN, 0.0},
                                    {0.0, N, 0.0},
                                    {W, C, E}};
    }
    for (int j = SacrificeID2; j < SacrificeID2+1; j++) {
        int i = SacrificeID1;
        
        double C = (-3.0/2.0)*Dphideta[i][j];
        double E = (1.0/2.0)*Dphidxi[i][j];
        double W = (-1.0/2.0)*Dphidxi[i][j];
        double N = (4.0/2.0)*Dphideta[i][j];
        double NN = (-1.0/2.0)*Dphideta[i][j];

        matrixCoefficients[i][j] = {{ 0.0, NN, 0.0},
                                    { 0.0, N, 0.0},
                                    { W, C, E}};
                                    
    }
    
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Eigen::Triplet<double>> tripletList;

    for (int j = 0; j < columns-1; j++) {
        int i = 1;
        int id = (i-1)*(columns - 1)+j;
        double NW = matrixCoefficients[i][j][0][0];
        double N = matrixCoefficients[i][j][0][1];
        double NE = matrixCoefficients[i][j][0][2];
        if (j == 0) {
            Pconstants[id] = - (NW*phi.F[i-1][columns-2] + N*phi.F[i-1][j] + NE*phi.F[i-1][j+1]);
        }
        else {
            Pconstants[id] = - (NW*phi.F[i-1][j-1] + N*phi.F[i-1][j] + NE*phi.F[i-1][j+1]);
        }
    }
    for (int j = 0; j < columns-1; j++) {
        int i = 1;
        for (int k = 1; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                tripletList.push_back(Eigen::Triplet<double>(matrixID(i,j,lines,columns), matrixID(i + k-1, j + l-1, lines, columns), matrixCoefficients[i][j][k][l]));
            }
        }
    }
    for (int i = 2; i < lines-1; i++) {
        for (int j = 0; j < columns-1; j++) {
            if (i == SacrificeID1 && j == SacrificeID2) {
                continue;
            }
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    tripletList.push_back(Eigen::Triplet<double>(matrixID(i,j,lines,columns), matrixID(i + k-1, j + l-1, lines, columns), matrixCoefficients[i][j][k][l]));
                }
            }
        }
    }
    std::cout << "rodou ate aq 6\n";
    for (int j = 0; j < columns-1; j++) {
        int i = lines -1;
        if (i == SacrificeID1 && j == SacrificeID2) {
            continue;
        }
        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                tripletList.push_back(Eigen::Triplet<double>(matrixID(i,j,lines,columns), matrixID(i + k-2, j + l-1, lines, columns), matrixCoefficients[i][j][k][l]));
            }
        }
    }
    
    for (int j = SacrificeID2; j < SacrificeID2+1; j++) {
        int i = SacrificeID1;
        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                tripletList.push_back(Eigen::Triplet<double>(matrixID(i,j,lines,columns), matrixID(lines-1 + k-2, 0 + l-1, lines, columns), matrixCoefficients[i][j][k][l]));
            }
        }
    }
    
    Eigen::SparseMatrix<double> M((lines-1)*(columns-1)+1,(lines-1)*(columns-1));
    M.setFromTriplets(tripletList.begin(), tripletList.end());
    M.makeCompressed(); 
    
    Eigen::VectorXd b((lines-1)*(columns-1)+1);
        
    for (int i = 0; i <(lines-1)*(columns-1)+1; i++) {
        b(i) = Pconstants[i];
    }
        
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    
    Eigen::SparseMatrix<double> MtM = M.transpose() * M;
    Eigen::VectorXd Mtb = M.transpose() * b;

    solver.compute(MtM);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return 1;
    }

        // 7. Solve the system
    Eigen::VectorXd solution = solver.solve(Mtb);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
        return 1;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << duration.count() << " segundos na iteracao \n";

    for (int i = 1; i < lines; i++) {
        for (int j = 0; j < columns-1; j++){
            phi.F[i][j] = solution[(i-1)*(columns - 1)+j];
        } 
    }
    for (int i = 1; i < lines; i++) {
        int j = columns-1;
        phi.F[i][j] = phi.F[i][0];
    }
    std::cout << "valor de phi no LE e: "<< phi.F[lines-1][(columns+1)/2] << '\n';
    
    for (int i = 1; i < lines; i++) {
        for (int j = 1; j < columns; j++) {
            u.F[i][j] = phi.Ddx(i,j);
            v.F[i][j] = phi.Ddy(i,j);
        }
        u.F[i][0] = u.F[i][columns-1];
        v.F[i][0] = v.F[i][columns-1];

    }
    std::cout << "valor de u no TE e: "<< u.F[lines-1][0] << '\n';
    std::cout << "valor de v no TE e: "<< v.F[lines-1][0] << '\n';
    // output dos dados dos vetores u.F, v.F e P.F para o arquivo .tec 
    // x.F;y.F;u.F;v.F;P.F
    tecplot.open("malha.vtk",std::fstream::out|std::fstream::trunc);
    if(tecplot.is_open()){
        tecplot <<"# vtk DataFile Version 3.0"<<std::endl<<"Campo 2D"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
        tecplot<<"\n";
        tecplot<<"DIMENSIONS "<<columns<<" "<< lines<<" 1"<<std::endl;
        tecplot<<"POINTS "<<columns*lines<<" float"<<std::endl;
        tecplot<<"\n";
        //Dados espaciais
        for (int i = 0; i < lines; i++) {
            for (int j = 0; j < columns; j++) {
                tecplot << x.F[i][j]<< " " << y.F[i][j]<< " "<<"0.0"<<"\n";
            }
        }
        //Dados escalares
        
        tecplot << "\nPOINT_DATA " << columns * lines << "\n";
        tecplot << "SCALARS Pressure float 1\n";
        tecplot << "LOOKUP_TABLE default\n";
        for (int i = 0; i < lines; i++){
            for (int j = 0; j < columns; j++){
                tecplot << phi.F[i][j] << "\n";
            }
        }
        // Velocidade (u, v)
        tecplot << "\nVECTORS Velocity float\n";
        for (int i = 0; i < lines; i++){
            for (int j = 0; j < columns; j++){
                tecplot << u.F[i][j] << " " << v.F[i][j] << " 0.0\n";
            }
        }

    }

    tecplot.close();
    std::cout << "terminou" << '\n';
    
    return 0;
}
       


bool checkConvergence(std::vector<std::vector<double>> M, std::vector<std::vector<double>> M0, int lin, int col) {
    double delta = 1.0e-6;
    
    for (int i = 1; i < lin; i++) {
        for (int j = 1; j < col-1; j++) {
            if (std::abs(M[i][j] - M0[i][j]) >= delta) {
                return false;
            }
        }
    }

    return true;
}

double NACAsup(int m, int p, int t, double c, double x) {
    
    double thicknees = c*0.05*t*(0.2969*std::sqrt(x/c) - 0.1260*x/c - 0.3516*std::pow(x/c, 2) + 0.2843*std::pow(x/c, 3) - 0.1036*std::pow(x/c, 4));
    
    if (x <= 0) {
        return 0;
    }
    else if (x <= p*c*0.1) {
        return 0.01*m*(0.2*p*x - std::pow(x,2)/c)/std::pow(0.1*p, 2) + thicknees*std::cos(std::atan(0.01*m*(0.2*p - 2*x/c)/std::pow(0.1*p, 2)));
    }
    else if (x <= c) {
        return c*0.01*m*(1 - 0.2*p + (0.2*p*x/c - std::pow(x/c,2)))/std::pow(1-0.1*p,2) + thicknees*std::cos(std::atan(0.01*m*(0.2*p - 2*x/c)/std::pow(1-0.1*p, 2)));
    }
    return 0;
}

double NACAinf(int m, int p, int t, double c, double x){

    double thicknees = c*0.05*t*(0.2969*std::sqrt(x/c) - 0.1260*x/c - 0.3516*std::pow(x/c, 2) + 0.2843*std::pow(x/c, 3) - 0.1036*std::pow(x/c, 4));
    
    if (x <= 0) {
        return 0;
    }
    else if (x <= p*c*0.1) {
        return c*0.01*m*(0.2*p*x/c - std::pow(x/c,2))/std::pow(0.1*p, 2) - thicknees*std::cos(std::atan(c*0.01*m*(0.2*p/c - 2*x/(c*c))/std::pow(0.1*p, 2)));
    }
    else if (x <= c) {
        return c*0.01*m*(1 - 0.2*p + (0.2*p*x/c - std::pow(x/c,2)))/std::pow(1-0.1*p,2) - thicknees*std::cos(std::atan(c*0.01*m*(0.2*p/c - 2*x/(c*c))/std::pow(1-0.1*p, 2))) ;
    }
    return 0;
}

double Dydx(double x, std::string pos, int m, int p, int t, double c){
    std::vector<double> n;
    double derivative;
    double h = 1.0e-9;
    
    if (pos == "NACAsup") {
        derivative = (NACAsup(m,p,t,c,x+h) - NACAsup(m,p,t,c,x-h))/(2*h);
    }
    else if (pos == "NACAinf") {
        derivative = (NACAinf(m,p,t,c,x+h) - NACAinf(m,p,t,c,x-h))/(2*h);
    }
    return derivative;
}

int matrixID(int i, int j, int lin, int col) {
    if (j == col -1) {
        j = 0;
    }
    if (j == col) {
        j = 1;
    }
    if (j == -1) {
        j = col -2;
    }
    return (i-1)*(col - 1)+j;
}