#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <omp.h>
#include <chrono>


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
                    Ddxi[i][j] = (F[i][j+1] - F[i][j-1])/2;
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

        void constructor (scalarField0 x0,scalarField0 y0){
            x=x0;
            y=y0;
            lines = x0.lines;
            columns = x0.columns;
            std::vector<double> temp; 
            for (int j = 0; j<columns; j++) {
                temp.push_back(0.0);
            } 
            for (int i = 0; i<lines; i++) {
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
            if (i==0) {
                return (F[i+2][j] - 4.0*F[i+1][j] + 3.0*F[i][j])/2.0;
            }
            if (i==lines-1) {
                return (- F[i-2][j] + 4.0*F[i-1][j] - 3.0*F[i][j])/2.0;
            }
            return (F[i-1][j] - F[i+1][j])/2.0;
        }
        double D2dxi2(int i,int j) {
            if (j == 0 || j == columns-1) {
                return (F[i][1] -2.0*F[i][j] + F[i][columns-2]);
            }
            return (F[i][j+1] -2.0*F[i][j] + F[i][j-1]);
        }
        double D2deta2(int i,int j) {
            if (i==0) {
                return (F[i+3][j] - 4.0*F[i+2][j] + 5.0*F[i+1][j] - 2.0*F[i][j]);
            }
            if (i==lines-1) {
                return (-F[i-3][j] + 4.0*F[i-2][j] - 5.0*F[i-1][j] + 2.0*F[i][j]);
            }
            return (F[i-1][j] -2.0*F[i][j] + F[i+1][j]);
        }
        double D2dxideta(int i,int j) {
            if (j == 0 || j == columns-1) {
                return (F[i-1][1] + F[i+1][columns-2] - F[i-1][columns-2] - F[i+1][1])/4.0;
            }
            return (F[i-1][j+1] + F[i+1][j-1] - F[i-1][j-1] - F[i+1][j+1])/4.0;
        }
};

class flowSim{
    public:
    scalarField rho, u, v, p, e, T;
    scalarField mu, k;
    scalarField TAUxx, TAUxy, TAUyy, qx, qy;

    scalarField U1, U2, U3, U4;
    scalarField E1, E2, E3, E4;
    scalarField F1, F2, F3, F4;

    double T0, mu0;
    double Pr;
    double Cp, Cv, R;
    
    void constructor(double T1, double mu1, double Pr0, double Cp0, double Cv0, scalarField0 x, scalarField0 y) {
        T0 = T1;
        mu0 = mu1;
        Pr = Pr0; 
        Cp = Cp0;
        Cv = Cv0;
        R = Cp - Cv;
        rho.constructor(x,y);
        u.constructor(x,y);
        v.constructor(x,y);
        p.constructor(x,y);
        e.constructor(x,y);
        T.constructor(x,y);
        TAUxx.constructor(x,y);
        TAUxy.constructor(x,y);
        TAUyy.constructor(x,y);
        qx.constructor(x,y);
        qy.constructor(x,y);
        mu.constructor(x,y);
        k.constructor(x,y);
        U1.constructor(x,y);
        U2.constructor(x,y);
        U3.constructor(x,y);
        U4.constructor(x,y);
        E1.constructor(x,y);
        E2.constructor(x,y);
        E3.constructor(x,y);
        E4.constructor(x,y);
        F1.constructor(x,y);
        F2.constructor(x,y);
        F3.constructor(x,y);
        F4.constructor(x,y);
    }
    void calculateUEF() {
        dynamicViscosity();
        calculateTAUxx();
        calculateTAUxy();
        calculateTAUyy();
        thermalConductivity();
        calculateQx();
        calculateQy();
        #pragma omp parallel for
        for (int i = 0; i < u.x.lines; i++) {
            for (int j = 0; j < u.x.columns; j++) {
                U1.F[i][j] = rho.F[i][j];
                U2.F[i][j] = rho.F[i][j]*u.F[i][j];
                U3.F[i][j] = rho.F[i][j]*v.F[i][j];
                U4.F[i][j] = rho.F[i][j]*(e.F[i][j] + (std::pow(u.F[i][j],2) + std::pow(v.F[i][j],2))/2.0);
                E1.F[i][j] = rho.F[i][j]*u.F[i][j];
                E2.F[i][j] = rho.F[i][j]*u.F[i][j]*u.F[i][j] + p.F[i][j] - TAUxx.F[i][j];
                E3.F[i][j] = rho.F[i][j]*u.F[i][j]*v.F[i][j] - TAUxy.F[i][j];
                E4.F[i][j] = (U4.F[i][j] + p.F[i][j])*u.F[i][j] + qx.F[i][j] - u.F[i][j]*TAUxx.F[i][j] - v.F[i][j]*TAUxy.F[i][j];
                F1.F[i][j] = rho.F[i][j]*v.F[i][j];
                F2.F[i][j] = rho.F[i][j]*u.F[i][j]*v.F[i][j] - TAUxy.F[i][j];
                F3.F[i][j] = rho.F[i][j]*v.F[i][j]*v.F[i][j] + p.F[i][j] - TAUyy.F[i][j];
                F4.F[i][j] = (U4.F[i][j] + p.F[i][j])*v.F[i][j] + qy.F[i][j] - u.F[i][j]*TAUxy.F[i][j] - v.F[i][j]*TAUyy.F[i][j];
            }
        }
    }
    
    void decodeU() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                rho.F[i][j] = U1.F[i][j];
                u.F[i][j] = U2.F[i][j]/U1.F[i][j];
                v.F[i][j] = U3.F[i][j]/U1.F[i][j];
                e.F[i][j] = U4.F[i][j]/U1.F[i][j] - (std::pow(u.F[i][j],2.0) + std::pow(v.F[i][j],2.0))/2.0;
                T.F[i][j] = e.F[i][j]/Cv;
                p.F[i][j] = rho.F[i][j]*R*T.F[i][j];
            }
        }
    }
    void calculateTAUxx() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                TAUxx.F[i][j] = 2.0*mu.F[i][j]*(-(u.Ddx(i,j) + v.Ddy(i,j))/3.0 + u.Ddx(i,j)); }}
    }
    void calculateTAUxy() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                TAUxy.F[i][j] = mu.F[i][j]*(u.Ddy(i,j) + v.Ddx(i,j)); }}
    }
    void calculateTAUyy() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                TAUyy.F[i][j] = 2.0*mu.F[i][j]*(-(u.Ddx(i,j) + v.Ddy(i,j))/3.0 + v.Ddy(i,j)); }}
    }
    void calculateQx() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                qx.F[i][j] = - k.F[i][j]*T.Ddx(i,j); }}
    }
    void calculateQy() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                qy.F[i][j] = - k.F[i][j]*T.Ddy(i,j); }}
    }
    void dynamicViscosity() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                mu.F[i][j] = mu0*std::pow(T.F[i][j]/T0 , 1.5)*(T0 + 110.0)/(T.F[i][j] + 110.0); }}
    }
    void thermalConductivity() {
        #pragma omp parallel for
        for (int i = 0; i < u.lines; i++) {
            for (int j = 0; j < u.columns; j++) {
                k.F[i][j] = mu.F[i][j]*Cp/Pr; }}
    }
};

class boundaryConditionsClass{
    public:
        double rho_infinity;
        double u_infinity;
        double v_infinity;
        double p_infinity;
        double T_wall;
        boundaryConditionsClass(double rho, double u, double v, double p, double T) {
            rho_infinity = rho;
            u_infinity = u;
            v_infinity = v;
            p_infinity = p;
            T_wall = T;
        }
        void apply(flowSim &flow) {
            for (int j = 0; j < flow.u.columns; j++) {
                int i = 0;
                flow.rho.F[i][j] = rho_infinity;
                flow.u.F[i][j] = u_infinity;
                flow.v.F[i][j] = v_infinity;
                flow.p.F[i][j] = p_infinity;
                flow.e.F[i][j] = flow.Cv*p_infinity/((flow.R)*rho_infinity);
                flow.T.F[i][j] = p_infinity/(flow.R*rho_infinity);
            }
            for (int j = 0; j < flow.u.columns; j++) {
                int i = flow.u.lines-1;
                flow.u.F[i][j] = 0.0;
                flow.v.F[i][j] = 0.0;
                flow.T.F[i][j] = T_wall;
                flow.e.F[i][j] = flow.Cv*T_wall;
            }
            #pragma omp parallel for
            for (int i = 0; i < flow.u.lines; i++) {
                for (int j = 0; j < flow.u.columns; j++) {
                    if (flow.T.F[i][j] < 0) {
                        flow.T.F[i][j] = 0.15;
                    }
                    if (flow.e.F[i][j] < 0) {
                        flow.e.F[i][j] = 0.15;
                    }
                    if (flow.p.F[i][j] < 0) {
                        flow.p.F[i][j] = 0.15;
                    }
                    if (flow.rho.F[i][j] < 0) {
                        flow.rho.F[i][j] = 0.15;
                    }
                }
            }
        }
};

void readGridFile(std::string fileName, scalarField0 &x, scalarField0 &y);
void writeVTKfile(flowSim &flow, std::string fileName);
double stepSize(flowSim &flow);
void macCormack(flowSim &flow, flowSim &predictedFlow, flowSim &correctedFlow, boundaryConditionsClass &boundaryConditions);
bool checkConvergence(flowSim &flow0, flowSim &flow, double delta);
bool deuNan(flowSim &flow);
void ondeDeuNan(flowSim &flow);


int main() {
    std::string gridFileName = "cilindro.grid";
    std::string outputFileName = "escoamento.vtk";

    double T0 = 288.16;
    double mu0 = 1.7894e-5;
    double Pr = 0.71; 
    double Cp = 1004.0;
    double Cv = 717.0;

    double rho_infinity = 1.225;
    double u_infinity = 1.0e1;
    double v_infinity = 0.0;
    double p_infinity = 1.01325e5;
    double T_wall = p_infinity/((Cp-Cv)*rho_infinity);

    double precision = 0.00000001;



    scalarField0 x,y;
    readGridFile(gridFileName, x,y);
    
    flowSim initialFlow, predictedFlow, correctedFlow;
    initialFlow.constructor(T0,mu0,Pr,Cp,Cv,x,y);
    predictedFlow.constructor(T0,mu0,Pr,Cp,Cv,x,y);
    correctedFlow.constructor(T0,mu0,Pr,Cp,Cv,x,y);

    boundaryConditionsClass boundCond(rho_infinity, u_infinity, v_infinity, p_infinity, T_wall);
    boundCond.apply(correctedFlow);
    for (int i = 1; i < initialFlow.u.lines-1; i++) {
        for (int j = 0; j < initialFlow.u.columns; j++) {
            correctedFlow.rho.F[i][j] = rho_infinity;
            correctedFlow.u.F[i][j] = u_infinity;
            correctedFlow.v.F[i][j] = v_infinity;
            correctedFlow.p.F[i][j] = p_infinity;
            correctedFlow.e.F[i][j] = correctedFlow.Cv*p_infinity/((correctedFlow.R)*rho_infinity);
            correctedFlow.T.F[i][j] = p_infinity/(correctedFlow.R*rho_infinity);
        }
    }
    for (int j = 0; j < initialFlow.u.columns; j++) {
        int i = initialFlow.u.lines-1;
        correctedFlow.rho.F[i][j] = rho_infinity;
        correctedFlow.p.F[i][j] = p_infinity;
        correctedFlow.e.F[i][j] = correctedFlow.Cv*p_infinity/((correctedFlow.R)*rho_infinity);
        correctedFlow.T.F[i][j] = p_infinity/(correctedFlow.R*rho_infinity);
    }
    
    int MAXiterations = 1000000;
    for (int iteration = 1; iteration <= MAXiterations; iteration++) {
        auto start = std::chrono::high_resolution_clock::now();
        
        initialFlow = correctedFlow;
        macCormack(initialFlow, predictedFlow, correctedFlow, boundCond);
        
        if (deuNan(predictedFlow)) {
            std::cout << "Deu nan predicted flow\n";
            ondeDeuNan(predictedFlow);
            writeVTKfile(initialFlow, "antesDedarRuim.vtk");
            writeVTKfile(predictedFlow, outputFileName);
            break;
        }

        if (deuNan(correctedFlow)) {
            std::cout << "Deu nan corrected flow\n";
            ondeDeuNan(correctedFlow);
            writeVTKfile(correctedFlow, outputFileName);
            break;
        }
        //writeVTKfile(correctedFlow, outputFileName);
        if (iteration%1000 == 1) {
            writeVTKfile(predictedFlow, outputFileName);
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << duration.count() << " milisegundos na iteracao "<<iteration<< '\n';
    }
    
    writeVTKfile(predictedFlow, outputFileName);
    return 0;
}

void macCormack(flowSim &flow, flowSim &predictedFlow, flowSim &correctedFlow, boundaryConditionsClass &boundaryConditions) {
    double dt = stepSize(flow);

    flow.calculateUEF();

    #pragma omp parallel for
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            predictedFlow.U1.F[i][j] = flow.U1.F[i][j] - (flow.E1.Ddx(i,j) + flow.F1.Ddy(i,j))*dt;
            predictedFlow.U2.F[i][j] = flow.U2.F[i][j] - (flow.E2.Ddx(i,j) + flow.F2.Ddy(i,j))*dt;
            predictedFlow.U3.F[i][j] = flow.U3.F[i][j] - (flow.E3.Ddx(i,j) + flow.F3.Ddy(i,j))*dt;
            predictedFlow.U4.F[i][j] = flow.U4.F[i][j] - (flow.E4.Ddx(i,j) + flow.F4.Ddy(i,j))*dt;
        }
    }

    predictedFlow.decodeU();
    boundaryConditions.apply(predictedFlow);
    predictedFlow.calculateUEF();

    #pragma omp parallel for
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            correctedFlow.U1.F[i][j] = flow.U1.F[i][j] -0.5*((flow.E1.Ddx(i,j) + flow.F1.Ddy(i,j)) + (predictedFlow.E1.Ddx(i,j) + predictedFlow.F1.Ddy(i,j)))*dt;
            correctedFlow.U2.F[i][j] = flow.U2.F[i][j] -0.5*((flow.E2.Ddx(i,j) + flow.F2.Ddy(i,j)) + (predictedFlow.E2.Ddx(i,j) + predictedFlow.F2.Ddy(i,j)))*dt;
            correctedFlow.U3.F[i][j] = flow.U3.F[i][j] -0.5*((flow.E3.Ddx(i,j) + flow.F3.Ddy(i,j)) + (predictedFlow.E3.Ddx(i,j) + predictedFlow.F3.Ddy(i,j)))*dt;
            correctedFlow.U4.F[i][j] = flow.U4.F[i][j] -0.5*((flow.E4.Ddx(i,j) + flow.F4.Ddy(i,j)) + (predictedFlow.E4.Ddx(i,j) + predictedFlow.F4.Ddy(i,j)))*dt;
        }
    }
    correctedFlow.decodeU();
    boundaryConditions.apply(correctedFlow);
}

double stepSize(flowSim &flow) {
    /*
    std::vector<std::vector<double>> dt(flow.u.lines, std::vector<double>(flow.u.columns, 0.0));
    for (int i = 0; i < flow.u.lines; i++) {
        for (int i = 0; i < flow.u.lines; i++) {

        }
    }*/
   return 1.0e-6;
}

bool checkConvergence(flowSim &flow0, flowSim &flow, double delta) {
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (abs(flow.u.F[i][j] - flow0.u.F[i][j]) > delta) {
                return false;
            }
        }
    }
    return true;
}

void readGridFile(std::string fileName, scalarField0 &x, scalarField0 &y) {
    std::ifstream malha;
    std::string linha;
    int lin;
    int col;

    // ler o valor de "I" e armazenar em "columns" e "J" em "lines" do arquivo
    malha.open(fileName, std::ios::in);
    if (malha.is_open()) {
        std::getline(malha, linha); 

        if (std::getline(malha, linha)) {
            
            std::regex pattern(R"(I=\s*(\d+).*J=\s*(\d+))");
            std::smatch match;

            if (std::regex_search(linha, match, pattern)) {
                col = std::stoi(match[1]);
                lin = std::stoi(match[2]);
                //i=columns;
                //j=lines;
            }
        }
    }
    std::cout << lin << '\n' << col << '\n';
    x.constructor(lin,col);
    y.constructor(lin,col);
    
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
            if (index2==x.columns) {
                index2=0;
                index1++;
                break;
            }
        }
        
        if (index1 == x.lines) {
            break;
        }   
    }
    x.calculateMetrics();
    y.calculateMetrics();
    std::cout << "malha lida" << '\n';
}

void writeVTKfile(flowSim &flow, std::string fileName) {
    std::ofstream tecplot;
    tecplot.open(fileName,std::fstream::out|std::fstream::trunc);
        if(tecplot.is_open()){
            tecplot <<"# vtk DataFile Version 3.0"<<std::endl<<"Campo 2D"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
            tecplot<<"\n";
            tecplot<<"DIMENSIONS "<<flow.u.columns<<" "<< flow.u.lines<<" 1"<<std::endl;
            tecplot<<"POINTS "<<flow.u.columns*flow.u.lines<<" float"<<std::endl;
            tecplot<<"\n";
            //Dados espaciais
            for (int i = 0; i < flow.u.lines; i++) {
                for (int j = 0; j < flow.u.columns; j++) {
                    tecplot << flow.u.x.F[i][j]<< " " << flow.u.y.F[i][j]<< " "<<"0.0"<<"\n";
                }
            }
            //Dados escalares
            tecplot << "\nPOINT_DATA " << flow.u.columns * flow.u.lines << "\n";
            tecplot << "SCALARS Pressure float 1\n";
            tecplot << "LOOKUP_TABLE default\n";
            for (int i = 0; i < flow.u.lines; i++){
                for (int j = 0; j < flow.u.columns; j++){
                    tecplot << flow.p.F[i][j] << "\n";
                }
            }
            tecplot << "SCALARS Density float 1\n";
            tecplot << "LOOKUP_TABLE default\n";
            for (int i = 0; i < flow.u.lines; i++){
                for (int j = 0; j < flow.u.columns; j++){
                    tecplot << flow.rho.F[i][j] << "\n";
                }
            }
            tecplot << "SCALARS Temperature float 1\n";
            tecplot << "LOOKUP_TABLE default\n";
            for (int i = 0; i < flow.u.lines; i++){
                for (int j = 0; j < flow.u.columns; j++){
                    tecplot << flow.T.F[i][j] << "\n";
                }
            }
            tecplot << "SCALARS Energy float 1\n";
            tecplot << "LOOKUP_TABLE default\n";
            for (int i = 0; i < flow.u.lines; i++){
                for (int j = 0; j < flow.u.columns; j++){
                    tecplot << flow.e.F[i][j] << "\n";
                }
            }
            // Velocidade (u, v)
            tecplot << "\nVECTORS Velocity float\n";
            for (int i = 0; i < flow.u.lines; i++){
                for (int j = 0; j < flow.u.columns; j++){
                    tecplot << flow.u.F[i][j] << " " << flow.v.F[i][j] << " 0.0\n";
                }
            }

        }
        tecplot.close();
}

bool deuNan(flowSim &flow) {
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.rho.F[i][j])) { std::cout << "Nan: rho i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.u.F[i][j])) { std::cout << "Nan: u i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.v.F[i][j])) { std::cout << "Nan: v i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.p.F[i][j])) { std::cout << "Nan: p i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.e.F[i][j])) { std::cout << "Nan: e i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.T.F[i][j])) { std::cout << "Nan: T i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.mu.F[i][j])) { std::cout << "Nan: mu i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.k.F[i][j])) { std::cout << "Nan: k i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.TAUxx.F[i][j])) { std::cout << "Nan: TAUxx i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.TAUxy.F[i][j])) { std::cout << "Nan: TAUxy i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.TAUyy.F[i][j])) { std::cout << "Nan: TAUyy i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.qx.F[i][j])) { std::cout << "Nan: qx i,j=" << i << "," << j <<'\n'; return true;}
            else if (std::isnan(flow.qy.F[i][j])) { std::cout << "Nan: qy i,j=" << i << "," << j <<'\n'; return true;}
        }
    }
    return false;
}

void ondeDeuNan(flowSim &flow) {
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.rho.F[i][j])) {std::cout << "Nan: rho i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.u.F[i][j])) {std::cout << "Nan: u i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.v.F[i][j])) {std::cout << "Nan: v i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.p.F[i][j])) {std::cout << "Nan: p i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.e.F[i][j])) {std::cout << "Nan: e i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.T.F[i][j])) {std::cout << "Nan: T i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.mu.F[i][j])) {std::cout << "Nan: mu i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.k.F[i][j])) {std::cout << "Nan: k i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.TAUxx.F[i][j])) {std::cout << "Nan: TAUxx i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.TAUxy.F[i][j])) {std::cout << "Nan: TAUxy i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.TAUyy.F[i][j])) {std::cout << "Nan: TAUyy i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.qx.F[i][j])) {std::cout << "Nan: qx i,j=" << i << "," << j <<'\n';}}}
    for (int i = 0; i < flow.u.lines; i++) {
        for (int j = 0; j < flow.u.columns; j++) {
            if (std::isnan(flow.qy.F[i][j])) {std::cout << "Nan: qy i,j=" << i << "," << j <<'\n';}}}
}
