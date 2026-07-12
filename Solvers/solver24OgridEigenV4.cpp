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
        //Norte
        double Ddx_N(int i,int j) {
            return (Ddxi(i,j)*y.Ddeta[i][j] - Ddeta_N(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy_N(int i,int j) {
            return (Ddeta_N(i,j)*x.Ddxi[i][j] - Ddxi(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2_N(int i, int j) {
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            return (Ddxi(i,j)*A + Ddeta_N(i,j)*B + D2dxideta_N(i,j)*C + D2dxi2(i,j)*D + D2deta2_N(i,j)*E)/std::pow(J[i][j],3); 
        }
        double D2dy2_N(int i, int j) {
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi(i,j)*A + Ddeta_N(i,j)*B + D2dxideta_N(i,j)*C + D2dxi2(i,j)*D + D2deta2_N(i,j)*E)/std::pow(J[i][j],3);
        }
        //Sul
        double Ddx_S(int i,int j) {
            return (Ddxi(i,j)*y.Ddeta[i][j] - Ddeta_S(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy_S(int i,int j) {
            return (Ddeta_S(i,j)*x.Ddxi[i][j] - Ddxi(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2_S(int i, int j) {
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            
            return (Ddxi(i,j)*A + Ddeta_S(i,j)*B + D2dxideta_S(i,j)*C + D2dxi2(i,j)*D + D2deta2_S(i,j)*E)/std::pow(J[i][j],3);    
        }
        double D2dy2_S(int i, int j) {
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi(i,j)*A + Ddeta_S(i,j)*B + D2dxideta_S(i,j)*C + D2dxi2(i,j)*D + D2deta2_S(i,j)*E)/std::pow(J[i][j],3);
        }
        //Leste
        double Ddx_E(int i,int j) {
            return (Ddxi_E(i,j)*y.Ddeta[i][j] - Ddeta(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy_E(int i,int j) {
            return (Ddeta(i,j)*x.Ddxi[i][j] - Ddxi_E(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2_E(int i, int j) {
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            
            return (Ddxi_E(i,j)*A + Ddeta(i,j)*B + D2dxideta_E(i,j)*C + D2dxi2_E(i,j)*D + D2deta2(i,j)*E)/std::pow(J[i][j],3);   
        }
        double D2dy2_E(int i, int j) {
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi_E(i,j)*A + Ddeta(i,j)*B + D2dxideta_E(i,j)*C + D2dxi2_E(i,j)*D + D2deta2(i,j)*E)/std::pow(J[i][j],3);
        }
        //Oeste
        double Ddx_W(int i,int j) {
            return (Ddxi_W(i,j)*y.Ddeta[i][j] - Ddeta(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy_W(int i,int j) {
            return (Ddeta(i,j)*x.Ddxi[i][j] - Ddxi_W(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2_W(int i, int j) {
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            
            return (Ddxi_W(i,j)*A + Ddeta(i,j)*B + D2dxideta_W(i,j)*C + D2dxi2_W(i,j)*D + D2deta2(i,j)*E)/std::pow(J[i][j],3);  
        }
        double D2dy2_W(int i, int j) {
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi_W(i,j)*A + Ddeta(i,j)*B + D2dxideta_W(i,j)*C + D2dxi2_W(i,j)*D + D2deta2(i,j)*E)/std::pow(J[i][j],3);
        }
        
        //Sudoeste
        double Ddx_SW(int i,int j) {
            return (Ddxi_W(i,j)*y.Ddeta[i][j] - Ddeta_S(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy_SW(int i,int j) {
            return (Ddeta_S(i,j)*x.Ddxi[i][j] - Ddxi_W(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2_SW(int i, int j) {
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            
            return (Ddxi_W(i,j)*A + Ddeta_S(i,j)*B + D2dxideta_SW(i,j)*C + D2dxi2_W(i,j)*D + D2deta2_S(i,j)*E)/std::pow(J[i][j],3);  
        }
        double D2dy2_SW(int i, int j) {
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi_W(i,j)*A + Ddeta_S(i,j)*B + D2dxideta_SW(i,j)*C + D2dxi2_W(i,j)*D + D2deta2_S(i,j)*E)/std::pow(J[i][j],3);
    
        }
        //Sudeste
        double Ddx_SE(int i,int j) {
            return (Ddxi_E(i,j)*y.Ddeta[i][j] - Ddeta_S(i,j)*y.Ddxi[i][j])/J[i][j];
        }
        double Ddy_SE(int i,int j) {
            return (Ddeta_S(i,j)*x.Ddxi[i][j] - Ddxi_E(i,j)*x.Ddeta[i][j])/J[i][j];
        }
        double D2dx2_SE(int i, int j) {
            double A,B,C,D,E;
            A = (J[i][j]*y.D2dxideta[i][j] - DJdxi[i][j]*y.Ddeta[i][j])*y.Ddeta[i][j] - (J[i][j]*y.D2deta2[i][j] - DJdeta[i][j]*y.Ddeta[i][j])*y.Ddxi[i][j];
            B = (DJdxi[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxi2[i][j])*y.Ddeta[i][j] - (DJdeta[i][j]*y.Ddxi[i][j] - J[i][j]*y.D2dxideta[i][j])*y.Ddxi[i][j];
            C = (-2.0*J[i][j]*y.Ddxi[i][j]*y.Ddeta[i][j]);
            D = J[i][j]*std::pow(y.Ddeta[i][j],2);
            E = J[i][j]*std::pow(y.Ddxi[i][j],2);
            
            return (Ddxi_E(i,j)*A + Ddeta_S(i,j)*B + D2dxideta_SE(i,j)*C + D2dxi2_E(i,j)*D + D2deta2_S(i,j)*E)/std::pow(J[i][j],3);   
        }
        double D2dy2_SE(int i, int j) {
            double A,B,C,D,E;
            A = (DJdeta[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2deta2[i][j])*x.Ddxi[i][j] - (DJdxi[i][j]*x.Ddeta[i][j] - J[i][j]*x.D2dxideta[i][j])*x.Ddeta[i][j];
            B = (J[i][j]*x.D2dxideta[i][j] - DJdeta[i][j]*x.Ddxi[i][j])*x.Ddxi[i][j] - (J[i][j]*x.D2dxi2[i][j] - DJdxi[i][j]*x.Ddxi[i][j])*x.Ddeta[i][j];
            C = (-2.0*J[i][j]*x.Ddxi[i][j]*x.Ddeta[i][j]);
            D = J[i][j]*std::pow(x.Ddeta[i][j],2);
            E = J[i][j]*std::pow(x.Ddxi[i][j],2);
            return (Ddxi_E(i,j)*A + Ddeta_S(i,j)*B + D2dxideta_SE(i,j)*C + D2dxi2_E(i,j)*D + D2deta2_S(i,j)*E)/std::pow(J[i][j],3);
        }
        double Ddxi(int i,int j) {
            /*if (j > 1 && j < columns -2) {
                return (-F[i][j+2] + 8.0*F[i][j+1] - 8.0*F[i][j-1] + F[i][j-2])/12.0;
            }*/
            
            return (F[i][j+1] - F[i][j-1])/2.0;
        }
        double Ddeta(int i,int j) {
            /*if (i > 1 && i < lines -2) {
                return (-F[i-2][j] + 8.0*F[i-1][j] - 8.0*F[i+1][j] + F[i+2][j])/12.0;
            }*/
            return (F[i-1][j] - F[i+1][j])/2.0;
        }
        double D2dxi2(int i,int j) {
            
            if (j > 1 && j < columns -2) {
                return (-F[i][j+2] + 16.0*F[i][j+1] - 30.0*F[i][j] + 16.0*F[i][j-1] - F[i][j-2])/12.0;
            }
            
            return (F[i][j+1] -2.0*F[i][j] + F[i][j-1]);
        }
        double D2deta2(int i,int j) {
            
            if (i > 1 && i < lines-2) {
                return (-F[i-2][j] + 16.0*F[i-1][j] - 30.0*F[i][j] + 16.0*F[i+1][j] - F[i+2][j])/12.0;
            }
                
            return (F[i-1][j] -2.0*F[i][j] + F[i+1][j]);
        }
        double D2dxideta(int i,int j) {
            return (F[i-1][j+1] + F[i+1][j-1] - F[i-1][j-1] - F[i+1][j+1])/4.0;
        }
        //Norte
        double Ddeta_N(int i,int j) {
            return (F[i+2][j] - 4.0*F[i+1][j] + 3.0*F[i][j])/2;
        }
        double D2deta2_N(int i,int j) {
            return (F[i+3][j] - 4.0*F[i+2][j] + 5.0*F[i+1][j] - 2.0*F[i][j]);
        }
        double D2dxideta_N(int i,int j) {
            return (Ddeta_N(i,j+1) - Ddeta_N(i,j-1))/2.0;
        }
        //Sul
        double Ddeta_S(int i,int j) {
            return (-F[i-2][j] + 4.0*F[i-1][j] - 3.0*F[i][j])/2.0;
        }
        double D2deta2_S(int i,int j) {
            return (-F[i-3][j] + 4.0*F[i-2][j] - 5.0*F[i-1][j] + 2.0*F[i][j]);
        }
        double D2dxideta_S(int i,int j) {
            return (Ddeta_S(i,j+1) - Ddeta_S(i,j-1))/2.0;
        }
        //Oeste
        double Ddxi_W(int i,int j) {
            return (F[i][j+1] - F[i][columns-2])/2.0;
        }
        double D2dxi2_W(int i,int j) {
            return (F[i][j+1] -2.0*F[i][j] + F[i][columns-2]);
        }
        double D2dxideta_W(int i,int j) {
            return (Ddxi_W(i-1,j) - Ddxi_W(i+1,j))/2.0;
        }
        //Leste
        double Ddxi_E(int i,int j) {
            return (F[i][j+1] - F[i][columns-2])/2.0;
        }
        double D2dxi2_E(int i,int j) {
            return (F[i][j+1] -2.0*F[i][j] + F[i][columns-2]);;
        }
        double D2dxideta_E(int i,int j) {
            return (Ddxi_E(i-1,j) - Ddxi_E(i+1,j))/2.0;
        }
        //Noroeste
        double D2dxideta_NW(int i,int j) {
            return (Ddeta_N(i,j+1) - Ddeta_N(i,columns-2))/2.0;
        }
        //Nordeste
        double D2dxideta_NE(int i,int j) {
            return (Ddeta_N(i,j+1) - Ddeta_N(i,columns-2))/2.0;
        }
        //Sudeste
        double D2dxideta_SE(int i,int j) {
            return (Ddeta_S(i,1) - Ddeta_S(i,j-1))/2.0;
        }
        //Sudoeste
        double D2dxideta_SW(int i,int j) {
            return (Ddeta_S(i,1) - Ddeta_S(i,j-1))/2.0;
        }
        /*
        //Oeste
        double Ddxi_W(int i,int j) {
            return (-F[i][j+2] + 4*F[i][j+1] - 3*F[i][j])/2;
        }
        double D2dxi2_W(int i,int j) {
            return (-F[i][j+3] + 4*F[i][j+2] - 5*F[i][j+1] + 2*F[i][j]);
        }
        double D2dxideta_W(int i,int j) {
            return (Ddxi_W(i-1,j) - Ddxi_W(i+1,j))/2;
        }
        //Leste
        double Ddxi_E(int i,int j) {
            return (F[i][j-2] - 4*F[i][j-1] + 3*F[i][j])/2;
        }
        double D2dxi2_E(int i,int j) {
            return (F[i][j-3] - 4*F[i][j-2] + 5*F[i][j-1] - 2*F[i][j]);
        }
        double D2dxideta_E(int i,int j) {
            return (Ddxi_E(i-1,j) - Ddxi_E(i+1,j))/2;
        }
        //Noroeste
        double D2dxideta_NW(int i,int j) {
            return (-Ddeta_N(i,j+2) + 4*Ddeta_N(i,j+1) - 3*Ddeta_N(i,j))/2;
        }
        //Nordeste
        double D2dxideta_NE(int i,int j) {
            return (Ddeta_N(i,j-2) - 4*Ddeta_N(i,j-1) + 3*Ddeta_N(i,j))/2;
        }
        //Sudeste
        double D2dxideta_SE(int i,int j) {
            return (Ddeta_S(i,j-2) - 4*Ddeta_S(i,j-1) + 3*Ddeta_S(i,j))/2;
        }
        //Sudoeste
        double D2dxideta_SW(int i,int j) {
            return (-Ddeta_S(i,j+2) + 4*Ddeta_S(i,j+1) - 3*Ddeta_S(i,j))/2;
        }*/
        
};

bool checkConvergence(std::vector<std::vector<double>> M, std::vector<std::vector<double>> M0, int lin, int col);


int main() {
    std::cout << "esta rodando\n";
    int lines;
    int columns;
    std::ifstream malha;
    std::ofstream tecplot;
    std::string linha;
    // ler o valor de "I" e armazenar em "columns" e "J" em "lines" do arquivo
    //malha.open("OgridNACA20011.grid", std::ios::in);//Arquivo 1 .grid
    malha.open("OgridNACA20011(medio).grid", std::ios::in);
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
    
    //
    x.calculateMetrics();
    y.calculateMetrics();

    
    
    double density = 1.2;
    double viscosity = 1.48e-5; //viscosidade cinematica
    double Vinfinity = 20;
    double Pinfinity = 1.0e5;

    scalarField u0(lines,columns,x,y);
    scalarField v0(lines,columns,x,y);
    scalarField P0(lines,columns,x,y);
    
    scalarField u(lines,columns,x,y);
    scalarField v(lines,columns,x,y);
    scalarField P(lines,columns,x,y);

    scalarField uTemp(lines,columns,x,y);
    scalarField vTemp(lines,columns,x,y);
    scalarField PTemp(lines,columns,x,y);

    for (int i = 0; i<lines; i++) {
        for(int j = 0; j<columns; j++) {
            if (std::isnan(u.J[i][j]) || std::isnan(u.DJdeta[i][j]) || std::isnan(u.DJdxi[i][j])) {
                std::cout << "deu nan i="<<i<<", j="<<j;
            }
        }
    }
    
    
    for (int j = 0; j < columns; j++) {
        u.F[0][j] = Vinfinity;
        P.F[0][j] = Pinfinity;
    }
    /*
    for (int iteration = 0; iteration < 50000; iteration++){
        for (int i = 1; i < lines-1; i++) {
            u.F[i][0] = (u.F[i+1][0] + u.F[i-1][0] + u.F[i][1] + u.F[i][columns -2])/4;
            
            P.F[i][0] = (P.F[i+1][0] + P.F[i-1][0] + P.F[i][1] + P.F[i][columns -2])/4;
            
            for(int j = 1; j < columns-1; j++) {
                u.F[i][j] = (u.F[i+1][j] + u.F[i-1][j] + u.F[i][j+1] + u.F[i][j-1])/4;
                P.F[i][j] = (P.F[i+1][j] + P.F[i-1][j] + P.F[i][j+1] + P.F[i][j-1])/4;
            }
            u.F[i][columns-1] = u.F[i][0];
            P.F[i][columns-1] = P.F[i][0];
        }
    }
        */
    
    for (int i = 1; i < lines-1; i++) {
        for (int j = 0; j< columns; j++) {
            u.F[i][j] = Vinfinity;
        }
    }
    
    for (int i = 1; i<lines; i++) {
        for (int j = 0; j< columns; j++) {   
            P.F[i][j] = Pinfinity;
        }
    }
    
    int iterations=100000;
    double dt=0.000001;

    std::vector<double> tempVector;
    std::vector<std::vector<double>> A,B,C,D,E;
    std::vector<std::vector<double>> DP2dxi2, DP2deta2, DP2dxideta, DPdxi, DPdeta;
    std::vector<std::vector<double>> P_N, P_NE, P_E, P_SE, P_S, P_SW, P_W, P_NW, P_C, P_NN, P_NNN;
    std::vector<double> tempVector2;
    std::vector<double> Pconstants;
    for (int i = 0; i < (lines-2)*(columns-1); i++) {
        tempVector.push_back(0.0);
        Pconstants.push_back(0.0);
    }
    for (int j = 0; j< columns; j++) {
        tempVector.push_back(0.0);
    }
    for (int i = 0; i < lines; i++) {
        A.push_back(tempVector);
        B.push_back(tempVector);
        C.push_back(tempVector);
        D.push_back(tempVector);
        E.push_back(tempVector);
        P_N.push_back(tempVector);
        P_NE.push_back(tempVector);
        P_E.push_back(tempVector);
        P_SE.push_back(tempVector);
        P_S.push_back(tempVector);
        P_SW.push_back(tempVector);
        P_W.push_back(tempVector);
        P_NW.push_back(tempVector);
        P_C.push_back(tempVector);
        P_NN.push_back(tempVector);
        P_NNN.push_back(tempVector);
        DP2dxi2.push_back(tempVector);
        DP2deta2.push_back(tempVector);
        DP2dxideta.push_back(tempVector);
        DPdxi.push_back(tempVector);
        DPdeta.push_back(tempVector);
    }
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j< columns; j++) {
            A[i][j] = std::pow(x.Ddeta[i][j],2.0) + std::pow(y.Ddeta[i][j],2.0);
            B[i][j] = std::pow(x.Ddxi[i][j],2.0) + std::pow(y.Ddxi[i][j],2.0);
            C[i][j] = x.Ddxi[i][j]*x.Ddeta[i][j] + y.Ddxi[i][j]*y.Ddeta[i][j];
            D[i][j] = x.D2dxideta[i][j]*x.Ddeta[i][j] - x.D2deta2[i][j]*x.Ddxi[i][j] - y.D2deta2[i][j]*y.Ddxi[i][j] + y.D2dxideta[i][j]*y.Ddeta[i][j];
            E[i][j] = x.D2dxideta[i][j]*x.Ddxi[i][j] - x.D2dxi2[i][j]*x.Ddeta[i][j] - y.D2dxi2[i][j]*y.Ddeta[i][j] + y.D2dxideta[i][j]*y.Ddxi[i][j];
        }
    }
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j< columns; j++) {
            DP2dxi2[i][j] = P.J[i][j]*A[i][j]/std::pow(P.J[i][j], 3.0);
            DP2deta2[i][j] = P.J[i][j]*B[i][j]/std::pow(P.J[i][j], 3.0);
            DP2dxideta[i][j] = - 2.0*P.J[i][j]*C[i][j]/std::pow(P.J[i][j], 3.0);
            DPdxi[i][j] = (P.J[i][j]*D[i][j] - P.DJdxi[i][j]*A[i][j] + P.DJdeta[i][j]*C[i][j])/std::pow(P.J[i][j], 3.0);
            DPdeta[i][j] = (P.J[i][j]*E[i][j] + P.DJdxi[i][j]*C[i][j] - P.DJdeta[i][j]*B[i][j])/std::pow(P.J[i][j], 3.0);
        }
    }
    /*
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j< columns; j++) {
            P_C[i][j] = -2.0*DP2dxi2[i][j] -2.0*DP2deta2[i][j] - (3.0/2.0)*DPdxi[i][j];
            P_E[i][j] = DP2dxi2[i][j] + (4.0/2.0)*DPdxi[i][j];
            P_W[i][j] = DP2dxi2[i][j];
            P_N[i][j] = DP2deta2[i][j] + DPdeta[i][j]/2.0;
            P_S[i][j] = DP2deta2[i][j] - DPdeta[i][j]/2.0;
            P_NE[i][j] = DP2dxideta[i][j]/4.0;
            P_NW[i][j] = - DP2dxideta[i][j]/4.0;
            P_SE[i][j] = - DP2dxideta[i][j]/4.0;
            P_SW[i][j] = DP2dxideta[i][j]/4.0;
            P_EE[i][j] = - DPdxi[i][j]/2.0;
        }
    }
        */
       
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j< columns; j++) {
            P_C[i][j] = -2.0*DP2dxi2[i][j] -2.0*DP2deta2[i][j];
            P_E[i][j] = DP2dxi2[i][j] + DPdxi[i][j]/2.0;
            P_W[i][j] = DP2dxi2[i][j] - DPdxi[i][j]/2.0;
            P_N[i][j] = DP2deta2[i][j] + DPdeta[i][j]/2.0;
            P_S[i][j] = DP2deta2[i][j] - DPdeta[i][j]/2.0;
            P_NE[i][j] = DP2dxideta[i][j]/4.0;
            P_NW[i][j] = - DP2dxideta[i][j]/4.0;
            P_SE[i][j] = - DP2dxideta[i][j]/4.0;
            P_SW[i][j] = DP2dxideta[i][j]/4.0;
        }
    }
    /*    
    for (int i = 3; i < lines; i++) {
        for (int j = 0; j< columns; j++) {
            P_C[i][j] = -2.0*DP2dxi2[i][j] - (20.0/12.0)*DP2deta2[i][j];
            P_E[i][j] = DP2dxi2[i][j] + DPdxi[i][j]/2.0;
            P_W[i][j] = DP2dxi2[i][j] - DPdxi[i][j]/2.0;
            P_N[i][j] = (6.0/12.0)*DP2deta2[i][j] + DPdeta[i][j]/2.0;
            P_NN[i][j] = (4.0/12.0)*DP2deta2[i][j];
            P_NNN[i][j] = (-1.0/12.0)*DP2deta2[i][j];
            P_S[i][j] = (11.0/12.0)*DP2deta2[i][j] - DPdeta[i][j]/2.0;
            
            P_NE[i][j] = DP2dxideta[i][j]/4.0;
            P_NW[i][j] = - DP2dxideta[i][j]/4.0;
            P_SE[i][j] = - DP2dxideta[i][j]/4.0;
            P_SW[i][j] = DP2dxideta[i][j]/4.0;
        }
    }*/
    

    for (int iteration = 0; iteration < iterations; iteration++) {
        
        P0.F = P.F;
        u0.F = u.F;
        v0.F = v.F;
        uTemp.F = u.F;
        vTemp.F = v.F;
        /*
        for (int i = 1; i < lines-1; i++) {
            uTemp.F[i][0] = u0.F[i][0] + dt*(-(2.0*u0.F[i][0]*u0.Ddx_W(i,0) + v0.F[i][0]*u0.Ddy_W(i,0)+u0.F[i][0]*u0.Ddy_W(i,0) + P0.Ddx(i,0)/density) + viscosity*(u0.D2dx2_W(i,0)+u0.D2dy2_W(i,0)));
            for (int j = 1; j < columns -1; j++) {
                uTemp.F[i][j] = u0.F[i][j] + dt*(-(2.0*u0.F[i][j]*u0.Ddx(i,j) + v0.F[i][j]*u0.Ddy(i,j)+u0.F[i][j]*v0.Ddy(i,j) + P0.Ddx(i,j)/density) + viscosity*(u0.D2dx2(i,j)+u0.D2dy2(i,j)));
                //+ P0.Ddx(i,j)/density
            }
            //uTemp.F[i][columns-1] = u0.F[i][columns-1] + dt*(-(u0.F[i][columns-1]*u0.Ddx_E(i,columns-1)*2 + v0.F[i][columns-1]*u0.Ddy_E(i,0) + u0.F[i][columns-1]*v0.Ddy_E(i,columns-1)) + viscosity*(u0.D2dx2_E(i,columns-1)+u0.D2dy2_E(i,columns-1)));
            uTemp.F[i][columns-1] = uTemp.F[i][0];
        }

        for (int i = 1; i < lines-1; i++) {
            vTemp.F[i][0] = v0.F[i][0] + dt*(-(u0.F[i][0]*v0.Ddx_W(i,0)+v0.F[i][0]*u0.Ddx_W(i,0) + 2.0*v0.F[i][0]*v0.Ddy_W(i,0) + P0.Ddy_W(i,0)/density) + viscosity*(v0.D2dx2_W(i,0)+v0.D2dy2_W(i,0)));
            for (int j = 1; j < columns -1; j++) {
                vTemp.F[i][j] = v0.F[i][j] + dt*(-(u0.F[i][j]*v0.Ddx(i,j)+v0.F[i][j]*u0.Ddx(i,j) + 2.0*v0.F[i][j]*v0.Ddy(i,j) + P0.Ddy(i,j)/density) + viscosity*(v0.D2dx2(i,j)+v0.D2dy2(i,j)));
                //+ P0.Ddy(i,j)/density
            }
            //vTemp.F[i][columns-1] = v0.F[i][columns-1] + dt*(-(v0.F[i][columns-1]*u0.Ddx_E(i,columns-1) + u0.F[i][columns-1]*v0.Ddx_E(i,columns-1) + v0.F[i][columns-1]*v0.Ddy_E(i,columns-1)*2 ) + viscosity*(v0.D2dx2_E(i,columns-1)+v0.D2dy2_E(i,columns-1)));
            vTemp.F[i][columns-1] = vTemp.F[i][0];
        }
        */
        for (int i = 1; i < lines-1; i++) {
            uTemp.F[i][0] = u0.F[i][0] + dt*(-(u0.F[i][0]*u0.Ddx_W(i,0) + v0.F[i][0]*u0.Ddy_W(i,0) + P0.Ddx(i,0)/density) + viscosity*(u0.D2dx2_W(i,0)+u0.D2dy2_W(i,0)));
            for (int j = 1; j < columns -1; j++) {
                uTemp.F[i][j] = u0.F[i][j] + dt*(-(u0.F[i][j]*u0.Ddx(i,j) + v0.F[i][j]*u0.Ddy(i,j)+ P0.Ddx(i,j)/density) + viscosity*(u0.D2dx2(i,j)+u0.D2dy2(i,j)));
                //+ P0.Ddx(i,j)/density
            }
            //uTemp.F[i][columns-1] = u0.F[i][columns-1] + dt*(-(u0.F[i][columns-1]*u0.Ddx_E(i,columns-1)*2 + v0.F[i][columns-1]*u0.Ddy_E(i,0) + u0.F[i][columns-1]*v0.Ddy_E(i,columns-1)) + viscosity*(u0.D2dx2_E(i,columns-1)+u0.D2dy2_E(i,columns-1)));
            uTemp.F[i][columns-1] = uTemp.F[i][0];
        }

        for (int i = 1; i < lines-1; i++) {
            vTemp.F[i][0] = v0.F[i][0] + dt*(-(u0.F[i][0]*v0.Ddx_W(i,0) + v0.F[i][0]*v0.Ddy_W(i,0) + P0.Ddy_W(i,0)/density) + viscosity*(v0.D2dx2_W(i,0)+v0.D2dy2_W(i,0)));
            for (int j = 1; j < columns -1; j++) {
                vTemp.F[i][j] = v0.F[i][j] + dt*(-(u0.F[i][j]*v0.Ddx(i,j) + v0.F[i][j]*v0.Ddy(i,j) + P0.Ddy(i,j)/density) + viscosity*(v0.D2dx2(i,j)+v0.D2dy2(i,j)));
                //+ P0.Ddy(i,j)/density
            }
            //vTemp.F[i][columns-1] = v0.F[i][columns-1] + dt*(-(v0.F[i][columns-1]*u0.Ddx_E(i,columns-1) + u0.F[i][columns-1]*v0.Ddx_E(i,columns-1) + v0.F[i][columns-1]*v0.Ddy_E(i,columns-1)*2 ) + viscosity*(v0.D2dx2_E(i,columns-1)+v0.D2dy2_E(i,columns-1)));
            vTemp.F[i][columns-1] = vTemp.F[i][0];
        }

        bool deuNan = false;
        for (int i = 0; i<lines; i++) {
            for(int j = 0; j<columns; j++) {
                if (std::isnan(uTemp.F[i][j]) || std::isnan(vTemp.F[i][j])) {
                    deuNan = true;
                }
            }
        }
        if (deuNan) {
            std::cout << "deu nan \n";
            break;
        }
        auto start = std::chrono::high_resolution_clock::now();
         
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        typedef Eigen::Triplet<double> TripletType;
        typedef Eigen::VectorXd VectorType;
        
        std::vector<TripletType> tripletList;
        for (int i = 1; i < lines-1; i++) {
            for (int j = 0; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                if (j == 0) {
                    Pconstants[id] = density*(uTemp.Ddx_W(i,j) + vTemp.Ddy_W(i,j))/dt;
                }
                else {
                    Pconstants[id] = density*(uTemp.Ddx(i,j) + vTemp.Ddy(i,j))/dt;
                }
            }
        }
        
        for (int i = 1; i < lines-2; i++) {
            for (int j = 0; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id, P_C[i][j]));
            }
        }
        for (int i = 1; i < lines-2; i++) {
            for (int j = 0; j < columns-2; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id+1, P_E[i][j]));
            }
        }
        /*
        for (int i = 1; i < lines-1; i++) {
            for (int j = 0; j < columns-3; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id+1, P_EE[i][j]));
            }
        }*/
        for (int i = 1; i < lines-2; i++) {
            for (int j = 1; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id-1, P_W[i][j]));
            }
        }
        for (int i = 2; i < lines-1; i++) {
            for (int j = 0; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id - (columns-1), P_N[i][j]));
            }
        }/*
        for (int i = 3; i < lines-1; i++) {
            for (int j = 0; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id - 2*(columns-1), P_NN[i][j]));
            }
        }
        for (int i = 4; i < lines-1; i++) {
            for (int j = 0; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id - 3*(columns-1), P_NNN[i][j]));
            }
        }*/
        for (int i = 1; i < lines-2; i++) {
            for (int j = 0; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id + (columns-1), P_S[i][j]));
            }
        }
        for (int i = 2; i < lines-1; i++) {
            for (int j = 0; j < columns-2; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id+1 - (columns-1), P_NE[i][j]));
            }
        }
        for (int i = 2; i < lines-1; i++) {
            for (int j = 1; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id-1 - (columns-1), P_NW[i][j]));
            }
        }
        for (int i = 1; i < lines-2; i++) {
            for (int j = 0; j < columns-2; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id+1 + (columns-1), P_SE[i][j]));
            }
        }
        for (int i = 1; i < lines-2; i++) {
            for (int j = 1; j < columns-1; j++) {
                int id = (i-1)*(columns - 1)+j;
                tripletList.push_back(TripletType(id, id-1 + (columns-1), P_SW[i][j]));
            }
        } // centro restricoes
        for (int j = 0; j < columns - 1; j++) {
            int i = lines-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id, P_C[i][j] + P_S[i][j]));
        } // leste restricoes
        for (int j = 0; j < columns - 2; j++) {
            int i = lines-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id+1, P_E[i][j] + P_SE[i][j]));
        }
        for (int i = 1; i < lines-2; i++) {
            int j = columns-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id+1-(columns-1), P_E[i][j]));
        } // leste leste restricoes
        /*
        for (int i = 1; i < lines-1; i++) {
            int j = columns-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id+2-(columns-1), P_EE[i][j]));
        }
        for (int i = 1; i < lines-1; i++) {
            int j = columns-3;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id+2-(columns-1), P_EE[i][j]));

        } */// oeste restricoes 
        for (int j = 1; j < columns - 1; j++) {
            int i = lines-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id-1, P_W[i][j] + P_SW[i][j]));
        }
        for (int i = 1; i < lines-2; i++) {
            int j = 0;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id-1+(columns-1), P_W[i][j]));
        } // nordeste restricoes
        for (int i = 2; i < lines-1; i++) {
            int j = columns-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id-(columns-1)+1-(columns-1), P_NE[i][j]));
        } //noroeste restricoes
        for (int i = 2; i < lines-1; i++) {
            int j = 0;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id-(columns-1)-1+(columns-1), P_NW[i][j]));
        } //sudeste restricoes
        for (int i = 1; i < lines-2; i++) {
            int j = columns-2;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id+(columns-1)+1-(columns-1), P_SE[i][j]));
        } //sudoeste restricoes
        for (int i = 1; i < lines-2; i++) {
            int j = 0;
            int id = (i-1)*(columns - 1)+j;
            tripletList.push_back(TripletType(id, id+(columns-1)-1+(columns-1), P_SW[i][j]));
        }
        int i = lines-2;
        int j = columns-2;
        int id = (i-1)*(columns - 1)+j;
        tripletList.push_back(TripletType(id, id+1-(columns-1), P_E[i][j] + P_SE[i][j]));
        j = 0;
        id = (i-1)*(columns - 1)+j;
        tripletList.push_back(TripletType(id, id-1+(columns-1), P_W[i][j] + P_SW[i][j]));
        
        SparseMatrixType M((lines-2)*(columns-1),(lines-2)*(columns-1));
        M.setFromTriplets(tripletList.begin(), tripletList.end());
        M.makeCompressed(); 
        
        VectorType b((lines-2)*(columns-1));
        
        for (int i = 0; i <(lines-2)*(columns-1); i++) {
            b(i) = Pconstants[i];
        }
        
        Eigen::SparseLU<SparseMatrixType> solver;
        solver.compute(M);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Decomposition failed!" << std::endl;
            return 1;
        }

        // 7. Solve the system
        VectorType solution = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solving failed!" << std::endl;
            return 1;
        }

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << duration.count() << " segundos na iteracao "<<iteration<< '\n';

        for (int i = 1; i < lines-1; i++) {
            for (int j = 0; j < columns-1; j++){
                PTemp.F[i][j] = solution[(i-1)*(columns - 1)+j];
            } 
        }
        for (int j = 0; j< columns-1; j++) {
            int i = lines-1;
            PTemp.F[i][j] = PTemp.F[i-1][j];
        }
        for (int i = 1; i < lines; i++) {
            int j = columns-1;
            PTemp.F[i][j] = PTemp.F[i][0];
        }
        for (int i = 0; i<lines; i++) {
            for(int j = 0; j<columns; j++) {
                if (std::isnan(P.F[i][j])) {
                    deuNan = true;
                }
            }
        }
        if (deuNan) {
            std::cout << "deu nan na pressao \n";
            break;
        }
        
        for (int i = 1; i < lines -1; i++) {
            u.F[i][0] = uTemp.F[i][0] - dt*PTemp.Ddx_W(i,0)/density;
            v.F[i][0] = vTemp.F[i][0] - dt*PTemp.Ddy_W(i,0)/density;
            u.F[i][columns-1] = u.F[i][0];
            v.F[i][columns-1] = v.F[i][0];
            for (int j = 1; j < columns -1; j++) {
                u.F[i][j] = uTemp.F[i][j] - dt*PTemp.Ddx(i,j)/density;
                v.F[i][j] = vTemp.F[i][j] - dt*PTemp.Ddy(i,j)/density;
            }
        }
        
        
        for (int i = 1; i < lines; i++) {
            for (int j = 0; j < columns; j++) {
                P.F[i][j] = PTemp.F[i][j] + P0.F[i][j];
            }
        }

        
        bool residualP = checkConvergence(P.F,P0.F,lines,columns);
        bool residualU = checkConvergence(u.F,u0.F,lines,columns);
        bool residualV = checkConvergence(v.F,v0.F,lines,columns);
        
        std::cout<<"iteracao "<<iteration<< " concluida\n" ;
        if (residualP && residualU && residualV) {
            std::cout << "foram realizadas " << iteration << " iteracoes" << '\n'; 
            break;
        }
        
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
                    tecplot << P.F[i][j] << "\n";
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
        
    }
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
                tecplot << P.F[i][j] << "\n";
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

