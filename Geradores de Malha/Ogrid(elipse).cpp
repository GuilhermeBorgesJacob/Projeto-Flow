// Autor: Guilherme Borges Jacob
// Data: outubro de 2025
// Descrição: este programa distribui pontos ao redor de um aerofólio NACA de quatro digitos de maneira
// estruturada a fim de criar uma malha em formato de "O" que pode ser utilizada para análises de CFD

#define _USE_MATH_DEFINES //utilizado para obter o valor de pi
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>

long double NACAsup(int m, int p, int t, long double c, long double x);
long double NACAinf(int m, int p, int t, long double c, long double x);
bool checkConvergence(std::vector<std::vector<long double>> M, std::vector<std::vector<long double>> M0, int lines, int columns);

int main() {
    //Numeros NACA (mptt)
    int m = 0;
    int p = 0;
    int t = 11;
    //
    
    long double c = 1.0;     //comprimento da corda
    
    int lines = 120;            //84,121           // (48,81) (56,93)
    int columns = 150;        //SEMPRE IMPAR  //quantidade de pontos do aerofolio
    
    long double maxD = 15.0/(2.5*t/100);   //raio do circulo externo

    int iterations = 100000; //maximo de iteracoes

    //********************************************************************************************
    //inicializacao das matrizes x (posicao dos pontos (i,j) em no eixo x)
    //e y (posicao dos pontos (i,j) em no eixo y)
    std::vector<std::vector<long double>> x;
    std::vector<std::vector<long double>> y;

    std::vector<long double> temp;
    for (int i = 0; i < columns; i++) {
        temp.push_back(0.0);
    }
    for (int i = 0; i < lines; i++) {
        x.push_back(temp);
        y.push_back(temp);
    }
    
    //*********************************CONDICOES DE CONTORNO*************************************
    
    //distribuição de pontos distantes do aerofolio (primeira linha das matrizes)
    for (int i = 0; i < columns; i++) {
        x[0][i] = c/2 + maxD*c*std::cos(-(2*M_PI*i)/(columns-1));
        y[0][i] = maxD*c*std::sin(-(2*M_PI*i)/(columns-1))*(2.5*t/100);
    }

    //curva 2 & curva 3 (laterais das matrizes)
    long double d = 0.1;
    for (int i = 0; i < lines; i++) {
        x[i][columns-1] = d*std::exp((lines-1-i)*std::log(((maxD-0.5)*c+d)/d)/(lines-1)) + c - d;
        x[i][0] = d*std::exp((lines-1-i)*std::log(((maxD-0.5)*c+d)/d)/(lines-1)) + c - d;
        y[i][columns-1] = 0;
        y[i][0] = 0;
    }
    
    //distribuição de pontos no aerofolio (ultima linha das matrizes)
    long double var;
    for (int i = 0; i < (columns+1)/2; i++) {
        //var = M_PI*(1.0 - 1.0*i/((columns+1)/2-1));
        var = M_PI*((columns-1)/2-i)/((columns-1)/2);
        var = (1.0 - std::cos(var))/2.0;
        /*
        var = 1.0*((columns-1)/2-i)/((columns-1)/2);
        if (var >= 0.5) {
            var = 0.5*c*(5*std::atan((2*var - 1)*(std::sqrt(5 + 2*std::sqrt(5))))/(2*M_PI) + 1);
        }
        else {
            var = c*2*std::pow(var, 2);
        }
        */
        x[lines-1][i] = c*var;
        y[lines-1][i] = NACAinf(m, p, t, c, var);
        
    }    
    for (int i = (columns+1)/2; i < columns; i++) {
        //var = M_PI*(i - (columns+1)/2)/(columns-1 - (columns+1)/2);
        var = M_PI*(i-(columns-1)/2)/((columns-1)/2);
        var = (1.0 - std::cos(var))/2.0;
        /*
        var = 1.0*(i-(columns-1)/2)/((columns-1)/2);
        if (var <= 0.5) {
            var = c*2*std::pow(var, 2);
        }
        else {
            var = 0.5*c*(5*std::atan((2*var - 1)*(std::sqrt(5 + 2*std::sqrt(5))))/(2*M_PI) + 1);
        }
        */

        x[lines-1][i] = c*var;
        y[lines-1][i] = NACAsup(m, p, t, c, var);
        
    }

    //*********************************ITERACOES SOBRE AS MATRIZES X e Y************************************
    std::vector<std::vector<long double>> x0;
    std::vector<std::vector<long double>> y0;
    for (int iteration = 1; iteration <= iterations; iteration++) { 
        x0 = x;
        y0 = y;
        #pragma omp parallel for
        for (int i = 1; i < lines-1; i++) {
            x[i][0] = (x[i][1] + x[i][columns-2] + x[i-1][0] + x[i+1][0])/4;
            //y[i][0] = (y[i][1] + y[i][columns-2] + y[i-1][0] + y[i+1][0])/4;
            for(int j = 1; j < columns-1; j++) {
                x[i][j] = (x[i][j+1] + x[i][j-1] + x[i-1][j] + x[i+1][j])/4;
                y[i][j] = (y[i][j+1] + y[i][j-1] + y[i-1][j] + y[i+1][j])/4;
            }
            x[i][columns-1] = x[i][0];
            //y[i][columns-1] = y[i][0];
        }
        

        if (iteration>3000) {
            if (checkConvergence(x,x0,lines,columns) && checkConvergence(y,y0,lines,columns)) {
                std::cout << "Foram realizadas " << iteration << " iteracoes.";
                break;
            }
        }
    }
        
    std::fstream arquivo;
    //================================================================================
    //Print dos dados .grid
    std::string nome = "OgridNACA2" + std::to_string(m) + std::to_string(p) + std::to_string(t) + ".grid";
    arquivo.open(nome,std::fstream::out|std::fstream::trunc);
    if(arquivo.is_open()){
        arquivo<<"VARIABLES= \"X\", \"Y\""<<std::endl;
        arquivo<<"ZONE T=\"DEFORMED\", I=         "<< columns <<", J=          "<<lines<<" , DATAPACKING=POINT"<<std::endl;
        for(int i=0;i<lines;i++){
            for(int j=0;j<columns;j++){
                arquivo<< x[i][j]<<"  "<<y[i][j]<<std::endl;
            }
        }
        
    }
    arquivo.close();
    //========================================================================================
    
    return 0;
}

long double NACAsup(int m, int p, int t, long double c, long double x) {
    
    long double thicknees = c*0.05*t*(0.2969*std::sqrt(x/c) - 0.1260*x/c - 0.3516*std::pow(x/c, 2) + 0.2843*std::pow(x/c, 3) - 0.1036*std::pow(x/c, 4));
    
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

long double NACAinf(int m, int p, int t, long double c, long double x){

    long double thicknees = c*0.05*t*(0.2969*std::sqrt(x/c) - 0.1260*x/c - 0.3516*std::pow(x/c, 2) + 0.2843*std::pow(x/c, 3) - 0.1036*std::pow(x/c, 4));
    
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

bool checkConvergence(std::vector<std::vector<long double>> M, std::vector<std::vector<long double>> M0, int lines, int columns) {
    long double delta = 1.0e-9;
    
    for (int i = 1; i < lines-1; i++) {
        for (int j = 1; j < columns-1; j++) {
            if (std::abs(M[i][j] - M0[i][j]) >= delta) {
                return false;
            }
        }
    }
    return true;
}
