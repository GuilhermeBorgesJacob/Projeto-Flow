#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define BOLD    "\033[1m"

double deltaY = 0.10;
double gradiente_valorExato = 1582.0;
double constante = 3.7373e-7 * 12; // 4.48476e-6
double tauW_valorExato = constante * gradiente_valorExato;

void exportarCSV(const vector<double>& gradientes,
                 const vector<double>& tauW,
                 const vector<string>& ordens) {
                     ofstream arquivo("resultados.csv");

    if (!arquivo.is_open()) {
        cerr << "Erro: Não foi possível criar o arquivo CSV!" << endl;
        return;
    }

    arquivo << "Ordem,Gradiente (du/dy),Tau_w (Pa),Erro (%)\n";

    for(size_t i = 0; i < gradientes.size(); i++){
        double erro = (1 - tauW[i]/tauW_valorExato) * 100;
        arquivo << ordens[i] << ","
                << fixed << setprecision(6) << gradientes[i] << ","
                << fixed << setprecision(6) << tauW[i] << ","
                << fixed << setprecision(4) << erro << "\n";
    }

    arquivo.close();
    cout << "\n Dados exportados para 'resultados.csv'" << endl;
}

void erro (double tauW, double tauW_valorExato){
    cout << "Erro: " << (1 - tauW/tauW_valorExato)*100 << "%" << endl;
}

int main(){
    vector<double> velocidades = {0, 150.54, 286.77, 410.03};

    double primeiraOrdem = (velocidades[1] - velocidades[0]) / deltaY;
    double segundaOrdem = (-3*velocidades[0] + 4*velocidades[1] - velocidades[2]) / (2*deltaY);
    double terceiraOrdem = (-11*velocidades[0] + 18*velocidades[1] - 9*velocidades[2] + 2*velocidades[3]) / (6 * deltaY);

    vector<double> gradientes = {primeiraOrdem, segundaOrdem, terceiraOrdem};
    vector<double> tauW;
    vector<string> ordens = {"1ª", "2ª", "3ª"};

    for(double gradiente : gradientes) {
        tauW.push_back(constante * gradiente);
    }

    cout << BOLD << CYAN << "╔══════════════════════════════════════════════════╗" << RESET << endl;
    cout << BOLD << CYAN << "║          ANÁLISE DE DIFERENÇAS FINITAS           ║" << RESET << endl;
    cout << BOLD << CYAN << "╚══════════════════════════════════════════════════╝" << RESET << endl;
    cout << "\n" << BOLD << "Valor Exato: " << RESET << tauW_valorExato << " lb/ft²\n" << endl;

    for(int i = 0; i < 3; i++){
        vector<string> simbolos = {"❶", "❷", "❸"};

        cout << BOLD << YELLOW << "\n┌──────────────────────────────────────────────────┐" << RESET << endl;
        cout << BOLD << YELLOW << "│ " << simbolos[i] << " DIFERENÇA DE " << ordens[i] << " ORDEM" << RESET << endl;
        cout << BOLD << YELLOW << "├──────────────────────────────────────────────────┤" << RESET << endl;
        cout << BOLD << YELLOW << "│ Gradiente du/dy: " << fixed << setprecision(2) << setw(12) << gradientes[i] << RESET << endl;
        cout << BOLD << YELLOW << "│ τ_w:             " << scientific << setprecision(4) << setw(12) << tauW[i] << " lb/ft²" << RESET << endl;
        cout << BOLD << YELLOW << "│ ";
        erro(tauW[i], tauW_valorExato);
        cout << BOLD << YELLOW << "└──────────────────────────────────────────────────┘" << RESET << endl;
    }
    exportarCSV(gradientes, tauW, ordens);
    return 0;
}
