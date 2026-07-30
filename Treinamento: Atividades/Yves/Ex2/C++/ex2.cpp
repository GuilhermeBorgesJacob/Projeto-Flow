#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

const string PASTA_SAIDA = "./Ex2/";

// PARÂMETROS FÍSICOS
double comprimento;      // Comprimento da barra (m)
double alpha;            // Difusividade térmica (m²/s)
double tempEsquerda;     // Temperatura na borda esquerda (°C)
double tempDireita;      // Temperatura na borda direita (°C)

// PARÂMETROS NUMÉRICOS (GRID)
int Nx;          // Número de pontos da grid (incluindo bordas)
double dx;       // Espaçamento espacial
double dt;       // Passo de tempo
int Nt;          // Número de passos de tempo
double fo;       // Número de Fourier

void salvarGrid(const vector<double>& T, int passo, double tempo) {
    string nome = PASTA_SAIDA + "grid_tempo_" + to_string(passo) + ".csv";
    ofstream arquivo(nome);
    if (!arquivo.is_open()) {
        cerr << "Erro ao criar arquivo!" << endl;
        return;
    }

    // Cabeçalho
    arquivo << "# Passo temporal: " << passo << endl;
    arquivo << "# Tempo: " << tempo << " s" << endl;
    arquivo << "# Nx: " << Nx << ", dx: " << dx << ", dt: " << dt << endl;
    arquivo << "# fo = alpha*dt/dx^2 = " << fo << endl;
    arquivo << "x,T\n";

    for (size_t i = 0; i < T.size(); i++) {
        double x = i * dx;
        arquivo << x << "," << T[i] << "\n";
    }

    arquivo.close();
}

void salvarEvolucao(const vector<vector<double>>& historico,
                    const vector<double>& tempos) {
    string nome = PASTA_SAIDA + "evolucao_completa" + ".csv";
    ofstream arquivo(nome);

    if (!arquivo.is_open()) {
        cerr << "Erro ao criar arquivo" << endl;
        return;
    }

    // Cabeçalho
    arquivo << "x";
    for (size_t n = 0; n < tempos.size(); n++) {
        arquivo << ",T_" << fixed << setprecision(2) << tempos[n] << "s";
    }
    arquivo << "\n";

    // Dados
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        arquivo << x;
        for (size_t n = 0; n < historico.size(); n++) {
            arquivo << "," << historico[n][i];
        }
        arquivo << "\n";
    }

    arquivo.close();
}

void salvarPerfilFinal(const vector<double>& T, double tempo_final) {
    string nome = PASTA_SAIDA + "perfil_final" + ".csv";
    ofstream arquivo(nome);

    if (!arquivo.is_open()) {
        cerr << "Erro ao criar arquivo!" << endl;
        return;
    }

    arquivo << "x,T\n";
    for (size_t i = 0; i < T.size(); i++) {
        double x = i * dx;
        arquivo << x << "," << T[i] << "\n";
    }

    arquivo.close();
    cout << "\n Dados exportados para 'resultados.csv'" << endl;
}

double maximo(const vector<double>& T) {
    double max_val = T[1];
    for (size_t i = 1; i < T.size(); i++) {
        if (T[i] > max_val) max_val = T[i];
    }
    return max_val;
}

int main() {
    cout << "Insira os parâmetros fisicos para a simulação:" << endl;
    cout << "  Comprimento da barra (m): ";
    cin >> comprimento;
    cout << "  Difusividade termica (m²/s): ";
    cin >> alpha;
    cout << "  Temperatura na borda esquerda (°C): ";
    cin >> tempEsquerda;
    cout << "  Temperatura na borda direita (°C): ";
    cin >> tempDireita;

    Nx = 21;                              // Número de pontos (incluindo bordas)
    dx = comprimento / (Nx - 1);          // Espaçamento espacial
    dt = 0.001;                           // Passo de tempo (pode ser ajustado)
    Nt = 2000;                            // Número de passos de tempo
    fo = alpha * dt / (dx * dx);           // Número de Fourier

    // ESTABILIDADE
    if (fo > 0.5) {
        cout << "\n ATENÇÃO: fo = " << fo << " > 0.5" << endl;
        cout << "The error will progressively become larger and will eventually cause the numerical marching solution to 'blow up' on the computer" << endl;
        return 1;
    } else {
        cout << "\nr = " << fo << " <= 0.5" << endl;
        cout << "The numerical solution will proceed in a stable manner";
    }

    // GRID
    vector<double> tempAtual(Nx, 0.0);       // Temperatura atual (tempo n)
    vector<double> tempFutura(Nx, 0.0);      // Temperatura futura (tempo n+1)

    // Condição inicial: T(x,0) = 0°C (barra fria)
    // Bordas têm as temperaturas fixas

    // CONDIÇÕES DE CONTORNO
    tempAtual[0] = tempEsquerda;       // x = 0
    tempAtual[Nx-1] = tempDireita;     // x = comprimento
    tempFutura[0] = tempEsquerda;
    tempFutura[Nx-1] = tempDireita;

    // MÉTODO EXPLÍCITO
    vector<vector<double>> historico;
    vector<double> tempos;

    historico.push_back(tempAtual);
    tempos.push_back(0.0);

    cout << "\n Evolução temporal" << endl;

    for (int n = 1; n <= Nt; n++) {
        double tempo_atual = n * dt;

        // EQUAÇÃO (4.36) DO ANDERSON
        // T_i^{n+1} = T_i^n + fo * (T_{i+1}^n - 2*T_i^n + T_{i-1}^n)
        for (int i = 1; i < Nx - 1; i++) {
            tempFutura[i] = tempAtual[i] + fo * (tempAtual[i+1] - 2*tempAtual[i] + tempAtual[i-1]);
        }

        tempAtual.swap(tempFutura);

        int frequencia_salvamento = 200;

        if ((n + 1) % frequencia_salvamento == 0 || n == Nt - 1) {
            historico.push_back(tempAtual);
            tempos.push_back(tempo_atual);

            salvarGrid(tempAtual, n + 1, tempo_atual);

            cout << "  t = " << fixed << setprecision(1) << tempo_atual
                 << " s | T_max = " << fixed << setprecision(2) << maximo(tempAtual) << " °C" << endl;
        }

        // Manter temperatura fixa nas bordas
        tempAtual[0] = tempEsquerda;
        tempAtual[Nx-1] = tempDireita;
        tempFutura[0] = tempEsquerda;
        tempFutura[Nx-1] = tempDireita;
    }

    salvarPerfilFinal(tempAtual, Nt * dt);
    salvarEvolucao(historico, tempos);

    return 0;
}
