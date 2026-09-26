#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

// Atualiza com Diferenças Finitas
void atualizarTemperatura(
    const vector<double>& T,
    vector<double>& Tnovo,
    double Fo)
{
    const size_t N = T.size();

    for(size_t i = 1; i < N - 1; ++i)
    {
        Tnovo[i] =
            T[i] +
            Fo *
            (T[i+1]
            -2.0*T[i]
            +T[i-1]);
    }

    // Condições de contorno

    Tnovo.front() = 200.0;
    Tnovo.back()  = 0.0;
}

//------------------------------------------------------------

int main()
{

    //----------------------------------------------------
    // Constantes
    //----------------------------------------------------

    constexpr double L = 2.0;
    constexpr int N = 101;

    constexpr double alpha = 0.01;

    constexpr double dt = 0.001;

    constexpr int passos = 100000;

    constexpr int salvarCada = 50;

    constexpr double dx = L/(N-1);

    constexpr double Fo = alpha*dt/(dx*dx);

    //----------------------------------------------------

    if(Fo > 0.5)
    {
        cout << "Metodo instavel\n";
        cout << "Fo = " << Fo << endl;
        return 0;
    }

    //----------------------------------------------------

    vector<double> T(N,20.0);
    vector<double> Tnovo(N);

    T.front() = 200.0;
    T.back() = 0.0;

    //----------------------------------------------------

    ofstream arquivo("temperaturas.csv");

    arquivo<<"Tempo";

    for(int i=0;i<N;i++)
        arquivo<<",x"<<i;

    arquivo<<"\n";

    //----------------------------------------------------

    cout<<fixed<<setprecision(2);

    //----------------------------------------------------

    for(int n=0;n<=passos;n++)
    {
    //============================================
         //---------------------------------------------
        // Salva apenas algumas iterações para otimizar
        //---------------------------------------------

        if(n % salvarCada == 0)
        {
            arquivo<<n*dt;

            for(double temp : T)
                arquivo<<","<<temp;

            arquivo<<"\n";

            cout
            << "Tempo = "
            << setw(6)
            << n*dt
            << " s\r";

            cout.flush();
        }

        //---------------------------------------------

        atualizarTemperatura(T,Tnovo,Fo);

        //---------------------------------------------

        swap(T,Tnovo);

    }

    arquivo.close();

    cout<<"\n\nSimulacao finalizada.\n";

    return 0;
}
