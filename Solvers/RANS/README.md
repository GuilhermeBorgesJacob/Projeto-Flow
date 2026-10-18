modelo k-e
Algumas ressalvas--- o modelo de turbulencia não está rodando sem erros ainda, apenas com a turbulência desligada
Estou tentando uma transição do MacCormack para AUSM mas por enquanto parte segue o passo do MacCormack e parte dos calculos seguem o passo AUSM
Parte dos comentários e correção de bugs feita provisoriamente por IA, mas mudarei isso logo
Por enquanto não é capaz de aceitar todas os tipos de grid mas quando eu arrumar as contas focarei nisso
Variáveis estão direto no código, tentarei colocar num arquivo separado depois que ele estiver funcionando bem
Muitas gambiarras provisórias feitas parte por mim e parte pela IA para evitar o erro NaN por hora
para compilar: g++ -O3 -march=native -fopenmp -std=c++17 solver_rans_ke_M03.cpp -o rans
para rodar: OMP_NUM_THREADS=4 ./rans cilindro_turb.grid 20000 0.0 0.1 1  #### esse ultimo "1" é para turbulência ligada, o que não está funcionando atualmente, "0" se for testar normal
