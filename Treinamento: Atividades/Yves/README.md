## Como executar

### Exercício 1: Exemplo 4.2 do Anderson

```bash
cd "Ex 1/C++"
g++ -Wall -std=c++17 gpda-ex1.cpp -o gpda-ex1
./gpda-ex1
```

Como eu tinha achado meio simples, decidi adicionar umas firulas no terminal kkkk

### Exercício 2: Solução numérica da Eq. 3.28

```bash
"Ex 2/C++"
g++ -Wall -std=c++17 ex2.cpp -o ex2
./ex2_interativo
``` 
**OBS.:** Para este exercício, decidi deixar com que o usuário colocasse uma input qualquer. Ou seja, os resultados variam de acordo com a entrada. Na pasta do exercício, é possível visualizar gráficos e fotos gerados para as seguintes entradas: temperatura esquerda = 100°C, temperatura direita = 0°C, difusividade térmica = 0.01m²/s e comprimento = 1m.

## Visualizar com Python
```bash
cd "Ex X/Python"
python -m venv venv
source venv/bin/activate
pip install pandas matplotlib numpy
python grafico.py
```

**OBS.:** alterar o primeira linha de comando (cd "Ex X/Python") de acordo com a pasta desejada. 
Problemas na execução do código podem ser por conta de falta de algum dos pacotes necessários (pandas, matplotlib ou numpy)

## Estrutura da pasta
```
├── Ex 1/ 
│ ├── C++/
│ │ └── gpda-ex1.cpp 
│ ├── Python/
│ │ └── grafico.py 
│ ├── analise_resultados.png 
│ └── resultados.csv 
│
└── Ex 2/ 
├── C++/
│ └── ex2_interativo.cpp
├── Python/
│ └── grafico.py 
├── evolucao_completa.csv 
└── perfil_final.csv 
```
