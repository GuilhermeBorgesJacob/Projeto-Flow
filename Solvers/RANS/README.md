# Solver 2D compressível RANS *k–ε* — MacCormack

Solver em C++ para escoamento compressível bidimensional sobre uma malha estruturada curvilínea, com opção de turbulência pelo modelo padrão **k–ε**. O código calcula as equações de Navier–Stokes/RANS usando diferenças direcionais do tipo MacCormack (preditor–corretor), propriedades viscosas efetivas e dissipação artificial. Descrição do código, equações e melhorias planejadas


> **cautela:** PArte desse readme foi escrito ou corrigido pela IA para ganhar tempo, principalmente a escrita das equações porque ia demorar muito fazendo essa parte. Ou seja, não encontrei nenhum erro mas é possível que haja. 

---

## 1. Equações governantes

O código representa um gás ideal caloricamente perfeito, com propriedades constantes `Cp`, `Cv` e `R`, viscosidade molecular variável pela lei de Sutherland e, quando habilitada, viscosidade turbulenta obtida por *k–ε*.

### 1.1 Conservação de massa

$$
\frac{\partial \rho}{\partial t}
+\frac{\partial (\rho u)}{\partial x}
+\frac{\partial (\rho v)}{\partial y}=0.
$$

Onde \(\rho\) é a densidade e \(u,v\) são as componentes da velocidade média.

### 1.2 Conservação da quantidade de movimento

Na direção \(x\):

$$
\frac{\partial (\rho u)}{\partial t}
+\frac{\partial (\rho u^2+p)}{\partial x}
+\frac{\partial (\rho uv)}{\partial y}
=\frac{\partial \tau_{xx}}{\partial x}
+\frac{\partial \tau_{xy}}{\partial y}.
$$

Na direção \(y\):

$$
\frac{\partial (\rho v)}{\partial t}
+\frac{\partial (\rho uv)}{\partial x}
+\frac{\partial (\rho v^2+p)}{\partial y}
=\frac{\partial \tau_{xy}}{\partial x}
+\frac{\partial \tau_{yy}}{\partial y}.
$$

O tensor de tensões usa a viscosidade efetiva \(\mu_e=\mu+\mu_t\) e inclui a parcela isotrópica turbulenta quando a turbulência está ligada:

$$
\tau_{xx}=\frac{2}{3}\mu_e\left(2\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}\right)-\frac{2}{3}\rho k,
$$

$$
\tau_{yy}=\frac{2}{3}\mu_e\left(2\frac{\partial v}{\partial y}-\frac{\partial u}{\partial x}\right)-\frac{2}{3}\rho k,
\qquad
\tau_{xy}=\mu_e\left(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\right).
$$

Quando `useTurbulence = 0`, a contribuição \(2\rho k/3\) é removida e \(\mu_t=0\).

### 1.3 Conservação da energia média

O estado conservativo de energia usado no código é

$$
U_4=\rho E=\rho\left(e+\frac{u^2+v^2}{2}\right),
\qquad e=C_vT.
$$

A equação é escrita em forma de fluxo:

$$
\frac{\partial (\rho E)}{\partial t}
+\frac{\partial[(\rho E+p)u]}{\partial x}
+\frac{\partial[(\rho E+p)v]}{\partial y}
=\frac{\partial}{\partial x}\left(u\tau_{xx}+v\tau_{xy}-q_x\right)
+\frac{\partial}{\partial y}\left(u\tau_{xy}+v\tau_{yy}-q_y\right),
$$

com fluxo térmico modelado por

$$
q_x=-\lambda_e\frac{\partial T}{\partial x},\qquad
q_y=-\lambda_e\frac{\partial T}{\partial y},
$$

$$
\lambda_e=C_p\left(\frac{\mu}{Pr}+\frac{\mu_t}{Pr_t}\right).
$$

**Convenção importante:** neste arquivo, \(U_4\) não inclui \(\rho k\); a energia cinética turbulenta é transportada separadamente por \(U_5=\rho k\). O fechamento energético e o acoplamento com a produção turbulenta precisam ser verificados por testes de conservação/benchmark antes de se considerar a implementação fisicamente validada. Isso foi feito para facilitar mas pode ser alterado

### 1.4 Modelo de turbulência padrão *k–ε*

A variável \(k\) é a energia cinética turbulenta por unidade de massa e \(\varepsilon\) sua taxa de dissipação. As equações implementadas são, esquematicamente,

$$
\frac{\partial(\rho k)}{\partial t}+\nabla\cdot(\rho\mathbf{u}k)
=\nabla\cdot\left[\left(\mu+\frac{\mu_t}{\sigma_k}\right)\nabla k\right]+P_k-\rho\varepsilon,
$$

$$
\frac{\partial(\rho\varepsilon)}{\partial t}+\nabla\cdot(\rho\mathbf{u}\varepsilon)
=\nabla\cdot\left[\left(\mu+\frac{\mu_t}{\sigma_\varepsilon}\right)\nabla\varepsilon\right]
+C_{\varepsilon1}\frac{\varepsilon}{k}P_k
-C_{\varepsilon2}\rho\frac{\varepsilon^2}{k}.
$$

A viscosidade turbulenta é

$$
\mu_t=\min\left(\rho C_\mu\frac{k^2}{\varepsilon},\;10^5\mu\right).
$$

As constantes padrão no código (alteráveis no main mas precisa recompilar depois) são

| Constante | Valor |
|---|---:|
| \(C_\mu\) | 0,09 |
| \(C_{\varepsilon1}\) | 1,44 |
| \(C_{\varepsilon2}\) | 1,92 |
| \(\sigma_k\) | 1,0 |
| \(\sigma_\varepsilon\) | 1,3 |
| \(Pr_t\) | 0,9 |
| \(\kappa\) | 0,41 |
| \(E\) (lei logarítmica) | 9,793 |

A produção é calculada a partir das derivadas de velocidade e da dilatação compressível; em seguida é limitada numericamente a

$$
0\leq P_k\leq 10\rho\varepsilon.
$$

Os termos de destruição das equações de \(k\) e \(\varepsilon\) recebem tratamento implícito pontual. Isso é uma escolha numérica de estabilização, não é necessariamente uma validação para o modelo com a turbulência ligada

### 1.5 Gás ideal e viscosidade molecular

$$
p=\rho RT,\qquad e=C_vT,\qquad \gamma=\frac{C_p}{C_v}.
$$

A viscosidade molecular é calculada pela lei de Sutherland:

$$
\mu(T)=\mu_0\left(\frac{T}{T_0}\right)^{3/2}\frac{T_0+S}{T+S},
$$

com \(T_0=288{,}16\,\mathrm{K}\), \(\mu_0=1{,}7894\times10^{-5}\,\mathrm{Pa\,s}\) e \(S=110\,\mathrm{K}\). A função limita a temperatura usada nesse cálculo a pelo menos 10 K para evitar problemas.

---

## 2. Método numérico e passo de tempo

### 2.1 Malha curvilínea e transformação de coordenadas

O código recebe coordenadas físicas \(x,y\) de uma malha estruturada e calcula o jacobiano aproximado

$$
J=x_\xi y_\eta-x_\eta y_\xi.
$$

As derivadas físicas são recuperadas pela regra da cadeia:

$$
\frac{\partial \phi}{\partial x}
=\frac{\phi_\xi y_\eta-\phi_\eta y_\xi}{J},
\qquad
\frac{\partial \phi}{\partial y}
=\frac{\phi_\eta x_\xi-\phi_\xi x_\eta}{J}.
$$

As diferenças computacionais pressupõem espaçamento uniforme no espaço computacional (incrementos unitários). A malha física pode ser curvilínea, mas a qualidade da transformação e a ausência de jacobianos próximos de zero são essenciais.

### 2.2 MacCormack preditor–corretor

O avanço convectivo-viscoso é baseado no método de MacCormack. Primeiro, constrói os fluxos físicos nas direções \(x\) e \(y\), com derivadas direcionais progressivas no preditor; depois atualiza os fluxos no estado predito e usa derivadas regressivas no corretor.(Não consegui discrever melhor)

Em notação resumida:

$$$
U^*=U^n-\Delta t_{i,j}\left(D_x^+ E(U^n)+D_y^+F(U^n)\right),
$$

$$
U^{n+1}=\frac12\left[U^n+U^*-\Delta t_{i,j}\left(D_x^-E(U^*)+D_y^-F(U^*)\right)\right].
$$

Aqui \(U\) representa o vetor das seis variáveis conservativas e \(E,F\) os fluxos nas direções cartesianas físicas. O código também atualiza as fontes de turbulência e acrescenta dissipação artificial depois do corretor.

### 2.3 Estimativa do passo de tempo local

O passo local é calculado por uma estimativa convectiva–difusiva:

$$
\Delta t_{i,j}=\frac{CFL}{
\frac{|u|+a}{\Delta s_1}+\frac{|v|+a}{\Delta s_2}
+2\nu_e\left(\frac{1}{\Delta s_1^2}+\frac{1}{\Delta s_2^2}\right)},
$$

onde

$$
a=\sqrt{\gamma RT},\qquad \nu_e=\frac{\mu+\mu_t}{\rho}.
$$

\(\Delta s_1\) e \(\Delta s_2\) são escalas de comprimento estimadas a partir das distâncias entre pontos vizinhos nas duas direções computacionais. O método calcula um \(\Delta t_{i,j}\) por célula e a rotina retorna também o menor valor para diagnóstico.

### 2.4 Soft start e dissipação artificial

O modelo ajusta alguns fatores com o tempo 

Nas primeiras 3000 iterações, o programa varia linearmente:

- o CFL de 0,01 até o CFL solicitado;
- o coeficiente de dissipação artificial de 0,05 até o valor solicitado.

A dissipação é adicionada como um laplaciano discreto das variáveis conservativas. Ela pode reduzir oscilações, mas não equivale a um fluxo de face conservativo com limitador de positividade. Portanto, seu efeito sobre os balanços globais deve ser medido, e o coeficiente precisa ser testado por sensibilidade.

### 2.5 Positividade e limitações do avanço

Após cada predição/correção, o código aplica pisos para densidade, temperatura/energia interna, \(k\) e \(\varepsilon\), além de substituir velocidades não finitas por zero. Isso evita algumas falhas numéricas imediatas, mas modifica o estado quando os pisos são acionados. Quantidade de células corrigidas não é reportada; assim, uma simulação pode continuar sem deixar claro o quanto a solução foi alterada pelas proteções.

---

## 3. Estrutura do programa

| Componente | Função |
|---|---|
| `scalarField0` | Armazena coordenadas e calcula derivadas computacionais da malha. |
| `gridMetrics` | Calcula métricas e jacobiano da transformação curvilínea. |
| `scalarField` | Armazena campos escalares e calcula derivadas físicas a partir das métricas. |
| `flowSim` | Armazena as variáveis primitivas, conservativas, tensões e fluxos; calcula propriedades viscosas. |
| `boundaryConditionsClass` | Impõe as condições de campo distante e de parede codificadas no programa. |
| `computeTimeStep` | Calcula o passo local por célula a partir de escalas convectivas e difusivas. |
| `macCormack` | Executa as etapas preditora e corretora, inclui fontes turbulentas e calcula resíduos. |
| `addArtificialDissipation` | Adiciona dissipação artificial às seis variáveis conservativas. |
| `readGridFile` | Lê a malha estruturada de um arquivo de coordenadas. |
| `writeVTKfile` | Exporta a geometria e os campos para visualização. |
| `main` | Lê argumentos, inicializa os campos, aplica o soft start e executa as iterações. |

O vetor conservativo é

$$
U=[\rho,\rho u,\rho v,\rho E,\rho k,\rho\varepsilon]^T.
$$


---

## 4. Condições de contorno implementadas

As condições de contorno são fixas na classe `boundaryConditionsClass`, e não são selecionadas por um arquivo de configuração (mudarei isso assim que possível):

- **Primeira linha (`i = 0`):** campo distante. \(\rho,u,v,p,T\) são definidos pelos valores livres codificados; \(k\) e \(\varepsilon\) dependem do sentido estimado do escoamento normal.
- **Última linha (`i = L-1`):** parede estacionária sem deslizamento, com \(u=v=0\), temperatura imposta `T_wall` e pressão extrapolada da linha interna com piso de 10 Pa. Para turbulência ligada, \(\varepsilon\) na linha interna é recalculado com uma expressão algébrica de parede.
- **Colunas:** as diferenças e o passo de tempo usam índices periódicos nas colunas; o programa foi estruturado para uma malha tipo O com fechamento periódico nessa direção.

A função de parede estima

$$
u_\tau=C_\mu^{1/4}\sqrt{k},\qquad y^+=\frac{\rho u_\tau y}{\mu}.
$$

Para \(y^+>11{,}225\), usa uma lei logarítmica para a tensão de parede; abaixo desse limite, usa uma expressão linear viscosa. A função é simplificada e depende da malha: com turbulência ativa, a primeira célula deve ser compatível com o tratamento de parede de alto Reynolds, em geral com \(y^+\) na faixa de parede logarítmica recomendada pelo modelo e pela aplicação.


---

## 5. Malha de entrada

O programa recebe o arquivo de malha como primeiro argumento. O leitor procura um cabeçalho que contenha `I = ...` e `J = ...`; a partir da linha seguinte, espera pares numéricos `x y` até preencher a malha declarada.

---

## 6. Compilação

O código usa C++ e OpenMP. Em Linux com `g++` instalado:

```bash
g++ -O3 -march=native -fopenmp -std=c++17 solver_rans_ke_M03.cpp -o rans_m03
```

---

## 7. Execução e parâmetros

O programa não lê `parametros.txt` ou outro arquivo de configuração ainda, mas isso será feito no futuro. As condições físicas e várias propriedades estão fixas no código-fonte. Os argumentos opcionais da linha de comando são posicionais:

```bash
OMP_NUM_THREADS=4 ./rans_m03 cilindro.grid 500000 0.002 0.4 1
```

| Posição | Argumento | Padrão | Significado |
|---:|---|---:|---|
| 1 | `gridFileName` | `cilindro.grid` | Caminho para o arquivo de malha. |
| 2 | `maxIterations` | `500000` | Número máximo de iterações. |
| 3 | `dissipationCoeff` | `0.002` | Coeficiente final da dissipação artificial; começa em 0,05 e decai durante o soft start. |
| 4 | `cfl` | `0.4` | CFL final usado no passo local; começa em 0,01 e aumenta por 3000 iterações. |
| 5 | `useTurbulence` | `1` | `1` habilita o caminho turbulento k–ε; `0` desliga a viscosidade e as fontes turbulentas. |

Exemplo de execução curta para depuração:

```bash
OMP_NUM_THREADS=2 ./rans_m03 cilindro.grid 100 0.002 0.1 0
```

Os valores físicos definidos em `main` incluem:

| Propriedade | Valor |
|---|---:|
| \(T_0\) | 288,16 K |
| \(\mu_0\) | \(1,7894\times10^{-5}\) Pa·s |
| \(Pr\) | 0,71 |
| \(C_p\) | 1004 J/(kg·K) |
| \(C_v\) | 717 J/(kg·K) |
| \(R=C_p-C_v\) | 287 J/(kg·K) |
| \(p_\infty\) | 101325 Pa |
| \(\rho_\infty\) | 1,225 kg/m³ |
| \(u_\infty\) | 102 m/s |
| \(v_\infty\) | 0 m/s |
| Intensidade turbulenta inicial | 1% |
| Razão inicial \(\mu_t/\mu\) usada para estimar \(\varepsilon_\infty\) | 10 |
| Temperatura da parede | igual à temperatura livre calculada |

Para mudar essas condições, altere os valores em `main` e recompile.

---

## 8. Saída e visualização no ParaView

O código grava arquivos VTK no diretório de execução:

```text
escoamento_00000001.vtk
escoamento_00000200.vtk
escoamento_00000400.vtk
...
```

A primeira gravação ocorre na iteração 1; depois, a cada `writeEvery` iterações (padrão 200). Os arquivos contêm geometria `STRUCTURED_GRID`, dados em pontos (`POINT_DATA`) e os seguintes campos:

- `Pressure` — pressão;
- `Density` — densidade;
- `Temperature` — temperatura;
- `Mach` — número de Mach calculado pelo código;
- `TKE` — energia cinética turbulenta \(k\);
- `Epsilon` — dissipação \(\varepsilon\);
- `Velocity` — vetor velocidade \((u,v,0)\).


## 9. Monitoramento e convergência

A cada `reportEvery` iterações (padrão 100), o programa imprime o CFL, o coeficiente de dissipação, seis normas de resíduos e o tempo computacional. Os resíduos são baseados na diferença entre o estado corrigido e o estado anterior dividida pelo passo local.

A normalização usa, para cada equação, o maior resíduo observado até aquele momento. Portanto, o primeiro valor normalizado tende a ser 1 e a queda dessas razões não é, isoladamente, uma prova suficiente de convergência. O código não possui critério automático de parada por convergência; ele termina após `maxIterations` ou em caso de erro fatal.


---

## 10. Limitações conhecidas e próximas mudanças

1. **Não usa AUSM+.** O método espacial é MacCormack com dissipação artificial; não há divisão de Mach/pressão AUSM nem fluxo convectivo de Riemann por face.
2. **Variáveis e parametros separados** As variáveis e parâmetros atualmente estão no próprio código e a cada alteração é necessáro recompilar. O código será modificado para ler um arquivo .txt separado ou algo assim
3. **Tratamento de parede simplificado.** A parede está fixada na última linha e a função de parede pressupõe um modelo de alto Reynolds; a implementação não é um modelo de baixo Reynolds.
4. **Adequação às malhas.** O código ainda não suporta todos os tipo de malha sem precisar alterar o código
5. **Pisos numéricos** Densidade, energia interna, \(k\) e \(\varepsilon\) podem ser corrigidos após etapas; ainda não é reportada o niumero de correções feitas dessa forma.
6. **Passo local para regime permanente.** Não é um método apropriado, sem alterações, para interpretar todas as iterações como um transiente físico, então tentarei adaptar para suportar regime transiente também
7. **Incluir validação** 
8. **Angulo de ataque** Não existe atualmente uma opção de ângulo de ataque, pressão de saída ou temperatura de parede via linha de comando.Isso será adicionado
9. **Volumes finitos** Tentarei talvez adequar para volumes finitos caso encontre alguma vantagem clara para isso

## É necessário validar ainda

---
