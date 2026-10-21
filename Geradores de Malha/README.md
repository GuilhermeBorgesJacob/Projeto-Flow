# Geradores de malha automáticos por equações diferencias
A principal vantagem deste método de geração de malha consiste na sua generalidade, podendo ser aplicado a partir do mesmo algorítimo para diversas geometrias.
Sua principal desvantem se dá pela dificuldade do controle da qualidade da malha perto de regiões "muito" concavas ou convexas, além disso, estes métodos são baseados em geometrias de malhas estruturadas em que por si só trazem limitações em geometrias mais complexas.

## Equações elípticas

As equações com comportamneto elíptico utilizadas para gerar as malhas são: 

$$
\frac{\partial^2 x}{\partial \xi^2} + \frac{\partial^2 x}{\partial \eta^2} = 0
$$

$$
\frac{\partial^2 y}{\partial \xi^2} + \frac{\partial^2 y}{\partial \eta^2} = 0
$$

Onde as coordenadas $\xi$ e $\eta$ são coordenadas computacionais.

### Solução numérica explicita

Aplicando diferenças finitas centrais em torno do ponto com índice (i,j) temos:

$$
\frac{x_{i+1,j} - 2x_{i,j} + x_{i-1,j}}{\Delta\xi^2} + \frac{x_{i,j+1} - 2x_{i,j} + x_{i,j-1}}{\Delta\eta^2} = 0
$$

$$
\frac{y_{i+1,j} - 2y_{i,j} + y_{i-1,j}}{\Delta\xi^2} + \frac{y_{i,j+1} - 2y_{i,j} + y_{i,j-1}}{\Delta\eta^2} = 0
$$

Como os valores de $\Delta\xi^2$ e $\Delta\eta^2$ são arbitrários assume-se que são iguais a um. Assim, isolando $x_{i,j}$ e $y_{i,j}$ as seguintes expressões são obtidas:

$$
x_{i,j} = \frac{x_{i+1,j} + x_{i-1,j} + x_{i,j+1} + x_{i,j-1}}{4}
$$

$$
y_{i,j} = \frac{y_{i+1,j} + y_{i-1,j} + y_{i,j+1} + y_{i,j-1}}{4}
$$

Quando aplicadas as equações acima para cada ponto do domínio um sistema linear de equações é formado. Para solucionar o sistema condições de contorno devem ser aplicadas em todo limite do domínio computacional. **É importante notar que o comportamento do sistema discreto condiz com o comportamento do sistema de equações originais, onde ambos necessitam de aplicação de condições de contorno em toda borda do domínio e ambos são lineares, logo a discretização escolhida é razoável**.
As principais técnicas numéricas utilizadas para a solução do sistema de equações descritas acima são a método de gauss-seidel e a técnica de relaxação.
