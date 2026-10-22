# Escoamento incompressível

As equações que descrevem este regime são a equação da continuade

$$
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0
$$

e as equações de conservação da quantidade de movimento:

$$
\frac{\partial u}{\partial t} + \frac{\partial u^2}{\partial x} + \frac{\partial uv}{\partial y} = - \frac{1}{\rho}\frac{\partial p}{\partial x} + \nu \bigg(  \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\bigg)
$$

$$
\frac{\partial v}{\partial t} + \frac{\partial uv}{\partial x} + \frac{\partial v^2}{\partial y} = - \frac{1}{\rho}\frac{\partial p}{\partial y} + \nu \bigg(  \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2}\bigg)
$$

Na forma vetorial:

$$
\vec{\nabla} \cdot \vec{V} = 0
$$

$$
\frac{\partial \vec{V}}{\partial t} + ( \vec{V} \cdot \vec{\nabla}) \vec{V} = - \frac{1}{\rho} \vec{\nabla} p + \nu \nabla^2 \vec{V}
$$

Três equações são suficientes para descrever o escoamento, pois apenas o campo de velocidades e pressão são incógnitas.

## A técnica de correção de pressão

A técnica consiste em um procedimento iterativo entre os campos de velocidade e de pressão. A partir de uma aproximação inicial para a pressão, a equação da quantidade de movimento pode ser resolvida para determinar o campo de velocidade. O campo de velocidade obtido não satisfaz a equação da continuidade, que exige divergência nula e, portanto, precisa ser corrigido. Como essa correção afeta o campo de pressão, define-se uma correção de pressão, obtida ao impor que a velocidade corrigida satisfaça a equação da continuidade. Isso resulta em uma equação de Poisson para a correção da pressão (HIRSH 1888).



