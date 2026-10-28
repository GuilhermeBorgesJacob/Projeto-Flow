$$
\frac{\partial \rho}{\partial t}  +\frac{\partial( \rho u)}{\partial x} +\frac{\partial( \rho v)}{\partial y} = 0
$$
$$
\frac{\partial (\rho u )}{\partial t} + \frac{\partial (\rho u^2 )}{\partial x} + \frac{\partial (\rho uv )}{\partial y}= -\frac{\partial p}{\partial x} + \frac{\partial \tau_{xx}}{\partial x} + \frac{\partial \tau_{yx}}{\partial y} + \rho f_x
$$
$$
\frac{\partial (\rho v )}{\partial t} + \frac{\partial (\rho uv )}{\partial x} + \frac{\partial (\rho v^2 )}{\partial y}= -\frac{\partial p}{\partial y} + \frac{\partial \tau_{xy}}{\partial x} + \frac{\partial \tau_{yy}}{\partial y} + \rho f_y
$$
$$
\frac{\partial E_t }{\partial t} + \frac{\partial (E_t u )}{\partial x} + \frac{\partial (E_t v)}{\partial y}= \rho \dot{q} + \frac{\partial}{\partial x} \bigg( k \frac{\partial T}{\partial x} \bigg) + \frac{\partial}{\partial y} \bigg( k \frac{\partial T}{\partial y} \bigg) - \frac{\partial (up)}{\partial x} - \frac{\partial(vp)}{\partial y} + \frac{\partial( u \tau_{xx})}{\partial x} + \frac{\partial( u \tau_{yx})}{\partial y} + \frac{\partial( v \tau_{xy})}{\partial x} + \frac{\partial( v \tau_{yy})}{\partial y} + \rho f \cdot V
$$
$$
\frac{\partial (\rho u )}{\partial t} + \frac{\partial (\rho u^2 +p-\tau_{xx})}{\partial x} + \frac{\partial (\rho uv - \tau_{yx})}{\partial y}=  \rho f_x
$$
$$
\frac{\partial (\rho v )}{\partial t} + \frac{\partial (\rho uv -\tau_{xy})}{\partial x} + \frac{\partial (\rho v^2 +p - \tau_{yy})}{\partial y}=  \rho f_y
$$
$$
\frac{\partial ( E_t)}{\partial t} + \frac{\partial ((E_t+p)u+q_x - u\tau_{xx} - v\tau_{xy})}{\partial x} + \frac{\partial ((E_t+p)v+q_y - u\tau_{xy} - v\tau_{yy})}{\partial y}= 0 
$$
$$
p = \rho R T
$$
$$
e = e(T,p)
$$
$$
e = c_v T
$$
$$
e = e(p,\rho)
$$
$$
e = \frac{1}{\gamma-1}\frac{p}{\rho}
$$
$$
\tau_{xy} = \tau_{yx} = \mu \bigg(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \bigg)
$$
$$
\tau_{xx} = \lambda(\nabla \cdot V) +2\mu \frac{\partial u}{\partial x}
$$
$$
\tau_{yy} = \lambda(\nabla \cdot V) +2\mu \frac{\partial v}{\partial y}
$$
$$
\lambda = -\frac{2}{3}\mu
$$
$$
\mu = \mu_0 \bigg(\frac{T}{T_0}\bigg)^{3/2}\frac{T_0+110}{T+110}
$$
$$
Pr =  0.71 = \frac{\mu c_p}{k}
$$
$$
\frac{\partial (\rho e)}{\partial t}+\frac{\partial(\rho e u)}{\partial x} + \frac{\partial (\rho ev)}{\partial y} = \rho \dot{q} + k\bigg(\frac{\partial^2 T}{\partial x^2} +\frac{\partial^2 T}{\partial y^2}\bigg) + \frac{\partial k}{\partial x}\frac{\partial T}{\partial x} +\frac{\partial k}{\partial y}\frac{\partial T}{\partial y} -p(\nabla \cdot V) + \lambda (\nabla \cdot V)^2 + \mu \bigg[2\bigg(\frac{\partial u}{\partial x}\bigg)^2 + 2\bigg(\frac{\partial v}{\partial y}\bigg)^2 + \bigg( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \bigg)^2\bigg]
$$
