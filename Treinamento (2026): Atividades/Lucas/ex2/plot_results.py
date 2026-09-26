import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.interpolate import griddata

plt.rcParams.update({
    "figure.dpi": 130,
    "font.size": 10,
    "axes.grid": True,
    "grid.alpha": 0.3,
})

R = 1.0

# 1) Campos 2D: Cp, u, v (numerico) sobre o dominio
#    Usamos a malha ESTRUTURADA (r,theta) -> pcolormesh, o que evita
#    artefatos de triangulacao (tricontourf conectaria pontos atraves
#    do "buraco" do cilindro de forma espuria).

field = pd.read_csv("field_reference.csv")

Nr = field["i"].max() + 1
Nt = field["j"].max() + 1

def to_grid(colname):
    piv = field.pivot(index="i", columns="j", values=colname)
    return piv.values  # shape (Nr, Nt)

r_grid = to_grid("r")[:, 0]           # (Nr,)
theta_grid = to_grid("theta")[0, :]   # (Nt,)
CP = to_grid("cp_num")
U = to_grid("u_num")
V = to_grid("v_num")

theta_closed = np.concatenate([theta_grid, [theta_grid[0] + 2*np.pi]])
CP_c = np.concatenate([CP, CP[:, [0]]], axis=1)
U_c = np.concatenate([U, U[:, [0]]], axis=1)
V_c = np.concatenate([V, V[:, [0]]], axis=1)

RR, TT = np.meshgrid(r_grid, theta_closed, indexing="ij")
X = RR * np.cos(TT)
Y = RR * np.sin(TT)

def make_pcolor(ax, Z, title, cmap="RdBu_r"):
    pc = ax.pcolormesh(X, Y, Z, shading="auto", cmap=cmap)
    ax.add_patch(Circle((0, 0), R, facecolor="white", edgecolor="k", zorder=5, linewidth=1.2))
    ax.set_aspect("equal")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.set_xlabel("x/R"); ax.set_ylabel("y/R")
    ax.set_title(title)
    return pc

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
pc0 = make_pcolor(axes[0], CP_c, "Coeficiente de pressão $C_p$ (numérico)")
plt.colorbar(pc0, ax=axes[0], shrink=0.85)
pc1 = make_pcolor(axes[1], U_c, "Velocidade $u$ (numérico)", cmap="coolwarm")
plt.colorbar(pc1, ax=axes[1], shrink=0.85)
pc2 = make_pcolor(axes[2], V_c, "Velocidade $v$ (numérico)", cmap="coolwarm")
plt.colorbar(pc2, ax=axes[2], shrink=0.85)
fig.suptitle("Escoamento potencial em torno de um cilindro — malha 120×240", y=1.03)
fig.tight_layout()
fig.savefig("fig1_campos_P_u_v.png", bbox_inches="tight")
plt.close(fig)

# 2) Comparação numérico x analítico ao longo de uma linha radial (theta=0)
#    e campo de vetores de velocidade
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# 2a) vetores de velocidade (subamostrados) sobre contorno de Cp
ax = axes[0]
pc = ax.pcolormesh(X, Y, CP_c, shading="auto", cmap="RdBu_r")
plt.colorbar(pc, ax=ax, shrink=0.85, label="$C_p$")
x_flat, y_flat = X[:, :-1].ravel(), Y[:, :-1].ravel()
u_flat, v_flat = U.ravel(), V.ravel()
xi = np.linspace(-6, 6, 34)
yi = np.linspace(-6, 6, 34)
XI, YI = np.meshgrid(xi, yi)
Ui = griddata((x_flat, y_flat), u_flat, (XI, YI), method="cubic")
Vi = griddata((x_flat, y_flat), v_flat, (XI, YI), method="cubic")
inside = (XI**2 + YI**2) < (1.15 * R)**2
Ui[inside] = np.nan; Vi[inside] = np.nan
ax.quiver(XI, YI, Ui, Vi, color="k", scale=25, width=0.003, alpha=0.8)
ax.add_patch(Circle((0, 0), R, facecolor="white", edgecolor="k", zorder=5, linewidth=1.2))
ax.set_aspect("equal"); ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
ax.set_xlabel("x/R"); ax.set_ylabel("y/R")
ax.set_title("Campo de velocidade (vetores) sobre $C_p$")

# 2b) Cp na superficie do cilindro: numerico x analitico
surf = pd.read_csv("surface_reference.csv")
surf = surf.sort_values("theta")
ax = axes[1]
theta_deg = np.degrees(surf["theta"])
ax.plot(theta_deg, surf["cp_exact"], "-", color="k", lw=2, label="Analítico: $C_p=1-4\\sin^2\\theta$")
ax.plot(theta_deg, surf["cp_num"], "o", color="tab:red", ms=3.5, label="Numérico (malha 120×240)")
ax.set_xlabel(r"$\theta$ (graus)")
ax.set_ylabel("$C_p$")
ax.set_title("Distribuição de $C_p$ na superfície do cilindro")
ax.legend()
ax.set_xlim(0, 360)

fig.tight_layout()
fig.savefig("fig2_vetores_e_Cp_superficie.png", bbox_inches="tight")
plt.close(fig)

# 3) Perfis radiais de u e v ao longo de theta = 0 e theta = 90
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for ax, theta_target, label in zip(axes, [0.0, np.pi/2], [r"$\theta=0°$ (linha de estagnação)", r"$\theta=90°$ (topo do cilindro)"]):
    sub = field[np.isclose(field["theta"], theta_target, atol=1e-6)].sort_values("r")
    if len(sub) == 0:
        # pega o mais proximo disponivel na malha
        thetas = field["theta"].unique()
        closest = thetas[np.argmin(np.abs(thetas - theta_target))]
        sub = field[field["theta"] == closest].sort_values("r")
    ax.plot(sub["r"], sub["u_exact"], "-", color="tab:blue", lw=2, label="$u$ analítico")
    ax.plot(sub["r"], sub["u_num"], "o", color="tab:blue", ms=3, mfc="none", label="$u$ numérico")
    ax.plot(sub["r"], sub["v_exact"], "-", color="tab:orange", lw=2, label="$v$ analítico")
    ax.plot(sub["r"], sub["v_num"], "s", color="tab:orange", ms=3, mfc="none", label="$v$ numérico")
    ax.set_xlabel("r/R")
    ax.set_ylabel("Velocidade")
    ax.set_title(label)
    ax.legend(fontsize=8)
    ax.set_xlim(1, 6)

fig.tight_layout()
fig.savefig("fig3_perfis_radiais_u_v.png", bbox_inches="tight")
plt.close(fig)

# 4) Estudo de convergência de malha (erro vs. h, log-log)
conv = pd.read_csv("convergence.csv")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.loglog(conv["h"], conv["L2_phi"], "o-", label=r"$L_2(\phi)$")
ax.loglog(conv["h"], conv["Linf_phi"], "s-", label=r"$L_\infty(\phi)$")
# linha de referencia de segunda ordem
h_ref = conv["h"].values
C = conv["L2_phi"].values[0] / h_ref[0]**2
ax.loglog(h_ref, C * h_ref**2, "k--", lw=1, label="Referência $O(h^2)$")
ax.set_xlabel("h = $\\Delta r$ (espaçamento radial)")
ax.set_ylabel("Erro")
ax.set_title("Convergência do potencial $\\phi$")
ax.legend()

ax = axes[1]
ax.loglog(conv["h"], conv["L2_cp_surf"], "o-", color="tab:red", label=r"$L_2(C_p)$ na superfície")
ax.loglog(conv["h"], conv["Linf_cp_surf"], "s-", color="tab:purple", label=r"$L_\infty(C_p)$ na superfície")
C2 = conv["L2_cp_surf"].values[0] / h_ref[0]**2
ax.loglog(h_ref, C2 * h_ref**2, "k--", lw=1, label="Referência $O(h^2)$")
ax.set_xlabel("h = $\\Delta r$ (espaçamento radial)")
ax.set_ylabel("Erro")
ax.set_title("Convergência de $C_p$ na superfície")
ax.legend()

fig.tight_layout()
fig.savefig("fig4_convergencia_malha.png", bbox_inches="tight")
plt.close(fig)

# tabela resumo em texto
print(conv.to_string(index=False))

# ordem de convergencia observada 
print("\nOrdem de convergência observada (L2 phi):")
for k in range(1, len(conv)):
    h1, h0 = conv["h"].iloc[k], conv["h"].iloc[k-1]
    e1, e0 = conv["L2_phi"].iloc[k], conv["L2_phi"].iloc[k-1]
    p = np.log(e1/e0) / np.log(h1/h0)
    print(f"  {conv['Nr'].iloc[k-1]:>4}->{conv['Nr'].iloc[k]:<4}  ordem p = {p:.2f}")

print("\nFiguras geradas: fig1_campos_P_u_v.png, fig2_vetores_e_Cp_superficie.png, "
      "fig3_perfis_radiais_u_v.png, fig4_convergencia_malha.png")
