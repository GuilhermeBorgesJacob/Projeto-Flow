#!/usr/bin/env python3
"""Figura Mach + linhas de corrente + Cp a partir de um .vts do solver.  Uso: python3 plot_vaca.py campo_final.vts saida.png"""
import sys, numpy as np, xml.etree.ElementTree as ET
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.interpolate import griddata
def ler(path):
    pc = ET.parse(path).getroot().find(".//Piece"); e = [int(v) for v in pc.get("Extent").split()]; NI, NJ = e[1] + 1, e[3] + 1
    P = np.array(pc.find("Points/DataArray").text.split(), float).reshape(NJ, NI, 3)
    d = {}
    for a in pc.find("PointData"):
        n = int(a.get("NumberOfComponents")); v = np.array(a.text.split(), float).reshape(NJ, NI, n); d[a.get("Name")] = v[..., 0] if n == 1 else v
    return P[..., 0], P[..., 1], d
X, Y, d = ler(sys.argv[1])
M, Cp, V = d["Mach"], d["Cp"], d["Velocity"]
body_x, body_y = X[0], Y[0]
fig, ax = plt.subplots(2, 1, figsize=(9.5, 10))
xl, yl = (-2.2, 5.0), (-2.4, 2.4)
for a, f, nome, cmap, lim in ((ax[0], M, "Numero de Mach", "turbo", (0, 0.8)), (ax[1], Cp, "Coeficiente de pressao Cp", "RdBu_r", (-1.5, 1.2))):
    pm = a.pcolormesh(X, Y, f, shading="gouraud", cmap=cmap, vmin=lim[0], vmax=lim[1], rasterized=True)
    a.fill(body_x, body_y, color="0.2", zorder=5); a.set_xlim(*xl); a.set_ylim(*yl); a.set_aspect("equal"); a.set_title(nome)
    fig.colorbar(pm, ax=a, fraction=0.025, pad=0.01)
# linhas de corrente (interpola em grade regular; mascara o corpo)
gx, gy = np.linspace(*xl, 400), np.linspace(*yl, 300); GX, GY = np.meshgrid(gx, gy)
pts = np.stack([X.ravel(), Y.ravel()], 1)
U_ = griddata(pts, V[..., 0].ravel(), (GX, GY), method="linear"); V_ = griddata(pts, V[..., 1].ravel(), (GX, GY), method="linear")
from matplotlib.path import Path
inside = Path(np.stack([body_x, body_y], 1)).contains_points(np.stack([GX.ravel(), GY.ravel()], 1)).reshape(GX.shape)
U_[inside] = np.nan; V_[inside] = np.nan
ax[0].streamplot(gx, gy, U_, V_, color="w", density=2.2, linewidth=0.5, arrowsize=0.6)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=100); print("salvo", sys.argv[2])
