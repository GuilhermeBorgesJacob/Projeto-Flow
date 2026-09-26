import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

dados = pd.read_csv('../evolucao_completa.csv')

x = dados['x'].values

colunas_t = [col for col in dados.columns if col != 'x']
tempos = [float(col.split('_')[1].replace('s', '')) for col in colunas_t]

fig, ax = plt.subplots(figsize=(12, 8))

cores = plt.cm.viridis(np.linspace(0, 1, len(colunas_t)))

for i, col in enumerate(colunas_t):
    ax.plot(x, dados[col], color=cores[i],
            label=f't = {tempos[i]:.2f}s', linewidth=2)

ax.set_xlabel('Posição x (m)', fontsize=12)
ax.set_ylabel('Temperatura T (°C)', fontsize=12)
ax.set_title('Evolução da Temperatura na Barra', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_ylim([-5, 105])

plt.tight_layout()
plt.savefig('evolucao_temperatura.png', dpi=150, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(figsize=(12, 8))

# Criar matriz de temperaturas
T_matrix = dados[colunas_t].values.T

# Plotar mapa de cores
im = ax.imshow(T_matrix, aspect='auto', origin='lower',
               extent=[x.min(), x.max(), min(tempos), max(tempos)],
               cmap='hot', interpolation='bilinear')

ax.set_xlabel('Posição x (m)', fontsize=12)
ax.set_ylabel('Tempo t (s)', fontsize=12)
ax.set_title('Mapa de Calor', fontsize=14)

cbar = plt.colorbar(im, ax=ax, label='Temperatura (°C)')

plt.tight_layout()
plt.savefig('mapa_calor.png', dpi=150, bbox_inches='tight')
plt.show()

perfil_final = pd.read_csv('../perfil_final.csv')

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(perfil_final['x'], perfil_final['T'], 'b-', linewidth=2, label='Solução Numérica')

ax.set_xlabel('Posição x (m)', fontsize=12)
ax.set_ylabel('Temperatura T (°C)', fontsize=12)
ax.set_title('Perfil Final de Temperatura', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig('perfil_final.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Temperatura final no centro: {perfil_final['T'].iloc[len(x)//2]:.2f} °C")
