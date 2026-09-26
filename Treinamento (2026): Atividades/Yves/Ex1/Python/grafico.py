import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VALOR_EXATO = 0.00709227

def main():
    if not os.path.exists('../resultados.csv'):
        print("Arquivo resultados.csv não encontrado!")
        print("Execute o código em C++ primeiro.")
        return

    dados = pd.read_csv('../resultados.csv')

    erros = (1 - dados['Tau_w (Pa)'] / VALOR_EXATO) * 100

    # ==================== GRAFICO 1 ====================
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))
    fig1.set_facecolor('white')
    fig1.suptitle('Análise de Diferenças Finitas', fontsize=14, fontweight='bold')

    # Resultados τ_w
    cores = ['#FF6B6B', '#4ECDC4', '#45B7D1']
    bars = ax1.bar(dados['Ordem'], dados['Tau_w (Pa)'], color=cores)
    ax1.axhline(y=VALOR_EXATO, color='red', linestyle='--', label='Exato')
    ax1.set_title('Resultados')
    ax1.set_ylabel('τ_w (Pa)')
    ax1.set_xlabel('Ordem')
    ax1.legend()
    ax1.set_facecolor('white')
    ax1.grid(True, alpha=0.3)

    # Erro
    colors = ['green' if abs(e) < 5 else 'orange' for e in dados["Erro (%)"]]
    ax2.bar(dados['Ordem'], dados['Erro (%)'], color=colors)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_title('Erros')
    ax2.set_ylabel('Erro (%)')
    ax2.set_xlabel('Ordem')
    ax2.set_facecolor('white')
    ax2.grid(True, alpha=0.3)

    # Relação de precisão
    ax3.scatter(dados['Gradiente (du/dy)'], dados['Tau_w (Pa)'], s=100, color='#2E86AB')
    ax3.set_xlabel('Gradiente (du/dy)')
    ax3.set_ylabel('τ_w')
    ax3.set_title('Relação de precisão')
    ax3.set_facecolor('white')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=VALOR_EXATO, color='red', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig('../analise_resultados.png', dpi=150, facecolor='white', transparent=False)
    plt.close()

    # ==================== GRÁFICO 2 ====================
    fig2, ax = plt.subplots(figsize=(10, 6))
    fig2.set_facecolor('white')
    ax.set_facecolor('white')

    ordens = [1, 2, 3]
    ax.plot(ordens, erros, 'o-', linewidth=2, markersize=10, color='#2E86AB')
    ax.fill_between(ordens, 0, erros, alpha=0.2, color='#2E86AB')

    ax.set_xlabel('Ordem da Diferença Finita', fontsize=12)
    ax.set_ylabel('Erro (%)', fontsize=12)
    ax.set_title('Convergência do Método', fontsize=14, fontweight='bold')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(ordens)

    for x, y in zip(ordens, dados['Erro (%)']):
        ax.annotate(f'{y:.4f}%', (x, y), xytext=(5, 5),
                   textcoords='offset points', fontsize=10)

    plt.tight_layout()
    plt.savefig('../convergencia_metodo.png', dpi=150, facecolor='white', transparent=False)
    plt.close()

if __name__ == "__main__":
    main()
