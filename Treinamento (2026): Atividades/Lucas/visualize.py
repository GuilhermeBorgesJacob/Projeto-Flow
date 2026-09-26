import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Lê o arquivo gerado

dados = pd.read_csv("temperaturas.csv")

# Primeira coluna = tempo
tempo = dados.iloc[:, 0].values

# outras colunas = temperaturas
temperaturas = dados.iloc[:, 1:].values

numero_nos = temperaturas.shape[1]

posicao = np.linspace(0, 2.0, numero_nos)

# Cria o mapa de calor

plt.figure(figsize=(12,6))

imagem = plt.imshow(
    temperaturas,
    aspect='auto',
    origin='lower',
    cmap='hot',
    extent=[
        posicao[0],
        posicao[-1],
        tempo[0],
        tempo[-1]
    ]
)

plt.colorbar(imagem, label="Temperatura (°C)")

plt.xlabel("Posição na barra (m)")
plt.ylabel("Tempo (s)")
plt.title("Condução de calor unidimensional")

plt.tight_layout()

plt.show()
