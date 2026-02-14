def mochila(valores, pesos, Wmax):
    n = len(valores)
    dp = [[0 for x in range(Wmax + 1)] for x in range(n + 1)]

    # Preenchendo a matriz dp de baixo para cima
    for i in range(1, n + 1):
        for w in range(Wmax + 1):
            if pesos[i - 1] <= w:
                dp[i][w] = max(valores[i - 1] + dp[i - 1][w - pesos[i - 1]],
                               dp[i - 1][w])
            else:
                dp[i][w] = dp[i - 1][w]

    # Mostrando a matriz dp
    print("Matriz de programação dinâmica (dp):")
    for linha in dp:
        print(linha)

    # Identificando quais itens foram colocados na mochila
    itens_selecionados = []
    w = Wmax
    for i in range(n, 0, -1):
        if dp[i][w] != dp[i - 1][w]:
            itens_selecionados.append(
                i)  # Armazenando o índice do item (1-based)
            w -= pesos[i - 1]

    itens_selecionados.reverse()  # Para manter a ordem original dos itens

    return dp[n][Wmax], itens_selecionados


valores = [5, 4, 7, 7]
pesos = [5, 6, 8, 4]
Wmax = 13

valor_maximo, itens_selecionados = mochila(valores, pesos, Wmax)

# Exibindo o resultado
print(f"\nValor máximo colocado na mochila: R$ {valor_maximo}")
print("Itens a serem colocados na mochila:", itens_selecionados)

# Escrevendo o resultado em um arquivo de texto
with open("resultado.txt", "w") as arquivo:
    arquivo.write(f"Valor máximo colocado na mochila: R$ {valor_maximo}\n")
    arquivo.write("Itens a serem colocados na mochila:\n")
    for item in itens_selecionados:
        arquivo.write(
            f" - Item {item} (Valor: R$ {valores[item - 1]}, Peso: {pesos[item - 1]} kg)\n"
        )
