def alinhamento_local_com_matriz(seq1, seq2, match=5, mismatch=-3, gap=-4):
    # Inicializar a matriz com zeros
    n = len(seq1)
    m = len(seq2)
    matriz = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    melhor_escore = 0
    melhor_posicao = (0, 0)

    # Preencher a matriz
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                escore = match
            else:
                escore = mismatch

            # Calcular o valor para cada célula da matriz
            matriz[i][j] = max(
                0,
                matriz[i - 1][j - 1] + escore,  # Diagonal (match/mismatch)
                matriz[i - 1][j] + gap,  # De cima (gap)
                matriz[i][j - 1] + gap)  # Da esquerda (gap)

            # O melhor escore
            if matriz[i][j] > melhor_escore:
                melhor_escore = matriz[i][j]
                melhor_posicao = (i, j)

    # Exibir a matriz de escores
    print("Matriz de Escores:")
    for linha in matriz:
        print(" ".join(f"{v:3}" for v in linha))

    # Reconstruir o alinhamento a partir da melhor posição
    alinhamento1, alinhamento2, alinhamento_simb = "", "", ""
    i, j = melhor_posicao

    while matriz[i][j] != 0:
        if seq1[i - 1] == seq2[j - 1]:
            alinhamento1 = seq1[i - 1] + alinhamento1
            alinhamento2 = seq2[j - 1] + alinhamento2
            alinhamento_simb = "|" + alinhamento_simb  # Match
            i -= 1
            j -= 1
        elif matriz[i - 1][j] + gap == matriz[i][j]:
            alinhamento1 = seq1[i - 1] + alinhamento1
            alinhamento2 = "-" + alinhamento2
            alinhamento_simb = " " + alinhamento_simb  # Gap
            i -= 1
        else:
            alinhamento1 = "-" + alinhamento1
            alinhamento2 = seq2[j - 1] + alinhamento2
            alinhamento_simb = " " + alinhamento_simb  # Gap
            j -= 1

    # Mostrar as respostas
    print("\nAlinhamento 1: ", alinhamento1)
    print("Alinhamento   : ", alinhamento_simb)
    print("Alinhamento 2: ", alinhamento2)
    print("Escore total: ", melhor_escore)

    # MOstrar os escores
    escores = ""
    for k in range(len(alinhamento1)):
        if alinhamento1[k] == alinhamento2[k]:
            escores += f"{match} "
        elif alinhamento1[k] == "-" or alinhamento2[k] == "-":
            escores += f"{gap} "
        else:
            escores += f"{mismatch} "

    print("Escores      : ", escores)


# Sequencas da atividade
seq1 = "GAATTCAGTTA"
seq2 = "GGATCGA"

alinhamento_local_com_matriz(seq1, seq2)
