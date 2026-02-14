from Bio import pairwise2
from Bio.Align import substitution_matrices

# Sequências dadas pelo professor
seq1 = ("GGGCAATATGAA---TTTAAA---GTAGAATACCAAATGATAGAAACAGACTGCCTGA-TTGATCATTTT"
        "GATTTTTTAAAGTGTGTA------TAAATTGCTGTTCCTTAATTTGATTA")
seq2 = ("GGGCAATATGAAATTTTTAAAGGAGTAGAATACGAAATGATAGACACAGACTGCCTGAATTGAGGATTTT"
        "GATTTCTTAAATTGTGTTTCTTTCTAAATTGCTGTTCCTTAATTTTATTA")

# Definição das matrizes de identidade e de escore
matriz_identidade = {
    'A': {'A': 1, 'T': 0, 'G': 0, 'C': 0, '-': 0},
    'T': {'A': 0, 'T': 1, 'G': 0, 'C': 0, '-': 0},
    'G': {'A': 0, 'T': 0, 'G': 1, 'C': 0, '-': 0},
    'C': {'A': 0, 'T': 0, 'G': 0, 'C': 1, '-': 0},
    '-': {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
}

matriz_escore = {
    ('A', 'A'): 5, ('A', 'C'): -1, ('A', 'G'): -2, ('A', 'T'): -1, ('A', '-'): -3,
    ('C', 'A'): -1, ('C', 'C'): 5, ('C', 'G'): -3, ('C', 'T'): -2, ('C', '-'): -4,
    ('G', 'A'): -2, ('G', 'C'): -3, ('G', 'G'): 5, ('G', 'T'): -2, ('G', '-'): -2,
    ('T', 'A'): -1, ('T', 'C'): -2, ('T', 'G'): -2, ('T', 'T'): 5, ('T', '-'): -1,
    ('-', 'A'): -3, ('-', 'C'): -4, ('-', 'G'): -2, ('-', 'T'): -1, ('-', '-'): 0
}

# Função para calcular o escore de identidade e escore total
def calcular_escores(seq1, seq2, matriz):
    escore_total = 0
    escores_individuais = []
    for a, b in zip(seq1, seq2):
        escore = matriz.get(a, {}).get(b, -1)
        escore_total += escore
        escores_individuais.append(escore)
    return escore_total, escores_individuais

# Cálculo dos escores
escore_total_identidade, escores_identidade = calcular_escores(seq1, seq2, matriz_identidade)
escore_total, escores = calcular_escores(seq1, seq2, matriz_escore)

# Cálculo porcentagem
percentual_identidade = (escores_identidade.count(1) / len(escores_identidade)) * 100

print(f"Escore total (Identidade): {escore_total_identidade}")
print(f"Escore total (Matriz de escores): {escore_total}")
print(f"Identidade (%): {percentual_identidade:.2f}%")

# Definição da matriz de substituição personalizada para pairwise2
substitution_matrix = substitution_matrices.Array(alphabet="ACGT-", dims=2)
for (i, char1) in enumerate("ACGT-"):
    for (j, char2) in enumerate("ACGT-"):
        substitution_matrix[i, j] = matriz_escore.get((char1, char2), -1)

# Configuração do alinhamento
alignments = pairwise2.align.globaldx(seq1, seq2, substitution_matrix, open_gap_score=-2, extend_gap_score=-0.5, one_alignment_only=True)

# Imprime os alinhamentos e escores
for alignment in alignments:
    print("\nAlinhamento:")
    print(pairwise2.format_alignment(*alignment))
    print(f"Escore do alinhamento: {alignment[2]}")

# Cálculo do escore médio
escore_medio = escore_total / len(seq1)

print(f"Escore médio: {escore_medio:.2f}")
