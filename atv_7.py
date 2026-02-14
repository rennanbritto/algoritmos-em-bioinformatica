from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# função para ler um arquivo FASTA e retornar uma lista de sequências
def ler_fasta(caminho_arquivo):
    sequencias = []
    for registro in SeqIO.parse(caminho_arquivo, "fasta"):
        sequencias.append(str(registro.seq)) #converte o objeto seq em uma string
    return sequencias


# função para encontrar a maior sobreposição entre duas sequências
def encontrar_sobreposicao(seq1, seq2, sobreposicao_min=3):
    sobreposicao_max = min(len(seq1), len(seq2))
    for i in range(sobreposicao_max, sobreposicao_min - 1, -1):
        if seq1[-i:] == seq2[:i]: # Compara o final da string com o início da outra
            return i
    return 0


# função para montar o contig a partir das sequências fornecidas
def montar_contig(sequencias):
    contig_montado = sequencias.pop(0) #retira a primeira sequência da lista
    while sequencias:
        encontrou_sobreposicao = False
        for i, seq in enumerate(sequencias):
            sobreposicao = encontrar_sobreposicao(contig_montado, seq)
            if sobreposicao > 0:
                contig_montado += seq[sobreposicao:] #junta a sobreposição no contig
                sequencias.pop(i)
                encontrou_sobreposicao = True
                break
        if not encontrou_sobreposicao:
            print("Nenhuma sobreposição encontrada. Saindo do loop.")
            break
    return contig_montado


# função para escrever a sequência do contig em um arquivo FASTA
def escrever_fasta(sequencia, caminho_arquivo, nome_contig="contig"):
    registro = SeqRecord(Seq(sequencia), id=nome_contig, description="")
    SeqIO.write(registro, caminho_arquivo, "fasta")


# função principal
if __name__ == "__main__":
    arquivo_entrada = "reads4.fasta"
    arquivo_saida = "contig.fasta"

    # lendo as sequências do arquivo FASTA de entrada
    sequencias = ler_fasta(arquivo_entrada)

    if sequencias:
        # montando o contig a partir das sequências lidas
        contig = montar_contig(sequencias)

        # escrevendo o contig montado em um arquivo FASTA de saída
        escrever_fasta(contig, arquivo_saida)

        print(f"Contig montado e salvo em {arquivo_saida}")
    else:
        print("Nenhuma sequência encontrada no arquivo.")
