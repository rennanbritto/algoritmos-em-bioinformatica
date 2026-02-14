from Bio import SeqIO
from Bio.Seq import Seq

# Dicionário de RNA para aminoácidos
rna_para_aa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


# Função para traduzir RNA para aminoácidos
def traduzir_rna_para_proteina(rna_seq):
  seq_proteina = ""
  for i in range(0, len(rna_seq), 3):
    codon = rna_seq[i:i + 3]
    if len(codon) == 3:
      seq_proteina += rna_para_aa.get(codon, "X")
  return seq_proteina


# Abrir arquivo FASTA com as sequências de DNA
arquivo_entrada = "dengue.fasta"
arquivo_saida_rna = "sequencias_rna.fasta"

with open(arquivo_saida_rna, "w") as saida_rna:
  # Ler sequências de DNA e transcrever para RNA
  for sequencia in SeqIO.parse(arquivo_entrada, "fasta"):
    seq_dna = sequencia.seq
    seq_rna = seq_dna.transcribe()

    # Escrever no arquivo de saída
    saida_rna.write(f">Sequencia de DNA: {sequencia.id}\n")
    saida_rna.write(f"{seq_dna}\n")
    saida_rna.write(f">Sequencia de RNA: {sequencia.id}\n")
    saida_rna.write(f"{seq_rna}\n")

    # Traduzir RNA para proteína em 6 frames
    for frame in range(3):
      seq_proteina = traduzir_rna_para_proteina(seq_rna[frame:])
      with open(f"proteina_frame_{sequencia.id}_{frame + 1}.fasta",
                "w") as arquivo_proteina:
        arquivo_proteina.write(f">Sequencia de Proteina (Frame {frame + 1})\n")
        arquivo_proteina.write(seq_proteina + "\n")

      seq_rna_reverso = seq_rna.reverse_complement()
      seq_proteina_reversa = traduzir_rna_para_proteina(
          seq_rna_reverso[frame:])
      with open(f"proteina_frame_reverso_{sequencia.id}_{frame + 1}.fasta",
                "w") as arquivo_proteina:
        arquivo_proteina.write(
            f">Sequencia de Proteina (Frame Reverso {frame + 1})\n")
        arquivo_proteina.write(seq_proteina_reversa + "\n")

print("Traducao para RNA e para proteina concluida.")
