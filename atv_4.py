from Bio import SeqIO

arquivo_entrada = "denguejunto.fasta"
arquivo_saida_proteina = "proteinas_dengue.fasta"

# Lista de códons de terminação
terminacao_codon = ["TAA", "TAG", "TGA"]

# Inicializando a variável de contagem de proteínas
num_proteina = 1

with open(arquivo_saida_proteina, "w") as saida_proteina:
  # Ler sequência do genoma da Dengue
  for sequencia in SeqIO.parse(arquivo_entrada, "fasta"):
    seq_dna = sequencia.seq

    # Traduzir DNA para RNA
    seq_rna = seq_dna.transcribe()

    # Iterar sobre os três frames de leitura
    for frame in range(3):
      frame_seq = seq_rna[frame:]

      # Inicializar a variável de contagem de proteínas para cada frame
      num_proteina_frame = num_proteina

      # Encontrar e traduzir regiões codificadoras de proteínas
      start_index = frame_seq.find("AUG")
      while start_index != -1:
        end_index = len(frame_seq)
        for term_codon in terminacao_codon:
          term_index = frame_seq.find(term_codon, start_index)
          if term_index != -1 and term_index < end_index:
            end_index = term_index

        # Traduzir a sequência do códon de início ao códon de término
        seq_proteina = frame_seq[start_index:end_index].translate(table=1)

        # Escrever no arquivo de saída, incluindo o número da proteína e o frame
        saida_proteina.write(
            f">{sequencia.description}, Frame {frame+1}, proteína {num_proteina_frame}\n"
        )
        saida_proteina.write(f"{seq_proteina}\n")

        # Procurar o próximo códon de início após o término encontrado
        start_index = frame_seq.find("AUG", end_index)

        # Adicionar o número da proteína para a próxima iteração
        num_proteina_frame += 1

    # Adiconar o número total de proteínas
    num_proteina = num_proteina_frame

print("Identificação e tradução de regiões codificadoras concluída.")
