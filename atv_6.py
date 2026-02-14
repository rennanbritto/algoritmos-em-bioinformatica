def extrair_segmentos(texto):
  segmentos = {}
  lines = texto.split("\n")
  for line in lines:
    if line.startswith("     Region"):
      info = line.split()
      if len(info) >= 6:  # Verifica se tem pelo menos 6 elementos na lista
        inicio, fim = map(int, info[1].split(".."))
        nome = info[5].replace('"', "")
        segmentos[nome] = (inicio, fim)
  return segmentos


def extrair_sequencia(texto, inicio, fim):
  sequencia = texto[inicio - 1:fim]
  sequencia = sequencia.replace(" ", "").replace("\n", "")
  return sequencia


def escrever_fasta(sequencia, nome_segmento):
  with open(f"{nome_segmento}.fasta", "w") as arquivo:
    arquivo.write(f">{nome_segmento}\n")
    for i in range(0, len(sequencia), 80):
      arquivo.write(sequencia[i:i + 80] + "\n")


# Ler o conteúdo do arquivo .pg
caminho_arquivo = "sequence.gp"
with open(caminho_arquivo, "r") as arquivo:
  texto_blast = arquivo.read()

# Extrair os segmentos e suas posições
segmentos = extrair_segmentos(texto_blast)

# Extrair e escrever as sequências dos segmentos em arquivos FASTA
for nome, (inicio, fim) in segmentos.items():
  sequencia = extrair_sequencia(texto_blast, inicio, fim)
  escrever_fasta(sequencia, nome)

print("Arquivos FASTA gerados com sucesso!")
