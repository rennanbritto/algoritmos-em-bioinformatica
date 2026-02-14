from Bio import SeqIO
import math
import matplotlib.pyplot as plt


def ler_arquivo_fasta(caminho_arquivo):
  sequencias = {}
  numero_sequencias = 0
  with open(caminho_arquivo, 'r') as arquivo:
    linha = arquivo.readline()
    while linha:
      if linha.startswith('>'):
        cabecalho = linha.strip()[1:]
        sequencias[cabecalho] = ''
        linha = arquivo.readline()
        while linha and not linha.startswith('>'):
          sequencias[cabecalho] += linha.strip()
          linha = arquivo.readline()
          numero_sequencias += 1
      else:
        linha = arquivo.readline()
  return sequencias, numero_sequencias


def calcular_contagem_nucleotideos(sequencia):
  count_A = sequencia.count("A")
  count_T = sequencia.count("T")
  count_C = sequencia.count("C")
  count_G = sequencia.count("G")
  total_nucleotideos = len(sequencia)
  return total_nucleotideos, count_A, count_T, count_C, count_G


def calcular_conteudo_GC(sequencia):
  total_nucleotideos, count_A, count_T, count_C, count_G = calcular_contagem_nucleotideos(
      sequencia)
  conteudo_GC = (count_G + count_C) / total_nucleotideos * 100
  return conteudo_GC


def calcular_temperatura_anelamento(conteudo_GC, total_nucleotideos):
  return 64.9 + 0.41 * conteudo_GC - (500 / total_nucleotideos)


def escrever_resultados(caminho_arquivo_saida, resultados):
  with open(caminho_arquivo_saida, "w") as arquivo_saida:
    arquivo_saida.write(f"Número total de sequências: {numero_sequencias}\n\n")
    for cabecalho, (sequencia, conteudo_GC,
                    temperatura_anelamento) in resultados.items():
      arquivo_saida.write(f"Sequência: {cabecalho}\n")
      arquivo_saida.write(f"Conteúdo GC: {conteudo_GC:.2f}%\n")
      arquivo_saida.write(
          f"Temperatura de Anelamento: {temperatura_anelamento:.2f} °C\n")

      # Adiciona a contagem de cada nucleotídeo
      arquivo_saida.write(f"Contagem de A: {sequencia.count('A')}\n")
      arquivo_saida.write(f"Contagem de T: {sequencia.count('T')}\n")
      arquivo_saida.write(f"Contagem de C: {sequencia.count('C')}\n")
      arquivo_saida.write(f"Contagem de G: {sequencia.count('G')}\n\n")

      # Adiciona a temperatura de melting ao dicionário de resultados
      resultados[cabecalho] = (sequencia, conteudo_GC, temperatura_anelamento)


def main(caminho_arquivo_entrada, caminho_arquivo_saida):
  sequencias = ler_arquivo_fasta(caminho_arquivo_entrada)
  resultados = {}

  # Cálculo da temperatura de melting para cada sequência
  for cabecalho, sequencia in sequencias.items():
    total_nucleotideos, count_A, count_T, count_C, count_G = calcular_contagem_nucleotideos(
        sequencia)
    conteudo_GC = calcular_conteudo_GC(sequencia)
    temperatura_anelamento = calcular_temperatura_anelamento(
        conteudo_GC, total_nucleotideos)
    resultados[cabecalho] = (sequencia, conteudo_GC, temperatura_anelamento)

  # Escreve os resultados em um arquivo
  escrever_resultados(caminho_arquivo_saida, resultados, numero_sequencias)

  # Criação do gráfico
  conteudos_GC = [conteudo_GC for _, (_, conteudo_GC, _) in resultados.items()]
  temperaturas_melting = [
      temperatura_anelamento
      for _, (_, _, temperatura_anelamento) in resultados.items()
  ]

  plt.scatter(conteudos_GC, temperaturas_melting)
  plt.xlabel('Conteúdo GC (%)')
  plt.ylabel('Temperatura de Melting (°C)')
  plt.title('Conteúdo GC vs. Temperatura de Melting')
  plt.grid(True)
  plt.savefig("grafico.png")  # Salva o gráfico como um arquivo de imagem
  plt.show()


if __name__ == "__main__":
  caminho_arquivo_entrada = 'deng.fasta'
  caminho_arquivo_saida = 'resultados.txt'
  main(caminho_arquivo_entrada, caminho_arquivo_saida)
