import os
import itertools
import nw
import numpy as np
import distancias
from tabulate import tabulate


def main():
    # Estos son los fasta que vamos a leer
    archivos = ['apiscoi.fasta', 'cucarachacoi.fasta', 'dromecoi.fasta', 'shrimpcoi.fasta']
    
    taxa = leer_archivos(archivos)

    taxa = alinear_secuencias(taxa)

    jukes_cantor = tabla_jk(taxa, archivos)

def leer_archivos(archivos):
    taxa = []
    for archivo in archivos:
        secuencia = ""
        with open(os.getcwd() + '/resources/' + archivo) as file:
            for linea in itertools.islice(file, 1, None):
                secuencia += linea.rstrip()
        taxa.append(secuencia)
    return taxa
    
def alinear_secuencias(taxas):
    #(como nw es alineación por pares y no múltiple este no es la alineación más exacta)
    #Escogimos alinearlos todos conforme el str más largo
    l = len(taxas)
    taxa_largo, index_longest = "", -1
    for i in range(l):
        if len(taxas[i]) > len(taxa_largo):
            taxa_largo = taxas[i]
            index_longest = i
    # for i in range(l):
    #     for j in range(l):
    #         if i != j:
    #             print("Alineando", i, j)
    #             alineacion = nw.needleman_wunsch(taxas[i], [j], 1, -1, -1)

    for secuencia in range(l):
        if secuencia > 0:
            print("Alineando con secuencia", secuencia)
            alineadas = nw.needleman_wunsch(taxas[index_longest], taxas[secuencia], 1, -1, -3)
            taxas[index_longest] = alineadas[0]
            taxas[secuencia] = alineadas[1]
            try:
                with open('pruebas/' + str(secuencia), 'w') as f:
                    f.write(taxas[secuencia])
            except FileNotFoundError:
                print("AAA")

    return taxas

def tabla_jk(taxas, archivos):
    tamanio = len(taxas)
    matrix = np.zeros((tamanio, tamanio), dtype = float)
    for i in range(tamanio):
        for j in range(tamanio):
            matrix[i][j] = distancias.jukes_cantor(taxas[i], taxas[j])
    return matrix


main()