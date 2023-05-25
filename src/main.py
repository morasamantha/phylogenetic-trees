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

    for i in range(1, l):
        print("Alineando", i)
        alineacion = nw.needleman_wunsch(taxas[i-1], taxas[i], 1, -1, -1)
        taxas[i] = alineacion[1]
        taxas[i-1] = alineacion[0]

    return taxas

def tabla_jk(taxas, archivos):
    tamanio = len(taxas)
    matrix = np.zeros((tamanio, tamanio), dtype = float)
    for i in range(tamanio):
        for j in range(tamanio):
            matrix[i][j] = distancias.jukes_cantor(taxas[i], taxas[j])
    return matrix


main()