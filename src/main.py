import os
import itertools
import nw
import numpy as np
import distancias
from tabulate import tabulate


def main():
    # Estos son los fasta que vamos a leer
    archivos = ['apiscoi.fasta', 'cucarachacoi.fasta', 'dromecoi.fasta', 'shrimpcoi.fasta']
    
    taxa_dic = leer_archivos(archivos)

    taxa_dic = alinear_secuencias(taxa_dic)

    jukes_cantor = tabla_jk(taxa_dic, archivos)

def leer_archivos(archivos):

    taxa = dict()
    for archivo in archivos:
        secuencia = ""
        with open(os.getcwd() + '/resources/' + archivo) as file:
            for linea in itertools.islice(file, 1, None):
                secuencia += linea.rstrip()

        taxa[archivo] = secuencia
    return taxa
    
def alinear_secuencias(taxas):
    #(como nw es alineación por pares y no múltiple este no es la alineación más exacta)
    #Escogimos alinearlos todos conforme el str más largo
    taxa_largo, key_taxa_largo = "", None
    for key, value in taxas.items():
        if len(value) > len(taxa_largo):
            taxa_largo = value
            key_taxa_largo = key

    for taxa, secuencia in taxas.items():
        if taxa != key_taxa_largo:
            print("Alineando las secuencias: ", key_taxa_largo, taxa)
            alineadas = nw.needleman_wunsch(taxas[key_taxa_largo], secuencia, 1, -1, -1)
            taxas[key_taxa_largo] = alineadas[0]
            taxas[taxa] = alineadas[1]
    return taxas
    # for taxa, secuencia in taxas.items():
    #     print("Alineando", taxa)
    #     alineadas = nw.needleman_wunsch(taxas['shrimpcoi.fasta'], secuencia, 1, -1, -1)
    #     taxas['shrimpcoi.fasta'] = alineadas[0]
    #     taxas[taxa] = alineadas[1]

def tabla_jk(taxas, archivos):
    tamanio = len(taxas)
    matrix = np.zeros((tamanio, tamanio), dtype = float)
    for i in range(tamanio):
        for j in range(tamanio):
            print("Comparando:", taxas.get())
            matrix[i][j] = distancias.jukes_cantor(taxas[archivos[i]], taxas[archivos[j]])
    return matrix


main()