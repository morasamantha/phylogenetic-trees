import os
import itertools
import nw
import msa
import numpy as np
import distancias
from tabulate import tabulate
import trees

def main():
    # Estos son los fasta que vamos a leer. DEBEN estar en la carpeta de resources
    archivos = ['apiscoi.fasta', 'cucarachacoi.fasta', 'dromecoi.fasta', 'shrimpcoi.fasta']
    
    taxas = leer_archivos(archivos)

    # Para guardar las tablas
    modelos_distancia = ['naive', 'jk', 'mk', 'tamura']
    tablas_pares, tablas_multiples = {},{} # de hecho son dicts pero weno xp

    # Para la alineación uno a uno
    for modelo in modelos_distancia: 
        ruta_p = '/tablas/nj/' + modelo  + '.txt'
        tabla_p = alinear_secuencias_pares(taxas, archivos, modelo, ruta_p)
        tabla_pares[modelo] = tabla_p
    
    # Para la alineación múltiple
    taxas_multiples = alinear_secuencias_multiples(taxas)

    for modelo in modelos_distancia:
        ruta_m = '/tablas/msa/' + modelo  + '.txt'
        tabla_m = tabla(taxas, modelo, archivos, ruta_m)
        tablas_multiples[modelo] = tabla_m

    # Árboles
    for modelo, tabla in tablas_pares.items():
        ruta = '/graphs/p/' + modelo  + '.txt'
        trees.neighbor_joining(tabla, archivos, ruta)

    for modelo, tabla in tablas_pares.items():
        ruta = '/graphs/m/' + modelo  + '.txt'
        trees.neighbor_joining(tabla, archivos, ruta)

def leer_archivos(archivos):
    taxa = []
    for archivo in archivos:
        secuencia = ""
        with open(os.getcwd() + '/resources/' + archivo) as file:
            for linea in itertools.islice(file, 1, None):
                secuencia += linea.rstrip()
        taxa.append(secuencia)
    return taxa
    
def alinear_secuencias_pares(taxas, nombres, modelo, ruta):
    tamanio = len(taxas)
    matrix = np.zeros((tamanio, tamanio), dtype = float)

    for i in range(tamanio):
        for j in range(i+1, tamanio):
            alineadas = nw.needleman_wunsch(taxas[i], taxas[j], 1, -1, -1)
            distancia = distancias.calcula_modelo(alineadas[0], alineadas[1], modelo)
            matrix[i][j] = distancia 
            matrix[j][i] = distancia

    with open(os.getcwd() + ruta, 'w') as new_file:
        new_file.write(tabulate(matrix, headers=nombres, tablefmt="pretty"))

    return matrix

def alinear_secuencias_multiples(taxas):
    return msa.msa(len(taxas), taxas)

def tabla(taxas, modelo, nombres, ruta):
    tamanio = len(taxas)
    matrix = np.zeros((tamanio, tamanio), dtype = float)

    for i in range(tamanio):
        for j in range(tamanio):
            if i != j: 
                matrix[i][j] = distancias.calcula_modelo(taxas[i], taxas[j], modelo)
    
    with open(os.getcwd() + ruta , 'w') as new_file:
        new_file.write(tabulate(taxas, headers=nombres, tablefmt="pretty"))
    
    return matrix

main()