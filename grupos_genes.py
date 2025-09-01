import os
import glob
import pandas as pd
import numpy as np
tops = []
tails = []
meds = []
# Recorremos las distintas réplicas
for i in [10,11,12]:
    # Este tocho para sacar los archivos de RIN max y min de cada rep
    max = 0
    min = 10
    for x in glob.glob(f"../Procesado/outputs/TS{i}*"):
        rin = float(x.split("_")[2])
        if rin >= max:
            fich_max = x
            max = rin
        if rin <= min:
            fich_min = x
            min = rin
    print(fich_max)
    print(fich_min)
    # Leemos los csv
    RINmin = pd.read_csv(fich_min, sep = "\t")
    RINmax = pd.read_csv(fich_max, sep = "\t")
    # Y ahora nos quedamos con el porcentaje de la isoforma más frecuente de cada gen
    # Para el de RIN mínimo
    RINmin["porc max"] = RINmin.groupby(by="associated_gene")["porcentaje"].transform("max")
    genes_RINmin = RINmin[["associated_gene","porc max"]].drop_duplicates()
    # Para el de RIN máximo
    RINmax["porc max"] = RINmax.groupby(by="associated_gene")["porcentaje"].transform("max")
    genes_RINmax = RINmax[["associated_gene","porc max"]].drop_duplicates()
    print("He llegado aquí")
    # Joineamos por nombre de gen
    fusion = pd.merge(genes_RINmin,genes_RINmax,on = ["associated_gene", "associated_gene"], how = "inner")
    print(len(fusion))
    # Sacamos la diferencia absoluta
    fusion["dif absoluta"] = abs(fusion["porc max_x"] - fusion["porc max_y"])
    fusion = fusion.sort_values(by="dif absoluta", ascending = False).reset_index(drop=True)
    # Nos quedamos con un top de genes 
    top1830 = fusion.head(650)
    print(f"El mínimo valor de dif absoluta es: {top1830["dif absoluta"].max()}")
    top_genes = top1830["associated_gene"]
    tops.append(top_genes)
    tail_genes = fusion.tail(90)["associated_gene"]
    tails.append(tail_genes)
    # Sacamos la mediana de las diferencias y vemos dónde está
    mediana = np.median(fusion["dif absoluta"])
    mediana_pos = fusion.index[fusion["dif absoluta"] == mediana]
    print(mediana)
    n = 400
    genes_medios = fusion[mediana_pos[0]-n:mediana_pos[0]+n]["associated_gene"]
    meds.append(genes_medios)
    
# Aquí vemos las intersecciones
print(os.listdir())
with open("output/top.txt","w") as f:
    print("Intersección de los top entre todos")
    inter = set(tops[0]) & set(tops[1]) & set(tops[2])
    print(inter)
    for i in inter:
        f.write(f'{i}\n')
    
with open("output/tail.txt","w") as f:
    print("Intersección de los ultimos entre todos")
    inter = set(tails[0]) & set(tails[1]) & set(tails[2])
    print(inter)
    for i in inter:
        f.write(f'{i}\n')
    
with open("output/med.txt","w") as f:
    print("Intersección de los medianos entre todos")
    inter = set(meds[0]) & set(meds[1]) & set(meds[2])
    print(inter)
    for i in inter:
        f.write(f'{i}\n')
