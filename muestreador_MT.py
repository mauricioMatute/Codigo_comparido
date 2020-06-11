#!/usr/local/bin/anaconda3/bin/python

import random
import csv
from Bio import SeqIO
import numpy as np                                                       

archivo="Caenorhabditis_elegans.WBcel235.dna.chromosome.MtDNA.fa" #Aquí escribe el nombre del archivo que quieres leer
iteraciones=5000 #Número de ventanas sobre las cuales calcular las variables, es decir, el número de observaciones en la base de datos.

#Esta rutina extrae la secuencia del archivo .fasta (suponemos que solo hay una por archivo) y calcula su longitud
def lee_sec(arch):
    fasta_sequences = SeqIO.parse(arch,'fasta')
    for seq_record in fasta_sequences:
        secuencia=str(seq_record.seq).upper()
        longi=len(seq_record)
    return longi, secuencia

#Esta rutina calcula el valor de 'iteraciones' observaciones de los sensores para el archivo dado, tomando ventanas de tamaño aleatorio (puedes cambiarlo a que sea un tamaño fijo más abajo) y arma una base en formato .csv 
def calcula_sensores(longitud, secuencia):
    respuesta=[i for i in range(11)]
    with open('datos_'+archivo+'.csv','w') as newFile: #Aquí puedes cambiar el nombre del archivo que se genera con la base de datos.
        newFileWriter=csv.writer(newFile)
        newFileWriter.writerow(['CG_content','CpG_ratio','dws','dry','dmk','LTP','HTP','VTP','ITP','tam_ventana','inicio_ventana'])
        for i in range(1,iteraciones):
            window=random.randint(100,1000) #Este es el tamaño de la ventana, ahorita es aleatorio uniforme entre 100 y 1000, pero podría ser fijo, solo cambiando esta línea.
            inicio=random.randint(1,longitud-window)
            CGcont=0
            CpGcont=0
            Nww=0
            Nss=0
            Nws=0
            Nsw=0
            Nw=0
            Ns=0
            Nrr=0
            Nyy=0
            Nry=0
            Nyr=0
            Nr=0
            Ny=0
            Nmm=0
            Nkk=0
            Nmk=0
            Nkm=0
            Nm=0
            Nk=0
            LTPcont=0
            HTPcont=0
            VTPcont=0
            ITPcont=0
            cAux='N'
            for j in range(inicio,inicio+window+1):
                c=secuencia[j]
                if c=='C':
                    CGcont+=1
                    Ns+=1
                    Ny+=1
                    Nm+=1
                    if cAux=='C':
                        Nss+=1
                        Nyy+=1
                        Nmm+=1
                        LTPcont+=1
                    elif cAux=='G':
                        Nss+=1
                        Nry+=1
                        Nkm+=1
                        HTPcont+=1
                    elif cAux=='A':
                        Nws+=1
                        Nry+=1
                        Nmm+=1
                        ITPcont+=1
                    elif cAux=='T':
                        Nws+=1
                        Nyr+=1
                        Nkm+=1
                        HTPcont+=1
                elif c=='G':
                    CGcont+=1
                    Ns+=1
                    Nr+=1
                    Nk+=1
                    if cAux=='C':
                        CpGcont+=1
                        Nss+=1
                        Nyr+=1
                        Nmk+=1
                        VTPcont+=1
                    elif cAux=='G':
                        Nss+=1
                        Nrr+=1
                        Nkk+=1
                        HTPcont+=1
                    elif cAux=='A':
                        Nws+=1
                        Nrr+=1
                        Nkk+=1
                        LTPcont+=1
                    elif cAux=='T':
                        Nws+=1
                        Nyr+=1
                        Nkk+=1
                        VTPcont+=1
                elif c=='A':
                    Nw+=1
                    Nr+=1
                    Nm+=1
                    if cAux=='C':
                        Nsw+=1
                        Nyr+=1
                        Nmm+=1
                        VTPcont+=1
                    elif cAux=='G':
                        Nsw+=1
                        Nrr+=1
                        Nkm+=1
                        HTPcont+=1
                    elif cAux=='A':
                        Nww+=1
                        Nrr+=1
                        Nmm+=1
                        LTPcont+=1
                    elif cAux=='T':
                        Nww+=1
                        Nyr+=1
                        Nkm+=1
                        VTPcont+=1
                elif c=='T':
                    Nw+=1
                    Ny+=1
                    Nk+=1
                    if cAux=='C':
                        Nsw+=1
                        Nyy+=1
                        Nmk+=1
                        LTPcont+=1
                    elif cAux=='G':
                        Nsw+=1
                        Nry+=1
                        Nkk+=1
                        ITPcont+=1
                    elif cAux=='A':
                        Nww+=1
                        Nry+=1
                        Nmk+=1
                        ITPcont+=1
                    elif cAux=='T':
                        Nww+=1
                        Nyy+=1
                        Nkk+=1
                        LTPcont+=1
                cAux=c
            CGcontent=CGcont/window
            if CGcont!=0:
                CpGratio=CpGcont/(((CGcont/2)**2)/window)
            else:
                CpGratio=0
            if Ns!=0 and Nw!=0:
                dws=(Nww*Nss-Nsw*Nws)/(Ns*Nw)
            else:
                dws=0
            if Nr!=0 and Ny!=0:
                dry=(Nrr*Nyy-Nyr*Nry)/(Nr*Ny)
            else:
                dry=0
            if Nm!=0 and Nk!=0:
                dmk=(Nmm*Nkk-Nkm*Nmk)/(Nm*Nk)
            else:
                dmk=0
            LTP=LTPcont/(window-1)
            HTP=HTPcont/(window-1)
            VTP=VTPcont/(window-1)
            ITP=ITPcont/(window-1)

            respuesta[0]=CGcontent
            respuesta[1]=CpGratio
            respuesta[2]=dws
            respuesta[3]=dry
            respuesta[4]=dmk
            respuesta[5]=LTP
            respuesta[6]=HTP
            respuesta[7]=VTP
            respuesta[8]=ITP
            respuesta[9]=window
            respuesta[10]=inicio

            newFileWriter.writerow(respuesta)

var1, var2=lee_sec(archivo)
calcula_sensores(var1, var2)
