# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:05:43 2016

@author: blowv6
"""
from Bio import AlignIO
import math
#TASK 1 
#read file
alignment = AlignIO.read("tRNA.stock","stockholm")

#TASK 2
#pi(x) in column_dist; H(i) in col_entropy
seq_num = len(alignment)
seq_len = len(alignment[0])
column_dist = {}
for col in range(0,seq_len):
    dist = {"A":0,"U":0,"G":0,"C":0,"-":0}
    for row in range(0,seq_num):
        dist[alignment[row][col]] += 1./seq_num
    column_dist[col] = dist


column_entropy = {}
for col in range(0,seq_len):
    h = 0;
    for val in column_dist[col].values():
        if val != 0:
            #-sigma(pi*log2(pi))
            h -= val*math.log2(val)
    column_entropy[col] = round(h,6)

#TASK3 p i,j(x,y)
nul = ["A","U","G","C","-"]
colpair_dist = {}
for i in range(0,seq_len):
    for j in range(i,seq_len):
        d = {}
        for nul1 in nul:
            for nul2 in nul:
                d[nul1+nul2] = 0
        for row in range(0,seq_num):
            d[alignment[row][i]+alignment[row][j]] += 1./seq_num
        colpair_dist[(i,j)] = d
                     
#TASK4 I(i,j)
colpair_mi = {}
for i in range(0,seq_len):
    for j in range(i,seq_len):
        mi = 0
        for nul1 in nul:
            for nul2 in nul:
                if colpair_dist[(i,j)][nul1+nul2] != 0:
                    mi += colpair_dist[(i,j)][nul1+nul2] * math.log2(\
                                   colpair_dist[(i,j)][nul1+nul2] \
                                                /(column_dist[i][nul1]\
                                                  *column_dist[j][nul2]))
        colpair_mi[(i,j)] = round(mi,6)
