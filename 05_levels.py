# -*- coding: utf-8 -*-
import itertools, json, re, sys, os
from collections import Counter
import operator
import random

umi_folder = sys.argv[1]
sampling_factor = float(sys.argv[2])
out_file = sys.argv[3]
out_file2 = sys.argv[4]
out_file3 = sys.argv[5]

def HammingDistance(seq1, seq2):
    return sum(itertools.imap(operator.ne, seq1, seq2))

def sampling_UMIs(UMIs, ratio):
    UMIs_sample = {}
    for UMI in UMIs:
        UMIs_sample[UMI] = sum([1 for x in range(UMIs[UMI]) if random.random()<=ratio])
    return {x:UMIs_sample[x] for x in UMIs_sample if UMIs_sample[x]>0}

def calculate_UMI_with_mismatch(UMIs):
    if len(UMIs.keys()) == 1:
        return [x for x in UMIs if UMIs[x]>0]
    UMIs = sorted(UMIs.items(), key=lambda k: k[1], reverse=True)
    UMI_info = {x[0]:x[1] for x in UMIs}
    umi_num = len(UMIs)
    if umi_num <= 10:
        for idx1 in range(0, umi_num-1):
            for idx2 in range(idx1+1, umi_num):
                umi_1 = UMIs[idx1][0]
                umi_2 = UMIs[idx2][0]
                if HammingDistance(umi_1, umi_2) <= 1:
                    UMI_info[umi_1] += UMI_info[umi_2]
                    UMI_info[umi_2] = 0
    return [x for x in UMI_info if UMI_info[x]>0]

# test = {'TCCGAT': 33, 'TCCCTA': 10, 'ACCGAT': 1, 'CACGGG': 7, 'ACCGAA':100}
# calculate_UMI_with_mismatch(test)
#umi_folder = "./reads_10x1"

paths = [os.path.join(umi_folder, x) for x in os.listdir(umi_folder) if re.match("UMI.*txt",x)]

genes = {}
UMIs_count = {}
readCounts = {}
barcodes = []


for filepath in paths:
    with open(filepath, 'r') as file:
        for line in file:
            infos = re.split('\t', line)
            barcode = infos[0]
            barcodes.append(barcode)
            umi_counts = json.loads(infos[1])
            for gene in umi_counts:
                gene_UMIs = sampling_UMIs(umi_counts[gene], sampling_factor)
                bcs = calculate_UMI_with_mismatch(gene_UMIs)
                reads = sum(gene_UMIs.values())
                for bc in bcs:
                    if bc in UMIs_count:
                        UMIs_count[bc] += 1
                    else:
                        UMIs_count[bc] = 1
                if gene in genes:
                    genes[gene][barcode] = len(bcs)
                    readCounts[gene][barcode] = reads
                else:
                    genes[gene] = {}
                    genes[gene][barcode] = len(bcs)
                    readCounts[gene] = {}
                    readCounts[gene][barcode] = reads

#Export Genes Which express in at least 10 samples
header = "\t".join(["genes"]+barcodes)
with open(out_file, 'w') as file:
    file.write(header+"\n")
    for gene in genes:
        umis = [gene]
        if len(genes[gene].keys())<= 0.02*len(barcodes) or sum(genes[gene].values())<=0.1*len(barcodes):
            continue
        for bc in barcodes:
            if bc in genes[gene]:
                umis.append(genes[gene][bc])
            else:
                umis.append(0)
        umis_str = "\t".join([str(x) for x in umis])
        file.write(umis_str+"\n")

with open(out_file2, 'w') as file:
    file.write(header+"\n")
    for gene in readCounts:
        umis = [gene]
        if sum(readCounts[gene].values())<=0.1*len(barcodes):
            continue
        for bc in barcodes:
            if bc in readCounts[gene]:
                umis.append(readCounts[gene][bc])
            else:
                umis.append(0)
        umis_str = "\t".join([str(x) for x in umis])
        file.write(umis_str+"\n")

with open(out_file3, 'w') as file:
    for bc in UMIs_count:
        if UMIs_count[bc]>=100:
            file.write(bc+"\t"+ str(UMIs_count[bc])+ "\n")
