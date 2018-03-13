# -*- coding: utf-8 -*-
import json, time, sys
import subprocess

reads_p1 = sys.argv[1]
reads_p2 = sys.argv[2]
barcode_files = sys.argv[3]
split_folder = sys.argv[4]
protocol = sys.argv[5]

with open(barcode_files, 'r') as infile:
    bc_infos = json.load(infile)

barcode_lists = set(bc_infos['barcode'].keys())

p1 = subprocess.Popen(["zcat", reads_p1], stdout = subprocess.PIPE, bufsize=100000000)
read1 = p1.stdout
p2 = subprocess.Popen(["zcat", reads_p2], stdout = subprocess.PIPE, bufsize=100000000)
read2 = p2.stdout

prefixs = []
outfiles = {}
outbuffers = {}
bases = ['A','T', 'C', 'G']
for base1 in bases:
    for base2 in bases:
        outfiles[base1+base2] = open(split_folder + '/' +base1+base2+ ".fq", 'wt')
        outbuffers[base1+base2] = []
        prefixs.append(base1+base2)

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0]!=x[1]])

def mutate_single_base(seq):
    mutated = []
    bases = ['A', 'T', 'C', 'G']
    for index in range(0, len(seq)):
        temp = list(seq)
        base_raw = temp[index]
        for base in bases:
            if base != base_raw:
                temp[index] = base
                mutated.append(''.join(temp))
    return mutated

def get_barcode(protocol, seq):
    if protocol == "10X":
        return seq[0:16]
    if protocol == "dropseq":
        return seq[0:12]
    if protocol == "indrop":
        w1 = "GAGTGATTGCTTGTGACGCCTT"
        if w1 in seq:
            w1_pos = seq.find(w1)
            if 7 < w1_pos < 12:
                return seq[0:w1_pos] + seq[w1_pos + 22:w1_pos + 22 + 8]
        else:
            for i in range(8, 12):
                w1_mutate = seq[i:i + 22]
                if HammingDistance(w1_mutate, w1) < 2:
                    return seq[0:i] + seq[i + 22:i + 22 + 8]
                    break
        return ""

def barcode_corrected(protocol, barcode):
    if barcode in barcode_lists:
        return barcode
    if protocol == "10X" or protocol == "indrop":
        for mut_bc in mutate_single_base(barcode):
            if mut_bc in barcode_lists:
                return mut_bc
        return ""
    if protocol == "dropseq":
        bc_dels = [barcode[0:11] + x for x in 'ATCG']
        for bc_del in bc_dels:
            if bc_del in barcode_lists:
                return bc_del
        return ""

index = 0
invalid = 0
starttime = time.time()
while True:
    index += 1
    header1 = read1.readline().strip()
    if not header1: break
    seq1 = read1.readline().strip()
    read1.readline()
    read1.readline()

    header2 = read2.readline().strip()
    seq2 = read2.readline().strip()
    read2.readline()
    qual2 = read2.readline().strip()

    barcode = get_barcode(protocol, seq1)
    barcode_corr = barcode_corrected(protocol, barcode)

    if barcode_corr == "":
        continue

    barcode_pre = barcode_corr[0:2]

    if barcode_pre[0] == 'N' or barcode_pre[1] == 'N':
        continue

    if protocol == "10X":
        UMI = seq1[16:26]
    if protocol == "dropseq":
        if bc_infos['barcode'][barcode_corr]['bc_del']:
            UMI = seq1[11:19]
        else:
            UMI = seq1[12:20]
    if protocol == "indrop":
        UMI = seq1[len(barcode)+22:len(barcode)+22+6]

    outbuffers[barcode_pre].append('_'.join(['@', barcode_corr, UMI, str(index)])+ "\n")
    outbuffers[barcode_pre].append(seq2[0:100]+"\n")
    outbuffers[barcode_pre].append("+\n")
    outbuffers[barcode_pre].append(qual2[0:100]+"\n")

    if index % 1000000 == 0:
        print("Processed", index, time.time()-starttime)
        starttime = time.time()
        for prefix in prefixs:
            outfiles[prefix].writelines(outbuffers[prefix])
            outbuffers[prefix] = []

for prefix in prefixs:
    outfiles[prefix].writelines(outbuffers[prefix])
    outbuffers[prefix] = []
    outfiles[prefix].close()
