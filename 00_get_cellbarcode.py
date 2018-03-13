# -*- coding: utf-8 -*-
import re, sys, time, json
import subprocess
import cPickle as pickle

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0]!=x[1]])

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
                    return seq[0:i] + seq[i + 22 : i + 22 + 8]
                    break
        return ""

configFile = sys.argv[1]
readFile = sys.argv[2]
output = sys.argv[3]
protocol = sys.argv[4]
min_reads = int(sys.argv[5])

#Prepare The Output
infos = {}

#Read Json
with open(configFile, 'r') as infile:
    configs = json.load(infile)

#Read Files
bc_counts = {}
initTime = time.time()
starttime = time.time()
process_reading = subprocess.Popen(["zcat", readFile], stdout = subprocess.PIPE, bufsize=100*1000000)
infile = process_reading.stdout
index = 0
while True:
    index += 1
    header = infile.readline().strip()
    if not header: break
    seq = infile.readline().strip()
    temp = infile.readline()
    temp = infile.readline()
    bc = get_barcode(protocol, seq)
    if bc == "":
        continue
    if bc in bc_counts:
        bc_counts[bc] += 1
    else:
        bc_counts[bc] = 1
    if index % 1000000 == 1:
        print("Processed", index-1, time.time() - starttime)
        starttime = time.time()

bc_counts = {k:v for k,v in bc_counts.items() if v>=min_reads}

infos['bc_counts'] = bc_counts
infos['total_reads'] = index
infos['filter_low_depth'] = sum(bc_counts.values())

with open(output, 'wb') as outfile:
    pickle.dump(infos, outfile)
