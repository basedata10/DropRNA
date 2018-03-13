# -*- coding: utf-8 -*-
import re, sys, time, json
import cPickle as pickle

#获取一个Barcode所有可能的突变序列
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

#对序列进行反向互补
def rev_comp(seq):
    ___tbl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(___tbl[s] for s in seq[::-1])

#计算两端序列之间的差异
def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0]!=x[1]])

#读入whitelist
def prepare_whitelist(configs, protocol):
    if protocol == "10X":
        bc_white = [set()]
        path = configs['10X']['whitelist']
        with open(path, 'r') as infile:
            for line in infile:
                bc_white[0].add(line.strip())
    if protocol == "indrop":
        bc_white = [set(), set()]
        path = configs['indrop']['whitelist']
        with open(path[0], 'r') as infile:
            for line in infile:
                bc_rev = rev_comp(line.strip())
                bc_white[0].add(bc_rev)
        with open(path[1], 'r') as infile:
            for line in infile:
                bc_rev = rev_comp(line.strip())
                bc_white[1].add(bc_rev)
    if protocol == "dropseq":
        bc_white = []
    return bc_white

def check_barcode_valid_whitelist(protocol, bc_white, barcode):
    if protocol == "10X":
        if barcode in bc_white[0]:
            return 1
        else:
            return 0
    if protocol == "dropseq":
        if barcode in ['TCAAAAGCAGTG']:
            return 0
        else:
            return 1
    if protocol == "indrop":
        if barcode[0:-8] in bc_white[0] and barcode[-8:] in bc_white[1]:
            return 1
        else:
            return 0

def get_top_items(datas, nums=1):
    return [x[0] for x in sorted(datas.items(), key = lambda (k, v): v, reverse = True)[0:nums]]

configFile = sys.argv[1]
barcode_count_file = sys.argv[2]
output = sys.argv[3]
protocol = sys.argv[4]
max_cells = int(sys.argv[5])
min_reads = int(sys.argv[6])

infos = {}
bc_counts = {}
with open(barcode_count_file, 'r') as infile:
    infos = pickle.load(infile)

if 'bc_counts' in infos:
    bc_counts = infos['bc_counts']
else:
    bc_counts = infos

def get_count(bc):
    if bc in bc_counts:
        return bc_counts[bc]
    else:
        return 0

#Read Json
with open(configFile, 'r') as infile:
    configs = json.load(infile)

#Read Barcode WhiteList
bc_white = prepare_whitelist(configs, protocol)

#Build Basic Stats
infos = {}
infos['process'] = {
    "readFile" : "",
    "configFile" : configFile,
    "output" : output,
    "protocol" : protocol,
    "max_cells" : max_cells,
    "min_reads" : min_reads
}

infos['stats'] = {}

#Stats All Barcode
bc_lists = bc_counts.keys()
infos['stats']['all'] = {"reads" : sum([bc_counts[x] for x in bc_lists]), "bc": len(bc_lists)}

#Filter By Min_Reads Step1 And Max Return N_max cells
bc_lists = get_top_items(bc_counts, int(max_cells))
bc_lists = [bc for bc in bc_lists if bc_counts[bc] >= min_reads * 0.2]
infos['stats']['abundent_bc_S1'] = {"reads" : sum([bc_counts[x] for x in bc_lists]), "bc": len(bc_lists)}

#Sort Barcode Lists By Abundance
bc_lists = sorted(bc_lists, key=lambda key:get_count(key), reverse = True)

#Estimate the Mismatch Ratio
def add_mismatch_infos(protocol, bc_lists):
    bc_infos = {}
    mutated_bcs = []

    for idx, bc in enumerate(bc_lists):
        raw = bc_counts[bc]
        if bc in mutated_bcs:
            continue

        bc_infos[bc] = {
            'raw': raw,
            'mis': 0
        }

        if protocol == '10X' or protocol == "indrop":
            for bc_mis in mutate_single_base(bc):
                bc_infos[bc]['mis'] += get_count(bc_mis)
                if bc_mis in bc_lists:
                    mutated_bcs.append(bc_mis)
            bc_infos[bc]['total'] = raw + bc_infos[bc]['mis']

        if protocol == 'dropseq':
            bc_del = [bc[0:11] + x for x in 'ATCG' if x != bc[-1]]
            del_counts = [get_count(x) for x in bc_del]
            for barcode in [x for x in mutate_single_base(bc) if x not in bc_del]:
                if get_count(barcode)>0:
                    bc_infos[bc]['mis'] += get_count(barcode)
                    mutated_bcs.append(barcode)
            if sum([1 for x in del_counts if x>=1000 and x>=0.3*raw and x<=3*raw]) == 3:
                bc_infos[bc]['bc_del'] = 1
                bc_infos[bc]['total'] = raw + sum(del_counts)
                mutated_bcs.extend(bc_del)
            else:
                bc_infos[bc]['bc_del'] = 0
                bc_infos[bc]['total'] = raw + bc_infos[bc]['mis']

    return bc_infos

bc_infos = add_mismatch_infos(protocol, bc_lists)

#Filter By Min_Reads Step2
bc_lists = [bc for bc, v in bc_infos.items() if v['total'] >= min_reads]
infos['stats']['abundent_bc_S2'] = {"reads" : sum([bc_infos[x]['total'] for x in bc_lists]), "bc": len(bc_lists)}

#Filter By Valid
bc_lists = [bc for bc in bc_lists if check_barcode_valid_whitelist(protocol, bc_white, bc)]
infos['stats']['valid_bc'] = {"reads" : sum([bc_infos[x]['total'] for x in bc_lists]), "bc": len(bc_lists)}

#Filter By Mismatch Ratio
#bc_lists = [k for k,v in bc_infos.items() if v['mis_count'] <= 0.1 * v['raw_count']]
#print(len(bc_lists), sum([bc_counts[x] for x in bc_lists]))
#infos['stats']['_bc'] = {"reads" : sum([bc_counts[x] for x in bc_lists]), "bc": len(bc_lists)}

#Write Mismatch Ratio
infos['barcode'] = {}
for bc in bc_lists:
    infos['barcode'][bc] = bc_infos[bc]

##Write The Files
with open(output,'w') as outputfile:
    json.dump(infos, outputfile)
