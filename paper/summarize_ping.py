# %% Initialization
import pickle
import textwrap
import sys
import math, statistics
import os
import collections
import multiprocessing as mp
from pprint import pprint, pformat
from itertools import repeat, groupby
from operator import itemgetter
import subprocess as sp
import pysam
import glob
import itertools
import re
import mappy
from tqdm import tqdm
sys.path=[".."]+sys.path
from common import  timeit, timing
import parasail
from Bio import Entrez
Entrez.email = "inumanag@uvic.ca"

from common import Region, Allele, Gene
from Bio import SeqIO

@timeit
def initialize(file="/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/kir.pickle"):
    from common import Region, Allele, Gene
    import time

    print(time.ctime(os.path.getmtime(file)))
    with open(file, "rb") as fo:
        (genes, _, _, minor_to_major) = pickle.load(fo)
    return genes, minor_to_major

@timeit
def calculate_maps(genes):
    for g in genes.values():
        for a in g.alleles.values():
            # Prepare translation maps
            a.idx_seq = a.seq
            wildtype_seq = a.gene.wildtype.seq
            prev, off = 0, 0
            pieces = []
            for pos, op in sorted(a.ops):
                if prev < pos:
                    pieces.append((prev, prev + off, wildtype_seq[prev:pos], ""))
                    prev = pos
                if op.startswith("ins"):
                    pieces.append((pos, pos + off, op[3:], (pos, op)))
                    off += len(op) - 3
                elif op.startswith("del"):
                    pieces.append((pos, pos + off, "", (pos, op)))
                    off -= len(op) - 3
                    prev += len(op) - 3
                else:
                    pieces.append((pos, pos + off, op[2], (pos, op)))
                    prev = pos + 1
            if prev < len(wildtype_seq):
                pieces.append((prev, prev + off, wildtype_seq[prev : len(wildtype_seq)], ""))

            # Mutation map (allele-indexed)
            # - tuple (mutation itself) for mutations within allele
            # - True for functional mutations that are not part of the allele
            # - None otherwise
            a.mutmap = {}
            other_func = {m[0]: m for m in a.gene.functional if m not in a.mutations}
            for wi, ai, s, m in pieces:
                if m and m[1].startswith("ins") and m in a.mutations:
                    for i in range(-1, len(s) + 1): a.mutmap[ai + i] = (0, m)
                if m and m[1].startswith("del") and m in a.mutations:
                    a.mutmap[ai] = a.mutmap[ai + 1] = (0, m)
                for j in range(len(s)):
                    if m and m in a.mutations: a.mutmap[ai + j] = (0, m)
                    elif wi + j in other_func: a.mutmap[ai + j] = (1, other_func[wi + j])

if __name__ == "__main__":
    genes, minor_to_major = initialize()
    calculate_maps(genes)

    # Check which alleles have mismatches with IPD_KIR sequences
    PMM, PGO, PGE = parasail.matrix_create("ACGT", 2, -8), 12, 2
    def get_cigar(a):
        l = re.split(r'([A-Z=])', a.cigar.decode.decode())[:-1]
        return [(int(l[i]), l[i+1]) for i in range(0, len(l), 2)]

    for g in genes.values():
        break
        orig = {}
        with pysam.FastxFile(f"../IPDKIR/fasta/{g.name}_gen.fasta") as fa:
            for r in fa:
                orig[r.comment.split()[0].split('*')[1]] = r.sequence
        for a in g.alleles.values():
            if a.name not in orig:  # exon-only sequences
                pass
            elif a.seq != orig[a.name]:
                assert len(a.seq) >= len(orig[a.name])
                aln = parasail.sg_qx_trace_scan_32(a.seq, orig[a.name], PGO, PGE, PMM)
                cigar = get_cigar(aln)
                rr = collections.defaultdict(int)
                if cigar[0][1] == 'I':
                    i = 0
                    for r in a.regions.values():
                        if i >= cigar[0][0]: break
                        rr[r.name] += 1
                        i += len(r.seq)
                if len(cigar) > 1 and cigar[-1][1] == 'I':
                    i = 0
                    for r in list(a.regions.values())[::-1]:
                        if i >= cigar[-1][0]: break
                        rr[r.name] += 1
                        i += len(r.seq)
                if len(cigar) == 3 and cigar[1][1] == '=' and set(rr) <= a.new_regions:
                    continue
                print(g.name, a.name, '\n  ',
                    [(r, len(a.regions[r].seq)) for r in a.new_regions],
                    any(r.partial for r in a.regions.values()),
                    len(a.seq), len(orig[a.name]), '\n  ',
                    list(rr.keys()), aln.cigar.decode.decode())

ids = [ 'NT_113949.2', 'NT_187636.1', 'NT_187637.1', 'NT_187638.1', 'NT_187639.1', 'NT_187640.1', 'NT_187641.1',
        'NT_187642.1', 'NT_187643.1', 'NT_187644.1', 'NT_187645.1', 'NT_187668.1', 'NT_187669.1', 'NT_187670.1',
        'NT_187671.1', 'NT_187672.1', 'NT_187673.1', 'NT_187674.1', 'NT_187675.1', 'NT_187676.1', 'NT_187677.1',
        'NT_187683.1', 'NT_187684.1', 'NT_187685.1', 'NT_187686.1', 'NT_187687.1', 'NT_187693.1',
        'NW_003571054.1', 'NW_003571055.2', 'NW_003571056.2', 'NW_003571057.2', 'NW_003571058.2', 'NW_003571059.2',
        'NW_003571060.1', 'NW_003571061.2', 'NW_016107300.1', 'NW_016107301.1', 'NW_016107302.1', 'NW_016107303.1',
        'NW_016107304.1', 'NW_016107305.1', 'NW_016107306.1', 'NW_016107307.1', 'NW_016107308.1', 'NW_016107309.1',
        'NW_016107310.1', 'NW_016107311.1', 'NW_016107312.1', 'NW_016107313.1', 'NW_016107314.1']

#%% HPRC nice
ground = {}
idx = [g for g in genes]
with open("/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/geny/paper/ground/hprc_gt_0801.csv") as f:  # based on HPRC-kir-annotations.extra.tsv
    for l in f:
        if l.startswith('ID'):
            idx = l.strip().split(',')[1:]
            continue
        l = l.strip().split(',')
        ground[l[0].split('.')[0]] = [[z[1:] if z and z[0] == '`' else z for z in i.split(';')] for i in l[1:]]
#%%
# pi = idx.index('KIR3DP1')  patches to HPRC-kir-annotations.extra.tsv
# ground['HG00741'][pi] = ['0060103', '007?','009?']
# ground['HG01358'][pi] = ['0079991', '0030202']  #7-snp
# ground['HG02257'][pi] = ['0090102', '0030203']
# ground['HG02572'][pi] = ['0150101', '0109991']  #10+14
# ground['HG02622'][pi] = ['0030202', '0090102', '0100102']  #9+snp!
# ground['HG02630'][pi] = ['0030216', '0030216', '0119991']  #011+004
# ground['HG02717'][pi] = ['0090102', '0100102', '0030202']  #9+4
# ground['HG03540'][pi] = ['017', '0100102']
#%%
ping = {}
# idx = []
# def parse(i):
#     i = i.split('*')[1]
#     return i[:3] if i != 'unresolved' else '???'
# with open("/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/PING_tuned_calls.csv") as f:  # based on HPRC-kir-annotations.extra.tsv
#     for l in f:
#         if l.startswith(','):
#             idx = l.strip().split(',')[1:]
#             continue
#         l = l.strip().split(',')
#         ping[l[0]] = {
#             idx[i]:
#             list(set([
#                 tuple(
#                     parse(xx)
#                     for xx in z.split('+')
#                     if not xx.endswith('*null')
#                 )
#                 for z in l[i + 1].split(' ')
#             ]))
#             for i in range(len(idx))
#         }

ping_pre = {}
idx = []

# def parse(ii):
#     if '*' not in ii:
#         return '???'
#     i = ii.split('*')[1]
#     gene = ii.split('*')[0]
#     return f'{gene}*{i[:3]}' if i != 'unresolved' else '???'
def parse(i):
    if '*' not in i:
        return '???'
    
    i = i.split('*')[1]
    return i[:3] if i != 'unresolved' else '???'

def parse_special(i):
    if '*' not in i:
        return '???'
    g = i.split('*')[0]
    i = i.split('*')[1]
    return f'{g}*{i[:3]}' if i != 'unresolved' else '???'

with open("/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/ping.csv") as f:
    for l in f:
        if l.startswith('"",'):
            idx = [x.strip('"') for x in l.strip().split(',')[1:]]
            continue
        l = l.strip().split(',')
        ping_pre[l[0].strip('"')] = {
            idx[i]: list(set([
                tuple(  
                    parse_special(xx)
                    for xx in z.split('+')
                    if not xx.endswith('*null')
                )
                for z in l[i + 1].strip('"').split(' ')
            ]))
            for i in range(len(idx)) # non last 3 colmns
        }
    combined = ['KIR2DL23','KIR2DS35','KIR3DL1S1','KIR2DL5']
    def remove_kir(tuple_list):
        return [
            tuple(
                element.split('*')[1] if '*' in element else element.replace('KIR', '')
                for element in t
            )
            for t in tuple_list
        ]

    def extract_genes(gene_tuple):
        
        # Extract gene names by splitting each string at '*' and taking the first part
        gene_names = [gene.split('*')[0] for gene in gene_tuple]
    
        # Remove duplicates by converting to a set, then back to a list
        unique_genes = list(set(gene_names))
        if len(unique_genes) == 1 and unique_genes[0] == '???':return []
        return unique_genes

    def parse_multi(ping,t,id):
        for r in t:
            if r != '???':
                gene = r.split('*')[0]
                allele = r.split('*')[1][:3]
                if not (allele,) in ping[id][gene]:
                    ping[id][gene].append((allele,))
                
        

    ping = {}
    for id in ping_pre:
        ping[id] = {}
        for gene in ping_pre[id]:
            if not gene in combined:
                ping[id][gene] = remove_kir(ping_pre[id][gene])
            else:
                if gene == 'KIR2DL23':
                    ping[id]['KIR2DL2'] = []
                    ping[id]['KIR2DL3'] = []
                    for t in ping_pre[id][gene]:
                        
                        genes_e = extract_genes(t)
                        if len(genes_e) == 1:
                            ping[id][genes_e[0]].append(remove_kir([t])[0])
                        else:
                            
                            parse_multi(ping, t,id)
                if gene == 'KIR2DS35':
                    ping[id]['KIR2DS3'] = []
                    ping[id]['KIR2DS5'] = []
                    for t in ping_pre[id][gene]:
                        
                        genes_e = extract_genes(t)
                        if len(genes_e) == 1:
                            ping[id][genes_e[0]].append(remove_kir([t])[0])
                        else:
                            parse_multi(ping, t,id)
                if gene == 'KIR3DL1S1':
                    ping[id]['KIR3DL1'] = []
                    ping[id]['KIR3DS1'] = []
                    for t in ping_pre[id][gene]:
                        genes_e = extract_genes(t)
                        if len(genes_e) == 1:
                            print(id, t)

                            ping[id][genes_e[0]].append(remove_kir([t])[0])
                        else:

                            parse_multi(ping, t,id)

                if gene == 'KIR2DL5':
                    ping[id]['KIR2DL5A'] = []
                    ping[id]['KIR2DL5B'] = []
                    for t in ping_pre[id][gene]:
                        genes_e = extract_genes(t)
                        if len(genes_e) == 1:
                            print(id, t)

                            ping[id][genes_e[0]].append(remove_kir([t])[0])
                        else:

                            parse_multi(ping, t,id)


#%%

#%%
t1k = collections.defaultdict(lambda: collections.defaultdict(list))
# parallel --bar -j13 './run-t1k -b {}  -c kiridx/kiridx_dna_coord.fa --od _x -t 8 --preset kir-wgs -f kiridx/kiridx_dna_seq.fa >_x/{/.}.log 2>&1' ::: ~/scratch/hprc_tmp/*.bam
def read_t1k(path):
    t1k = collections.defaultdict(list)
    with open(path) as f:
        for l in f:
            l = l.split('\t')
            # gene_name num_diff_alleles allele_1 abundance_1 quality_1 allele_2 abundance_2 quality_2
            if int(l[1]) > 0 and float(l[4]) > 0:
                t1k[l[0]].append('/'.join(x.split('*')[1][:3] for x in l[2].split(',')))
            if int(l[1]) > 1 and float(l[7]) > 0:
                t1k[l[0]].append('/'.join(x.split('*')[1][:3] for x in l[5].split(',')))
    return t1k
#%% 
# base_dir = "/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/out_old"
# base_dir = '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/T1K_manuscript_evaluation/T1K_1.0.5'
base_dir = '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/out'
for subfolder in glob.glob(os.path.join(base_dir, "*.final")):
    id = os.path.basename(subfolder).split('.')[0]
    genotype_file = os.path.join(subfolder, f"T1K_{id}_genotype.tsv")
    if os.path.isfile(genotype_file):
        t1k[id] = read_t1k(genotype_file)

def read_geny(path):
    res = collections.defaultdict(list)
    with open(path) as f:
        for l in f:
            lp = l.strip().split()
            if l.startswith('[ilp] KIR'):
                g, a = lp[1], lp[2]
                res[g].append(a[:3])
    return res
results = []
parents = collections.defaultdict(list)
geny_ran = set()
for dirr in [
    # '/project/shared/inumanag-kir/aldy-kir-new-1215/aldy-kir/HPRC_ZC',
    # '/home/inumanag/scratch/kir/aldy-kir/hprcr_cn2',
    '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/results_validation'
    # '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/HPRC_window_150_new_scc_02'
    # '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/Hyper_Test_HPRC/9'
    # '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/hprc'
]:
    res = collections.defaultdict(lambda: collections.defaultdict(list))
    for fid in sorted(glob.glob(f'{dirr}/*.log')):
        if '_' in fid.split('/')[-1]: continue
        id = fid.replace('.final', '').split('.')[0].split('/')[-1].split('_')[0]
        geny_ran.add(id)
        res[id] = read_geny(fid)
    results.append(res)

for id in sorted(ground):
    if id not in geny_ran: continue
    for gi, g in enumerate(genes):
        gn1 = sorted(x for x in results[0][id][g])
        # gn2 = sorted(x for x in results[1][id][g])
        # if gn1 != gn2:
            # print(id, g, gn1, '->', gn2, [x[:3] for x in ground[id][list(genes).index(g)]])

def get_correct(ref, dat, exc='!'):
    dgrn = [x[:3] for x in sorted(ref)]
    ddat = ['' for _ in range(len(dgrn))]
    doff = []
    for x in dat:
        if x:
            x = x[:3]
            if any(x == z and ddat[(ii := zi)] == '' for zi, z in enumerate(dgrn)):
                ddat[ii] = x
            else:
                doff.append(x)

    tp, fp, fn = 0, 0, 0
    for di, d in enumerate(ddat):
        if not d:
            if doff:
                ddat[di] = exc + doff[0]
                doff = doff[1:]
                fp += 1
            else:
                fn += 1
                ddat[di] = exc
        else:
            tp += 1
    if doff:
        fp += len(doff)
        ddat.append(exc * len(doff) + ";".join(doff))
    return dgrn, ddat, [tp, fp, fn]

def summarize(ground, results, t1k, key=None):
    bad_geny, bad_t1k = collections.defaultdict(list), collections.defaultdict(list)
    stat_geny, stat_t1k = collections.defaultdict(list), collections.defaultdict(list)

    matrix = collections.defaultdict(lambda: collections.defaultdict(str))
    total = collections.defaultdict(int)
    for id in ground:
        # if id not in geny_ran: continue
        for gi, g in enumerate(genes):
            grn = sorted(x[:3] for x in ground[id][gi] if x)
            total[g] += len(grn)
            gny = sorted(x     for x in results[id][g])
            tok = t1k[id].get(g, [])

            dr, dg, statg = get_correct(grn, gny)
            bad_geny[g].append(sum(x.count('!') for x in dg))

            dt, statt = None, None
            best = 99999999
            if tok and isinstance(tok[0], tuple): # multi-solutions
                for t in tok:
                    _, ddt, sstatt = get_correct(grn, t, '#')
                    if (s := sum(x.count('#') for x in ddt)) < best:
                        best = s
                        dt, statt = ddt, sstatt
            else:
                _, dt, statt = get_correct(grn, tok, '#')
                best = sum(x.count('#') for x in dt)
            bad_t1k[g].append(best)
            stat_geny[g].append(statg)
            stat_t1k[g].append(statt)

            for i in range(max(len(dr), len(dg), len(dt))):
                matrix[id, i][g, 'ref'] = dr[i] if i < len(dr) else ''
                matrix[id, i][g, 'gen'] = dg[i] if i < len(dg) else ''
                matrix[id, i][g, 't1k'] = dt[i] if i < len(dt) else ''

    for id in sorted(ground, key=key):
        mm = max(i for ix, i in matrix if ix == id) + 1
        for i in range(mm):
            x = ['' if i else id]
            x += [ matrix[id, i][g, t] for g in genes for t in 'ref gen t1k'.split() ]
            print(','.join(x))
    print('BAD:', bg:=sum(sum(x) for x in bad_geny.values()), bt:=sum(sum(x) for x in bad_t1k.values()))

    fmt = '{:12} & {:3}   & {:3} {:3} & {:3} {:6.1%}% & {:3} {:6.1%}% & {:3} {:5.2f}   & {:3} {:3} & {:3} {:6.1%}% & {:3} {:6.1%}% & {:3} {:5.2f} \\\\'
    print('{:12} & {:>3}   & {:>7} & {:>11} & {:>11} & {:>9}   & {:>7} & {:>11} & {:>11} & {:>9} \\\\\n\\hline'.format(
        'Gene', 'TOT', 'GX', 'GP', 'GF', 'GR', 'TX', 'TP', 'TR', 'TF'
    ))
    def format_line(title, tot, geny, t1k):
        gx = sum(geny[0])
        gacc = 1 - gx / tot

        tp = sum(tp for [tp, fp, fn] in geny[1])
        fp = sum(fp for [tp, fp, fn] in geny[1])
        fn = sum(fn for [tp, fp, fn] in geny[1])
        gprc = tp / (tp + fp)
        grec = tp / (tp + fn)
        gf1  = 2 * tp / (2 * tp + fp + fn)

        tx = sum(t1k[0])
        tacc = 1 - tx / tot
        tp = sum(tp for [tp, fp, fn] in t1k[1])
        fp = sum(fp for [tp, fp, fn] in t1k[1])
        fn = sum(fn for [tp, fp, fn] in t1k[1])
        tprc = tp / (tp + fp) if tp + fp else 0.
        trec = tp / (tp + fn) if tp + fp else 0.
        tf1  = 2 * tp / (2 * tp + fp + fn)

        bf = lambda x,y,op=float.__ge__: ("\\bf",x) if op(x, y) else ("", x)
        return fmt.format(title, tot,
                          *bf(gx, tx, int.__le__), *bf(gprc, tprc), *bf(grec, trec), *bf(gf1, tf1),
                          *bf(tx, gx, int.__le__), *bf(tprc, gprc), *bf(trec, grec), *bf(tf1, gf1)).replace('%%', '\\%')
    for g in genes:
        print(format_line('\\it ' + g, total[g], (bad_geny[g], stat_geny[g]), (bad_t1k[g], stat_t1k[g])))
    print('\\hline')
    print(format_line('All', sum(total.values()),
                      ([i for gg in bad_geny.values() for i in gg], [i for gg in stat_geny.values() for i in gg]),
                      ([i for gg in bad_t1k.values() for i in gg], [i for gg in stat_t1k.values() for i in gg])))

summarize(ground, results[0], ping)
#%%