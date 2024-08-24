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
#%%
@timeit
def initialize(file="kir.pickle"):
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
for path in glob.glob("../T1K/_x/*_genotype.tsv"):
    id = os.path.basename(path)[4:11]
    t1k[id] = read_t1k(path)

def read_geny(path):
    res = collections.defaultdict(list)
    with open(path) as f:
        for l in f:
            lp = l.strip().split()
            if l.startswith('[ilp] KIR'):
                g, a = lp[1], lp[2]
                res[g].append(a[:3])
    return res
# results = []
# parents = collections.defaultdict(list)
# geny_ran = set()
# for dirr in [
#     '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/hprc_ER_5'
#     # '/home/inumanag/scratch/kir/aldy-kir/hprcr_cn2',
# ]:
#     res = collections.defaultdict(lambda: collections.defaultdict(list))
#     for fid in sorted(glob.glob(f'{dirr}/*.log')):
#         if '_' in fid.split('/')[-1]: continue
#         id = fid.replace('.final', '').split('.')[0].split('/')[-1].split('_')[0]
#         geny_ran.add(id)
#         res[id] = read_geny(fid)
#     results.append(res)

# for id in sorted(ground):
#     if id not in geny_ran: continue
#     for gi, g in enumerate(genes):
#         gn1 = sorted(x for x in results[0][id][g])
#         gn2 = sorted(x for x in results[1][id][g])
#         if gn1 != gn2:
#             print(id, g, gn1, '->', gn2, [x[:3] for x in ground[id][list(genes).index(g)]])

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

# summarize(ground, results[1], t1k)
# summarize(ground, results[1], ping)

#%% Simulations nice
ground = {}
idx = []
with open("/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/sim-ground.csv") as f:
    for l in f:
        if l.startswith('ID'):
            idx = l.strip().split(',')[1:]
            continue
        l = l.strip().split(',')
        ground[l[0]] = [[z[1:] if z and z[0] == '`' else z for z in i.split(';')] for i in l[1:]]

for path in sorted(glob.glob("new_sim/*.tsv")):
    ground_new = collections.defaultdict(list)
    id = path.split('/')[-1][:-4]
    with open(path) as f:
        pa = ''
        # sample	haplotype	gene	allele	start	end	new variants	new functional_variants	sequence
        for l in f:
            if l.startswith('sample'): continue
            _, _, _, ga, st, ed, _, _, seq, *_ = l.strip().split('\t')
            if '*' not in ga: continue
            g, a = ga.split('*')
            ground_new[list(genes).index(g)].append(a[:3])
    for gi, g in enumerate(genes):
        if g == 'KIR3DP1': continue
        ng = sorted(i[:3] for i in ground[id][gi] if i)
        ground_new[gi].sort()
        if ng != ground_new[gi]:
            print(id, g, ng, '->', ground_new[gi])


#%%
t1k = collections.defaultdict(lambda: collections.defaultdict(list))
# parallel --bar -j24 'bash run_t1k_sim.sh {} >_sim/{/.}.spl.log 2>&1' ::: /project/shared/inumanag-kir/hprc/sim/*.?.fa
for path in glob.glob("../T1K/_sim2/*.fq/*_genotype.tsv"):
    id = path.split('/')[-2][:-6]
    t1k[id] = read_t1k(path)
directory = '/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/sim_ER_5'
for filename in os.listdir(directory):
    if '.fa' in filename:
        # Create the new filename by removing '.fa'
        new_filename = filename.replace('.fa', '')
        
        # Construct full file paths
        old_file = os.path.join(directory, filename)
        new_file = os.path.join(directory, new_filename)
        
        # Rename the file
        os.rename(old_file, new_file)
        print(f"Renamed: {filename} -> {new_filename}")

res = collections.defaultdict(lambda: collections.defaultdict(list))
geny_ran = set()
# for fid in sorted(glob.glob('/home/inumanag/scratch/kir/aldy-kir/hprcr/sim/*.log')):
for fid in sorted(glob.glob('/home/qinghuiz/inumanag-kir/aldy-kir-fast-geny/sim_ER_5/*.log')):
    id = fid.split('/')[-1][:-4]
    if 'extract' in id: continue
    geny_ran.add(id)
    res[id] = read_geny(fid)

for s in res:
    if '_N' not in s: continue
    l, r = s.split('_N')
    r = 'N' + r
    ground[s] = [list(filter(lambda x: x, ground[l][i] + ground[r][i])) for i in range(len(ground[l]))]

summarize(ground, res, t1k, key=lambda x: ('_N' in x, x))


#%%