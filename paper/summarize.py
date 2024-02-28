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
import gzip
import mappy
from tqdm import tqdm
from common import bisect_left, timeit, timing
import parasail
from common import Region, Allele, Gene
from Bio import SeqIO

@timeit
def initialize(file="database.geny"):
    from common import Region, Allele, Gene
    import time

    print(time.ctime(os.path.getmtime(file)))
    with gzip.open(file, "rb") as fo:
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
with open("../hprc-ground-new.csv") as f:  # based on HPRC-kir-annotations.extra.tsv
    for l in f:
        if l.startswith('ID'):
            idx = l.strip().split(',')[1:]
            continue
        l = l.strip().split(',')
        ground[l[0].split('.')[0]] = [[z[1:] if z and z[0] == '`' else z for z in i.split(';')] for i in l[1:]]
# pi = idx.index('KIR3DP1')  patches to HPRC-kir-annotations.extra.tsv
# ground['HG00741'][pi] = ['0060103', '007?','009?']
# ground['HG01358'][pi] = ['0079991', '0030202']  #7-snp
# ground['HG02257'][pi] = ['0090102', '0030203']
# ground['HG02572'][pi] = ['0150101', '0109991']  #10+14
# ground['HG02622'][pi] = ['0030202', '0090102', '0100102']  #9+snp!
# ground['HG02630'][pi] = ['0030216', '0030216', '0119991']  #011+004
# ground['HG02717'][pi] = ['0090102', '0100102', '0030202']  #9+4
# ground['HG03540'][pi] = ['017', '0100102']

ping = {}
idx = []
def parse(i):
    i = i.split('*')[1]
    return i[:3] if i != 'unresolved' else '???'
with open("../hprc-ping-tuned.csv") as f:  # based on HPRC-kir-annotations.extra.tsv
    for l in f:
        if l.startswith(','):
            idx = l.strip().split(',')[1:]
            continue
        l = l.strip().split(',')
        ping[l[0]] = {
            idx[i]:
            list(set([
                tuple(
                    parse(xx)
                    for xx in z.split('+')
                    if not xx.endswith('*null')
                )
                for z in l[i + 1].split(' ')
            ]))
            for i in range(len(idx))
        }

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
results = []
parents = collections.defaultdict(list)
geny_ran = set()
for dirr in [
    '/project/shared/inumanag-kir/aldy-kir-new-1215/aldy-kir/HPRC_ZC',
    '/home/inumanag/scratch/kir/aldy-kir/hprcr_cn2',
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
        gn2 = sorted(x for x in results[1][id][g])
        if gn1 != gn2:
            print(id, g, gn1, '->', gn2, [x[:3] for x in ground[id][list(genes).index(g)]])

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

summarize(ground, results[1], t1k)
# summarize(ground, results[1], ping)

#%% Simulations nice
ground = {}
idx = []
with open("../sim-ground.csv") as f:
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

res = collections.defaultdict(lambda: collections.defaultdict(list))
geny_ran = set()
# for fid in sorted(glob.glob('/home/inumanag/scratch/kir/aldy-kir/hprcr/sim/*.log')):
for fid in sorted(glob.glob('/home/inumanag/scratch/kir/aldy-kir/hprcr/sim2/*.log')):
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



#%% Mendelian concordance
results = collections.defaultdict(lambda: collections.defaultdict(list))
parents = collections.defaultdict(list)
t1k = collections.defaultdict(lambda: collections.defaultdict(list))
for path in glob.glob("T1K/_x/*_genotype.tsv"):
    id = os.path.basename(path)[4:11]
    t1k[id] = read_t1k(path)

geny = collections.defaultdict(lambda: collections.defaultdict(list))
for fid in sorted(glob.glob(f'/project/shared/aldy-kir/hprcr/*.log')):
    id = fid.split('.')[0].split('/')[-1]
    if '_' in id:
        ch, id = id.split('_')
        parents[ch].append(id)

for fid in sorted(glob.glob(f'/home/inumanag/scratch/kir/aldy-kir/hprcr_cn2/*.log')):
    id = fid.split('/')[-1][:7]
    geny[id] = read_geny(fid)

def enumerate_alleles(x):
    x = collections.Counter([i[:3] for i in x])
    return [(i, j) for i, c in sorted(x.items()) for j in range(c)]
def mendel(results):
    z = 0
    for cid, [mid, fid] in parents.items():
        for g in genes:
            cg, mfg = set(enumerate_alleles(results[cid][g])), set(enumerate_alleles(results[mid][g] + results[fid][g]))
            if not (cg <= mfg):
                z += 1
                print(cid, g, 'C', cid, results[cid][g])
                print('               ', 'M', mid, results[mid][g])
                print('               ', 'F', fid, results[fid][g])
        print()
    print('bad hits', z)
mendel(geny)
mendel(t1k)

#%% Runtimes
def parse_elog(id, path):
    total = 0
    mem = 0
    fields = collections.defaultdict(float)
    with open(fid) as f:
        for l in f:
            if 'took' in l:
                f = l.split()[1]
                l = l[l.index('took'):].split()
                if f not in 'correct em em_init get_scc':
                    total += float(l[1][:-1])
                fields[f] += float(l[1][:-1])
                if f == 'ilp_solve': mem = float(l[4].replace(',',''))
    def sec2min(s):
        s=int(s)
        return f"{s//60:3}:{s%60:02d}"
    print(f"{id:40}: {sec2min(total)}: {sec2min(fields['minimap2'])} {sec2min(fields['filtering'])} {sec2min(fields['em_scc'])} {sec2min(fields['ilp_solve'])}   {mem=:6,.0f}")
for fid in sorted(glob.glob(f'/home/inumanag/scratch/kir/aldy-kir/hprcr_cn2/*.elog')):
    # if '_' in fid: continue
    id = fid.split('.')[0].split('/')[-1].split('_')[-1]
    parse_elog(id, path)
for fid in sorted(glob.glob(f'/scratch/inumanag/kir/aldy-kir/hprcr/sim/*.elog')):
    id = fid.split('/')[-1][:-5]
    if 'extract' in id: continue
    parse_elog(id, path)


#%% Problematic reasons
res = collections.defaultdict(lambda: collections.defaultdict(list))
reasons = [0, 0, 0]
for fid in sorted(glob.glob(
    # '/home/inumanag/scratch/kir/aldy-kir/hprcr_cn2/*.log'
    '/home/inumanag/scratch/kir/aldy-kir/hprcr/sim/*.log'
)):
    id = fid.split('/')[-1][:-4]
    # id = fid.replace('.final.cram', '.bam').split('.')[0].split('/')[-1].split('_')[0]
    if 'extract' in id or id not in ground: continue

    expected = ground[id]

    filter = collections.defaultdict(set)
    em = collections.defaultdict(set)
    ilp = collections.defaultdict(set)
    with open(fid) as f:
        for l in f:
            lp = l.split()
            if l.startswith('[filter] K'):
                filter[lp[1]].add(lp[2][:3])
            if l.startswith('[em] K'):
                em[lp[1]].add(lp[2][:3])
            if l.startswith('[ilp] K'):
                ilp[lp[1]].add(lp[2][:3])
    print(id, '->', '; '.join(
        g + "*" + a
        for gi, g in enumerate(genes)
        if (a := ','.join(i[:3] for i in ground[id][gi] if i))
    ))
    for gi, g in enumerate(genes):
        for ex in ground[id][gi]:
            if not ex: continue
            ex = ex[:3]
            if ex in ilp[g]: continue
            elif ex in em[g]:
                print(f'  {g:9} ilp: ', ex)
                reasons[0] += 1
            elif ex in filter[g]:
                print(f'  {g:9} em:  ', ex)
                reasons[1] += 1
            else:
                print(f'  {g:9} filt:', ex)
                reasons[2] += 1
    print()

print(reasons)

#%% Find and validate KIR annotations from the GenBank
def gb_walk(id, ost=0, st=0, en=-1, rc=False):
    from Bio import Entrez, SeqIO
    Entrez.email = "inumanag@uvic.ca"

    # Load GenBank & FASTA data
    if not os.path.exists(f"refs/{id}.gb"):
        with Entrez.efetch(db='nucleotide', id=id, rettype='gb') as h, open(f"refs/{id}.gb", "w") as fo:
            fo.write(h.read())
    gb = SeqIO.read(f"refs/{id}.gb", "genbank")
    try:
        seq = str(gb.seq)
    except:
        if not os.path.exists(f"refs/{id}.fa"):
            with Entrez.efetch(db='nucleotide', id=id, rettype='fasta') as h, open(f"refs/{id}.fa", "w") as fo:
                fo.write(h.read())
        with pysam.FastaFile(f"refs/{id}.fa") as fa:
            seq = str(fa.fetch(id))
    contigs = gb.annotations.get('contig', '-')
    if not [f for f in gb.features if f.type == 'gene' or f.type == 'CDS']:
        if not os.path.exists(f"refs/{id}_parts.gb"):
            with Entrez.efetch(db='nucleotide', id=id, rettype='gb', style='withparts') as h, open(f"refs/{id}_parts.gb", "w") as fo:
                fo.write(h.read())
        gb = SeqIO.read(f"refs/{id}_parts.gb", "genbank")
    if en == -1: en = len(seq)

    info = []
    for f in gb.features:
        if f.type == 'gene' or f.type == 'CDS':
            gn = f.qualifiers.get('gene', [''])[0]
            an = f.qualifiers.get('allele', [None])[0]
            # Look in the NOTE section
            if not gn.startswith('KIR'):
                no = f.qualifiers.get('note', [''])[0] + f.qualifiers.get('product', [''])[0]
                l = set(re.findall(r'(KIR\w+)', no))
                if len(l) == 1:
                    gn = list(l)[0]
                elif len(l) > 1:
                    print("-> {} {}", id, l)
            assert not an or gn.startswith('KIR')

            s, e = int(f.location.start), int(f.location.end)
            if e < st or s > en: continue  # not in the included region

            if rc: ns, ne = ost + (en - e), ost + (en - s)
            else: ns, ne = ost + s, ost + e
            loci = (s, e), (ns, ne)
            if gn.startswith('KIR'):
                if an: an = an.split('*')[1]
                if gn == 'KIR3DX1': gn = {'KIR3DP1', 'KIR3DL1', 'KIR3DS1'}
                elif gn == 'KIR2DL5': gn = {'KIR2DL5A', 'KIR2DL5B'}
                elif gn == 'KIR3DL1S1': gn = {'KIR3DL1'}
                else: gn = {gn}
                info.append((id, gn, {an} if an else set(), *loci))

    if contigs != '-':
        ost = 0
        for m in re.findall(r'((complement\()?([0-9A-Z_]{2}\d+\.\d):(\d+)..(\d+))|(gap\((\d+)\))', contigs):
            if m[5]: # gap
                ost += int(m[6])
            else:
                m = m[1:]
                rc, s, e = m[0] != '', int(m[2]), int(m[3])
                info += gb_walk(m[1], ost, s, e, rc)[0]
                ost += e - s
    return info, seq

if __name__ == "__main__":
    gb_regions = {}
    gb_seqs = {}
    for id in ids:
        kir_regions, gb_seq = gb_walk(id)
        new_regions = []
        overlap = lambda a, b: max(0, min(a[1], b[1]) - max(a[0], b[0]))
        kir_regions = sorted(kir_regions, key=lambda x: x[4])
        for gi, (_, gns, ans, _, (nst, nen)) in enumerate(kir_regions):
            if new_regions and (new_regions[-1][0] & gns) and overlap((nst, nen), new_regions[-1][-1]):
                new_regions[-1][0].update(gns)
                new_regions[-1][-1][0] = min(new_regions[-1][-1][0], nst)
                new_regions[-1][-1][1] = max(new_regions[-1][-1][1], nen)
                if ans: new_regions[-1][1].update(ans)
            else:
                new_regions.append((gns, ans, [nst, nen]))
            if new_regions[-1][-1][0] < 0 or new_regions[-1][-1][1] < 0 or new_regions[-1][-1][0] >= new_regions[-1][-1][1]:
                print(id, gns, ans, new_regions[-1][-1])
                # break
        gb_regions[id] = new_regions
        gb_seqs[id] = gb_seq
    print('regions done')
with open("_x", "w") as fo:
    for id in gb_seqs:
        print(id, gb_seqs[id], file=fo)
        for g, _, [s, e] in gb_regions[id]:
            print(id, '|'.join(g), s, e, file=fo)

def align(rs, s):
    aln = parasail.sg_qx_trace_scan_32(rs, s, PGO, PGE, PMM)

    seqr = mappy.revcomp(s)
    alnr = parasail.sg_qx_trace_scan_32(rs, seqr, PGO, PGE, PMM)
    rc = False
    if alnr.score > aln.score:
        aln, s = alnr, seqr
        rc = True

    cigar = get_cigar(aln)
    cstr = aln.cigar.decode.decode()
    return s, cigar, cstr, aln.score, rc

def find_closest_match(data):
    idx, gns, als, st, ed = data

    if st >= ed:
        print(f"[error] {idx} {st}-{ed} seq empty")
        return idx, []

    best = {}
    for gn in gns:
        if gn not in genes:
            if gn.startswith('KIR'): print(f"[error] {idx} {gn} unknown")
            continue

        wildtype_seq = genes[gn].wildtype.seq
        seq = gb_seqs[id][st:ed]
        ss, cigar, cigar_str, _, rc = align(wildtype_seq, seq)
        print(gns, ed-st, len(wildtype_seq))
        if cigar[0][1] == 'I' or cigar[-1][1] == 'I':
            pass
            # if rc: cigar = cigar[::-1]
            # if cigar[0][1] == 'I':
            #     st = max(0, st - cigar[0][0])
            # if len(cigar) > 1 and cigar[-1][1] == 'I':
            #     ed += cigar[-1][0]
            # seq = gb_seqs[id][st:ed]
            # seq, cigar, cigar_str, _, _ = align(wildtype_seq, seq)
        else:
            seq = ss

        ops = []
        wi, ri = 0, 0
        for sz, op in cigar:
            if op == "D":
                ops.append((wi, "ins" + seq[ri : ri + sz])); ri += sz
            elif op == "I":
                ops.append((wi, "del" + wildtype_seq[wi : wi + sz])); wi += sz
            else:
                if op == "X":
                    for i in range(sz):
                        ops.append((wi + i, f"{wildtype_seq[wi + i]}>{seq[ri + i]}"))
                ri += sz; wi += sz
        ops = set(ops)

        func_pos = {m[0]: m for m in genes[gn].functional}
        for a in genes[gn].alleles.values():
            common = a.ops & ops
            diff = a.ops ^ ops
            miss_func = set(m for m in a.ops if m not in ops if m in genes[gn].functional)
            for i in miss_func:  # KIR3DP1 issue...
                if i[0] == 335 and i[1].startswith('del'):
                    miss_func.remove(i)
                    break
            new_func = sorted(set(f'{m[0]}.{m[1]}#{func_pos[m[0]][0]}.{func_pos[m[0]][1]}'
                              for m in ops if m[0] in func_pos if m not in a.func))
            check_anyways = any(a.name.startswith(ax) for ax in als)
            if not miss_func or check_anyways:
                new = [check_anyways, gn, a.name, st, ed, diff, common, new_func, miss_func]
                key = gn, minor_to_major[gn, a.name]
                if key not in best or len(best[key][5]) > len(diff):
                    best[key] = new
    best = sorted(best.values(), key=lambda x: (not x[0], len(x[-1]), len(x[-2]), len(x[5])))

    for bi, [c, _, _, _, _, _, _, n, m] in enumerate(best):
        if not c and (m or n):
            extra = int(bi == 0)
            best = best[:bi + extra]
            break
    for b in best:
        [_, gn, an, _, _, _, _, _, _] = b
        _, _, cigar_str, new_score, _ = align(genes[gn].alleles[an].seq, seq)
        b.append(new_score)
        b.append(cigar_str)
    return idx, best

with mp.Pool(64) as pool, timing("alignment"):
    data = [
        ((id, gi), gns, als, st, ed)
        for id, rs in gb_regions.items()
        if id.startswith('NT_113949')
        for gi, [gns, als, [st, ed]] in enumerate(rs)
    ]
    res = pool.map(find_closest_match, data)

prev = None
for (id, gi), best in sorted(res, key=lambda x: (x[0][0], str(gb_regions[x[0][0]][x[0][1]][0]))):
    if id != prev:
        print('\n')
        prev = id
    [gns, als, _] = gb_regions[id][gi]

    try:
        old = {gn: '|'.join(i[:3] for i in ground[id][idx.index(gn)]) for gn in gns}
    except:
        old = {}
    print(f'-', f"{id:9} {gns} {als}: {old}")

    h = ' ' * 1
    found_any = False
    for [match_name, gn, an, st, ed, diff, common, new_func, miss_func, aln_score, aln_cigar] in best:
        match = not new_func and not miss_func
        old_match = an[:3] in old[gn]
        if old_match:
            found_any = True
        # if found_any:
            # break
        print(h, '✅' if match_name and match else '⚠️',
                f'{gn:10} {an:10} ({minor_to_major[gn, an][:3]}) {len(common)=:3} {len(diff)=:3} {st:10,}-{ed:10,} {aln_score=:10} ... ',
                '😄' if old_match else '😡',
                aln_cigar)
        if not match:
            # _, _, ass, new_score = align(genes[gn].alleles[name].seq, seq)
            new_func = sorted(new_func)
            miss_func = sorted(miss_func)
            print(h, '  ', f"{new_func=} {miss_func=}")

#%% Pad GenBank sequences
def aln(id, gn, qs):
    rs = genes[gn].wildtype.seq
    aln = parasail.sg_trace_scan_32(qs, rs, PGO, PGE, PMM)  # gaps at qs not penalized
    rs_rc = mappy.revcomp(rs)
    aln_rc = parasail.sg_trace_scan_32(qs, rs_rc, PGO, PGE, PMM)
    rc = False
    r_st, r_en = 0, len(rs)
    q_st, q_en = 0, len(qs)
    # D: deletion in WILDTYPE[R] / insert in GENBANK[Q]
    # I: deletion in GENBANK[Q] / insert in WILDTYPE[R]
    rc = False
    if aln_rc.score > aln.score:
        rc = True
        cigar = get_cigar(aln_rc)
        # print('   -', aln_rc.cigar.decode.decode())
        i = 0
        while cigar[i][1] in 'ID':
            if cigar[i][1] == 'D':
                r_en -= cigar[i][0]
            else:
                q_st += cigar[i][0]
            i += 1
        j = -1
        while cigar[j][1] in 'ID':
            if cigar[j][1] == 'D':
                r_st += cigar[j][0]
            else:
                q_en -= cigar[j][0]
            j -= 1
    else:
        cigar = get_cigar(aln)
        # print('   +', aln.cigar.decode.decode())
        i = 0
        while cigar[i][1] in 'ID':
            if cigar[i][1] == 'D':
                r_st += cigar[i][0]
            else:
                q_st += cigar[i][0]
            i += 1
        i = -1
        while cigar[i][1] in 'ID':
            if cigar[i][1] == 'D':
                r_en -= cigar[i][0]
            else:
                q_en -= cigar[i][0]
            i -= 1
    rg1 = genes[gn].wildtype.region(r_st)[0]
    rg2 = genes[gn].wildtype.region(r_en-1)[0]
    if rg1 and rg2 and (not rg1.name.startswith('utr') or not rg2.name.startswith('utr')) and r_st == 0:
        # missing functional parts
        print (f'ID {id}\n   GB: {q_st}-{q_en} {q_en-q_st}/{len(qs)}' +
               f' \n   WT: {r_st}-{r_en} {r_en-r_st}/{len(rs)} {"-" if rc else "+"} {rg1.name} {rg2.name}')
        if rc: return mappy.revcomp(rs[r_en:]), mappy.revcomp(rs[:r_st]), q_st, q_en, rc
        else: return rs[:r_st], rs[r_en:], q_st, q_en, rc
    return None

gb_seqs_new = {i: j for i, j in gb_seqs.items()}
seqs = gb_seqs
for id, r in gb_regions.items():
    mini = min([(i, rx[2]) for i, rx in enumerate(r)], key=lambda x: x[1][0])
    maxi = max([(i, rx[2]) for i, rx in enumerate(r)], key=lambda x: x[1][1])
    pre, post = '', ''
    if mini[1][0] < 30_000 - (mini[1][1] - mini[1][0]):
        rr = r[mini[0]]
        assert len(rr[0]) == 1
        qs = seqs[id][:mini[1][1] + 20_000]
        if (ar := aln(id, list(rr[0])[0], qs)):
            pr, po, qs, qe, rc = ar
            if qs == 0 and pr:
                pre = pr
                print(f'   {id} pre OK: {len(pre)}')
    if len(seqs[id]) - maxi[1][1] < 30_000 - (maxi[1][1] - maxi[1][0]):
        rr = r[maxi[0]]
        assert len(rr[0]) == 1
        qs = seqs[id][max(0, maxi[1][0] - 20_000):]
        if not qs:
            pass # print(maxi, len(seqs[id]))
        elif (ar := aln(id, list(rr[0])[0], qs)):
            pr, po, _, qe, rc = ar
            if qe == len(qs) and po:
                post = po
                print(f'   {id} post OK: {len(post)}')
    gb_seqs_new[id] = pre + gb_seqs[id] + post


#%% Exon mapping [TODO]
import mappy
id = 'NT_113949.2'
r_seq = gb_seqs[id]
idx = mappy.Aligner(seq=r_seq, preset='sr', k=10)

st = 0
for r in genes['KIR3DL2'].wildtype.regions.values():
    if r.name.startswith('e'):
        for h in idx.map(r.seq):
            print(f'{r.name}: hit {h.r_st}-{h.r_en}/{h.strand}')
            if h.strand != 1: continue

            ops = []
            r_start = h.r_st
            q_start = h.q_st  # q: read, r: reference
            for size, op in h.cigar:
                if op == 2:  # Deletion
                    ops.append((st + q_start, "ins" + r.seq[q_start:q_start + size]))
                    q_start += size
                elif op == 1:  # Insertion
                    ops.append((st + q_start, "del" + r_seq[r_start:r_start + size]))
                    r_start += size
                elif op == 4:  # Soft-clip
                    r_start += size
                elif op in [0, 7, 8]:  # M, X and =
                    for i in range(size):
                        if r.seq[q_start + i] != r_seq[r_start + i]:
                            ops.append((st + q_start + i, r.seq[q_start + i], r_seq[r_start + i]))
                    q_start += size
                    r_start += size
            print(' ', ops)
    st += len(r.seq)

#%% Simulations
if __name__ == "__main__":
    with pysam.FastaFile("/project/shared/inumanag-kir/refs/chr19.fa") as fa:
        seq = str(fa.fetch("chr19"))
        off = 1000
        lo = seq[54724235-off:54724235]
        ro = seq[54867216:54867216+off]
    ranges = {}
    for id in ids:
        # if id != 'NT_187636.1': continue
        print(id)
        ranges[id] = {}
        # with pysam.FastaFile(f"refs/{id}.fa") as fa:
        seq = lo + gb_seqs_new[id] + ro
        with open(f"sim_new/{id}.fa", "w") as fo, open(f"sim_new/{id}.sam", "w") as fs:
            print('@SQ', f'SN:{id}', f'LN:{len(seq)}', sep='\t', file=fs)
            rl = 100
            dist = 100
            for i in range(0, len(seq) - dist - 2 * rl, 10):
                ovl = {ranges[id].get(j, '') for j in range(i, i+rl+dist+rl)}
                ovl = '_'.join(z for z in sorted(ovl) if z)
                print(f">r_{i}/1 {i}_{ovl}", file=fo)
                print((s := seq[i:i + rl].upper()), file=fo)
                print(f">r_{i}/2 {i + rl + dist}_{ovl}", file=fo)
                print((s := mappy.revcomp(seq[i + rl + dist:i + rl + dist + rl].upper())), file=fo)

                print(f"r_{i}", 0,  id, i+1,         250, '100M', '=', i+rl+dist+1, -1, s, '*', sep='\t', file=fs)
                print(f"r_{i}", 16, id, i+rl+dist+1, 250, '100M', '=', i+1,         -1, s, '*', sep='\t', file=fs)

    pairs = [   ('NT_113949.2', 'NT_187645.1'),
                ('NT_187637.1', 'NT_187693.1'),
                ('NT_187638.1', 'NT_187645.1'),
                ('NT_187639.1', 'NT_187644.1'),
                ('NT_187641.1', 'NW_016107310.1'),
                ('NT_187670.1', 'NT_187677.1'),
                ('NT_187673.1', 'NW_016107306.1'),
                ('NT_187676.1', 'NT_187670.1'),
                ('NT_187676.1', 'NW_016107303.1'),
                ('NT_187677.1', 'NT_187686.1'),
                ('NT_187683.1', 'NT_187639.1'),
                ('NT_187683.1', 'NW_003571058.2'),
                ('NW_003571054.1', 'NT_187645.1'),
                ('NW_003571055.2', 'NT_187671.1'),
                ('NW_003571056.2', 'NT_187685.1'),
                ('NW_003571056.2', 'NW_003571055.2'),
                ('NW_016107304.1', 'NT_187684.1'),
                ('NW_016107307.1', 'NT_113949.2'),
                ('NW_016107307.1', 'NT_187676.1'),
                ('NW_016107313.1', 'NT_187637.1'),
                ('NW_016107314.1', 'NW_016107301.1'),  ]
    for p1, p2 in pairs:
        print(p1, p2)
        i = 0
        with open(f"sim_new/{p1}_{p2}.fa", "w") as fo:
            with open(f"sim_new/{p2}.fa") as f:
                for li, l in enumerate(f):
                    print(l.strip(), file=fo)
            with open(f"sim_new/{p1}.fa") as f:
                for li, l in enumerate(f):
                    if li % 2 == 1: print(l.strip(), file=fo)
                    else: print(l.replace(">r_", ">q_").strip(), file=fo)

#%%
for i in glob.glob('/home/inumanag/scratch/kir/aldy-kir/hprcr/sim/*.log'):
    if '_N' in i and 'extract' not in i:
        print(i.split('/')[-1][:-4])



#%% Check Mike's data

def check(sample, gn, an, seq):
    wildtype_seq = genes[gn].wildtype.seq
    seq, cigar, cigar_str, aln_score = align(wildtype_seq, seq)

    # Calculate changes w.r.t. wildtype allele
    ops = []
    wi, ri = 0, 0
    for sz, op in cigar:
        # Calculate differences between the wildtype and the current allele
        if op == "D":
            ops.append((wi, "ins" + seq[ri : ri + sz]))
            ri += sz
        elif op == "I":
            ops.append((wi, "del" + wildtype_seq[wi : wi + sz]))
            wi += sz
        else:
            if op == "X":
                for i in range(sz):
                    ops.append((wi + i, f"{wildtype_seq[wi + i]}>{seq[ri + i]}"))
            ri += sz
            wi += sz

    ops = set(ops)
    best = []
    func_pos = {m[0]: m for m in genes[gn].functional}  # position of functional mutat
    real = None
    for a in genes[gn].alleles.values():
        common = len(a.ops & ops)
        diff = len(a.ops ^ ops)
        miss_func = set(m for m in a.ops if m not in ops if m in genes[gn].functional)
        for i in miss_func:
            if i[0] == 335 and i[1].startswith('del'):
                miss_func.remove(i)
                break
        new_func = sorted(set(f'{m[0]}.{m[1]}#{func_pos[m[0]][0]}.{func_pos[m[0]][1]}'
                              for m in ops if m[0] in func_pos if m not in a.func))
        if not miss_func:
            best.append((a.name, new_func, common, diff, []))
        # if a.name.startswith('003') and not real:
        if a.name == an:
            real = (a.name, new_func, common, diff, miss_func)
    # we prefer those W/O new_func, max. number of common muts, and min number of missing muts
    best.sort(key=lambda x: (len(x[1]), -x[2], x[3]))

    print(f'-', f"{sample:9} {gn:8} {an:10}", '=>', f"{aln_score=:10} {cigar_str=}")
    def check(name, new_func, common, diff, miss_func, force=False):
        match = not new_func and not miss_func
        if not force and not match: return False, False
        match_name = name.startswith(an[:3])
        print(' '*20, '✅' if match_name and match else '⚠️',
              f'{name:10} ({minor_to_major[gn, name][:3]}) {common=:3} {diff=:3} ... ')
        if not match_name or force:
            _, _, ass, new_score = align(genes[gn].alleles[name].seq, seq)
            new_func = sorted(new_func)
            miss_func = sorted(miss_func)
            print(' '*20, '  ', f"{new_score=} {new_func=} {miss_func=} {ass=}")
        return match, match_name

    for b in best:
        m, mn = check(*b,force=True)
        # if mn: real = False
    if real:
        print(' '*20, '---------')
        check(*real, force=True)
        # an = '' if an == '???' else an
        # if name.startswith(an): return name + ("" if not new_func else " ?")
        # else: return f'{name}?{an}'

with open('HPRC-kir-annotations.extra.tsv') as f:
    pa = ''
    for l in f:
        asm, ctg, _, ga,  st, ed, _, _, seq, *_ = l.strip().split('\t')
        if '*' not in ga: continue
        g, a = ga.split('*')
        if g == 'KIR3DP1' and 'HG03540' in asm:
            check(asm[:7], g, a, seq)
        # if pa != asm: print()
        pa = asm





#%%
if __name__ == "__main__":
    ground = {}
    with open("../_exp/ground.csv") as f:
        for l in f:
            if l.startswith('ID'): continue
            l = l.strip().split(',')
            ground[l[0]] = l[1:]

ids = [ 'NT_113949.2', 'NT_187636.1', 'NT_187637.1', 'NT_187638.1', 'NT_187639.1', 'NT_187640.1', 'NT_187641.1',
        'NT_187642.1', 'NT_187643.1', 'NT_187644.1', 'NT_187645.1', 'NT_187668.1', 'NT_187669.1', 'NT_187670.1',
        'NT_187671.1', 'NT_187672.1', 'NT_187673.1', 'NT_187674.1', 'NT_187675.1', 'NT_187676.1', 'NT_187677.1',
        'NT_187683.1', 'NT_187684.1', 'NT_187685.1', 'NT_187686.1', 'NT_187687.1', 'NT_187693.1',
        'NW_003571054.1', 'NW_003571055.2', 'NW_003571056.2', 'NW_003571057.2', 'NW_003571058.2', 'NW_003571059.2',
        'NW_003571060.1', 'NW_003571061.2', 'NW_016107300.1', 'NW_016107301.1', 'NW_016107302.1', 'NW_016107303.1',
        'NW_016107304.1', 'NW_016107305.1', 'NW_016107306.1', 'NW_016107307.1', 'NW_016107308.1', 'NW_016107309.1',
        'NW_016107310.1', 'NW_016107311.1', 'NW_016107312.1', 'NW_016107313.1', 'NW_016107314.1']

overlap = lambda a, b: max(0, min(a[1], b[1]) - max(a[0], b[0]))
regions = collections.defaultdict(lambda: collections.defaultdict(list))
seqs = {}
def deep_search(id, oid, depth=0, ost=0, st=0, en=-1, rc=False):
    if not os.path.exists(f"../refs/{id}.gb"):
        with Entrez.efetch(db='nucleotide', id=id, rettype='gb') as h, open(f"../refs/{id}.gb", "w") as fo:
            fo.write(h.read())
    gb = SeqIO.read(f"../refs/{id}.gb", "genbank")
    try:
        seqs[id] = seq = str(gb.seq)
    except:
        with pysam.FastaFile(f"../refs/{id}.fa") as fa:
            seqs[id] = seq = str(fa.fetch(id))
    if en == -1: en=len(seq)
    def coo(f):
        s, e = int(f.location.start), int(f.location.end)
        if rc: ns, ne = en - e, en - s
        else: ns, ne = st + s, st + e
        return (s, e), (ost + ns, ost + ne)
    info = []
    gbx = gb
    if not [f for f in gb.features if f.type == 'gene' or f.type == 'CDS']:
        if not os.path.exists(f"../refs/{id}_parts.gb"):
            with Entrez.efetch(db='nucleotide', id=id, rettype='gb', style='withparts') as h, open(f"../refs/{id}_parts.gb", "w") as fo:
                fo.write(h.read())
        gbx = SeqIO.read(f"../refs/{id}_parts.gb", "genbank")
    for f in gbx.features:
        if f.type == 'gene' or f.type == 'CDS':
            gn = f.qualifiers.get('gene', [''])[0]
            an = f.qualifiers.get('allele', [''])[0]
            if not gn.startswith('KIR') and an: raise 'wooo'
            if gn.startswith('KIR') and an: info.append((an, coo(f)))
            elif gn.startswith('KIR'): info.append((gn + '*???', coo(f)))
    if info:
        info = sorted(set(info))
        results = collections.defaultdict(list)
        for gi, (ga, ((st, en), (nst, nen))) in enumerate(info):
            gn, an = ga.split('*')
            regions[oid][gn].append((nst, nen))
            regions[id][gn].append((st, en))
    if (c := gb.annotations.get('contig', '')) == '': return
    ost = 0
    for m in re.findall(r'((complement\()?([0-9A-Z_]{2}\d+\.\d):(\d+)..(\d+))|(gap\((\d+)\))', c):
        if m[5]:
            ost += int(m[6])
        else:
            m = m[1:]
            rc, s, e = m[0] != '', int(m[2]), int(m[3])
            l = e-s
            deep_search(m[1], oid, depth+2, ost, s, e, rc)
            ost += l

if __name__ == "__main__":
    for id in ids:
        deep_search(id, id)
    for iv in regions.values():
        for g in iv:
            l = sorted(iv[g])
            nl = []
            for i in l:
                if nl and overlap(nl[-1], i): nl[-1] = min(i[0], nl[-1][0]), max(i[1], nl[-1][1])
                else: nl.append(i)
            iv[g] = nl

#%%
def check_single(data):
    gb, a = data

    db = a.seq #['0030101'].seq
    if gb == "": return a.name, -1, ['ERROR']
    aln = parasail.sg_qx_trace_scan_32(db, gb, PGO, PGE, PMM)
    if aln.score < -1000:
        seqr = mappy.revcomp(gb)
        alnr = parasail.sg_qx_trace_scan_32(db, seqr, PGO, PGE, PMM)
        if alnr.score > aln.score: aln, gb = alnr, seqr
    cigar = get_cigar(aln)
    ops = []
    wi, ri = 0, 0
    for sz, op in cigar:
        # Calculate differences between the wildtype and the current allele
        if op == "D":
            ops.append((wi, "ins" + gb[ri : ri + sz]))
            ri += sz
        elif op == "I":
            ops.append((wi, "del" + db[wi : wi + sz]))
            wi += sz
        else:
            if op == "X":
                for i in range(sz):
                    ops.append((wi + i, f"{db[wi + i]}>{gb[ri + i]}"))
            ri += sz
            wi += sz
    if ops[0][1].startswith('del'): ops = ops[1:]
    if ops[-1][1].startswith('del'): ops = ops[:-1]
    return a.name, aln.score, [(pos, op if len(op) < 20 else f'{op[:20]}...',
                                (*a.mutmap[pos][1], *(['🔴'] if a.mutmap[pos][1] in a.gene.functional else [])))
                               for pos, op in ops if pos in a.mutmap]
def check(id, g):
    for st, ed in regions[id][g]:
        gb = seqs[id][st:ed]
        with mp.Pool(32) as pool:
            res = pool.map(check_single, [(gb, a) for a in genes[g].alleles.values()])
            max_s = max(res, key=lambda x: x[1])
            print(id, *max_s)

if __name__ == "__main__":
    ground = {}
    with open("../_exp/ground.csv") as f:
        for l in f:
            if l.startswith('ID'): continue
            l = l.strip().split(',')
            ground[l[0]] = l[1:]
    for id in ground:
        if id != "NT_187683.1": continue
        for gi, g in enumerate(genes):
            if '?' in ground[id][gi]:
                print('->', id, g, ground[id][gi])
                check(id, g)

#%%


#%%
