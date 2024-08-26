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
from common import bisect_left, timeit, timing
import parasail
from Bio import Entrez
Entrez.email = "inumanag@uvic.ca"

from common import Region, Allele, Gene
from Bio import SeqIO

# Parasail parameters
PMM, PGO, PGE = parasail.matrix_create("ACGT", 2, -8), 12, 2

# Parse CIGAR from Parasail alignment
def get_cigar(a):
    l = re.split(r'([A-Z=])', a.cigar.decode.decode())[:-1]
    return [(int(l[i]), l[i+1]) for i in range(0, len(l), 2)]

@timeit
def initialize(file="../kir.pickle"):
    from common import Region, Allele, Gene
    import time

    print(time.ctime(os.path.getmtime(file)))
    with open(file, "rb") as fo:
        (genes, _, _, minor_to_major) = pickle.load(fo)
    return genes, minor_to_major

if __name__ == "__main__":
    genes, minor_to_major = initialize()


#%% Simulations: search ground-truth alleles
ids = [ 'NT_113949.2', 'NT_187636.1', 'NT_187637.1', 'NT_187638.1', 'NT_187639.1', 'NT_187640.1', 'NT_187641.1',
        'NT_187642.1', 'NT_187643.1', 'NT_187644.1', 'NT_187645.1', 'NT_187668.1', 'NT_187669.1', 'NT_187670.1',
        'NT_187671.1', 'NT_187672.1', 'NT_187673.1', 'NT_187674.1', 'NT_187675.1', 'NT_187676.1', 'NT_187677.1',
        'NT_187683.1', 'NT_187684.1', 'NT_187685.1', 'NT_187686.1', 'NT_187687.1', 'NT_187693.1',
        'NW_003571054.1', 'NW_003571055.2', 'NW_003571056.2', 'NW_003571057.2', 'NW_003571058.2', 'NW_003571059.2',
        'NW_003571060.1', 'NW_003571061.2', 'NW_016107300.1', 'NW_016107301.1', 'NW_016107302.1', 'NW_016107303.1',
        'NW_016107304.1', 'NW_016107305.1', 'NW_016107306.1', 'NW_016107307.1', 'NW_016107308.1', 'NW_016107309.1',
        'NW_016107310.1', 'NW_016107311.1', 'NW_016107312.1', 'NW_016107313.1', 'NW_016107314.1']

# Find star-alleles of similations & create simple simulations
def find_closest(gn, an, seq, st=0, en=0, depth=0):
    """
    gn: gene name
    an: allele name (from annotation)
    st, en: location of annotation within seq (for debug purposes only)
    """

    if 'DP1' in gn and gn.endswith('DP1'): gn = gn[:7]  # fix for bad DP1 annotations
    if gn not in genes: return  # not a KIR gene

    rseq = genes[gn].wildtype.seq  # take wildtype
    aln = parasail.sg_qx_trace_scan_32(rseq, seq, PGO, PGE, PMM)  # do semiglobal alignment
    if aln.score < -1000:  # if score is too bad, try reverse complement if it is better
        seqr = mappy.revcomp(seq)
        alnr = parasail.sg_qx_trace_scan_32(rseq, seqr, PGO, PGE, PMM)
        if alnr.score > aln.score:
            aln, seq = alnr, seqr

    cigar = get_cigar(aln)
    ops = []
    wi, ri = 0, 0
    for sz, op in cigar:
        # Calculate differences between the wildtype and the current allele
        if op == "D":
            ops.append((wi, "ins" + seq[ri : ri + sz]))
            ri += sz
        elif op == "I":
            ops.append((wi, "del" + rseq[wi : wi + sz]))
            wi += sz
        else:
            if op == "X":
                for i in range(sz):
                    ops.append((wi + i, f"{rseq[wi + i]}>{seq[ri + i]}"))
            ri += sz
            wi += sz
    ops = set(ops)
    best = []
    func_muts = {m[0]: m for m in genes[gn].functional}
    for a in genes[gn].alleles.values():
        common = len(a.ops & ops)  # common mutations
        diff = len(a.ops ^ ops)  # unique mutations
        miss_func = set(m for m in a.ops if m not in ops if m in genes[gn].functional)  # missing functional mutations
        new_func = sorted(set((m, func_muts[m[0]]) for m in ops if m[0] in func_muts if m not in a.func))  # new functional mutations
        if not miss_func:  # consider all alleles whose functional are NOT missed
            best.append((a.name, new_func, common, diff))
    best.sort(key=lambda x: (len(x[1]), -x[2], x[3]))  # select best based on lowest number of new mutations

    if not best:  # no hit found!
        raise 'no hit found'
    for allele, new_func_muts, common_muts, diff_muts in best[:1]:
        print('  ' * depth, f'-', f"{gn:8} {an:10} ({st:8}-{en:8})", '=>',
              f'{allele:10} {common_muts:3} {diff_muts:3} {aln.score:6}',
              '' if allele.startswith(an) and not new_func_muts else f'CHECK! {sorted(new_func_muts)}')
        an = '' if an == '???' else an
        if allele.startswith(an): return allele + ("" if not new_func_muts else " ?")
        else: return f'{allele}?{an}'

def deep_search(id, depth=0, ost=0, st=0, en=-1, rc=False, ranges={}):
    """
    id: GenBank ID
    depth: depth of recursion (some GenBank records reference others and we need to analyze them)
    ost: start of the original sequence
    st: start of the sequence
    en: end of the sequence
    rc: is sequence reference-complemented?
    ranges: locations that each KIR gene covers (for debugging)
    """

    # Load GenBank data
    if not os.path.exists(f"../refs/{id}.gb"):
        with Entrez.efetch(db='nucleotide', id=id, rettype='gb') as h, open(f"../refs/{id}.gb", "w") as fo:
            fo.write(h.read())

    # Load sequence from GenBank or from FASTA (some GB files have no seq field)
    gb = SeqIO.read(f"../refs/{id}.gb", "genbank")
    try:
        seq = str(gb.seq)
    except:
        with pysam.FastaFile(f"../refs/{id}.fa") as fa:
            seq = str(fa.fetch(id))

    if en == -1: en=len(seq)
    print('  ' * depth, f'{id}: {st:,}-{en:,} ({en-st=:,}; {rc=}; {ost=:,})')

    def location(f):
        s, e = int(f.location.start), int(f.location.end)
        if rc: ns, ne = en - e, en - s
        else: ns, ne = st + s, st + e
        return (s, e), (ns, ne)  # return locations on the current sub-sequence and the original sequence
    info = []
    # If GenBank file has no gene information, load extended withparts file and check if it has it
    gbx = gb
    if not [f for f in gb.features if f.type == 'gene' or f.type == 'CDS']:
        if not os.path.exists(f"../refs/{id}_parts.gb"):
            with Entrez.efetch(db='nucleotide', id=id, rettype='gb', style='withparts') as h, open(f"../refs/{id}_parts.gb", "w") as fo:
                fo.write(h.read())
        gbx = SeqIO.read(f"../refs/{id}_parts.gb", "genbank")
    # Analyze GenBank data and find KIR locations and alleles
    # TODO: there are some genes that are not annotated properly; check them
    for f in gbx.features:
        if f.type == 'gene' or f.type == 'CDS':
            gn = f.qualifiers.get('gene', [''])[0]
            an = f.qualifiers.get('allele', [''])[0]
            if not gn.startswith('KIR') and an:
                raise 'not expected'
            if gn.startswith('KIR') and an:
                info.append((an, location(f)))
            elif gn.startswith('KIR'):
                info.append((gn + '*???', location(f)))

    if info:  # found annotations
        overlap = lambda a, b: max(0, min(a[1], b[1]) - max(a[0], b[0]))  # check if two intervals overlap
        info = sorted(set(info))
        results = collections.defaultdict(list)
        for gi, (ga, ((st, en), (nst, nen))) in enumerate(info):
            gn, an = ga.split('*')
            for i in range(nst, nen): ranges[i] = gn  # populate ranges
            # check if current annotations overlaps with the previous
            # (often annotations are messed up and overlap)
            overwrite = gi and info[gi - 1][0].split('*')[0] == gn and overlap((st, en), (info[gi - 1][1][0]))
            # find the best allele for this range
            r = find_closest(gn, an, str(seq[st:en]), nst, nen, depth+2)
            if overwrite: results[gn][-1] = r
            else: results[gn].append(r)
        # print best hits
        print('  ' * (depth+2), '@@', ','.join([id] + [';'.join(results[g]) for g in genes]))

    # analyze subsequences; we assume only depth-1 sequences (always the case, it seems)
    c = gb.annotations.get('contig', '-')
    if c == '-': return
    ost_new = 0  # TODO: fix for larger depths
    for m in re.findall(r'((complement\()?([0-9A-Z_]{2}\d+\.\d):(\d+)..(\d+))|(gap\((\d+)\))', c):
        if m[5]:  # gap, ignore
            print(' ' * (depth+4), f'GAP_{m[6]}')
            ost_new += int(m[6])
        else:  # hit, analyze nested sequence
            m = m[1:]
            rc, s, e = m[0] != '', int(m[2]), int(m[3])
            l = e - s
            deep_search(m[1], depth+2, ost_new, s, e, rc, ranges)
            ost_new += l

if __name__ == "__main__":
    ranges = collections.defaultdict(dict)
    for id in ids:
        deep_search(id, ranges=ranges[id])

#%% Simple paired-end perfect simulations
if __name__ == "__main__":
    with pysam.FastaFile("../../inumanag-kir/refs/chr19.fa") as fa:
        seq = str(fa.fetch("chr19"))
        off = 1000
        lo = seq[54724235-off:54724235]
        ro = seq[54867216:54867216+off]
    for id in ids:
        print(id)
        with pysam.FastaFile(f"../refs/{id}.fa") as fa:
            seq = lo + str(fa.fetch(id)) + ro
        # Create both FQs and SAMs
        with open(f"../sim/{id}.fa", "w") as fo, open(f"../sim/{id}.sam", "w") as fs:
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

#%% Check Mike's data
with open('scripts/HPRC-kir-annotations.tsv') as f:
    for l in f:
        asm, ctg, st, ed, g, a, seq = l.strip().split('\t')
        find_closest(asm[:7], g, a, seq)

#%% Check T1K and Geny on HPRC
ground = {}
with open("scripts/hprc-ground.csv") as f:
    for l in f:
        if l.startswith('ID'): continue
        l = l.strip().split(',')
        ground[l[0].split('.')[0]] = [[z[1:] if z and z[0] == '`' else z for z in i.split(';')] for i in l[1:]]

t1k = collections.defaultdict(lambda: collections.defaultdict(list))
for path in glob.glob("T1K/_x/*_genotype.tsv"):
    id = os.path.basename(path)[4:11]
    with open(path) as f:
        for l in f:
            l = l.split('\t')
            if int(l[1]) > 0:
                t1k[id][l[0]].append('/'.join(x.split('*')[1][:3] for x in l[2].split(',')))
            if int(l[1]) > 1:
                t1k[id][l[0]].append('/'.join(x.split('*')[1][:3] for x in l[5].split(',')))

results = collections.defaultdict(lambda: collections.defaultdict(list))
parents = collections.defaultdict(list)
for fid in glob.glob('aldy-kir/_logs/*.log'):
    id = fid.split('.')[0].split('/')[-1]
    cid, id = id.split('_') if '_' in id else (id, id)
    if cid != id: parents[cid].append(id)
    d_ilp = collections.defaultdict(list)
    with open(fid) as f:
        for l in f:
            lp = l.strip().split()
            if '-> copy number:' in l:
                g, a = lp[0].split('*')
                results[id][g] += [a[:3]] * int(float(lp[-1]))

def enumerate_alleles(x):
    x = collections.Counter([i[:3] for i in x])
    return [(i, j) for i, c in sorted(x.items()) for j in range(c)]
def mendel(results):
    z = 0
    for cid, (mid, fid) in parents.items():
        for g in genes:
            cg, mfg = set(enumerate_alleles(results[cid][g])), set(enumerate_alleles(results[mid][g] + results[fid][g]))
            if not (cg <= mfg):
                z += 1
                print('->', g, cid, results[cid][g])
                print('          ', mid, results[mid][g])
                print('          ', fid, results[fid][g])
    print('bad hits', z)
mendel(results)
mendel(t1k)

def get_correct(ref, dat):
    nb, nm = 0, 0
    dgrn = [x[:3] for x in sorted(ref)]
    for x in dat:
        if x:
            if x[:3] in dgrn:
                dgrn.remove(x[:3])
            elif any((z := y[:3]) in dgrn for y in x.split('/')[1:]):
                nm += 1
                dgrn.remove(z)
            else:
                nb += 1
    if nb <= len(dgrn):
        nb += len(dgrn) - nb
    return nb, nm

Nt, Gb, Tb, Tm = 0, 0, 0, 0
for id in ground:
    row = []
    print(f'{"":14}', end=' ')
    for _, g in enumerate(genes): print(f'| {g:<8}', end=' ')
    print()

    for gi, g in enumerate(genes):
        grn = sorted(x[:3] for x in ground[id][gi] if x)
        gny = sorted(x     for x in results[id][g])
        tok = sorted(x[:3] for x in t1k[id][g] if x)

        Nt += len(grn)
        e = ['+'.join(grn), '+'.join(gny), '+'.join(tok)]

        if g == 'KIR3DP1':
            row.append(e+ ['-', '-'])
            continue
        nb, _ = get_correct(grn, gny)
        Gb += nb
        e.append('#' * nb) # geny
        nb, nm = get_correct(grn, t1k[id][g])
        Tb += nb
        Tm += nm
        e.append(('!' * nb) + ('?' * nm)) #  t1k
        row.append(e)

    head = [id, 'geny', 't1k', 'g_e', 't_e']
    for ri in range(5):
        print(f'{head[ri]:14}', end=' ')
        for gi, g in enumerate(genes):
            print(f'| {row[gi][ri]:<8}', end=' ')
        print()
    print('-' * 200)
print(Nt, Gb, Gb/Nt, Tb, Tb/Nt, Tm)

#%% Check T1K and Geny on simulations
ground = {}
with open("scripts/sim-ground.csv") as f:
    for l in f:
        if l.startswith('ID'): continue
        l = l.strip().split(',')
        ground[l[0].split('.')[0]] = [[z[1:] if z and z[0] == '`' else z for z in i.split(';')] for i in l[1:]]

t1k = collections.defaultdict(lambda: collections.defaultdict(list))
for path in glob.glob("T1K/_sim/*.fq/T1K_genotype.tsv"):
    id = path.split('/')[2][:-5]
    with open(path) as f:
        for l in f:
            l = l.split('\t')
            if int(l[1]) > 0:
                t1k[id][l[0]].append('/'.join(x.split('*')[1][:3] for x in l[2].split(',')))
            if int(l[1]) > 1:
                t1k[id][l[0]].append('/'.join(x.split('*')[1][:3] for x in l[5].split(',')))

results = collections.defaultdict(lambda: collections.defaultdict(list))
for fid in glob.glob('aldy-kir/_logs/sim/*.log'):
    id = fid.split('/')[-1][:-6]
    d_ilp = collections.defaultdict(list)
    with open(fid) as f:
        for l in f:
            lp = l.strip().split()
            if '-> copy number:' in l:
                g, a = lp[0].split('*')
                results[id][g] += [a[:3]] * int(float(lp[-1]))

Nt, Gb, Tb, Tm = 0, 0, 0, 0
Ng, Gg, Tg = 0, 0, 0
for id in ground:
    row = []
    print(f'{"":14}', end=' ')
    for _, g in enumerate(genes): print(f'| {g:<8}', end=' ')
    print()

    for gi, g in enumerate(genes):
        grn = sorted(x[:3] for x in ground[id][gi] if x)
        gny = sorted(x     for x in results[id][g])
        tok = sorted(x[:3] for x in t1k[id][g] if x)

        Nt += len(grn)
        e = ['+'.join(grn), '+'.join(gny), '+'.join(tok)]

        if g == 'KIR3DP1':
            row.append(e+ ['-', '-'])
            continue
        nb, _ = get_correct(grn, gny)
        Gb += nb
        if nb: Gg += 1
        e.append('#' * nb) # geny
        nb, nm = get_correct(grn, t1k[id][g])
        Tb += nb
        Tm += nm
        if nb: Tg += 1
        e.append(('!' * nb) + ('?' * nm)) #  t1k
        row.append(e)

        if grn or gny or tok: Ng += 1

    head = [id, 'geny', 't1k', 'g_e', 't_e']
    for ri in range(5):
        print(f'{head[ri]:14}', end=' ')
        for gi, g in enumerate(genes):
            print(f'| {row[gi][ri]:<8}', end=' ')
        print()
    print('-' * 200)
print(f'{Nt} ({Ng})', f'{Gb} ({Gg}) {Tb} {Tm} ({Tg})')
















#%%
#%%


ground = {}
with open("../_exp/ground.csv") as f:
    for l in f:
        if l.startswith('ID'): continue
        l = l.strip().split(',')
        ground[l[0]] = [[z for z in i.split(';')] for i in l[1:]]
for id in ids:
    d_filt = collections.defaultdict(set)
    d_em = collections.defaultdict(set)
    d_ilp = collections.defaultdict(list)
    try:
        with open(f"../sim/full/{id}.fq.log_em") as f:
            for l in f:
                lp = l.split()
                if l.startswith('[filter] KIR'):
                    g, a = lp[1:3]
                    d_filt[g].add(a[:3])
                if l.startswith('[em] KIR'):
                    g, a = lp[1:3]
                    d_em[g].add(a[:3])
                if l.startswith('[ilp] KIR'):
                    g, a = lp[1:3]
                    d_ilp[g].append(a[:3])
    except:
        pass
    # print(d_filt)
    # print(d_em)
    # print(d_ilp)
    # print(ground[id])
    bad = 0
    total = 0
    dx = []
    print(f'{"":14}', end=' ')
    for _, g in enumerate(genes):
        print(f'| {g:<8}', end=' ')
    print()
    for gi, g in enumerate(genes):
        total += 1
        gg = [x[:3] for x in sorted(ground[id][gi]) if x]
        i = sorted(d_ilp[g])

        e = ['+'.join(gg), '+'.join(i)]
        e.append('+'.join('% '[int(a in d_filt[g])] for a in gg))
        e.append('+'.join('$ '[int(a in   d_em[g])] for a in gg))
        e.append('+'.join('# '[int(ga == a)]        for ga, a in itertools.zip_longest(gg, i, fillvalue='')))
        dx.append(e)

    head = [id, 'solved', 'filter', 'em', 'ilp']
    for ri in range(5):
        print(f'{head[ri]:14}', end=' ')
        for gi, g in enumerate(genes):
            print(f'| {dx[gi][ri]:<8}', end=' ')
        print()
    print('-' * 200)
    #     if not all(ga[:3] == ia for ga, ia in itertools.zip_longest(gg, i, fillvalue='')):
    #         # print(id, g, gg, i)
    #         bad += 1
    # print(f'=> {id:12} {bad:2} {total:2}')



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
ground['NT_187683.1']
# %%


#%% ground hprc

#%%
t1k

#%%

    bad = 0
    total = 0
    dx = []
    print(f'{"":14}', end=' ')
    for _, g in enumerate(genes):
        print(f'| {g:<8}', end=' ')
    print()
    for gi, g in enumerate(genes):
        total += 1
        gg = [x[:3] for x in sorted(ground[id][gi]) if x]
        i = sorted(d_ilp[g])

        e = ['+'.join(gg), '+'.join(i)]
        e.append('+'.join('% '[int(a in d_filt[g])] for a in gg))
        e.append('+'.join('$ '[int(a in   d_em[g])] for a in gg))
        e.append('+'.join('# '[int(ga == a)]        for ga, a in itertools.zip_longest(gg, i, fillvalue='')))
        dx.append(e)

    head = [id, 'solved', 'filter', 'em', 'ilp']
    for ri in range(5):
        print(f'{head[ri]:14}', end=' ')
        for gi, g in enumerate(genes):
            print(f'| {dx[gi][ri]:<8}', end=' ')
        print()
    print('-' * 200)



#%%
    # print(d_filt)
    # print(d_em)
    # print(d_ilp)
    # print(ground[id])
    bad = 0
    total = 0
    dx = []
    print(f'{"":14}', end=' ')
    for _, g in enumerate(genes):
        print(f'| {g:<8}', end=' ')
    print()
    for gi, g in enumerate(genes):
        total += 1
        gg = [x[:3] for x in sorted(ground[id][gi]) if x]
        i = sorted(d_ilp[g])

        e = ['+'.join(gg), '+'.join(i)]
        e.append('+'.join('% '[int(a in d_filt[g])] for a in gg))
        e.append('+'.join('$ '[int(a in   d_em[g])] for a in gg))
        e.append((wq := '+'.join('# '[int(ga == a)]        for ga, a in itertools.zip_longest(gg, i, fillvalue=''))))
        dx.append(e)

    head = [id, 'solved', 'filter', 'em', 'ilp']
    for ri in range(5):
        print(f'{head[ri]:14}', end=' ')
        for gi, g in enumerate(genes):
            print(f'| {dx[gi][ri]:<8}', end=' ')
        print()
    print('-' * 200)
    #     if not all(ga[:3] == ia for ga, ia in itertools.zip_longest(gg, i, fillvalue='')):
    #         # print(id, g, gg, i)
    #         bad += 1
    # print(f'=> {id:12} {bad:2} {total:2}')


