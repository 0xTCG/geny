# %% Imports & Initialization
import pickle
import os
import sys
import re
import numpy as np
import multiprocessing as mp
import subprocess as sp
import collections
import pysam
import mappy
import gzip
import logging
import argparse

from math import ceil
log = logging.info

def parse_args():
    parser = argparse.ArgumentParser(
        description="Geny: a tool for genotyping KIR genes",
    )
    parser.add_argument('file', help="input BAM/CRAM/FASTA/FASTQ file", metavar='file')

    parser.add_argument('-V', '--version', action="version", version="Geny v1.0 (Feb 2024)")
    parser.add_argument('-l', '--log', default=None, help="log file location")
    parser.add_argument('-c', '--coverage', type=int, default=0, help="sample coverage")
    parser.add_argument('-t', '--threads', type=int, default=8, help="number of threads to use")
    return parser.parse_args()

if __name__ == "__main__":
    # Use sys.argv when called from command line
    try:
        get_ipython().__class__.__name__
        NOTEBOOK = True
        log("[!!!] running in notebook")
        args = dotdict({'threads': 16, 'coverage': 20})
        logging.basicConfig(level=logging.INFO)
    except:
        NOTEBOOK = False
        args = parse_args()
        if args.log:
            logging.basicConfig(filename=args.log, filemode='w', encoding='utf-8', level=logging.INFO)
        else:
            logging.basicConfig(level=logging.INFO)
        log("[!!!] running as script")
    log(f"[!!!] using {args.threads} threads")


from common import timeit, timing, dotdict

class Database:
    def __init__(self, path):
        with gzip.open(path, "rb") as fo:
            (self.genes, _, _, self.minor_to_major) = pickle.load(fo)

    def alleles(self):
        for g in self.genes.values():
            yield from g.alleles.values()

if __name__ == "__main__":
    with timing("initializaiton"):
        db = Database("database.geny")


#%% Read BAM reads
class Hit:
    pass

class Sample:
    def __init__(self, path):
        self.path = path
        self.reads = []

    class Read:
        def __init__(self, name, pair, seq, comment=None, qual=None):
            self.name = name
            self.pair = pair
            self.seq = seq
            self.comment = comment
            self.alignments = {}

@timeit
def get_bam_reads(path: str) -> Sample:
    f = Sample(path)
    log(f'[sample] bam: {path}')
    regions = [
        "chr2:142731000-142734000",
        "chr12:98943000-98946000",
        "chr19:12733000-12736000",
        "chr19:41300000-41400000",
        "chr19:46000000-47000000",
        "chr19:52000000-56000000",
    ]
    with pysam.AlignmentFile(path) as bam:
        for contig in bam.header.references:
            if contig.startswith("chr19_"):
                regions.append(contig)
    for r in regions:
        with pysam.AlignmentFile(path) as bam:
            for read in bam.fetch(region=r):
                f.reads.append(Sample.Read(
                    read.query_name,
                    read.is_read2,
                    read.query_sequence,
                    read.query_qualities
                ))
    f.read_len = len(f.reads[0].seq)
    f.reads.sort(key=lambda r: (r.name, r.pair))
    return f

@timeit
def get_fq_reads(path: str) -> Sample:
    f = Sample(path)
    log(f'[sample] fasta: {path}')
    with pysam.FastxFile(path) as fq:
        for read in fq:
            name, pair = read.name.split("/")
            f.reads.append(Sample.Read(
                name,
                pair == "2",
                read.sequence,
                read.comment,
                '', # read.quality
            ))
    f.read_len = len(f.reads[0].seq)
    f.reads.sort(key=lambda r: (r.name, r.pair))
    return f

@timeit
def calculate_coverage(path: str, cn_region=("chr22", 19941772, 19969975)):
    cn = collections.defaultdict(int)
    with pysam.AlignmentFile(path) as bam:
        for read in bam.fetch(region=f"{cn_region[0]}:{cn_region[1]}-{cn_region[2]}"):
            if read.reference_end is None: continue
            a = (read.reference_start, read.reference_end)
            b = (cn_region[1], cn_region[2])
            if a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]:
                start = read.reference_start
                if read.cigartuples is None: continue
                if read.is_supplementary: continue
                for op, size in read.cigartuples:
                    if op in [0, 7, 8, 2]:
                        for i in range(size):
                            cn[start + i] += 1
                        start += size
    cov = sum(cn.values()) / (cn_region[2] - cn_region[1]) / 2
    return round(cov)

if __name__ == "__main__":
    path = args.file
    if path.endswith('.fa') or path.endswith('.fq') or path.endswith('.fastq'):
        sample = get_fq_reads(path)
        cov = args.coverage
    else:
        sample = get_bam_reads(path)
        cov = calculate_coverage(path)

    if args.coverage:
        cov = args.coverage
    if not cov:
        raise ValueError("coverage must be above zero")

    sample.expected_coverage = cov
    sample.min_coverage = 3
    log(f"[read] {os.path.basename(sample.path)}: {len(sample.reads):,} reads loaded ({sample.expected_coverage=:.1f}x)")

#%% Map reads
def align_minimap(file, database, MAX_NM, threads):
    alignments = {g: {a: [] for a in db.genes[g].alleles} for g in db.genes}
    if not os.path.exists(database):
        with open(database, "w") as fo:
            for a in db.alleles():
                print(f">{a.gene.name}.{a.name}", file=fo)
                print(a.seq, file=fo)
    with timing("minimap2"):
        cmd = [
            "minimap2",
            "-x", "sr", "--secondary=yes",  # short-read preset
            "-c",  # calculate CIGAR
            "-P", "--dual=no",  # do all-to-all mapping
            "-t", str(threads),
            database, file,
        ]
        p = sp.Popen(cmd, stderr=sp.DEVNULL, stdout=sp.PIPE, env=os.environ)
        total = 0
        print("Mapping reads: ", end='')
        for li, l in enumerate(iter(p.stdout.readline, "")):
            if li % 5_000_000 == 0: print(f"{li:,}", end='...')
            if not l: break

            rn, rl, qs, qe, st, ref, _, rs, re, _, _, _, nm, *_, cg = l.decode().split("\t")
            nm = int(nm[5:])
            rs, re, rl = int(rs), int(re), int(rl)
            if nm + rl - (re - rs) > MAX_NM: continue

            h = Hit()
            h.enabled = 0
            h.rid = int(rn)
            h.reversed = st != "+"
            h.cost = nm
            h.st, h.en = rs, re
            h.read_st, h.read_en = int(qs), int(qe)
            h.cigar = cg[5:].strip()
            h.ops = None
            g, a = ref.split(".")
            alignments[g][a].append(h)
            sample.reads[h.rid].alignments.setdefault(g, {}).setdefault(a, []).append(h)
            total += 1
        p.wait()
        print()
        assert p.returncode == 0
        for a in db.alleles():
            alignments[a.gene.name][a.name].sort(key=lambda h: h.st)
    log(f"[minimap2] {total:,}/{li:,} mappings handled")
    # os.unlink(file)  # QINGHUI: I keep this for debugging later
    return alignments

if __name__ == "__main__":
    fa = f"{sample.path}.extract.fa"
    with timing("prepare FASTA for mapping"):
        with open(fa, "w") as fo:
            for ri, r in enumerate(sample.reads):
                print(f">{ri}", file=fo)
                print(f"{r.seq}", file=fo)
    sample.alignments = align_minimap(fa, "kirdb.fa", MAX_NM=5, threads=args.threads)

#%% Validation
def validate_hit(r, seq, h, a):
    """Checks if a read maps to allele mutations correctly"""
    def get_cigar(a):
        l = re.split(r"([A-Z=])", a)[:-1]
        return [(int(l[i]), l[i + 1]) for i in range(0, len(l), 2)]

    r = mappy.revcomp(r) if h.reversed else r
    start, s_start = h.st, h.read_st
    cover = 0
    for size, op in get_cigar(h.cigar):
        if op == "D":
            for i in range(size):
                if a.mutmap.get(start + i, (2, 0))[0] < 2:
                    return (False, 0)
            start += size
        elif op == "I":
            if a.mutmap.get(start, (2, 0))[0] < 2:
                return (False, 0)
            s_start += size
        elif op == "S":
            s_start += size
        elif op == "M":
            for i in range(size):
                if (mm := a.mutmap.get(start + i, (2, 0)))[0] < 2:
                    if seq[start + i] != r[s_start + i]:
                        return (False, 0)
                    else:
                        mmi, mm = mm
                        if mmi == 0:
                            cover |= (0b001 if mm in a.func else 0b100)
                        else:
                            cover |= 0b011
            start += size
            s_start += size
    # 000: does not cover any mutation
    # 001: covers functional mutation
    # 010: covers functional shadow (_)
    # 100: covers minor mutation
    return (True, cover)

def filter_allele(a):
    span = [0 for _ in range(len(a.seq) + 2)]
    a.enabled = True
    a.max_span = 0

    for h in sample.alignments[a.gene.name][a.name]:
        if not h.enabled: continue
        a.max_span = max(h.en - h.st, a.max_span)
        for j in range(h.st, h.en):
            span[j] += 1

    minor_covered = {m: 0 for m in a.minor}
    alen = len(a.seq)

    for i, (mmi, mm) in a.mutmap.items():
        if mmi >= 2: continue
        if mmi == 1 or mm in a.func:
            till_end = max(0, alen - i)
            if till_end >= a.max_span and span[i] < sample.min_coverage:
                a.enabled = False
        elif (mmi == 0) and mm in a.minor:
            minor_covered[mm] = 1

    a.uncovered = sum(1 for i in span if i < sample.min_coverage)
    a.minor_uncovered = sum(1 for i in minor_covered.values() if not i) if minor_covered else 0
    if minor_covered:
        if a.minor_uncovered == len(minor_covered):
            a.enabled = False
        else:
            a.minor_miss = a.minor_uncovered / len(minor_covered)
    else:
        a.minor_miss = 0

def parse_allele_hits(a, MAX_IS = 2_000):
    for h in sample.alignments[a.gene.name][a.name]:
        h.enabled, h.cover = validate_hit(sample.reads[h.rid].seq, a.idx_seq, h, a)
    # Paired ends processing
    for h in sample.alignments[a.gene.name][a.name]:
        if not h.enabled: continue
        if h.rid and sample.reads[h.rid - 1].name == sample.reads[h.rid].name:
            rp = sample.reads[h.rid - 1]
        elif h.rid + 1 < len(sample.reads) and sample.reads[h.rid + 1].name == sample.reads[h.rid].name:
            rp = sample.reads[h.rid + 1]
        else:
            continue
        for hp in rp.alignments.get(a.gene.name, {}).get(a.name, []):
            if not hp.enabled: continue
            if abs(h.st - hp.st) < MAX_IS: break
        else:
            # did not find a pair. maybe this is at the edge of the allele? if so, keep it!
            if not (h.st < MAX_IS or h.en + MAX_IS > len(a.idx_seq)):
                h.enabled = False

def parse_gene_alignments(g):
    for a in g.alleles.values():
        parse_allele_hits(a)
    for a in g.alleles.values():
        filter_allele(a)
    return g.name, g.alleles, {an: [(h.enabled, h.cover) for h in hits] for an, hits in sample.alignments[g.name].items()}

if __name__ == "__main__":
    with mp.Pool(args.threads) as pool, timing("filtering"):
      res = pool.map(parse_gene_alignments, db.genes.values())
      for gn, ga, gal in res:
          db.genes[gn].alleles = ga
          for an, hits in sample.alignments[gn].items():
              for hi, h in enumerate(hits):
                  h.enabled, h.cover = gal[an][hi]

    for g in db.genes.values():
        ax = {}
        for a in g.alleles.values():
            if a.enabled:
                ax.setdefault(db.minor_to_major[g.name, a.name][:3], []).append((a.name, a.uncovered/len(a.idx_seq)))
        for a, aa in ax.items():
            al = '; '.join(f"{x[0]}: {x[1]:.1%}" for x in sorted(aa, key=lambda x: x[1]))
            log(f"[filter] {g.name} {a} => {al}")
    log(f"[filter] valid_alleles={sum(1 for a in db.alleles() if a.enabled):,}")

#%% EM & SCC
@timeit
def em(sample, iteration, ROUNDS=1, EPSILON=1e-4, CUTOFF=1e-3, err = 0.01):
    valid_alleles = [a for a in db.alleles() if a.enabled]
    valid_reads = [r for r in sample.reads
                   if any(h.enabled for a in valid_alleles for h in r.alignments.get(a.gene.name, {}).get(a.name, []))]
    log(f"[em] input: {len(valid_reads)=:,} {len(valid_alleles)=:,}")
    matrix_count = np.zeros((len(valid_reads), len(valid_alleles))) # mapping count
    matrix_err = err*np.ones((len(valid_reads), len(valid_alleles))) # error
    matrix_match = np.zeros((len(valid_reads), len(valid_alleles))) # bases of match
    matrix_mismatch = np.zeros((len(valid_reads), len(valid_alleles))) # bases of mismatch
    matrix_ones = np.ones((len(valid_reads), len(valid_alleles)))

    with timing("em_init"):
        for ai, a in enumerate(valid_alleles):
            total = len({i for h in sample.alignments[a.gene.name][a.name] if h.enabled for i in range(h.st, h.en)})
            for ri, r in enumerate(valid_reads):
                if cnt := sum(1 for h in r.alignments.get(a.gene.name, {}).get(a.name, []) if h.enabled):
                    matrix_count[ri][ai] = cnt / total #(scaling)
                    matrix_match[ri][ai] = cnt*len(r.seq)-sum(h.cost for h in r.alignments.get(a.gene.name, {}).get(a.name, []) if h.enabled)
                    matrix_mismatch[ri][ai] = sum(h.cost for h in r.alignments.get(a.gene.name, {}).get(a.name, []) if h.enabled)
    sample.err = err
    np.random.seed(0)
    best_LL, best_phi = -1e100, np.random.rand(len(valid_alleles))
    for _ in range(ROUNDS):
        phi = np.ones(len(valid_alleles))  # uniform initialization
        tot = np.sum(phi)
        for j in range(len(valid_alleles)):
            phi[j] /= tot  # initialization
        LL = 0
        for _ in range(iteration, -1, -1):
            # match
            m1 = np.power((matrix_ones-matrix_err),matrix_match)
            # mismatch
            m2 = np.power(matrix_err,matrix_mismatch)
            matrix = np.multiply(m1,m2)
            matrix = np.multiply(matrix, matrix_count)
            P_A_R = phi * matrix
            P_A_R += 1e-300
            new_LL = np.sum(np.log(np.sum(P_A_R, axis=1)))  # log-likelihood
            if abs(LL - new_LL) < EPSILON: break
            LL = new_LL

            P_A_R_norm = P_A_R.sum(axis=1).reshape(len(valid_reads), 1)
            P_A_R /= P_A_R_norm
            phi = np.sum(P_A_R, axis=0) / len(valid_reads)
            #update penalty
            aaa = np.sum(np.multiply(matrix_match,P_A_R))
            bbb = np.sum(np.multiply(matrix_mismatch,P_A_R))
            r = 1 + aaa / bbb
            new_err = 1 / r
            sample.err = new_err
            matrix_err = new_err*np.ones((len(valid_reads), len(valid_alleles))) # error
        if LL > best_LL:
            best_LL, best_phi = LL, phi

    indices = np.argsort(best_phi)
    em_results = {(valid_alleles[i].gene.name, valid_alleles[i].name): phi[i] for i in indices}
    for g in db.genes.values():
        for a in g.alleles.values():
            if not a.enabled: continue
            a.old_enabled = a.enabled
            a.enabled = (c := em_results.get((a.gene.name, a.name), 0)) > CUTOFF
            if a.enabled:
                log(f'[em] {a.gene.name} {a.name} => {c}')

@timeit
def get_scc(sample, depth, max_depth, top = 0.96):
    from itertools import groupby
    if depth >= max_depth: return
    def ranges(i):
        for a, b in groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
            b = list(b)
            yield b[0][1], b[-1][1]
    
    '''
    given start, end position on allele, get the
    hits whose start position fail within [start, end]
    '''
    def get_hits_in_range(g,a,start,end):
        alns = sample.alignments[g][a]
        for hi, h in enumerate(alns):
            if h.enabled and h.st >= start and h.st <= end:
                yield h

    with timing('correct hits'):
        for g in db.genes.values():
            g.new_funcs = set(m[0] for a in g.alleles.values() if a.enabled for m in a.func)
            for a in db.alleles():
                if not a.enabled: continue
                a.new_funcs = {pos for pos, (mmi, mm) in a.mutmap.items()
                               if (mmi == 0 and mm in a.func) or (mmi == 1 and mm[0] in g.new_funcs)}
                if not a.new_funcs: # 001 allele, enable all reads...
                    if all(h.enabled and (h.cover & 0b010) for h in sample.alignments[a.gene.name][a.name]):
                        continue
                else:
                    for h in sample.alignments[a.gene.name][a.name]:
                        if h.enabled and (h.cover & 0b010) and not any(h.st <= p < h.en for p in a.new_funcs):
                            h.cover &= 0b100
                            continue
    if depth == 0:
        for r in sample.reads:
            r.ok = False
            for gn, g in r.alignments.items():
                for an, hits in g.items():
                    if not db.genes[gn].alleles[an].enabled: continue
                    if any(h.enabled and (h.cover & 0b001) for h in hits):
                        r.ok = True
                        break
                if r.ok: break

    if depth == 0:
        FLANK = 150 # extend each window to left and right by FLANK
    else:
        FLANK = 0
    EDGE = 50
    sample.sccs = {}
    for g in db.genes.values():
        for a in g.alleles.values():

            if not a.enabled: continue
            scc = set()
            hits = sample.alignments[a.gene.name][a.name]

            ki = 0
            for hi, h in enumerate(hits):
                if not h.enabled or not sample.reads[h.rid].ok: continue
                for p in range(h.st, h.en):
                    scc.add(p)

            scc = list(ranges(sorted(list(scc))))

            extended = set() # added flank 
            for r in scc:
                for p in range(max(r[0]-FLANK,EDGE), min(r[1]+FLANK,len(a.seq)-EDGE)):
                    extended.add(p)
            
            sample.sccs[a.gene.name,a.name] = list(ranges(sorted(list(extended))))

    # get landmarks
    DIST = 50
    sample.landmarks = {}
    for g in db.genes.values():
        for a in g.alleles.values():
            if not a.enabled: continue
            sample.landmarks[a.gene.name, a.name] = []
            for scc in sample.sccs[a.gene.name, a.name]:
                l = list(range(scc[0], scc[1]+1, DIST))
                # print(l)
                for k in l:
                    sample.landmarks[a.gene.name, a.name].append(k)
    
    def get_hits(g, a, p):
        alns = sample.alignments[g][a]
        for hi, h in enumerate(alns):
            if h.enabled and h.st <= p < h.en:
                yield hi, h

    for g in db.genes.values():
        for a in g.alleles.values():
            if not a.enabled: continue
            landmarks = sample.landmarks[a.gene.name, a.name]
            for l in landmarks:
                for hi, h in get_hits(a.gene.name, a.name, l):
                    sample.reads[h.rid].ok = True
            
    captured = 0
    not_captured = 0
    for r in sample.reads:
        if r.ok:
            for gn, g in r.alignments.items():
                for an, hits in g.items():
                    if not db.genes[gn].alleles[an].enabled: continue
                    landmarks = sample.landmarks.get((gn, an), [])
                    for h in hits:
                        if not h.enabled:
                            continue
                        # check if there is one landmark capture the hit
                        can_capture = False
                        for landmark in landmarks:
                            if landmark >= h.st and landmark <= h.en:
                                can_capture = True
                        if not can_capture:
                            not_captured += 1
                        else:
                            captured += 1
    perc = captured / (captured + not_captured)
    print(f'percentage of captured is {perc}')
    if perc >= top: return
    get_scc(sample, depth+1, max_depth)

if __name__ == "__main__":
    with timing("em_scc"):
        em(sample, 1_000)
        for g in db.genes.values():
            used = set()
            for a in g.alleles.values():
                if not a.enabled: continue
                ma = db.minor_to_major[g.name, a.name]
                if ma in used:
                    a.old_enabled, a.enabled = a.enabled, False
                else:
                    used.add(ma)
        get_scc(sample)

#%%
@timeit
def ilp_solve(sample):
    from math import ceil
    from gurobipy import GRB, Model, quicksum, Env

    # env = Env(params = params)
    env = Env(empty=True)
    # env.setParam("OutputFlag", 0)
    env.start()
    
    def get_hits(g, a, p):
        alns = sample.alignments[g][a]
        for hi, h in enumerate(alns):
            if h.enabled and sample.reads[h.rid].ok and h.st <= p < h.en:
                # print(g, a, h.cost)
                yield hi, h

    valid_alleles = sorted((a.gene.name, a.name, 0) for a in db.alleles() if a.enabled)
    # print(f'[ilp] {len(valid_alleles)=:,}')
    valid_reads = set(i for i, r in enumerate(sample.reads) if r.ok)
    print(f'[ilp] {len(valid_reads)=:,}')

    # run the ILP and get result
    V = {}  # assignment variable (read -> allele)
    A = {}  # allele chosing variable (1 if an allele is chosen)
    D = {}
    model = Model('geny', env=env)
    model.setParam("Threads", THREADS)
    model.setParam("Presolve", 2)  # agressive presolve
    model.setParam("PreSparsify", 2)  # matrix sparsification
    model.setParam("NonConvex", 2)
    model.setParam('MIPGap', 0.05)
    # model.setParam('PoolSolutions', 50)
    model.setParam('TimeLimit',180*60)
    # Best objective 9.060000000000e+02, best bound 9.060000000000e+02, gap 0.0000%
    for i, (gn, an, _) in enumerate(valid_alleles[::]):  # allele
        max_cn = min(sum(1 for h in sample.alignments[gn][an]
                           if h.enabled if sample.reads[h.rid].ok
                           if h.st <= l < h.en)
                     for l in sorted(sample.landmarks.get((gn, an), [])))
        max_cn = ceil(max_cn / sample.expected_coverage) + 1
        for cni in range(1, min(max_cn,2)):
            valid_alleles.append((gn, an, cni))
    valid_alleles.sort()
    for i in range(len(valid_alleles)):
        A[i] = model.addVar(vtype=GRB.BINARY, name=f'S{i}_{cni}')
    for i in range(1, len(valid_alleles)):
        if valid_alleles[i - 1][:2] == valid_alleles[i][:2]:
            model.addConstr(A[i] <= A[i - 1])
    for j in valid_reads:
        D[j] = model.addVar(vtype=GRB.BINARY, name=f'D{j}')

    # add read assignment variable
    mapped_alleles_positions = {j: set() for j in valid_reads}
    mapped_reads_positions = {i: set() for i in range(len(valid_alleles))}
    for i, (gene, allele, _) in enumerate(valid_alleles):
        for landmark in sample.landmarks.get((gene, allele), []):
            for hi, h in get_hits(gene, allele, landmark):
                V[i, hi, h.rid] = model.addVar(vtype=GRB.BINARY, name=f'V{i}_{hi}_{h.rid}')
                mapped_alleles_positions[h.rid].add((i, hi))
                mapped_reads_positions[i].add((h.rid, hi))

    for j in valid_reads:  # reads
        model.addConstr(D[j] + quicksum(V[i, hi, j] for i, hi in mapped_alleles_positions[j]) <= 1)
        model.addConstr(D[j] + quicksum(V[i, hi, j] for i, hi in mapped_alleles_positions[j]) >= 1)
    for j in valid_reads:  # read
        for i, hi in mapped_alleles_positions[j]:  # gene
            model.addConstr(V[i, hi, j] <= A[i])
    for i in range(len(valid_alleles)):  # allele
        model.addConstr(A[i] <= quicksum(V[i, hi, j] for j, hi in mapped_reads_positions[i]))

    max_land = max(len(l) for l in sample.landmarks.values())
    Z = {}
    Z_costs = []
    DEL = {}
    for i, (gene, allele, _) in enumerate(valid_alleles):
        landmarks = sample.landmarks.get((gene, allele), [])
        avg = []
        for pos in landmarks:
            Z[i, pos] = model.addVar(lb=0, vtype=GRB.CONTINUOUS)
            actual_coverage_on_lm = quicksum(V[i, hi, h.rid] for hi, h in get_hits(gene, allele, pos))
            actual_ferf_coverage_on_lm = quicksum(V[i, hi, h.rid] for hi, h in get_hits(gene, allele, pos) if h.cost == 0)

            e = A[i] * (sample.expected_coverage) - actual_coverage_on_lm

            model.addConstr(Z[i, pos] + e >= 0)
            model.addConstr(Z[i, pos] - e >= 0)
            Z_costs.append(Z[i, pos] / sample.expected_coverage)
            avg.append(actual_coverage_on_lm)

        g_low = 0.5
        g_high = 1.5
        model.addConstr(quicksum(c for c in avg)/len(landmarks) >= g_low*A[i]*sample.expected_coverage)
        model.addConstr(quicksum(c for c in avg)/len(landmarks) <= g_high*A[i]*sample.expected_coverage)
    A_costs = []
    allele_selection_cost = 0
    for i, _ in enumerate(valid_alleles):
        A_costs.append(A[i] * allele_selection_cost)

    # constraint the percentage of dropped reads
    D_costs = []
    read_drop_cost = 0.08
    NUM_DROPPED = quicksum(D[j] for j in valid_reads)
    MAX_DROP_PERC = 0.20
    model.addConstr(NUM_DROPPED/len(valid_reads) <= MAX_DROP_PERC)


    for j in valid_reads: # read
        D_costs.append(D[j] * read_drop_cost)

    NM_costs = []
    nm_cost = 0.05
    penalty = 4 # it has to be this, not 2.7
    max_cost = 0
    prob = ()
    gn_prob = ''
    an_prob = ''
    hi_prob = 0
    for (i, hi, _), v in V.items():
        gn, an, _ = valid_alleles[i]
        sa = sample.alignments[gn][an]
        # print(gn, an, hi, sa[hi].cost, '&')
        # NM_costs.append(nm_cost * (np.exp(sa[hi].cost)-1) * v)
        NM_costs.append(nm_cost * sa[hi].cost * v)
        # NM_costs.append(nm_cost * (pow(penalty, sa[hi].cost)-1) * v)
        max_cost = max(max_cost, sa[hi].cost)
        gn_prob = gn
        an_prob = an
        hi_prob = hi

    print(f'max_cost is {max_cost} {gn_prob} {an_prob} {hi_prob}')

    Z_cost = quicksum(Z_costs)
    A_cost = quicksum(A_costs)
    D_cost = quicksum(D_costs)
    NM_cost = quicksum(NM_costs)
    model.setObjective(Z_cost + A_cost + D_cost + NM_cost)

    model.optimize()
    print(f'[ilp] {model.objVal=:.1f}')
    print(f'[ilp] {Z_cost.getValue()=:.1f}')
    print(f'[ilp] {A_cost.getValue()=:.1f}')
    print(f'[ilp] {D_cost.getValue()=:.1f}')
    print(f'[ilp] {NM_cost.getValue()=:.1f}')

    sol_A = []
    sol_R = collections.defaultdict(list)
    for (i, hi, j), v in V.items():
        if round(v.x) > 0:
            gn, an, _ = valid_alleles[i]
            sol_R[gn, an].append(j)
    for i, (gene, allele, cn) in enumerate(valid_alleles):
        if round(A[i].x) > 0:
            sol_A.append((gene, allele))
            num_reads_total_assigned = set()
            for j in valid_reads:
                for ii, hi in mapped_alleles_positions[j]:
                    if ii == i and (i, hi, j) in V and round(V[i, hi, j].X) > 0:
                        num_reads_total_assigned.add(j)
            print(f'[ilp] {gene} {allele} => {len(num_reads_total_assigned)}')
    dropped = [j for j in valid_reads if round(D[j].x) > 0]
    print(f'[ilp] {len(dropped)=:,}')

    for ai, (gn, an, _) in enumerate(valid_alleles):
        break
        if round(A[ai].x) == 0: continue
        a = db.genes[gn].alleles[an]
        key_landmarks = sorted(a.new_funcs)
        print('Ⓜ️', f'{a.gene.name:8} {a.name:8}', key_landmarks)
        for l in sorted(sample.landmarks.get((a.gene.name, a.name), [])):
            i = 0
            if l in key_landmarks: i = 2 + int(a.mutmap[l][0] == 0)

            # total read coverage
            cov = [(hi, h)
                    for hi, h in enumerate(sample.alignments[a.gene.name][a.name])
                    if h.enabled if sample.reads[h.rid].ok
                    if h.st <= l < h.en]
            cov_sol = [(hi, h.rid)
                        for hi, h in enumerate(sample.alignments[a.gene.name][a.name])
                        if h.rid in sol_R[a.gene.name, a.name] if h.st <= l < h.en]
            cov_cross = []
            diff = abs(len(cov_sol) - sample.expected_coverage)
            cross_gene = {h.rid for _, h in cov
                        if any(hx.enabled
                                for gn, ga in sample.reads[h.rid].alignments.items()
                                if gn != a.gene.name
                                for an, hl in ga.items()
                                if db.genes[gn].alleles[an].enabled
                                for hx in hl)}
            print(f'   {a.gene.name:8} {a.name:8}', f'{l:6}', '⚪️🟠🔴'[i],
                f'{len(cov):5,}', f'({len({h.rid for _, h in cov}):5,})',
                f'cross={len(cross_gene):5}',
                #f'{cov_exp:5,}', f'{cov_bad:5,}', f'{len(cov_fo):5,}',
                "|",
                f'{len(cov_sol):5,}',
                #f'{len(cov_sol_exp):5,}',
                f'{diff:5}' if diff else '',
                a.mutmap.get(l, '')
                )
            if diff and False:
                for hi, h in cov:
                    print(f'     {h.rid:6} {h.st:5} {h.en:5} {h.cigar:10} {sample.reads[h.rid].comment:20}:', end=' ')
                    if x := [ga for ga, r in sol_R.items() if h.rid in r]:
                        print(*x[0])
                    elif h.rid in dropped:
                        print("-")
                    else:
                        print([(valid_alleles[i][0], valid_alleles[i][1], round(x.x))
                            for (i, hi, j), x in V.items()
                            if j == h.rid if round(x.x) > 0])
        print()

    # for d in dropped:
    #     print(f'd: {sample.reads[d].name, sample.reads[d].comment}')
    return sol_A, sol_R, dropped


answer = ilp_solve(sample)
