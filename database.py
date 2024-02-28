#%% imports
import itertools
import os
import sys
import re
import copy
import pickle
import collections
import subprocess as sp

import pysam
import yaml
import parasail

from Bio import Entrez
from Bio import SearchIO, SeqIO
from Bio.Seq import Seq
from aldy.gene import Gene as AldyGene

from common import Region, Gene

# Goodies
class dotdict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def powerset(iterable):
    from itertools import chain, combinations
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


# YAML formatting
class literal(str):
    pass


def literal_presenter(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")


yaml.add_representer(literal, literal_presenter)
yaml.Dumper.ignore_aliases = lambda *_: True

# BioPython
Entrez.email = "inumanag@uvic.ca"

# Aldy
import logbook
sh = logbook.StderrHandler( format_string="{record.message}", level=logbook.DEBUG )
sh.push_application()

# Parasail parameters
PMM, PGO, PGE = parasail.matrix_create("ACGT", 2, -8), 12, 2

# Parse CIGAR from Parasail alignment
def get_cigar(a):
    l = re.split(r'([A-Z=])', a.cigar.decode.decode())[:-1]
    return [(int(l[i]), l[i+1]) for i in range(0, len(l), 2)]

def allele_align(a):
    """Compare allele to the wildtype allele and calculate differences"""

    a.new_regions = set()
    regions, ops = {}, []
    w_pos = 0
    for wr in a.gene.wildtype.regions.values():
        if (
            (wr.name not in a.regions and not wr.name.startswith("e")) and
            not (a.gene.name == 'KIR3DL3' and a.name == '005' and wr.name == 'i1')  # special case, intron deletion
        ):  # fill only missing introns / UTRs
            # A region is not present in the allele (e.g., intron or equal exon).
            # Copy it from the wildtype allele.
            a.new_regions.add(wr.name)
            regions[wr.name] = wr
            w_pos += len(wr.seq)
            continue

        if wr.name not in a.regions:
            print(f'  => {a.gene.name}.{a.name}: missing {wr.name}')
        r = regions[wr.name] = copy.copy(
            a.regions.setdefault(wr.name, Region(a.gene.name, wr.name, "", False, False))
        )
        a.regions[wr.name].partial = False
        # Align the region to the wildtype allele region.
        wi, ri = 0, 0
        if r.partial:
            # Partial exons must be aligned in infix (HW) mode to remove side gaps
            aln = parasail.sg_qx_trace_scan_32(wr.seq, r.seq, PGO, PGE, PMM)
            cigar = get_cigar(aln)
            if cigar[0][1] == 'I':
                ri = wi = cigar[0][0]
                cigar = cigar[1:]
            if cigar[-1][1] == 'I':
                cigar.pop()
            # For partial exons: pad the sequence with the missing content
            r.seq = wr.seq[:wi] + r.seq
        else:
            if not r.seq:
                cigar = [(len(wr.seq), 'I')]
                if wr.name.startswith("e"):
                    a.gene.functional[w_pos + wi, "del" + wr.seq] = f'{wr.name}_DEL'
            else:
                aln = parasail.nw_trace_scan_32(wr.seq, r.seq, PGO, PGE, PMM)
                cigar = get_cigar(aln)

        for sz, op in cigar:
            # Calculate differences between the wildtype and the current allele
            if op == "D":
                ops.append((w_pos + wi, "ins" + r.seq[ri : ri + sz]))
                ri += sz
            elif op == "I":
                ops.append((w_pos + wi, "del" + wr.seq[wi : wi + sz]))
                wi += sz
            else:
                if op == "X":
                    for i in range(sz):
                        ops.append(( w_pos + wi + i, f"{wr.seq[wi + i]}>{r.seq[ri + i]}"))
                ri += sz
                wi += sz
        if wi < len(wr.seq):
            # For partial exons: pad the sequence with the missing content
            assert r.partial, (a.gene.name, a.name)
            r.seq += wr.seq[wi:]
        w_pos += len(wr.seq)
    a.old = a.regions  # Keep the original (old) regions
    a.regions, a.ops = regions, ops

def generate_full_kir_seq(fa, genes):
    """
    Generate the complete KIR locus with all 17 KIR genes.
    Assumes that refs/chr19.fa (hg38's chr19) exists.
    Gene order taken from: https://onlinelibrary.wiley.com/doi/full/10.1111/imm.12847
    """

    # GenBank IDs that contain missing KIR genes
    ids = ["GU182347.1", "NW_003571055.2"]
    seqs = {}
    with pysam.FastaFile(fa) as fa:
        seqs["chr19"] = str(fa.fetch("chr19"))
    for id in ids:
        with Entrez.efetch(db="nucleotide", id=id, rettype="gb") as h, open(
            f"refs/{id}.gb", "w"
        ) as fo:
            fo.write(h.read())
        with Entrez.efetch(db="nucleotide", id=id, rettype="fasta") as h, open(
            f"refs/{id}.fa", "w"
        ) as fo:
            fo.write(h.read())
        with pysam.FastaFile(f"refs/{id}.fa") as fa:
            seqs[id] = str(fa.fetch(id))

    # First KIR gene cluster (GU182347.1)
    id = ids[0]
    gb = SeqIO.read(f"refs/{id}.gb", "genbank")
    locs = []
    for f in gb.features:
        if not (
            f.type == "CDS" or (f.type == "gene" and "DP" in f.qualifiers["gene"][0])
        ):
            continue
        if not (g := f.qualifiers["gene"][0]).startswith("KIR"):
            continue
        g, a = f.qualifiers["allele"][0].split("*")

        ex = [(int(p.start), int(p.end)) for p in f.location.parts]
        # BLAT the database wildtype gene to the GenBank sequence to find the exact coordinates of the match
        with open("q.fa", "w") as fo:
            print(f">{g}", file=fo)
            print(genes[g].seq, file=fo)
        sp.check_call(
            [
                "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/blat/3.5/bin/blat",
                "-t=dna",
                f"refs/{id}.fa",
                "-q=dna",
                "q.fa",
                f"o_{g}.psl",
            ],
            stdout=sp.DEVNULL,
        )

        res = SearchIO.read(f"o_{g}.psl", "blat-psl")
        hit = res.hits[0]
        hsp = hit.hsps[1 if (g, a) == ("KIR2DS3", "002") else 0]
        # Good matches are thise that end with OK (i.e. cover the whole gene). All of them are.
        print(
            f"{g:8} {a:10} ex={len(ex):2} => {hit.id}:{hsp.hit_start:6}-{hsp.hit_end:6} vs {res.id:8}:{hsp.query_start:6}-{hsp.query_end:6}",
            "OK"
            if hsp.query_start == 0 and hsp.query_end == len(genes[g].seq)
            else "  ",
            abs((hsp.hit_end - hsp.hit_start) - (hsp.query_end - hsp.query_start)),
        )
        # KIR3DP1 alignment is screwed up. Fix it manually.
        if g == "KIR3DP1":
            st = 10677  # also 1800 bp prefix force-inserted as it is not mapped
        else:
            st = hsp.hit_start
        # Add the spacer region.
        if locs:
            locs.append([locs[-1][0] + "_POST", id, locs[-1][-1], st])
        locs.append([g, id, st, hsp.hit_end])
        os.unlink(f"o_{g}.psl")

    # Add the missing genes that are present in the hg38.
    #   KIR3DL3 54724235-54736633 R (already added)
    # > KIR2DL3 54738278-54753053 R
    # > KIR3DL1 54816234-54830779 R
    # > KIR2DS4 54832498-54848567 R
    #   KIR3DL2 54850208-54867216 R (already added)
    locs[8] = ("KIR2DL3", "chr19", 54738278, 54753053)
    locs[9][0] = "KIR2DL3_POST"
    locs[18:] = [
        ["KIR3DL1", "chr19", 54816234, 54830779],
        ["KIR3DL1_POST", "chr19", 54830779, 54832676],
    ] + locs[18:]
    locs[28:] = [
        ["KIR2DS4", "chr19", 54832498, 54848567],
        ["KIR2DS4_POST", "chr19", 54848567, 54850208],
    ] + locs[28:]

    # Add missing genes from the second KIR cluster (NW_003571055.2)
    id = ids[1]
    gb = SeqIO.read(f"refs/{id}.gb", "genbank")
    locs2 = []
    for f in gb.features:  # same as above
        if not (
            f.type == "CDS" or (f.type == "gene" and "DP" in f.qualifiers["gene"][0])
        ):
            continue
        if not (g := f.qualifiers["gene"][0]).startswith("KIR"):
            continue

        if g not in ["KIR2DS5", "KIR2DS1"]:
            continue  # only these two are of interest
        ex = [(int(p.start), int(p.end)) for p in f.location.parts]
        with open("q.fa", "w") as fo:
            print(f">{g}", file=fo)
            print(genes[g].seq, file=fo)
        sp.check_call(
            [
                "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/blat/3.5/bin/blat",
                "-t=dna",
                f"refs/{id}.fa",
                "-q=dna",
                "q.fa",
                "o.psl",
            ],
            stdout=sp.DEVNULL,
        )

        res = SearchIO.read("o.psl", "blat-psl")
        hit = res.hits[0]
        hsp = hit.hsps[0]
        print(
            f"{g:8} {a:10} ex={len(ex):2} => {hit.id}:{hsp.hit_start:6}-{hsp.hit_end:6} vs {res.id:8}:{hsp.query_start:6}-{hsp.query_end:6}",
            "OK"
            if hsp.query_start == 0 and hsp.query_end == len(genes[g].seq)
            else "  ",
            abs((hsp.hit_end - hsp.hit_start) - (hsp.query_end - hsp.query_start)),
        )
        st = hsp.hit_start
        if locs2:
            locs2.append([locs2[-1][0] + "_POST", id, locs2[-1][-1], st])
        locs2.append([g, id, hsp.hit_start, hsp.hit_end])
        os.unlink("o.psl")
        os.unlink("q.fa")
    locs[26:] = locs2[:-1] + locs[26:]

    # Generate the FINAL KIR cluster that contains ALL genes.
    # This cluster is then inserted in chr19:54724235-54867216 (old KIR cluster)
    # and the new chr19kir is generated.
    kir_seq = []
    for gene, ref, st, ed in locs:
        if "_POST" not in gene:
            genes[gene].ref_start = 54724235 + sum(len(s) for s in kir_seq)
            kir_seq.append(genes[gene].seq.upper())
        else:
            kir_seq.append(seqs[ref][st:ed].upper())
    patched_seq = seqs["chr19"][:54724235] + "".join(kir_seq) + seqs["chr19"][54867216:]
    with open("refs/chr19kir.fa", "w") as fo:
        print(">chr19kir", file=fo)
        print(patched_seq, file=fo)
    sp.check_call(["samtools", "faidx", "refs/chr19kir.fa"])


# Arguments: <kir.dat> <chr19.fa>
if __name__ == "__main__":  # This guard is needed to prevent multipricessing
    # from executing this in each process
    os.system("mkdir -p refs defs")
    genes = {}

    # Read kir.dat
    # for record in SeqIO.parse("kir.dat", "imgt"):
    for record in SeqIO.parse(sys.argv[1], "imgt"):
        g, a = record.description.split(",")[0].split("*")
        if m := re.search(r'identical to ([0-9A-Z\*]+)', record.description):
            print('Ignoring', g, a, record.description)
        elif m := re.search(r'renamed (KIR[0-9A-Z\*]+)', record.description):
            print('Ignoring', g, a, record.description)
        elif m := re.search(r'error', record.description):
            print('Ignoring', g, a, record.description)
        else:
            genes.setdefault(g, Gene(g, a)).add(a, record)
    # Manually set wildtype allele for the following genes to match
    # the complete sequence.
    genes["KIR2DL4"]._wildtype = "0010201"
    genes["KIR2DP1"]._wildtype = "0010201"
    genes["KIR2DS1"]._wildtype = "0020101"
    genes["KIR2DS3"]._wildtype = "0010301"
    genes["KIR2DS5"]._wildtype = "0020101"
    genes["KIR3DS1"]._wildtype = "0130101"
    for g in genes.values():
        # Make sure that our exons indeed procude the advertised protein.
        p = "".join(
            r.seq for r in g.wildtype.regions.values() if r.is_exon and not r.pseudo
        )
        p = str(Seq(p).translate()).replace("*", "X")
        if g.wildtype.protein and p != g.wildtype.protein:
            # This should only happen for KIR3DS1; however, this protein *is* identical except
            # that is it a few bases longer (and thus we don't care).
            print("=> Protein mismatch", g.name, g.wildtype.name , p, g.wildtype.protein)
        g.wildtype.protein = p

        # Align alleles to the wildtype
        for a in g.alleles.values():
            allele_align(a)
        # Find functional mutations
        for a in g.alleles.values():
            a.parse_mutations()
        print(f"{g.name}: {len(g.seq)=}; {len(g.alleles)=}; {len(g.functional)=}; {len(g.mutations)=}")

    # Generate chr19kir.fa
    generate_full_kir_seq(sys.argv[2], genes)
    span = list(itertools.chain(*[[g.ref_start, g.ref_start + len(g.seq)] for g in genes.values()]))
    span = min(span), max(span)
    print("KIR locus span", span)

    # Check is our reference really correct
    with pysam.FastaFile(f"refs/chr19kir.fa") as fa:
        seq = fa.fetch("chr19kir")
    for g in genes.values():
        r, s = seq[g.ref_start : g.ref_start + len(g.seq)], g.seq
        print(g.name, r == s)
    # os.system("bwa index refs/chr19kir.fa")

    # Assign names to each KIR locus location (e.g., 54825723 -> KIR2DL1:e1)
    regions, st, prev_g = {}, 0, ""
    for g in sorted(genes.values(), key=lambda x: x.ref_start):
        while st and st < g.ref_start:
            # This is a spacer region
            regions[st] = f"{prev_g}:_:{g.name}"
            st += 1
        st = g.ref_start
        prev_g = g.name
        for r in g.wildtype.regions.values():
            for i in range(len(r.seq)):
                regions[st] = g.name + ":" + r.name
                st += 1

    # Calculate keystones (single mutations that solely define allele)
    for g in genes.values():
        seen = {}
        muts = collections.defaultdict(set)
        for a in g.alleles.values():
            mm = frozenset(m for m in a.ops if m in g.functional)
            if mm in seen: continue
            seen[mm] = a.name
            for ps in mm:
                muts[ps].add(a.name)

        for a in g.alleles.values():
            mm = frozenset(m for m in a.ops if m in g.functional)
            major = seen[mm]
            a.keystones = set(
                m for m in mm if len(muts[m]) == 1
            )
            # if a.keystones: print('KEY', g.name, a.name, a.keystones)

    substitutes = {
        'KIR3DP1': {
            (335, 'delGGGGATGGAGATCTGGGCCCAGAGGTGGAGATATAGGCCTGGAGGTGGAGTTATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGATATATGGGCCTGGAGATGGAGTGATGGGCCTAGAAGTGGAGATCTGGGTCTGGAGTGGAGATATGGGCCTGGAGGTGGAGATATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGAGATAGGAACCTGGAGGGGAGATATGAGCCTGGAGTGAAGATATTGGCCTGGGATGGAGATATGGGCCTGGAGTGGAGACATGGGCCTGGAGGTGGAGATATGGGCCTGGAGGTGGAGACATGGGCCTAGAGGTGGATATCTGGGCCTGGAGTGGACATATGGGCCTAGGATGGAGATATGGGCCTGGGTGTGGAGATATGGGCTTGGGGTGGAGATATGGGCCTGGATTGGAGATATGGGTCTAGGGTGGAAATATTGGCCTGGAGTGGAGATATGGGCCTGGAGTGGAGATATGGGCTTGGGGTGGGGATAGGGGCCTGGGGTGCGGATATGGGCCTGCAGGCTGGGTCTCTACACAGCCGACAGCCCTGTTCTTGGGTGCAGGCTGGCACTGAGGGTGAGTTTCCCTTCAGCCCAGCAAGGGCCTGGCTACCAAGACTCACAGCCCAGTGGGGGCAGCAAGGGAGTCCTGGTTTGCCTGCAGATGGATGGTCCATCATGATCTTTCTTTCCAG'): None,
            (1090, 'delGTGAGTCCTTCTCCAAACCTTCGGGTGTCATCTCCCCACATAAGAGGATTTTCCTGAAACAGGAGGGAAGCCCGGTGGGGGATTTTCTTATAAACAAGGATGAGGAGACCCTGGGGTGCTCAGCCCACAGTTCCGACCTTGCCCTCCCCAGCCTTCCTTTCCCTTGGCTGAGTCAGGTTCTGTGGGAACCCGGGAGGGTAGACTGGGGTCCTCCAAGCTGGGCTGTGCGGCTGGGATGTGGTGTCACTGGCAGAGGAAGGGAGCAAAGCAGTGCTAGGAACAGCAGGCCTCTGAGGACAAAGGTGTAACTCACACCCTCCAGCGTTTCCATGACGGTAGGGGCTGCAGTGTGGCTGCTGTCATTCTACCTCAGAGGTGGGGGAACCCCAGCCAGGGCCCTGACCTTCCAAATCCTCTGTTGGGGGCTCAGTTGTGTATTGTGGTTCACACATTGGCTGATATTCCATTCACAAAGAACATGCCCTCGACTCCATGTCTATTTGTGTTGTTTTATGTGAGTAATCTTGCAGGATTAAAATCTAGTAGGAGTCCCTTACTCAGCACTTGCTCAAAGTTCTCAGCTGACACTTTTGTTGTAGAGAGACGCCAAGTCTATGCGGGGTGGGTCCTTCCTGTAGCCCTGGGCACCCAGGTGTGGTAGGAGCCTTAGAAAGTGGAAATGGGAGAATCTTCTGACACGTGGAGGGAGGGGCGGCTC'): None,
            (1054, 'delGGTTCTTCTTGCTGCAGGGGGCCTGGACACATGAGG'): (335, 'delGGGGATGGAGATCTGGGCCCAGAGGTGGAGATATAGGCCTGGAGGTGGAGTTATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGATATATGGGCCTGGAGATGGAGTGATGGGCCTAGAAGTGGAGATCTGGGTCTGGAGTGGAGATATGGGCCTGGAGGTGGAGATATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGAGATAGGAACCTGGAGGGGAGATATGAGCCTGGAGTGAAGATATTGGCCTGGGATGGAGATATGGGCCTGGAGTGGAGACATGGGCCTGGAGGTGGAGATATGGGCCTGGAGGTGGAGACATGGGCCTAGAGGTGGATATCTGGGCCTGGAGTGGACATATGGGCCTAGGATGGAGATATGGGCCTGGGTGTGGAGATATGGGCTTGGGGTGGAGATATGGGCCTGGATTGGAGATATGGGTCTAGGGTGGAAATATTGGCCTGGAGTGGAGATATGGGCCTGGAGTGGAGATATGGGCTTGGGGTGGGGATAGGGGCCTGGGGTGCGGATATGGGCCTGCAGGCTGGGTCTCTACACAGCCGACAGCCCTGTTCTTGGGTGCAGGCTGGCACTGAGGGTGAGTTTCCCTTCAGCCCAGCAAGGGCCTGGCTACCAAGACTCACAGCCCAGTGGGGGCAGCAAGGGAGTCCTGGTTTGCCTGCAGATGGATGGTCCATCATGATCTTTCTTTCCAGGGTTCTTCTTGCTGCAGGGGGCCTGGACACATGAGGGTGAGTCCTTCTCCAAACCTTCGGGTGTCATCTCCCCACATAAGAGGATTTTCCTGAAACAGGAGGGAAGCCCGGTGGGGGATTTTCTTATAAACAAGGATGAGGAGACCCTGGGGTGCTCAGCCCACAGTTCCGACCTTGCCCTCCCCAGCCTTCCTTTCCCTTGGCTGAGTCAGGTTCTGTGGGAACCCGGGAGGGTAGACTGGGGTCCTCCAAGCTGGGCTGTGCGGCTGGGATGTGGTGTCACTGGCAGAGGAAGGGAGCAAAGCAGTGCTAGGAACAGCAGGCCTCTGAGGACAAAGGTGTAACTCACACCCTCCAGCGTTTCCATGACGGTAGGGGCTGCAGTGTGGCTGCTGTCATTCTACCTCAGAGGTGGGGGAACCCCAGCCAGGGCCCTGACCTTCCAAATCCTCTGTTGGGGGCTCAGTTGTGTATTGTGGTTCACACATTGGCTGATATTCCATTCACAAAGAACATGCCCTCGACTCCATGTCTATTTGTGTTGTTTTATGTGAGTAATCTTGCAGGATTAAAATCTAGTAGGAGTCCCTTACTCAGCACTTGCTCAAAGTTCTCAGCTGACACTTTTGTTGTAGAGAGACGCCAAGTCTATGCGGGGTGGGTCCTTCCTGTAGCCCTGGGCACCCAGGTGTGGTAGGAGCCTTAGAAAGTGGAAATGGGAGAATCTTCTGACACGTGGAGGGAGGGGCGGCTC')
        }
    }
    remove_func = {  # Fake functional mutations that we need to remove
    # TODO: 2DS4: 7011 delCCCGGAGCTCCTATGACATGTA
        'KIR2DL3': {(14224, 'insAGATCCAAAGTTGTCTCCTGCCCA')},
        'KIR3DL1': {(14039, 'insCGAGCACCACAGTCAGGTCT')},

        'KIR2DL4': { }
    }

    extra_func = {  # Extra functional mutations that we need to include
        'KIR3DP1': {
            (335, 'delGGGGATGGAGATCTGGGCCCAGAGGTGGAGATATAGGCCTGGAGGTGGAGTTATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGATATATGGGCCTGGAGATGGAGTGATGGGCCTAGAAGTGGAGATCTGGGTCTGGAGTGGAGATATGGGCCTGGAGGTGGAGATATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGAGATAGGAACCTGGAGGGGAGATATGAGCCTGGAGTGAAGATATTGGCCTGGGATGGAGATATGGGCCTGGAGTGGAGACATGGGCCTGGAGGTGGAGATATGGGCCTGGAGGTGGAGACATGGGCCTAGAGGTGGATATCTGGGCCTGGAGTGGACATATGGGCCTAGGATGGAGATATGGGCCTGGGTGTGGAGATATGGGCTTGGGGTGGAGATATGGGCCTGGATTGGAGATATGGGTCTAGGGTGGAAATATTGGCCTGGAGTGGAGATATGGGCCTGGAGTGGAGATATGGGCTTGGGGTGGGGATAGGGGCCTGGGGTGCGGATATGGGCCTGCAGGCTGGGTCTCTACACAGCCGACAGCCCTGTTCTTGGGTGCAGGCTGGCACTGAGGGTGAGTTTCCCTTCAGCCCAGCAAGGGCCTGGCTACCAAGACTCACAGCCCAGTGGGGGCAGCAAGGGAGTCCTGGTTTGCCTGCAGATGGATGGTCCATCATGATCTTTCTTTCCAGGGTTCTTCTTGCTGCAGGGGGCCTGGACACATGAGGGTGAGTCCTTCTCCAAACCTTCGGGTGTCATCTCCCCACATAAGAGGATTTTCCTGAAACAGGAGGGAAGCCCGGTGGGGGATTTTCTTATAAACAAGGATGAGGAGACCCTGGGGTGCTCAGCCCACAGTTCCGACCTTGCCCTCCCCAGCCTTCCTTTCCCTTGGCTGAGTCAGGTTCTGTGGGAACCCGGGAGGGTAGACTGGGGTCCTCCAAGCTGGGCTGTGCGGCTGGGATGTGGTGTCACTGGCAGAGGAAGGGAGCAAAGCAGTGCTAGGAACAGCAGGCCTCTGAGGACAAAGGTGTAACTCACACCCTCCAGCGTTTCCATGACGGTAGGGGCTGCAGTGTGGCTGCTGTCATTCTACCTCAGAGGTGGGGGAACCCCAGCCAGGGCCCTGACCTTCCAAATCCTCTGTTGGGGGCTCAGTTGTGTATTGTGGTTCACACATTGGCTGATATTCCATTCACAAAGAACATGCCCTCGACTCCATGTCTATTTGTGTTGTTTTATGTGAGTAATCTTGCAGGATTAAAATCTAGTAGGAGTCCCTTACTCAGCACTTGCTCAAAGTTCTCAGCTGACACTTTTGTTGTAGAGAGACGCCAAGTCTATGCGGGGTGGGTCCTTCCTGTAGCCCTGGGCACCCAGGTGTGGTAGGAGCCTTAGAAAGTGGAAATGGGAGAATCTTCTGACACGTGGAGGGAGGGGCGGCTC'),
            (2114, 'C>T'),
            # (2101, 'A>C'),
            # (1888, 'T>C'),
        }
    }
    for g in genes.values():
        for o in extra_func.get(g.name, set()):
            g.functional[o] = 'CUSTOM'
        for o in remove_func.get(g.name, set()):
            if o in g.functional:
                del g.functional[o]
        for o, s in substitutes.get(g.name, {}).items():
            if o in g.functional:
                g.functional[s] = g.functional[o]
                del g.functional[o]
        muts = {}
        mut_to_allele = {}
        snp_mut_to_allele = {}

        for a in g.alleles.values():
            a.enabled = False
            a.ops = set(mx for o in a.ops if (mx := substitutes.get(g.name, {}).get(o, o)))
            for o in a.ops:
                muts.setdefault(o, set()).add(a.name)
            a.func = set(o for o in a.ops if o in g.functional)
            a.mstd = set(o for o in a.ops if o not in g.functional if o[1][1] == '>')
            a.mext = set(o for o in a.ops if o not in g.functional if o not in a.mstd)
        for a in sorted(g.alleles.values(), key=lambda a: len(a.mext)):
            a.duplicate = None
            # if (x := minor_to_major.get((g.name, a.name), '???'))[:3] != a.name[:3]:
            #     print(f'warn: {g.name} {a.name} has bad major {x}')
            key_mut = tuple(sorted(a.ops))
            if key_mut in mut_to_allele:
                print(f'skip: {g.name} {a.name} is complete duplicate of {mut_to_allele[key_mut]}')
                a.enabled = False
                a.duplicate = mut_to_allele[key_mut]
                continue
            else:
                mut_to_allele[key_mut] = a.name

            key_snp_mut = tuple(sorted(a.func | a.mstd))
            while key_snp_mut in snp_mut_to_allele:
                ca, pa = g.alleles[a.name].mext, g.alleles[snp_mut_to_allele[key_snp_mut]].mext
                if uniq := set(ca) - set(pa):
                    m = sorted(uniq, key=lambda x: len(x[1]))[0]
                    # print('select', g.name, a.name, m)
                    g.alleles[a.name].mext.remove(m)
                    key_snp_mut = tuple(sorted(key_snp_mut + (m, )))
                else:
                    print(f'skip: {g.name} {a.name} is not indel-distinct from {snp_mut_to_allele[key_snp_mut]}', ca, pa)
                    break
            snp_mut_to_allele[key_snp_mut] = a.name

            a.enabled = True
            a.mutations = set(key_snp_mut)
            a.minor = set(o for o in key_snp_mut if o not in a.func)
            a.has_major_indel = any(o[1][1] != '>' for o in a.func)
            a.has_minor_indel = any(o[1][1] != '>' for o in a.minor)
            a.should_remap = False ## a.has_major_indel or a.has_minor_indel
        g.mutations = set(m for a in g.alleles.values() if a.enabled for m in a.func | a.minor)


    # Generate Aldy YAMLs
    os.system("mkdir -p defs")
    minor_to_major = {}
    for g in genes:
        y = genes[g].yaml()
        with open(f"defs/{g.lower()}.yml", "w") as fo:
            yaml.dump(y, fo, sort_keys=False, default_flow_style=None)

        gene = AldyGene(name=g, path=f"defs/{g.lower()}.yml", genome="hg38")
        gene.madict = (
            {  # set of all mutations: useful for debugging simulations later on
                mi: (ma, mi, set(m.get_minor_mutations(mi)))
                for ma, m in gene.alleles.items()
                for mi in m.minors
                if mi
            }
        )
        for a in gene.alleles:
            for mi in gene.alleles[a].minors:
                minor_to_major[g, mi] = a

    # Generate allele translation maps
    for g in genes.values():
        for a in g.alleles.values():
            # Prepare translation maps
            a.idx_seq = a.seq  # same as .seq, but faster as .seq is function call (property)
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

            # Mutation map for each mutation (allele-indexed)
            # - (0, m) for mutations within allele
            # - (1, m) for functional mutations that are not part of the allele
            # - (2, m) for other gene mutations that are not part of the allele
            a.mutmap = {}
            other_func = {m[0]: m for m in a.gene.functional if m not in a.mutations}
            other_silent = {m[0]: m for m in g.mutations if m not in a.mutations if m not in other_func}
            for wi, ai, s, m in pieces:
                if m and m[1].startswith("ins") and m in a.mutations:
                    for i in range(-1, len(s) + 1): a.mutmap[ai + i] = (0, m)
                if m and m[1].startswith("del") and m in a.mutations:
                    a.mutmap[ai] = a.mutmap[ai + 1] = (0, m)
                for j in range(len(s)):
                    if m and m in a.mutations: a.mutmap[ai + j] = (0, m)
                    elif wi + j in other_func: a.mutmap[ai + j] = (1, other_func[wi + j])
                    elif wi + j in other_silent: a.mutmap[ai + j] = (2, other_silent[wi + j])

    with open("kir.pickle", "wb") as fo:
        pickle.dump((genes, span, regions, minor_to_major), fo)

    for g in genes.values():
        for a in g.alleles.values():
            if a.name[:3] != minor_to_major[g.name, a.name][:3]:
                print('DUPLICATE', g.name, a.name, minor_to_major[g.name, a.name])

    # print("Substitutes & etc.")
    # import pprint, mappy
    # gx = {}
    # for g in genes.values():
    #     # Uncomment to print functional and analyze muts to remove
    #     # pprint.pprint((g.name, sorted(g.functional)))
    #     gx[g.name] = {}
    #     seq = str(g.wildtype.seq)
    #     idx = mappy.Aligner(seq=seq, preset="sr", k=12)
    #     for a in g.alleles.values():
    #         for p, o in sorted(a.func|a.minor):
    #             if o[1]=='>': continue
    #             r = seq[p - 40:p + 60]
    #             if o.startswith('ins'):
    #                 r = r[:40] + o[3:] + r[40:]
    #             if o.startswith('del'):
    #                 r = r[:40] + r[40 + len(o) - 3:]
    #             for h in idx.map(r):
    #                 if h.strand != 1 or h.q_st != 0: continue
    #                 start, s_start = h.r_st, 0
    #                 ops = []
    #                 for size, op in h.cigar:
    #                     if op == 2:
    #                         ops.append((start, "del" + seq[start: start + size]))
    #                         start += size
    #                     elif op == 1:  # Insertion
    #                         ops.append((start, "ins" + r[s_start: s_start + size]))
    #                         s_start += size
    #                     elif op == 4:  # Soft-clip
    #                         s_start += size
    #                     elif op in [0, 7, 8]:  # M, X and =
    #                         start += size
    #                         s_start += size
    #                 gx[g.name].setdefault((p,o),set()).update(ops)
    # pprint.pprint({g: dict(sorted((a, b) for a, b in x.items() if {a}!=b)) for g, x in gx.items()})
