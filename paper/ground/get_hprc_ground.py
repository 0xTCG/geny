#%%
import os
import csv

# List of KIR genes
kir_genes = [
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
    "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
    "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1"
]

# List of IDs
ids = [
    "HG00438", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741",
    "HG01071", "HG01106", "HG01109", "HG01175", "HG01243", "HG01258",
    "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978",
    "HG02055", "HG02080", "HG02145", "HG02148", "HG02257", "HG02572",
    "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886",
    "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540",
    "HG03579", "NA18906", "NA19240", "NA20129"
]

# Base directory
base_dir = "/project/shared/inumanag-kir/aldy-kir/push/aldy-kir/paper/ground/bakir-annotations"

def process_file(file_path, id, haplotype):
    gene_alleles = {gene: [] for gene in kir_genes}
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return gene_alleles
    
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[0] == id and parts[1] == haplotype and parts[2] in kir_genes:
                gene = parts[2]
                allele = parts[3].split('*')[1]
                gene_alleles[gene].append(allele)
    
    return gene_alleles

results = {}
#%%
for id in ids:
    maternal_file = os.path.join(base_dir, f"{id}.maternal.f1_assembly_v2.tsv")
    paternal_file = os.path.join(base_dir, f"{id}.paternal.f1_assembly_v2.tsv")
    
    maternal_alleles = process_file(maternal_file, id, "maternal")
    paternal_alleles = process_file(paternal_file, id, "paternal")
    
    combined_alleles = {}
    for gene in kir_genes:
        combined = maternal_alleles[gene]+(paternal_alleles[gene])
        combined_alleles[gene] = ';'.join(f'`{allele}' for allele in combined) if combined else ''
    
    results[id] = combined_alleles

#%%
# Write results to CSV
output_file = "hprc_gt_0801.csv"
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    header = ['ID'] + kir_genes
    writer.writerow(header)
    
    for id in ids:
        row = [f"{id}.bam"] + [results[id][gene] for gene in kir_genes]
        writer.writerow(row)

print(f"Results have been written to {output_file}")