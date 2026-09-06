"""Generate synthetic test dataset for end-to-end pipeline verification.

Creates:
  - test_data/synthetic_bgc.fasta   (full genome with BGC for input)
  - test_data/synthetic_vector.fasta (capture vector with EcoRI site)
  - test_data/synthetic_genome.fasta (full genome including BGC + flanks)
  - test_data/synthetic_bgc.gbk     (GenBank annotation with 'cluster' feature)
"""
import os
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


random.seed(42)


def generate_gc_rich_region(length, gc_target=0.70):
    """Generate a GC-rich random sequence (mimicking Streptomyces)."""
    bases = "GC"
    result = []
    gc_count = 0
    for _ in range(length):
        if gc_count / max(len(result), 1) < gc_target:
            result.append(random.choice(bases))
            gc_count += 1
        else:
            result.append(random.choice("AT"))
    return "".join(result)


def insert_pam_sites(seq, interval=70, strand="both"):
    """Insert NGG PAM sites at regular intervals for testing."""
    seq_list = list(seq)
    positions = []
    for pos in range(20, len(seq_list) - 3, interval):
        if random.random() < 0.7:
            # Insert NGG
            seq_list[pos:pos+3] = list("AGG")
            positions.append(pos)
            # Shift subsequent positions
            for i in range(len(positions)):
                if positions[i] > pos:
                    positions[i] += 3
    return "".join(seq_list), positions


def generate_synthetic_bgc(size=3000):
    """Generate a synthetic BGC with gene-like structure and PAM sites."""
    # Build a BGC with alternating high-GC genes and intergenic regions
    parts = []
    # Gene 1: high GC, contains a PKS module
    gene1 = generate_gc_rich_region(500, 0.72)
    # Add a KS-like motif
    gene1 = gene1[:100] + "GGCAGCGGCGGCGGCAGCGG" + gene1[120:]
    # Ensure PAM sites
    gene1, _ = insert_pam_sites(gene1, 150)
    parts.append(gene1)

    # Intergenic region
    inter1 = generate_gc_rich_region(200, 0.68)
    parts.append(inter1)

    # Gene 2: AT-Less cassette
    gene2 = generate_gc_rich_region(400, 0.70)
    gene2 = gene2[:50] + "ATGTTTGGCCGCCTTCATCG" + gene2[70:]
    gene2, _ = insert_pam_sites(gene2, 120)
    parts.append(gene2)

    # Intergenic
    inter2 = generate_gc_rich_region(150, 0.65)
    parts.append(inter2)

    # Gene 3
    gene3 = generate_gc_rich_region(450, 0.71)
    gene3, _ = insert_pam_sites(gene3, 140)
    parts.append(gene3)

    # Intergenic
    inter3 = generate_gc_rich_region(100, 0.69)
    parts.append(inter3)

    # Gene 4 (Acyl-CoA dehydrogenase-like)
    gene4 = generate_gc_rich_region(400, 0.70)
    gene4, _ = insert_pam_sites(gene4, 130)
    parts.append(gene4)

    # Intergenic
    inter4 = generate_gc_rich_region(100, 0.66)
    parts.append(inter4)

    # Gene 5 (oxidase)
    gene5 = generate_gc_rich_region(350, 0.71)
    gene5, _ = insert_pam_sites(gene5, 110)
    parts.append(gene5)

    # Intergenic / terminator region
    inter5 = generate_gc_rich_region(80, 0.65)
    parts.append(inter5)

    # Gene 6 (regulatory)
    gene6 = generate_gc_rich_region(200, 0.70)
    gene6, _ = insert_pam_sites(gene6, 100)
    parts.append(gene6)

    bgc_seq = "".join(parts)
    # Pad to exact size if needed
    while len(bgc_seq) < size:
        bgc_seq += random.choice("GC")
    bgc_seq = bgc_seq[:size]
    return bgc_seq


def generate_flanking_region(size=1000, prefix="upstream"):
    """Generate flanking genomic sequence with some PAM sites."""
    seq = generate_gc_rich_region(size, 0.69)
    seq, _ = insert_pam_sites(seq, 100)
    return seq


def generate_vector(size=5500):
    """Generate a synthetic capture vector with a unique EcoRI site.

    The vector has:
    - ampicillin resistance gene
    - Origin of replication
    - Multiple cloning site with EcoRI site
    """
    # Build a circular vector
    # Backbone (high-copy origin)
    origin = generate_gc_rich_region(1500, 0.55)

    # Ampicillin resistance gene
    amp_r = generate_gc_rich_region(800, 0.52)

    # Multiple cloning site with EcoRI site (GAATTC)
    mcs = "GAATTCGCGGCCGCTCTAGAGGATCCACTAGTCGAGCCATGG"

    # Yeast replication origin (optional element)
    yeast_origin = generate_gc_rich_region(500, 0.45)

    # Additional backbone
    backbone2 = generate_gc_rich_region(2600, 0.50)

    vector_seq = origin + amp_r + mcs + yeast_origin + backbone2
    if len(vector_seq) > size:
        vector_seq = vector_seq[:size]
    elif len(vector_seq) < size:
        vector_seq += generate_gc_rich_region(size - len(vector_seq), 0.50)
    return vector_seq


def write_fasta(seq, path, seq_id, description=""):
    rec = SeqRecord(Seq(seq.upper()), id=seq_id, description=description)
    with open(path, "w") as fh:
        SeqIO.write(rec, fh, "fasta")


def write_genbank(seq, path, seq_id, features):
    rec = SeqRecord(Seq(seq.upper()), id=seq_id, description="Synthetic BGC")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["genbank_division"] = "BCT"
    rec.annotations["date"] = "06-SEP-2026"
    rec.annotations["source"] = "Synthetic test data"
    rec.features = features
    # Add a source feature
    src = SeqFeature(FeatureLocation(0, len(seq)), type="source",
                     qualifiers={"organism": "Streptomyces synthetic",
                                 "strain": "test strain"})
    rec.features.insert(0, src)
    with open(path, "w") as fh:
        SeqIO.write(rec, fh, "genbank")


def main():
    test_data_dir = Path(__file__).parent / "test_data"
    test_data_dir.mkdir(exist_ok=True)

    # Generate parts
    upstream = generate_flanking_region(1000, "upstream")
    bgc = generate_synthetic_bgc(3000)
    downstream = generate_flanking_region(1000, "downstream")
    genome = upstream + bgc + downstream
    vector = generate_vector(5500)

    # Calculate BGC coordinates in the genome
    bgc_start = 1000  # 0-based
    bgc_end = 1000 + len(bgc)

    # Add CDS features within the BGC
    features = []
    gene_starts = [100, 550, 1050, 1500, 1900, 2400]
    gene_ends = [500, 1000, 1450, 1850, 2350, 2900]
    gene_names = ["PKS_KS", "PKS_AT", "PKS_PC", "dehydrogenase", "oxidase", "regulator"]
    for i, (gs, ge, name) in enumerate(zip(gene_starts, gene_ends, gene_names)):
        feat = SeqFeature(FeatureLocation(bgc_start + gs, bgc_start + ge, strand=1),
                          type="CDS",
                          qualifiers={"gene": name, "product": f"{name} protein",
                                     "protein_id": f"SYN_{i+1}"})
        features.append(feat)

    # Add the cluster feature
    cluster_feat = SeqFeature(FeatureLocation(bgc_start, bgc_end, strand=1),
                              type="cluster",
                              qualifiers={"label": "typeI_PKS_BGC",
                                         "product": "Type I PKS"})
    features.append(cluster_feat)

    # Write files
    write_fasta(genome, test_data_dir / "synthetic_genome.fasta",
                "synthetic_genome", "Synthetic genome containing BGC")
    write_fasta(genome, test_data_dir / "synthetic_bgc.fasta",
                "synthetic_bgc", "Full sequence containing BGC (same as genome)")
    write_fasta(vector, test_data_dir / "synthetic_vector.fasta",
                "pSYN_CATCH", "Synthetic capture vector")
    write_genbank(genome, test_data_dir / "synthetic_bgc.gbk",
                  "synthetic_bgc", features)

    print(f"Generated synthetic test data in {test_data_dir}/")
    print(f"  Genome/BGC: {len(genome)} bp, BGC at {bgc_start}-{bgc_end} ({len(bgc)} bp)")
    print(f"  Vector: {len(vector)} bp")
    print(f"  GC content: genome={sum(1 for b in genome if b in 'GC')/len(genome)*100:.1f}%, "
          f"vector={sum(1 for b in vector if b in 'GC')/len(vector)*100:.1f}%")


if __name__ == "__main__":
    main()
