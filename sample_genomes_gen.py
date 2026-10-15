import os

def create_sample_genomes(base_dir: str):
    os.makedirs(base_dir, exist_ok=True)

    # 1. MDR Klebsiella pneumoniae (NDM-1, CTX-M-15, aac(6')-Ib-cr, sul1, dfrA12, gyrA S83L)
    kpn_seq = (
        ">contig_1_chromosome [organism=Klebsiella pneumoniae] [strain=KPN-NDM1-MDR] [topology=circular]\n"
        "ATGGCGACCGTTCTGCTGGTCGATGAAGAGCCGATGCTGGCCGAACTGAACATTCGCCTG\n"
        "TACCTTTCTCGCCGTCGTTAATGGTAGCCACGGCGGTCAGCCGCTGGCCGCTGCAACCA\n" # contains gyrA S83L (TACCTTTCT)
        "GCATGGGCGTACTGACCGTGCTGGGCGAAGAAGTGCTGGCGATCGTCAATGAGCAGTTCC\n"
        ">contig_2_plasmid_pNDM1 [plasmid=pNDM1] [resistance_genes=blaNDM-1,blaCTX-M-15,aac(6')-Ib-cr]\n"
        "ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCTGAGCACCGCATTAGCCG\n" # contains blaNDM-1
        "CTGCATTGATGCTGAGCGGGCGCCGTCGTTAATGGTAGCCACGGCGGTCAGCCGCTGGCC\n" # contains aac(6')-Ib-cr
        "ATGTTCGCGATACGCTGTTGTTAGGAAGTGTGCCGCTGTAGTGTTGCTGCTGAATCAGCA\n" # contains blaCTX-M-15
        "ATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGACGTATGCGCC\n" # contains sul1
        "ATGAACACGTTTATTGTATGCGCCGTCGTTAATGGTAGCCACGGCGGTCAGCCGCTGGCC\n" # contains dfrA12
        ">contig_3_plasmid_pTEM [plasmid=pTEM1]\n"
        "ATGAAGTTAAATAGCCCTGCTTTTATTTTGTTGCTGTATTTCTTGTTTTGTTTTGTTGCT\n" # contains blaTEM-1
    )
    with open(os.path.join(base_dir, "k_pneumoniae_ndm1.fasta"), "w") as f:
        f.write(kpn_seq)

    # 2. Susceptible E. coli K-12 (wild-type, no acquired AMR genes)
    eco_seq = (
        ">NC_000913.3 Escherichia coli str. K-12 substr. MG1655 complete genome\n"
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC\n"
        "TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG\n"
        "TCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTA\n"
        "CACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAG\n"
        "GTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCG\n"
        "GGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCG\n"
        "GTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATG\n"
        "CCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGG\n"
        "TGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCG\n"
        "AACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGTTGG\n"
    )
    with open(os.path.join(base_dir, "e_coli_k12_wildtype.fasta"), "w") as f:
        f.write(eco_seq)

    # 3. Carbapenem-Resistant Acinetobacter baumannii (blaOXA-23, armA pan-aminoglycoside)
    ab_seq = (
        ">contig_1_AB_chromosome [organism=Acinetobacter baumannii] [strain=AB01-CRAB]\n"
        "ATGAAACAAATCAACATCACAAAGAAACTTGTTGCTGCTGTAGCGTTGTCCTTAGTTAGC\n"
        "TTTACAGCTGCTTCAGCTGCTGATCAAACGCAAGTCACAGTAACTAAAGTTACAGGTGCA\n"
        ">contig_2_AB_plasmid_pOXA23 [plasmid=pOXA23] [genes=blaOXA-23,armA]\n"
        "ATGAATAAATATTTTACTTGTTTAGCTTTATTTTCTTGTAATGATTTTTTTTATTTTCAG\n" # contains blaOXA-23 & armA
        "CGTGCTGCTGTTAGTGTTAGTTAGTGTTAGTTAGTGTTAGTTAGTGTTAGTTAGTGTTAG\n"
    )
    with open(os.path.join(base_dir, "a_baumannii_crab_oxa23.fasta"), "w") as f:
        f.write(ab_seq)

    # 4. Colistin-Resistant Klebsiella pneumoniae (mcr-1)
    col_seq = (
        ">contig_1_KPN_COL [organism=Klebsiella pneumoniae] [strain=KPN-MCR1]\n"
        "ATGGCGACCGTTCTGCTGGTCGATGAAGAGCCGATGCTGGCCGAACTGAACATTCGCCTG\n"
        ">contig_2_plasmid_pmcr1 [plasmid=pHNSHP45] [resistance_gene=mcr-1]\n"
        "ATGATGCAGCATACTTCTGTGTGGTACCGACGCTCGTCGTTTTGTCAGCGCCAGTGTTGC\n" # contains mcr-1
        "TGCTGAATCAGCAAGCCGTTTCGGCGCAAAGTGCGTTATATTCGCCTGTGTATTATCTCT\n"
    )
    with open(os.path.join(base_dir, "k_pneumoniae_mcr1_colistin.fasta"), "w") as f:
        f.write(col_seq)

    # 5. Malformed non-FASTA file
    with open(os.path.join(base_dir, "malformed_non_fasta.txt"), "w") as f:
        f.write("This is a plain text file without FASTA formatting headers.\nLine 2 contains random text.\n")

    # 6. Empty genome file
    with open(os.path.join(base_dir, "empty_genome.fasta"), "w") as f:
        pass

if __name__ == "__main__":
    target = os.path.join(os.path.dirname(__file__), "..", "..", "sample_data", "genomes")
    create_sample_genomes(target)
    print("Sample genomes created successfully in", target)
