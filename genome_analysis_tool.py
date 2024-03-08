from Bio import SeqIO
from Bio.SeqUtils import molecular_weight, GC123
import matplotlib.pyplot as plt
import argparse
from collections import Counter



parser = argparse.ArgumentParser(description= "Input a fasta file for analysis")

parser.add_argument('--file', metavar='FILE', type=str, help="Genome data file")
parser.add_argument('--type', metavar='type', type=str, help='fasta or gb',
                  choices=['fasta','gb'])
parser.add_argument('--name', metavar='NAME', type=str, help="Name of the entity being analysed")
parser.add_argument('--genetics', metavar='DNA or RNA', type=str, choices=['dna','rna'], help="Indicate whether this is a DNA or RNA sample  " )
parser.add_argument('--sense', metavar='+/-', type=str, choices=['+','-'], help="This indicates whether the RNA is a positive or negative")
parser.add_argument('--protein_length', metavar="Protein length", default=50, type=int,help="Minimum amino acid sequence that codes for protein")

args = parser.parse_args()

if args.type == 'fasta' or args.type == 'gb':
    data = SeqIO.read(args.file,args.type)

    if args.genetics.lower() == 'dna':

        DNA = data.seq

        # DNA analysis
        genome_length = len(DNA)
        dna_weight = molecular_weight(DNA)
        gc_content = GC123(DNA)[0]

        nucleotide_count = {
            'A': DNA.count('A'),
            'T': DNA.count('T'),
            'G': DNA.count('G'),
            'C': DNA.count('C')
        }

        # Transcription and Translation
        mRNA = DNA.transcribe()
        amino_acid_seq = mRNA.translate()
        total_amino_acids = len(amino_acid_seq)
        individual_amino_acid_count = Counter(amino_acid_seq)
        proteins_seq = amino_acid_seq.split('*')
        total_potential_proteins = len(proteins_seq)

        for protein in proteins_seq[:]:
            if len(protein) < args.protein_length:
                proteins_seq.remove(protein)

        total_useful_proteins = len(proteins_seq)
        proteins_sorted = sorted(proteins_seq, key=len)


        #presentation
        print(f"Name: {args.name}, {args.genetics.upper()}\n")
        print(f"Accession No: {data.name}\n")
        print(f"Description: {data.description}\n")
        print(f"Genome size: {genome_length}\n")
        print(f"Molecular weight: {dna_weight}\n")
        print(f"GC content: {gc_content}\n")
        print(f"Nucleotide count: {nucleotide_count}\n")
        print(f"Total number of amino acids: {total_amino_acids}\n")
        print(f"Amino acid count: {individual_amino_acid_count}\n")
        print(f"Protein sequences: {total_potential_proteins}\n")
        print(f"Proteins with amino acid > {args.protein_length}: {total_useful_proteins}\n")

        for i in proteins_sorted[:]:
            print(f"[*]Protein: {i}")
            print(f"[#]Length: {len(i)} amino acids \n")


        # Graphs
        # width = 0.5
        # plt.bar(nucleotide_count.keys(), nucleotide_count.values(), width, color=['b','r','g','c'])
        # plt.xlabel("Nucleotides")
        # plt.ylabel("Frequency")
        # plt.title(f" Frequency of Nucleotides in {args.name} Genome")
        # plt.plot()
        # plt.show()
        #
        #
        # # Plot a graph of the amino acid distribution
        # width = 0.5
        # plt.bar(individual_amino_acid_count.keys(), individual_amino_acid_count.values(), width,color= ['b','r','g','c'])
        # plt.xlabel("Amino acids")
        # plt.ylabel("Frequency")
        # plt.title("Frequency of Amino acids")
        # plt.plot()
        # plt.show()


    if args.genetics.lower() == 'rna':
        RNA = data.seq
        rna_genome_length = len(RNA)
        #rna_weight = molecular_weight(RNA)
        rna_gc_content = GC123(RNA)[0]

        rna_nucleotide_count = {
            'A': RNA.count('A'),
            'T': RNA.count('T'),
            'G': RNA.count('G'),
            'C': RNA.count('C')
        }

        rna_amino_acid_seq = RNA.translate()
        rna_total_amino_acids = len(rna_amino_acid_seq)
        rna_individual_amino_acid_count = Counter(rna_amino_acid_seq)
        rna_proteins = rna_amino_acid_seq.split('*')
        rna_total_potential_proteins = len(rna_proteins)

        for protein in rna_proteins[:]:
            if len(protein) < 50:
                rna_proteins.remove(protein)

        rna_total_useful_proteins = len(rna_proteins)
        rna_proteins_sorted = sorted(rna_proteins, key=len)

        print(f"Name: {args.name}, {args.sense}{args.genetics.upper()}\n")
        print(f"Genome size: {rna_genome_length}\n")
       # print(f"Molecular weight: {rna_weight}\n")
        print(f"GC content: {rna_gc_content}\n")
        print(f"Nucleotide count: {rna_nucleotide_count}\n")
        print(f"Total number of amino acids: {rna_total_amino_acids}\n")
        print(f"Amino acid count: {rna_individual_amino_acid_count}\n")
        print(f"Total potential proteins: {rna_total_potential_proteins}\n")
        print(f"Total useful proteins: {rna_total_useful_proteins}\n")



