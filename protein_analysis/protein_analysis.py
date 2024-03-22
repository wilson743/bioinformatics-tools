from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
from Bio.SeqUtils import molecular_weight, GC123
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
from collections import Counter
import argparse


parser = argparse.ArgumentParser(description="Protein analysis tool")

parser.add_argument('--seq',metavar='amino sequence', type=str)
parser.add_argument('--name', metavar='name of protein', type=str)
parser.add_argument('--pH',metavar='Specify a pH is relevant for analysis', type=int)


args = parser.parse_args()

amino_acid_sequence = ProteinAnalysis(args.seq)

#presentation
print(f"Name: {args.name}\n")
print(f"Length: {amino_acid_sequence.length}\n")
print(f"Amino acid count: {amino_acid_sequence.count_amino_acids()}\n")
print(f"Amino acid percentage: {(amino_acid_sequence.get_amino_acids_percent())}\n")
print(f"Molecular Weight: {amino_acid_sequence.molecular_weight()} Da\n")
print(f"Isoelectric Point: {amino_acid_sequence.isoelectric_point()}\n")
print(f"Instability Index: {amino_acid_sequence.instability_index()}\n")
print(f"Aromacity: {amino_acid_sequence.aromaticity()}\n")
print(f"Gravy: {amino_acid_sequence.gravy()}\n")
print(f"Charge at pH {args.pH}: {amino_acid_sequence.charge_at_pH(args.pH)}\n")


