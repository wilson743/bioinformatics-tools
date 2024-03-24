import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight, GC123
# import matplotlib.pyplot as plt
from collections import Counter
from datetime import datetime
# from pathlib import Path
import pandas as pd
import os
import tempfile

fileName = ''

now = datetime.now()
st.title("Genome Analysis Tool")
st.subheader(f"""
    Report Section:  {now.strftime("%d/%m/%Y %H:%M:%S")}
""")
st.sidebar.subheader("Parameter Tuning Section")
uploaded_file = st.sidebar.file_uploader("Upload a sequence file")
if uploaded_file:
    temp_dir = tempfile.mkdtemp()
    fileName = os.path.join(temp_dir, uploaded_file.name)
    with open(fileName, "wb") as f:
         f.write(uploaded_file.getvalue())


entity_name = st.sidebar.text_input("Enter descriptive name for entity being studied")
filetype = st.sidebar.selectbox("Select file type",['Fasta','gb'])
genotype = st.sidebar.selectbox("Select genotype",['DNA','RNA'])
protein_length = st.sidebar.slider("Minimum protein length",min_value=0, max_value=100,step=1)
# report = st.sidebar.checkbox("Generate report")
analyse = st.sidebar.button('Begin Analysis',help="Click to begin the anaylsis")

def action(filetype,fileName,genotype,entity_name,protein_length):

    #
    # st.write(fileName)
    # st.write(str(p))
    if filetype == 'Fasta' or filetype == 'gb':
        data = SeqIO.read(fileName,filetype.lower())

        if genotype.upper() == 'DNA':
            DNA = data.seq

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
                if len(protein) < protein_length:
                    proteins_seq.remove(protein)

            total_useful_proteins = len(proteins_seq)
            proteins_sorted = sorted(proteins_seq, key=len)

            #presentation
            df = pd.DataFrame.from_dict(nucleotide_count,orient='index')
            aa_df = pd.DataFrame.from_dict(individual_amino_acid_count,orient='index')

            with st.expander("General result"):
                st.text(f"Name: {entity_name}")
                st.text(f"Accession No: {data.name}\n")
                st.text(f"Description: {data.description}\n")
                st.text(f"Genome size: {genome_length}\n")
                st.text(f"Molecular weight: {dna_weight}\n")
                st.text(f"GC content: {gc_content}\n")
                st.text(f"Nucleotide count: {nucleotide_count}\n")
                st.text(f"Total number of amino acids: {total_amino_acids}\n")
                st.text(f"Amino acid count: {individual_amino_acid_count}\n")
                st.text(f"Protein sequences: {total_potential_proteins}\n")
                st.text(f"Proteins with amino acid > {protein_length}: {total_useful_proteins}\n")

            with st.expander("View tables"):
                nucleotide_col, amino_acid_col = st.columns(2)
                with nucleotide_col:
                    st.table(df)
                with amino_acid_col:
                    st.table(aa_df)



            with st.expander("View graphs"):
                # st.write(df)
                st.pyplot(df.plot.barh().figure)
                st.pyplot(aa_df.plot.barh().figure)


            with st.expander("View protein sequences"):
                id =1
                for i in proteins_sorted[:]:
                    st.write(f"[*]Protein {id}, {len(i)} amino acids {i}")
                    id += 1


if analyse:
    action(filetype,fileName,genotype,entity_name,int(protein_length))
