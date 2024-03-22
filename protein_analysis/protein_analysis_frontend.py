import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from stmol import showmol
import py3Dmol
import requests


def render_protein(pdb):
    st.success("Prediction done successfully")
    view = py3Dmol.view()
    view.addModel(pdb,'pdb')
    view.setStyle({'cartoon':{'color':'spectrum'}})
    view.setBackgroundColor('white')
    view.zoomTo()
    view.zoom(2,800)
    view.spin(True)
    showmol(view, height=500, width=650)


header = {
    'Content-Type': 'application/x-www-form-urlencoded',
}

st.title("Protein Analysis Tool")
st.subheader("Report Section")
st.sidebar.subheader("Parameter Tuning Section")
st.sidebar.text("""
    This program uses the ESM Fold algorithm 
    to predict the structure of a protein, 
    given the amino acid sequence
""")


st.sidebar.warning("Structural prediction requires internet ")
job_name = st.sidebar.text_input("Enter job name")
sequence = st.sidebar.text_area("Enter your protein sequence",height=150)
# options = st.sidebar._multiselect("Properties of interest",['Length','Molecular weight','Isoelectric point', 'Instability index','Aromacity'
#                                                             'Gravy'])

# pH = st.sidebar.slider("Select your desired pH",min_value=0,max_value=14)
# algorithm = st.sidebar.selectbox('Structural Prediction Algorithm',['EMS Fold', 'Alphafold'])
analyse = st.sidebar.button("Analyse")
amino_acid_sequence = ''


if analyse:
    amino_acid_sequence = ProteinAnalysis(sequence)
    df = pd.DataFrame.from_dict(amino_acid_sequence.count_amino_acids(), orient='index')
    try:
        with st.expander("General result"):
            st.text(f"Name: {job_name}\n")
            st.text(f"Length: {amino_acid_sequence.length}\n")
            # st.write(f"Amino acid count: {amino_acid_sequence.count_amino_acids()}\n")
            # st.write(f"Amino acid percentage: {(amino_acid_sequence.get_amino_acids_percent())}\n")
            st.text(f"Molecular Weight: {amino_acid_sequence.molecular_weight()} Da\n")
            st.text(f"Isoelectric Point: {amino_acid_sequence.isoelectric_point()}\n")
            st.text(f"Instability Index: {amino_acid_sequence.instability_index()}\n")
            st.text(f"Aromacity: {amino_acid_sequence.aromaticity()}\n")
            st.text(f"Gravy: {amino_acid_sequence.gravy()}\n")
            # st.text(f"Charge at pH {pH}: {amino_acid_sequence.charge_at_pH(pH)}\n")
            st.table(df)


    except Exception:
        pass


    with st.expander("View Protein Structure"):
        st.markdown("### Predicted Protein Structure")
        c = 1
        while c == 1:
            with st.spinner("Please wait"):
                try:
                    response = requests.post(url='https://api.esmatlas.com/foldSequence/v1/pdb',headers=header, data=sequence, verify=False)
                    pdb_data = response.content.decode('utf-8')
                    render_protein(pdb_data)
                    st.download_button(label="Download Structure", data=pdb_data, file_name=f'{job_name}_structure.pdb', mime='text/plain')
                except Exception:
                    st.error('Unable to connect to server, Make sure you are connected to the internet')
            c = 0



