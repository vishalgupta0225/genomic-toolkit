import streamlit as st
import pandas as pd
from dna_toolkit import DNA_Toolkit # Importing your class

# 1. Page Configuration
st.set_page_config(
    page_title="Genomic Toolkit",
    page_icon="🧬",
    layout="centered"
)

# 2. Title and Description
st.title("🧬 Genomic Analysis Toolkit")
st.write("A Python-based tool for analyzing DNA sequences, calculating molecular weights, and finding enzymatic cut sites.")

# 3. Sidebar Input
st.sidebar.header("Input Sequence")
sequence_input = st.sidebar.text_area("Paste DNA Sequence (FASTA or Raw):", height=200, value="ATGCGTGAATTCACCGTTTAG")

# Clean the input (remove headers if FASTA, remove newlines)
if ">" in sequence_input:
    sequence_input = "".join(sequence_input.split("\n")[1:]) # Remove first line
clean_seq = sequence_input.upper().replace("\n", "").replace(" ", "")

# 4. Main Application Logic
if st.button("Analyze Sequence"):
    if clean_seq:
        try:
            # Initialize your toolkit
            tool = DNA_Toolkit(clean_seq)
            
            # --- SECTION A: Basic Info ---
            st.subheader("1. General Statistics")
            col1, col2, col3 = st.columns(3)
            col1.metric("Length", f"{len(clean_seq)} bp")
            col2.metric("GC Content", f"{tool.gc_content():.2f}%")
            col3.metric("Protein Weight", f"{tool.protein_weight()} Da")
            
            # --- SECTION B: Nucleotide Composition (Chart) ---
            st.subheader("2. Nucleotide Composition")
            counts = tool.count_nucleotides()
            chart_data = pd.DataFrame.from_dict(counts, orient='index', columns=['Count'])
            st.bar_chart(chart_data)

            # --- SECTION C: Sequence Outputs ---
            st.subheader("3. Molecular Products")
            
            with st.expander("View RNA Sequence"):
                st.code(tool.transcribe(), language='text')
                
            with st.expander("View Protein Sequence"):
                st.code(tool.translate(), language='text')
                
            with st.expander("View Reverse Complement"):
                st.code(tool.reverse_complement(), language='text')

            # --- SECTION D: Enzymatic Analysis ---
            st.subheader("4. Restriction Map (EcoRI)")
            sites = tool.find_enzyme_sites("GAATTC")
            if isinstance(sites, list):
                st.success(f"EcoRI sites found at indices: {sites}")
            else:
                st.info("No EcoRI sites found.")

        except ValueError as e:
            st.error(f"Error: {e}")
    else:
        st.warning("Please enter a sequence.")