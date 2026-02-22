import streamlit as st
import pandas as pd
from dna_toolkit import DNA_Toolkit

# 1. Page Configuration
st.set_page_config(
    page_title="Genomic Toolkit",
    page_icon="🧬",
    layout="centered"
)

st.title("🧬 Genomic Analysis Toolkit")
st.write("Select a feature from the sidebar to analyze your sequence.")

# ==========================================
# SIDEBAR: INPUT & MENU
# ==========================================
st.sidebar.header("1. Input Sequence")
seq_type = st.sidebar.radio(
    "Sequence Type:",
    ("mRNA / cDNA (Coding)", "Genomic DNA (Introns)")
)

sequence_input = st.sidebar.text_area("Paste DNA Sequence:", height=150, value="ATGCGTGAATTCACCGTTTAG")

# Clean the input sequence
clean_seq = sequence_input.upper().replace("\n", "").replace(" ", "")
if ">" in sequence_input: # Basic FASTA handling
    clean_seq = "".join(sequence_input.split("\n")[1:]).upper().replace(" ", "")

st.sidebar.divider()

st.sidebar.header("2. Choose Feature")
# This radio button acts as our navigation menu to separate the features
feature_choice = st.sidebar.radio(
    "Select the tool you want to use:",
    [
        "📊 Basic Stats & GC Content", 
        "🔥 Melting Temperature (Tm)", 
        "🔄 Transcription & Translation", 
        "✂️ Restriction Site Mapping",
        "🧬 Primer Evaluation Tool"
    ]
)

# ==========================================
# MAIN PAGE: SEPARATED FEATURES
# ==========================================
if clean_seq:
    try:
        # Load the sequence into our biology engine
        tool = DNA_Toolkit(clean_seq)
        
        st.subheader(feature_choice)

        # --- FEATURE 1: BASIC STATS ---
        if feature_choice == "📊 Basic Stats & GC Content":
            st.write("Calculate the length, GC percentage, and nucleotide composition of your sequence.")
            if st.button("Calculate Stats"):
                col1, col2 = st.columns(2)
                col1.metric("Sequence Length", f"{len(clean_seq)} bp")
                col2.metric("GC Content", f"{tool.gc_content():.2f}%")
                
                st.write("**Nucleotide Composition:**")
                counts = tool.count_nucleotides()
                chart_data = pd.DataFrame.from_dict(counts, orient='index', columns=['Count'])
                st.bar_chart(chart_data)

        # --- FEATURE 2: MELTING TEMP ---
        elif feature_choice == "🔥 Melting Temperature (Tm)":
            st.write("Calculate the theoretical melting temperature of your sequence.")
            if st.button("Calculate Tm"):
                st.metric("Melting Temperature", f"{tool.melting_temperature()} °C")
                st.info("Note: This uses the Wallace rule for short sequences (<14 bp) and a standard GC% formula for longer sequences.")

        # --- FEATURE 3: TRANSCRIPTION & TRANSLATION ---
        elif feature_choice == "🔄 Transcription & Translation":
            st.write("Convert your DNA into RNA, find the reverse complement, and translate it into a protein.")
            if st.button("Convert Sequence"):
                st.markdown("**RNA Sequence:**")
                st.code(tool.transcribe(), language='text')
                
                st.markdown("**Reverse Complement:**")
                st.code(tool.reverse_complement(), language='text')

                if seq_type == "mRNA / cDNA (Coding)":
                    st.markdown("**Protein Sequence:**")
                    st.code(tool.translate(), language='text')
                    st.metric("Estimated Protein Weight", f"{tool.protein_weight()} Da")
                else:
                    st.warning("💡 Translation to protein is disabled for Genomic DNA because it contains non-coding introns. Change the sequence type in the sidebar to enable.")

        # --- FEATURE 4: RESTRICTION SITES ---
        elif feature_choice == "✂️ Restriction Site Mapping":
            st.write("Scan your sequence for common restriction enzyme recognition sites.")
            if st.button("Find Cut Sites"):
                sites = tool.find_all_restriction_sites()
                
                if sites:
                    for enzyme, positions in sites.items():
                        st.success(f"**{enzyme}** site(s) found at position(s): {positions}")
                else:
                    st.info("No common restriction sites found in this sequence.")

        # --- FEATURE 5: PRIMER EVALUATION ---
        elif feature_choice == "🧬 Primer Evaluation Tool":
            st.write("Check if your sequence meets the standard design rules for a good PCR primer.")
            
            if st.button("Evaluate Primer"):
                eval_results = tool.evaluate_primer()
                
                st.subheader("Primer Quality Report")
                
                # Loop through the results and display them nicely
                for criteria, data in eval_results.items():
                    value = data["Value"]
                    status = data["Status"]
                    
                    if status == "Good":
                        st.success(f"**{criteria}:** {value} (✅ Good)")
                    else:
                        st.warning(f"**{criteria}:** {value} (⚠️ {status})")            

    except ValueError as e:
        st.error(f"Error: {e}")
else:
    st.info("👈 Please paste a valid DNA sequence in the sidebar to begin.")
