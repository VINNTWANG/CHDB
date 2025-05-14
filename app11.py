#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CHD1-9 ClinVar & HPO Explorer: Interactive Positional Analysis (Streamlit Version)
Modifications:
- v3.8.9: Removed HPO phenotype frequency bar chart from Overview, keeping only the table.
- v3.8.8: Fixed AttributeError for fig in comparison plots with empty data. Enlarged data source annotation. Made "Show Detailed Plots" checkbox more prominent.
- v3.8.7: Replaced use_column_width with use_container_width. Adjusted legend positions. Added domain details expander. Addressed domain track x-axis label overlap.
- v3.8.6: Adjusted sidebar logo padding. Refined HPO tab layout (replaced single-bar charts with metrics/text). Ensured lollipop plot color consistency.
- v3.8.5: Fixed Sankey diagram link color alpha issue by converting hex+alpha to rgba string.
- v3.8.4: Enhanced table styling for a "Nature-inspired" look. Optimized HPO tab with a new top phenotypes table.
- v3.8.3: Implemented more aggressive hiding of x-axis labels in binned plots for very dense scenarios (num_bins > 20). UI styling remains consistent.
- v3.8.2: Integrated binned plot x-axis label improvements. Enhanced UI styling for a cleaner, more modern "Salk-inspired" aesthetic.
- v3.8.1: Restored all previously defined plotting functions that were accidentally omitted.
- v3.8: Attempt at 'Salk Institute' inspired aesthetic improvements; refined binned plot X-axis further.
# ... (all previous modifications) ...
"""

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import json
import re
import os
from collections import Counter
from datetime import datetime
from wordcloud import WordCloud
import networkx as nx
from PIL import Image

# ────────────────────────── PAGE CONFIGURATION (MUST BE FIRST STREAMLIT COMMAND) ─────────────────
try:
    logo_path_for_icon = "logo.png"
    if os.path.exists(logo_path_for_icon):
        page_icon_img = Image.open(logo_path_for_icon)
    else:
        page_icon_img = "🧬" # Fallback icon
except Exception:
    page_icon_img = "🧬"

st.set_page_config(
    layout="wide",
    page_title="CHD Gene Family Variant & HPO Explorer",
    page_icon=page_icon_img
)

# --- Custom CSS for "Salk-inspired" & "Nature-inspired" aesthetics ---
st.markdown("""
<style>
    /* General Body & Font */
    body {
        font-family: 'Roboto', 'Helvetica Neue', 'Arial', sans-serif;
        font-size: 15px;
        color: #4A4A4A;
        background-color: #FDFDFD;
        line-height: 1.6;
    }

    /* Headings - Clear hierarchy */
    .stTitle, h1 {
        color: #333;
        font-weight: 600;
        letter-spacing: -0.5px;
        margin-bottom: 0.8em;
    }
    .stHeader, h2 {
        color: #3A3A3A;
        font-weight: 500;
        margin-top: 2.5em;
        margin-bottom: 1em;
        padding-bottom: 0.4em;
        border-bottom: 1px solid #EAEAEA;
    }
    .stSubheader, h3 {
        color: #404040;
        font-weight: 500;
        margin-top: 2em;
        margin-bottom: 0.8em;
        letter-spacing: -0.25px;
    }
    h5 { /* Used for plot group titles */
        color: #505050;
        font-weight: 500;
        margin-top: 1.8rem;
        margin-bottom: 0.6rem;
        padding-bottom: 0.3em;
        border-bottom: 1px dotted #D0D0D0;
    }

    /* Sidebar Styling */
    .stSidebar > div:first-child {
        background-color: #F8F9FA;
        padding: 15px 25px 25px 25px; /* Reduced top padding */
        border-right: 1px solid #E0E0E0;
    }
    .stSidebar .stImage { 
        margin-top: 0; 
        margin-bottom: 1rem; 
        text-align: center;
    }
    .stSidebar .stRadio > label,
    .stSidebar .stSlider > label,
    .stSidebar .stMultiSelect > label {
        font-size: 1.0em;
        font-weight: 500;
        color: #333;
        margin-bottom: 0.6rem !important;
    }
    .stSidebar .stRadio div[role="radiogroup"] > label {
        font-size: 0.95em;
        transition: color 0.2s ease-in-out;
    }
    .stSidebar .stRadio div[role="radiogroup"] > label:hover {
        color: #007AFF;
    }
    .stSidebar .stRadio div[role="radiogroup"] > label input:checked + div {
        color: #007AFF; 
        font-weight: 500;
    }
    .stSidebar .stRadio div[role="radiogroup"] > label input:checked + div::before {
        border-color: #007AFF !important;
        background-color: #007AFF !important;
    }
    .stSidebar .stSlider > div[data-baseweb="slider"] > div:nth-child(2) {
        background: #D0D0D0;
    }
    .stSidebar .stSlider > div[data-baseweb="slider"] > div:nth-child(3) {
        background: #007AFF;
    }
    .stSidebar .stSlider > div[data-baseweb="slider"] > div:nth-child(4) {
        border: 2px solid #007AFF;
        background: white;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    .stSidebar .stMultiSelect [data-baseweb="tag"] {
        background-color: #E5F1FF;
        color: #0056b3;
        border-radius: 4px;
        font-size: 0.9em;
    }
    .stSidebar .stMultiSelect [data-baseweb="select"] > div:first-child {
        border-radius: 6px;
        border: 1px solid #D0D0D0;
    }
    .stSidebar .stMultiSelect [data-baseweb="select"] > div:first-child:hover {
        border-color: #B0B0B0;
    }
    .stSidebar .stMultiSelect [data-baseweb="select"] > div[aria-expanded="true"] {
        border-color: #007AFF;
        box-shadow: 0 0 0 2px rgba(0, 122, 255, 0.2);
    }
    .sidebar-footer {
        font-size: 0.9em;
        color: #555;
        padding-top: 20px;
        border-top: 1px solid #D0D0D0;
        margin-top: 25px;
    }
    .sidebar-footer b {
        color: #222;
    }
    .sidebar-footer a {
        color: #007AFF;
        text-decoration: none;
    }
    .sidebar-footer a:hover {
        text-decoration: underline;
        color: #0056b3;
    }

    /* Tabs Styling */
    div[data-baseweb="tab-list"] button[data-baseweb="tab"] {
        font-size: 1.05em;
        padding: 10px 18px;
        font-weight: 500;
        color: #555;
        border-bottom: 2px solid transparent;
        transition: color 0.2s ease-in-out, border-color 0.2s ease-in-out;
    }
    div[data-baseweb="tab-list"] button[data-baseweb="tab"][aria-selected="true"] {
        color: #007AFF;
        border-bottom: 2px solid #007AFF;
        font-weight: 600;
    }

    /* Input Widgets (like search) */
    .stTextInput > div > div > input, .stNumberInput > div > div > input {
        border-radius: 6px;
        border: 1px solid #D0D0D0;
        padding: 8px 10px;
        transition: border-color 0.2s ease-in-out, box-shadow 0.2s ease-in-out;
    }
    .stTextInput > div > div > input:focus, .stNumberInput > div > div > input:focus {
        border-color: #007AFF;
        box-shadow: 0 0 0 2px rgba(0, 122, 255, 0.2);
    }

    /* Button Styling */
    .stButton > button {
        border-radius: 6px;
        padding: 0.4em 1em;
        font-weight: 500;
        transition: background-color 0.2s ease-in-out, border-color 0.2s ease-in-out, color 0.2s ease-in-out;
    }
    .stButton > button:not([kind="secondaryLink"]):not([kind="icon"]) {
        background-color: #007AFF;
        color: white;
        border: 1px solid #007AFF;
    }
    .stButton > button:not([kind="secondaryLink"]):not([kind="icon"]):hover {
        background-color: #0056b3;
        border-color: #0056b3;
    }
    .stButton > button:not([kind="secondaryLink"]):not([kind="icon"]):active {
        background-color: #004085;
        border-color: #004085;
    }
    .stButton > button:not([kind="secondaryLink"]):not([kind="icon"]):focus {
        box-shadow: 0 0 0 2px rgba(0, 122, 255, 0.3);
    }
    .stButton > button[kind="secondaryLink"] {
        border: 1px solid #007AFF !important;
        color: #007AFF !important;
        background-color: transparent !important;
    }
    .stButton > button[kind="secondaryLink"]:hover {
        background-color: #E5F1FF !important;
        color: #0056b3 !important;
        border-color: #0056b3 !important;
    }

    /* Expander (for search results) */
    .stExpander {
        border: 1px solid #EAEAEA;
        border-radius: 8px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.03);
        transition: box-shadow 0.2s ease-in-out;
    }
    .stExpander:hover {
        box-shadow: 0 2px 6px rgba(0,0,0,0.05);
    }
    .stExpander header {
        font-weight: 500;
        font-size: 1.0em;
    }

    /* Dataframes -- "Nature Inspired" */
    .stDataFrame {
        border: 1px solid #D0D0D0; 
        border-radius: 4px; 
        box-shadow: none; 
        overflow: hidden;
        font-size: 0.9em; 
    }
    .stDataFrame thead th {
        background-color: #F8F9FA; 
        color: #222; 
        font-weight: 600; 
        border-bottom: 1.5px solid #B0B0B0; 
        border-top: none;
        border-left: none;
        border-right: none;
        text-align: left;
        padding: 10px 12px;
        white-space: nowrap; 
    }
    .stDataFrame tbody td {
        border-bottom: 1px solid #EAEAEA; 
        border-top: none;
        border-left: none; 
        border-right: none; 
        padding: 9px 12px; 
        color: #333;
        vertical-align: middle; 
    }
    .stDataFrame tbody tr:nth-of-type(even) {
        background-color: transparent; 
    }
    .stDataFrame tbody tr:last-child td {
        border-bottom: none; 
    }
    .stDataFrame tbody tr:hover {
        background-color: #EFF6FF; 
    }
    .stDataFrame .ag-cell {
        border-left: none !important;
        border-right: none !important;
        border-top: none !important; 
    }
     .stDataFrame .ag-header-cell {
        border-left: none !important;
        border-right: none !important;
    }


    /* Alert Boxes */
    .stAlert {
        border-radius: 6px;
        border-left-width: 4px;
        padding: 12px 18px;
        font-size: 0.95em;
    }
    .stAlert[data-testid="stNotification"] {
         box-shadow: 0 2px 10px rgba(0,0,0,0.08);
    }

    /* Markdown links in main content */
    .main a {
        color: #007AFF;
        text-decoration: none;
    }
    .main a:hover {
        text-decoration: underline;
        color: #0056b3;
    }
</style>
""", unsafe_allow_html=True)


# ────────────────────────── CONFIGURATION ──────────────────────────
ALL_GENES = [f"CHD{i}" for i in range(1, 10)]
DEFAULT_GENE_LENGTH_FALLBACK = 1000
LOGO_PATH = "logo.png"

GENE_IDENTIFIERS = {
    "CHD1": {"uniprot": "O14646", "refseq_nm": "NM_001270.4"},
    "CHD2": {"uniprot": "O14647", "refseq_nm": "NM_001271.4"},
    "CHD3": {"uniprot": "Q12873", "refseq_nm": "NM_001030059.3"},
    "CHD4": {"uniprot": "Q14839", "refseq_nm": "NM_001273.3"},
    "CHD5": {"uniprot": "Q8TDI0", "refseq_nm": "NM_015264.5"},
    "CHD6": {"uniprot": "Q8TD26", "refseq_nm": "NM_020919.4"},
    "CHD7": {"uniprot": "Q9P2D1", "refseq_nm": "NM_017780.4"},
    "CHD8": {"uniprot": "Q9HCK8", "refseq_nm": "NM_001170629.2"},
    "CHD9": {"uniprot": "Q3L8U1", "refseq_nm": "NM_152283.4"},
}

SIGNIFICANCE_PALETTE = {"Pathogenic": "#D81B60", "Likely pathogenic": "#F06292", "VUS": "#1F78B4", "Benign": "#33A02C", "Likely benign": "#B2DF8A", "Conflicting": "#A6CEE3", "Unknown": "#B3B3B3", "Other": "#FB9A99"}
VARIANT_TYPE_PALETTE = { "SNV": "#4C78A8", "Deletion": "#F58518", "Duplication": "#E45756", "Insertion": "#72B7B2", "Copy number gain": "#54A24B", "Copy number loss": "#B279A2", "Other": "#9DA9C7", "Indel": "#88B04B" }
PHENOTYPE_COLORS = px.colors.qualitative.Pastel1
GENE_COLORS = px.colors.qualitative.Set3
SIGNIFICANCE_MARKERS = {"Pathogenic": "circle", "Likely pathogenic": "triangle-up", "VUS": "square", "Benign": "diamond", "Likely benign": "triangle-down", "Conflicting": "cross", "Unknown": "x", "Other": "star"}
DOMAIN_COLORS = ["#1B9E77", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#D95F02", "#7570B3", "#A0A0A0", "#636EFA", "#EF553B", "#00CC96"] 

PLOT_FONT_FAMILY = "Roboto, Helvetica Neue, Arial, sans-serif"; PLOT_FONT_SIZE_MAIN_TITLE = 17; PLOT_FONT_SIZE_AXIS_TITLE = 12; PLOT_FONT_SIZE_AXIS_TICKS = 9; PLOT_FONT_SIZE_LEGEND_TITLE = 11; PLOT_FONT_SIZE_LEGEND_ITEMS = 10
DOMAIN_TRACK_PAPER_HEIGHT = 0.09; MAIN_PLOT_Y_DOMAIN_BOTTOM = DOMAIN_TRACK_PAPER_HEIGHT + 0.03; 
FIG_HEIGHT = 550; BINNED_PLOT_HEIGHT = 420; GENE_EXPLORER_OTHER_PLOT_HEIGHT = 420; COMPARISON_PLOT_HEIGHT = 480; COMPLEX_PLOT_HEIGHT = 650
PATHOGENIC_SIGNIFICANCES = ["Pathogenic", "Likely pathogenic"]; DEFAULT_TOP_N_PHENOTYPES = 10
TOP_N_PHENOTYPES_GENE_EXPLORER = 7
HPO_PLOT_HEIGHT = 400 
HPO_FILE_DIR = "chd_phenotype_data"
BIN_SIZE_POSITIONAL = 75

# ────────────────────────── HELPER FUNCTIONS ──────────────────────────
def standardize_variant_type(variant_type_str):
    if pd.isna(variant_type_str) or not isinstance(variant_type_str, str): return "Other"
    vt = variant_type_str.lower().strip();
    if not vt: return "Other"
    if any(term in vt for term in ["single nucleotide variant", "snv", "substitution"]): return "SNV"
    if "indel" in vt: return "Indel"
    if "deletion" in vt: return "Deletion"
    if "duplication" in vt: return "Duplication"
    if "insertion" in vt: return "Insertion"
    if "copy number loss" in vt: return "Copy number loss"
    if "copy number gain" in vt: return "Copy number gain"
    if "microsatellite" in vt: return "Microsatellite"
    return "Other"

def classify_clinical_significance(sig_str):
    if pd.isna(sig_str) or not isinstance(sig_str, str): return "Unknown"
    s = sig_str.lower().strip();
    if not s: return "Unknown"
    if "pathogenic/likely pathogenic" in s: return "Pathogenic"
    if "pathogenic" in s and "likely pathogenic" not in s and "benign" not in s : return "Pathogenic"
    if "likely pathogenic" in s and "benign" not in s: return "Likely pathogenic"
    if "benign/likely benign" in s: return "Benign"
    if "benign" in s and "likely benign" not in s and "pathogenic" not in s: return "Benign"
    if "likely benign" in s and "pathogenic" not in s : return "Likely benign"
    if "uncertain significance" in s or "vus" in s: return "VUS"
    if "conflicting interpretations" in s or "conflicting data" in s: return "Conflicting"
    if any(term in s for term in ["drug response", "pharmacogenomic", "risk factor", "association", "protective", "affects"]): return "Other"
    return "Unknown"

def extract_protein_position(name_str):
    if pd.isna(name_str) or not isinstance(name_str, str): return None
    match = re.search(r'p\.(?:[A-Z][a-z]{2})?(\d+)(?:(?:[A-Z][a-z]{2}|\*|=)|(?:ext(?:[A-Z][a-z]{2}|\*)\d+)|(?:fs(?:\*?\d+)?))?', str(name_str))
    return int(match.group(1)) if match else None

def clean_phenotype_list(pheno_str):
    if pd.isna(pheno_str): return "N/A"
    phenos = str(pheno_str).split('|'); expanded_phenos = []
    for p in phenos: expanded_phenos.extend(p.split(';'))
    phenos = [re.sub(r'\[PMID:\d+\]|\(HP:HP:\d+\)', '', p).strip().lower() for p in expanded_phenos]
    phenos = [p for p in phenos if p and p not in ["not specified", "not provided", "n/a", "none", "unspecified"]]
    unique_phenos = sorted(list(dict.fromkeys(phenos)))
    return '; '.join(unique_phenos) if unique_phenos else "N/A"

def simplify_origin(origin_str):
    if pd.isna(origin_str): return "Unknown"
    origin_str = str(origin_str).lower()
    if "de novo" in origin_str: return "De novo";
    if "inherited" in origin_str: return "Inherited";
    if "somatic" in origin_str: return "Somatic";
    if "germline" in origin_str: return "Germline"
    if "unknown" in origin_str or origin_str == ".": return "Unknown"
    return "Other/Mixed"

def hex_to_rgba(hex_color, alpha=0.6): 
    """Converts a hex color string to an rgba string."""
    hex_color = hex_color.lstrip('#')
    hlen = len(hex_color)
    if hlen != 6: 
        return f'rgba(204,204,204,{alpha})' 
    try:
        rgb = tuple(int(hex_color[i:i + 2], 16) for i in range(0, 6, 2))
    except ValueError:
        return f'rgba(204,204,204,{alpha})' 
    return f'rgba({rgb[0]},{rgb[1]},{rgb[2]},{alpha})'

# ────────────────────────── GENE-SPECIFIC DATA LOADING & PREPROCESSING ──────────────────────────
@st.cache_data
def load_gene_data(gene_name):
    base_path = os.path.dirname(os.path.abspath(__file__))
    clinvar_file_path = os.path.join(base_path, f"{gene_name}_clinvar.txt"); domains_file_path = os.path.join(base_path, f"{gene_name}_domains.json")
    df_clinvar_data = pd.DataFrame(); domains_data_list = []; error_messages = []
    # Dummy file creation logic
    if not os.path.exists(clinvar_file_path):
        dummy_max_pos_variant = np.random.randint(500, 2800)
        dummy_clinvar_data = (f"Name\tVariationID\tClinicalSignificance\tPhenotypeList\tOrigin\tAssembly\tType\tPositionVCF\n"
            f"p.Arg123Cys\t{gene_name}1\tPathogenic\t{gene_name} syndrome|Coloboma;CHARGE syndrome\tde novo\tGRCh38\tSNV\t100\n"
            f"p.Ser{dummy_max_pos_variant // 2}Thr\t{gene_name}2\tPathogenic\t{gene_name} syndrome;Intellectual disability\tinherited\tGRCh38\tSNV\t{dummy_max_pos_variant // 2}\n"
            f"p.Gly{dummy_max_pos_variant}Ter\t{gene_name}3\tLikely pathogenic\tColoboma;Other phenotype\tunknown\tGRCh38\tDeletion\t{dummy_max_pos_variant}\n"
            f"p.Val50Phe\t{gene_name}4\tVUS\tTest Pheno\t.\tGRCh38\tSNV\t50\n")
        try:
            with open(clinvar_file_path, 'w', encoding='utf-8') as f: f.write(dummy_clinvar_data)
        except Exception as e: error_messages.append(f"Could not create dummy ClinVar file for {gene_name}: {e}")
    if not os.path.exists(domains_file_path):
        dummy_domain_end = np.random.randint(800, 2900)
        dummy_domains_data = [{"name": f"{gene_name} Dom1", "start": 50, "end": max(100, dummy_domain_end // 3)}, {"name": f"{gene_name} Dom2", "start": max(150, dummy_domain_end // 2), "end": dummy_domain_end}]
        try:
            with open(domains_file_path, 'w', encoding='utf-8') as f: json.dump(dummy_domains_data, f)
        except Exception as e: error_messages.append(f"Could not create dummy domain file for {gene_name}: {e}")
    
    try: df_clinvar_data = pd.read_csv(clinvar_file_path, sep="\t", low_memory=False, encoding="utf-8", na_filter=True)
    except Exception as e: error_messages.append(f"ERROR reading {clinvar_file_path}: {e}"); return pd.DataFrame(columns=['Name', 'Position', 'StandardType', 'SignificanceClass', 'PhenotypeListClean', 'VariationID', 'Gene', 'OriginSimple']), [], DEFAULT_GENE_LENGTH_FALLBACK, error_messages
    df_clinvar_data = df_clinvar_data[df_clinvar_data["Assembly"] == "GRCh38"].copy()
    df_clinvar_data.loc[:, "Position"] = df_clinvar_data["Name"].apply(extract_protein_position)
    df_clinvar_data.loc[:, "StandardType"] = df_clinvar_data["Type"].apply(standardize_variant_type)
    df_clinvar_data.loc[:, "SignificanceClass"] = df_clinvar_data["ClinicalSignificance"].apply(classify_clinical_significance)
    df_clinvar_data.loc[:, "PhenotypeListClean"] = df_clinvar_data["PhenotypeList"].apply(clean_phenotype_list)
    df_clinvar_data.loc[:, "Gene"] = gene_name
    df_clinvar_data.loc[:, "OriginSimple"] = df_clinvar_data["Origin"].apply(simplify_origin) if 'Origin' in df_clinvar_data.columns else "Unknown"
    df_positions = df_clinvar_data.dropna(subset=["Position"]).copy()
    if not df_positions.empty: df_positions.loc[:, "Position"] = df_positions["Position"].astype(int)
    for sig_class in df_clinvar_data["SignificanceClass"].dropna().unique():
        if sig_class not in SIGNIFICANCE_PALETTE: SIGNIFICANCE_PALETTE[sig_class] = "#DDDDDD"
        if sig_class not in SIGNIFICANCE_MARKERS: SIGNIFICANCE_MARKERS[sig_class] = "circle"
    for type_class in df_clinvar_data["StandardType"].dropna().unique():
        if type_class not in VARIANT_TYPE_PALETTE: VARIANT_TYPE_PALETTE[type_class] = "#CCCCCC"
    try:
        with open(domains_file_path, "r", encoding="utf-8") as f: domains_data_raw = json.load(f)
        for i, d_item in enumerate(domains_data_raw):
            if isinstance(d_item, dict) and all(k in d_item for k in ['start', 'end']) and isinstance(d_item['start'], (int, float)) and isinstance(d_item['end'], (int, float)) and d_item['start'] < d_item['end']:
                d_item['name'] = d_item.get('name', f'Domain {i+1}'); domains_data_list.append(d_item)
        domains_data_list.sort(key=lambda x: x['start'])
    except Exception: pass
    max_pos_data = df_positions['Position'].max() if not df_positions.empty and 'Position' in df_positions.columns else 0
    max_pos_domain = max([d['end'] for d in domains_data_list], default=0) if domains_data_list else 0
    true_max_coord = max(max_pos_data, max_pos_domain); calculated_gene_length = int(true_max_coord) if true_max_coord > 0 else DEFAULT_GENE_LENGTH_FALLBACK; calculated_gene_length = max(calculated_gene_length, 10)
    df_cols_to_keep = ['Name', 'Position', 'StandardType', 'SignificanceClass', 'PhenotypeListClean', 'VariationID', 'Gene', 'OriginSimple']
    existing_cols_for_reduction = [col for col in df_cols_to_keep if col in df_clinvar_data.columns]
    df_gene_variants_processed = df_clinvar_data[existing_cols_for_reduction].copy() if not df_clinvar_data.empty else pd.DataFrame(columns=existing_cols_for_reduction)
    return df_gene_variants_processed, domains_data_list, calculated_gene_length, error_messages

@st.cache_data
def load_all_chd_data_for_comparison():
    all_variants_list = []; summary_list = []; errors_dict = {}
    for gene_name in ALL_GENES:
        df_gene, _, gene_len, errors = load_gene_data(gene_name)
        if errors: errors_dict[gene_name] = errors
        if not df_gene.empty: all_variants_list.append(df_gene)
        num_total = len(df_gene); num_patho = len(df_gene[df_gene['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]) if not df_gene.empty else 0
        uniprot_id = GENE_IDENTIFIERS.get(gene_name, {}).get("uniprot", "N/A"); refseq_nm = GENE_IDENTIFIERS.get(gene_name, {}).get("refseq_nm", "N/A")
        summary_list.append({"Gene": gene_name, "UniProt ID": uniprot_id, "RefSeq NM": refseq_nm, "Protein Length (aa)": gene_len, "Total Variants": num_total, "Pathogenic/Likely Path. Variants": num_patho, "Data Loading Status": "Errors" if errors else "OK"})
    combined_df = pd.concat(all_variants_list, ignore_index=True) if all_variants_list else pd.DataFrame()
    comparison_summary_df = pd.DataFrame(summary_list)
    desired_cols_order = ["Gene", "UniProt ID", "RefSeq NM", "Protein Length (aa)", "Total Variants", "Pathogenic/Likely Path. Variants", "Data Loading Status"]
    existing_cols_in_order = [col for col in desired_cols_order if col in comparison_summary_df.columns]
    comparison_summary_df = comparison_summary_df[existing_cols_in_order]
    return combined_df, comparison_summary_df, errors_dict

@st.cache_data
def load_hpo_data(gene_name):
    file_path = os.path.join(HPO_FILE_DIR, f"{gene_name}_annotations.json")
    if os.path.exists(file_path):
        try:
            with open(file_path, "r", encoding='utf-8') as f: data = json.load(f)
            phenotypes = data.get("phenotypes", []); diseases = data.get("diseases", [])
            df_phenotypes = pd.DataFrame(phenotypes) if phenotypes else pd.DataFrame(columns=['id', 'name'])
            df_diseases = pd.DataFrame(diseases) if diseases else pd.DataFrame(columns=['id', 'name', 'mondoId'])
            return df_phenotypes, df_diseases, None
        except Exception as e: return pd.DataFrame(), pd.DataFrame(), f"Error loading HPO data for {gene_name}: {e}"
    else: return pd.DataFrame(), pd.DataFrame(), f"HPO file not found for {gene_name}: {file_path}"

def classify_hpo_phenotype_category(phenotype_name):
    phenotype_name = str(phenotype_name).lower()
    if any(keyword in phenotype_name for keyword in ["autism", "intellectual", "speech", "seizure", "hypotonia", "developmental", "motor", "behavioral", "neuro"]): return "Neurological/Developmental"
    elif any(keyword in phenotype_name for keyword in ["macrocephaly", "midface", "frontal", "eyebrow", "eyelashes", "palpebral", "chin", "craniofacial", "facial", "head", "skull"]): return "Craniofacial/Skeletal"
    elif any(keyword in phenotype_name for keyword in ["cardiac", "heart", "artery", "vascular", "septal defect"]): return "Cardiovascular"
    elif any(keyword in phenotype_name for keyword in ["eye", "vision", "optic", "retina", "ocular", "iris", "coloboma"]): return "Ocular"
    elif any(keyword in phenotype_name for keyword in ["ear", "hearing", "auditory"]): return "Auditory"
    elif any(keyword in phenotype_name for keyword in ["immunodeficiency", "allergy", "immune"]): return "Immune"
    elif any(keyword in phenotype_name for keyword in ["growth", "stature", "short", "tall", "failure to thrive"]): return "Growth"
    elif any(keyword in phenotype_name for keyword in ["genital", "gonadal", "renal", "urinary"]): return "Genitourinary"
    else: return "Other/Systemic"

# ─────────────────── PLOTLY FIGURE CREATION FUNCTIONS (Gene Explorer Tab - ClinVar) ───────────────────
# (All Gene Explorer plotting functions from v3.8.8 are here)
def add_domain_legend_traces(fig, domains_data_local):
    current_domain_colors = DOMAIN_COLORS
    for i, domain in enumerate(domains_data_local):
        color = current_domain_colors[i % len(current_domain_colors)]; name = domain.get("name", f"Domain {i+1}")
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',marker=dict(color=color, size=10, symbol='square'),name=name,legendgroup="domains", legendgrouptitle_text="Protein Domains"))

def create_domain_shapes_paper(domains_data_local, gene_length_local):
    current_domain_colors = DOMAIN_COLORS
    shapes = [dict(type="rect", xref="x", yref="paper",x0=0, y0=0, x1=gene_length_local, y1=DOMAIN_TRACK_PAPER_HEIGHT,fillcolor="#F0F2F6", line_width=0.5, line_color="#D1D5DB", layer="below", opacity=1)]
    if not domains_data_local: return shapes
    for i, domain in enumerate(domains_data_local):
        color = current_domain_colors[i % len(current_domain_colors)]
        start, end = domain["start"], domain["end"]
        name = domain.get("name", f"Domain {i+1}")
        
        if end - start > 0: 
            shapes.append(dict(
                type="rect", xref="x", yref="paper",
                x0=start, y0=DOMAIN_TRACK_PAPER_HEIGHT * 0.1,
                x1=end,   y1=DOMAIN_TRACK_PAPER_HEIGHT * 0.9,
                fillcolor=color, line_color='#374151', line_width=1, 
                layer="below", opacity=0.9,
                name=name, 
            ))
    return shapes

def _setup_single_plot_with_domain_space(fig_title, y_axis_title, y_axis_range, gene_length_local, plot_height=FIG_HEIGHT, show_xaxis_title=True):
    xaxis_title_text = "Protein Position (aa)" if show_xaxis_title else ""
    fig = make_subplots(rows=1, cols=1) 

    fig.update_layout(
        title_text=fig_title, title_font_size=PLOT_FONT_SIZE_MAIN_TITLE, font_family=PLOT_FONT_FAMILY,
        plot_bgcolor='white',
        xaxis=dict(
            title_text=xaxis_title_text, 
            range=[-gene_length_local * 0.02, gene_length_local * 1.02],
            gridcolor='rgba(200,200,200,0.15)', zeroline=False,
            tickfont_size=PLOT_FONT_SIZE_AXIS_TICKS,
            title_font_size=PLOT_FONT_SIZE_AXIS_TITLE,
            title_standoff=20, 
            automargin=True 
        ),
        yaxis=dict(
            title_text=y_axis_title, domain=[MAIN_PLOT_Y_DOMAIN_BOTTOM, 1.0], range=y_axis_range,
            gridcolor='rgba(200,200,200,0.15)', zeroline=False,
            tickfont_size=PLOT_FONT_SIZE_AXIS_TICKS, title_font_size=PLOT_FONT_SIZE_AXIS_TITLE, title_standoff=10
        ),
        height=plot_height,
        margin=dict(t=70, b=110 if show_xaxis_title else 90, l=80, r=170), 
        legend=dict(
            x=1.02, y=1, xanchor='left', yanchor='top', 
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor="#E0E0E0", borderwidth=1, 
            font_size=PLOT_FONT_SIZE_LEGEND_ITEMS,
            title_font_size=PLOT_FONT_SIZE_LEGEND_TITLE, 
            tracegroupgap=10
        )
    )
    return fig

def create_binned_stacked_variant_plot(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name, bin_size=BIN_SIZE_POSITIONAL):
    plot_title = f"{current_gene_name} ClinVar Variant Overview"
    data_source_annotation = dict(
        text="Data displayed here is from ClinVar's 2025年5月6日 release.",
        align='left', showarrow=False, xref='paper', yref='paper', x=0.01, 
        y= -0.35 if domains_data_local else -0.25, 
        font=dict(size=10, color="#555")
    )

    if df_positions_filtered.empty:
        fig = _setup_single_plot_with_domain_space(
            f"{plot_title} (No Data)", "ClinVar Variants", [0, 5], gene_length_local,
            plot_height=BINNED_PLOT_HEIGHT, show_xaxis_title=not domains_data_local
        )
        fig.update_xaxes(showticklabels=False) 
        if domains_data_local:
             domain_shapes_empty = create_domain_shapes_paper(domains_data_local, gene_length_local)
             fig.update_layout(shapes=domain_shapes_empty)
             fig.add_annotation(text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
                                x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1), 
                                xanchor='center', yanchor='top',
                                font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY))
        fig.add_annotation(**data_source_annotation)
        return fig

    df_copy = df_positions_filtered.copy()
    df_copy['Position'] = pd.to_numeric(df_copy['Position'], errors='coerce')
    df_copy.dropna(subset=['Position'], inplace=True)
    
    if df_copy.empty: 
        fig = _setup_single_plot_with_domain_space(
            f"{plot_title} (No Valid Positions)", "ClinVar Variants", [0, 5], gene_length_local,
            plot_height=BINNED_PLOT_HEIGHT, show_xaxis_title=not domains_data_local
        )
        fig.update_xaxes(showticklabels=False)
        if domains_data_local:
             domain_shapes_empty = create_domain_shapes_paper(domains_data_local, gene_length_local)
             fig.update_layout(shapes=domain_shapes_empty)
             fig.add_annotation(text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
                                x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1),
                                xanchor='center', yanchor='top',
                                font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY))
        fig.add_annotation(**data_source_annotation)
        return fig

    df_copy['Position'] = df_copy['Position'].astype(int)
    df_copy = df_copy[df_copy['Position'] >= 0]

    if df_copy.empty: 
        fig = _setup_single_plot_with_domain_space(
            f"{plot_title} (No Valid Non-Negative Positions)", "ClinVar Variants", [0, 5], gene_length_local,
            plot_height=BINNED_PLOT_HEIGHT, show_xaxis_title=not domains_data_local
        )
        fig.update_xaxes(showticklabels=False)
        if domains_data_local:
             domain_shapes_empty = create_domain_shapes_paper(domains_data_local, gene_length_local)
             fig.update_layout(shapes=domain_shapes_empty)
             fig.add_annotation(text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
                                x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1),
                                xanchor='center', yanchor='top',
                                font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY))
        fig.add_annotation(**data_source_annotation)
        return fig

    max_pos_data = df_copy['Position'].max() if not df_copy.empty else 0
    upper_bin_limit = max(gene_length_local, max_pos_data, bin_size if bin_size > 0 else 1)
    
    if upper_bin_limit == 0 and bin_size == 0: 
        bins = np.array([0, 1]) 
    elif bin_size == 0: 
         bins = np.arange(0, upper_bin_limit + 1, 1) 
         bin_size = 1 
    else:
        bins = np.arange(0, upper_bin_limit + bin_size, bin_size)

    if len(bins) < 2:
        if upper_bin_limit > 0:
            bins = np.array([0, float(max(upper_bin_limit, bin_size, 1))])
        else: 
            bins = np.array([0, float(max(bin_size, 1))])
    
    bin_labels_str = [f"{int(bins[i])}-{int(bins[i+1]-1)}" for i in range(len(bins)-1)]
    
    bin_starts = bins[:-1]
    bin_ends = bins[1:]
    bin_centers = (bin_starts + bin_ends) / 2.0

    if not bin_labels_str: 
        fig = _setup_single_plot_with_domain_space(
            f"{plot_title} (Cannot Create Bins)", "ClinVar Variants", [0, 5], gene_length_local,
            plot_height=BINNED_PLOT_HEIGHT, show_xaxis_title=not domains_data_local
        )
        fig.update_xaxes(showticklabels=False)
        if domains_data_local:
             domain_shapes_empty = create_domain_shapes_paper(domains_data_local, gene_length_local)
             fig.update_layout(shapes=domain_shapes_empty)
             fig.add_annotation(text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
                                x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1),
                                xanchor='center', yanchor='top',
                                font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY))
        fig.add_annotation(**data_source_annotation)
        return fig

    df_copy['Bin_Interval_Str'] = pd.cut(df_copy['Position'],
                                         bins=bins,
                                         labels=bin_labels_str,
                                         right=False, 
                                         include_lowest=True, 
                                         duplicates='drop')

    binned_counts = df_copy.groupby(['Bin_Interval_Str', 'SignificanceClass'], observed=False).size().unstack(fill_value=0)
    all_sig_classes_ordered = list(SIGNIFICANCE_PALETTE.keys())
    for sig_class in all_sig_classes_ordered:
        if sig_class not in binned_counts.columns:
            binned_counts[sig_class] = 0
    binned_counts = binned_counts[all_sig_classes_ordered]
    binned_counts_reindexed = binned_counts.reindex(bin_labels_str, fill_value=0).fillna(0)

    max_y_val_overall = binned_counts_reindexed.sum(axis=1).max() if not binned_counts_reindexed.empty else 0
    total_variants_in_plot = len(df_copy)

    fig = _setup_single_plot_with_domain_space(
        plot_title,
        f"ClinVar Variants (Total in View: {total_variants_in_plot})",
        [0, max_y_val_overall * 1.15 if max_y_val_overall > 0 else 5],
        gene_length_local,
        plot_height=BINNED_PLOT_HEIGHT,
        show_xaxis_title=not domains_data_local 
    )
    fig.update_layout(barmode='stack') 

    variant_types_per_bin_series = df_copy.groupby('Bin_Interval_Str', observed=False)['StandardType'].apply(lambda x: Counter(x) if not x.empty else Counter())
    variant_types_per_bin = {label: variant_types_per_bin_series.get(label, Counter()) for label in bin_labels_str}

    x_numerical_coords_for_bars = bin_centers 

    for sig_class in binned_counts_reindexed.columns:
        if binned_counts_reindexed[sig_class].sum() == 0:
            continue
        
        hover_texts = []
        for i, bin_label_for_hover in enumerate(bin_labels_str):
            current_bin_all_sig_counts_at_label = binned_counts_reindexed.loc[bin_label_for_hover]
            hover_content = [f"<b>Bin: {bin_label_for_hover} aa</b>"]
            hover_content.append("<u>Significance Counts:</u>")
            total_in_this_exact_bin = 0
            for s_class_iter, count_val_iter in current_bin_all_sig_counts_at_label.items():
                if count_val_iter > 0:
                    hover_content.append(f"  {s_class_iter}: {count_val_iter}")
                total_in_this_exact_bin += count_val_iter
            hover_content.append(f"<b>Total in Bin: {total_in_this_exact_bin}</b>")

            type_counts_in_this_bin = variant_types_per_bin.get(bin_label_for_hover, Counter())
            if isinstance(type_counts_in_this_bin, Counter) and type_counts_in_this_bin:
                hover_content.append("<u>Variant Type Counts (in bin):</u>")
                for v_type, v_count in type_counts_in_this_bin.items():
                    if v_count > 0: hover_content.append(f"  {v_type}: {v_count}")
            hover_texts.append("<br>".join(hover_content))

        fig.add_trace(go.Bar(
            x=x_numerical_coords_for_bars, 
            y=binned_counts_reindexed[sig_class].values,
            name=sig_class,
            marker_color=SIGNIFICANCE_PALETTE.get(sig_class, "#CCCCCC"),
            customdata=hover_texts, 
            hovertemplate="%{customdata}<extra></extra>",
            legendgroup="significance",
            legendgrouptitle_text="Variant Significance",
            width=bin_size * 0.98 
        ))

    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local)
    fig.update_layout(shapes=domain_shapes)
    
    num_bins_total_for_ticks = len(bin_labels_str)
    tick_values_to_use = list(x_numerical_coords_for_bars) 
    tick_texts_to_use = list(bin_labels_str) 

    fig.update_xaxes(showticklabels=False)

    if num_bins_total_for_ticks > 0:
        max_bins_show_all_labels = 10
        max_bins_show_subset_labels = 20
        target_num_ticks_for_subset = 5

        final_tickvals = []
        final_ticktext = []
        tick_font_size = PLOT_FONT_SIZE_AXIS_TICKS
        tick_angle = 0
        show_labels_flag = False

        if num_bins_total_for_ticks <= max_bins_show_all_labels:
            final_tickvals = tick_values_to_use
            final_ticktext = tick_texts_to_use
            tick_font_size = PLOT_FONT_SIZE_AXIS_TICKS -1
            tick_angle = 30 if num_bins_total_for_ticks > 5 else 0
            show_labels_flag = True
        elif num_bins_total_for_ticks <= max_bins_show_subset_labels:
            num_actual_ticks = min(target_num_ticks_for_subset, num_bins_total_for_ticks)
            if num_actual_ticks < 2 and num_bins_total_for_ticks >= 2: num_actual_ticks = 2
            elif num_bins_total_for_ticks == 1: num_actual_ticks = 1
            
            indices = []
            if num_actual_ticks == 1 and num_bins_total_for_ticks == 1:
                indices = [0]
            elif num_bins_total_for_ticks > 0:
                indices = np.round(np.linspace(0, num_bins_total_for_ticks - 1, num_actual_ticks)).astype(int)
                indices = sorted(list(set(indices)))
            
            final_tickvals = [tick_values_to_use[i] for i in indices if i < len(tick_values_to_use)]
            final_ticktext = [tick_texts_to_use[i] for i in indices if i < len(tick_texts_to_use)]
            
            if final_tickvals:
                tick_font_size = PLOT_FONT_SIZE_AXIS_TICKS - 2
                tick_angle = 45
                show_labels_flag = True
        else: 
            show_labels_flag = False 

        if show_labels_flag and final_tickvals:
            fig.update_xaxes(
                showticklabels=True,
                tickmode='array',
                tickvals=final_tickvals,
                ticktext=final_ticktext,
                tickfont={'size': tick_font_size},
                tickangle=tick_angle
            )
    if domains_data_local:
        fig.add_annotation(
            text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
            x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1), 
            xanchor='center', yanchor='top',
            font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY)
        )
    
    fig.add_annotation(**data_source_annotation)
    return fig

def create_interactive_lollipop(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name):
    if df_positions_filtered.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Lollipop (No Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
    
    significance_hierarchy = list(SIGNIFICANCE_MARKERS.keys())
    df_copy = df_positions_filtered.copy()
    df_copy['SignificanceOrdered'] = pd.Categorical(df_copy['SignificanceClass'], categories=significance_hierarchy, ordered=True)
    
    lollipop_agg_data = df_copy.groupby('Position', as_index=False).agg(
        TotalCount=('VariationID', 'nunique'),
        DominantSignificance=('SignificanceOrdered', 'min'),
        VariantTypes=('StandardType', lambda x: ', '.join(sorted(list(x.unique())))),
        VariantNames=('Name', lambda x: '<br>'.join(sorted(list(x.unique()))[:5]) + (f'<br>...and {len(x.unique())-5} more' if len(x.unique()) > 5 else '')),
        Phenotypes=('PhenotypeListClean', lambda x: '<br>'.join(sorted(list(x.unique()))[:3]) + (f'<br>...and {len(x.unique())-3} more' if len(x.unique()) > 3 else ''))
    )
    lollipop_agg_data['DominantSignificance'] = lollipop_agg_data['DominantSignificance'].astype(str)

    if lollipop_agg_data.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Lollipop (No Aggregated Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])

    max_total_count = lollipop_agg_data['TotalCount'].max() if not lollipop_agg_data.empty else 1
    y_min = -max_total_count * 0.05 if max_total_count > 1 else -0.1
    y_max = max_total_count * 1.15 if max_total_count > 0 else 1.15
    
    fig = _setup_single_plot_with_domain_space(
        f"{current_gene_name} Variant Lollipop Plot", 
        "Unique Variants per Position", 
        [y_min, y_max], 
        gene_length_local,
        show_xaxis_title=not domains_data_local 
    )
    fig.update_layout(yaxis_zeroline=True, yaxis_zerolinecolor='rgba(200,200,200,0.5)')

    for _, row in lollipop_agg_data.iterrows():
        fig.add_trace(go.Scatter(
            x=[row['Position'], row['Position']], 
            y=[0, row['TotalCount']], 
            mode='lines', 
            line=dict(color='#BCCCDC', width=1.5), 
            hoverinfo='none', 
            showlegend=False
        ))

    min_px_size, max_px_size = 6, 20 
    scaled_plotly_sizes = pd.Series(min_px_size, index=lollipop_agg_data.index)
    if max_total_count > 0:
        log_counts = np.log1p(lollipop_agg_data['TotalCount'])
        log_max_count = np.log1p(max_total_count)
        if log_max_count > 0 : 
            scaled_plotly_sizes = min_px_size + (max_px_size - min_px_size) * (log_counts / log_max_count)
        else: 
             scaled_plotly_sizes = pd.Series(min_px_size if max_total_count == 0 else (min_px_size + max_px_size)/2 , index=lollipop_agg_data.index)

    for sig_class in significance_hierarchy:
        subset = lollipop_agg_data[lollipop_agg_data['DominantSignificance'] == sig_class]
        if not subset.empty:
            fig.add_trace(go.Scatter(
                x=subset['Position'], y=subset['TotalCount'], mode='markers', 
                marker=dict(
                    color=SIGNIFICANCE_PALETTE.get(sig_class, '#CCCCCC'), 
                    size=scaled_plotly_sizes.loc[subset.index], 
                    symbol=SIGNIFICANCE_MARKERS.get(sig_class, 'circle'), 
                    line=dict(width=1, color='#333333') 
                ), 
                name=sig_class, legendgroup="significance", legendgrouptitle_text="Variant Significance", 
                customdata=subset[['Position', 'TotalCount', 'DominantSignificance', 'VariantTypes', 'VariantNames', 'Phenotypes']], 
                hovertemplate=(
                    "<b>Position:</b> %{customdata[0]}<br>"
                    "<b>Variant Count:</b> %{customdata[1]}<br>"
                    "<b>Dominant Significance:</b> %{customdata[2]}<br>"
                    "<b>Variant Types:</b> %{customdata[3]}<br>"
                    "<b>Names:</b><br>%{customdata[4]}<br>"
                    "<b>Phenotypes:</b><br>%{customdata[5]}"
                    "<extra></extra>"
                )
            ))
            
    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local)
    fig.update_layout(shapes=domain_shapes)
    
    if domains_data_local: 
        add_domain_legend_traces(fig, domains_data_local)
        fig.add_annotation(
            text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
            x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1), 
            xanchor='center', yanchor='top',
            font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY)
        )
    return fig

def create_interactive_waterfall(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name):
    if df_positions_filtered.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Waterfall (No Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
    agg_for_stacking = df_positions_filtered.groupby(["Position", "SignificanceClass"]).size().unstack(fill_value=0)
    sig_order_plot = [s for s in reversed(list(SIGNIFICANCE_PALETTE.keys())) if s in agg_for_stacking.columns] + [s for s in agg_for_stacking.columns if s not in SIGNIFICANCE_PALETTE.keys()]
    agg_for_stacking = agg_for_stacking.reindex(columns=sig_order_plot, fill_value=0); max_y_val_overall = agg_for_stacking.sum(axis=1).max() if not agg_for_stacking.empty else 0
    
    fig = _setup_single_plot_with_domain_space(
        f"{current_gene_name} Variant Waterfall Plot", 
        "Cumulative Variant Count", 
        [0, max_y_val_overall * 1.1 if max_y_val_overall > 0 else 1], 
        gene_length_local,
        show_xaxis_title=not domains_data_local 
    )
    fig.update_layout(barmode='stack')
    bottom_val = pd.Series(0, index=agg_for_stacking.index, dtype=float)
    for significance_class in sig_order_plot:
        if significance_class not in agg_for_stacking.columns: continue
        counts_this_sig = agg_for_stacking[significance_class]; active_positions = counts_this_sig[counts_this_sig > 0].index; active_counts_for_sig = counts_this_sig[counts_this_sig > 0]
        if not active_counts_for_sig.empty:
            current_bottoms = bottom_val.loc[active_positions]; hover_texts_for_sig = []
            for pos_idx, pos in enumerate(active_positions):
                variants_at_pos_sig = df_positions_filtered[(df_positions_filtered['Position'] == pos) & (df_positions_filtered['SignificanceClass'] == significance_class)]; variant_info_list = [f"- {row['Name']} ({row['StandardType']}): {row['PhenotypeListClean'][:80]}{'...' if len(row['PhenotypeListClean']) > 80 else ''}" for _, row in variants_at_pos_sig.head(5).iterrows()]; variant_info_str = "<br>".join(variant_info_list)
                if len(variants_at_pos_sig) > 5: variant_info_str += f"<br>...and {len(variants_at_pos_sig)-5} more."
                hover_texts_for_sig.append(f"<b>Position:</b> {pos}<br><b>Significance:</b> {significance_class}<br><b>Count (this segment):</b> {active_counts_for_sig.iloc[pos_idx]}<br><b>Total at position:</b> {agg_for_stacking.loc[pos].sum()}<br><b>Variants:</b><br>{variant_info_str}<extra></extra>")
            fig.add_trace(go.Bar(x=active_positions, y=active_counts_for_sig, base=current_bottoms,name=significance_class, marker_color=SIGNIFICANCE_PALETTE.get(significance_class, "#CCCCCC"),width=max(1.0, gene_length_local / 500), opacity=0.9, customdata=hover_texts_for_sig, hovertemplate="%{customdata}",legendgroup="significance", legendgrouptitle_text="Variant Significance"))
        bottom_val = bottom_val.add(counts_this_sig, fill_value=0)
    
    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local); 
    fig.update_layout(shapes=domain_shapes)
    if domains_data_local: 
        add_domain_legend_traces(fig, domains_data_local)
        fig.add_annotation(
            text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
            x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1),
            xanchor='center', yanchor='top',
            font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY)
        )
    return fig

def create_interactive_density(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name):
    data_source_annotation_density = dict(
        text="Data displayed here is from ClinVar's 2025年5月6日 release.",
        align='left', showarrow=False, xref='paper', yref='paper', x=0.01, 
        y= -0.18 if domains_data_local else -0.15, 
        font=dict(size=10, color="#555")
    )
    if df_positions_filtered.empty:
        fig = go.Figure().update_layout(title_text=f"{current_gene_name} Individual Variant Plot (No Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
        fig.add_annotation(**data_source_annotation_density)
        return fig

    fig = _setup_single_plot_with_domain_space(
        f"{current_gene_name} Individual Variant Plot", 
        "", [-0.5, 0.5], 
        gene_length_local, 
        plot_height=FIG_HEIGHT,
        show_xaxis_title=not domains_data_local 
    )
    fig.update_layout(yaxis_showgrid=False, yaxis_showticklabels=False, yaxis_zeroline=False)
    df_copy = df_positions_filtered.copy()
    df_copy['y_jitter'] = np.random.uniform(-0.4, 0.4, size=len(df_copy))
    significance_plot_order = list(SIGNIFICANCE_PALETTE.keys())
    for sig_class in significance_plot_order:
        group_df = df_copy[df_copy['SignificanceClass'] == sig_class]
        if group_df.empty: continue
        fig.add_trace(go.Scatter(x=group_df['Position'], y=group_df['y_jitter'], mode='markers',marker=dict(color=SIGNIFICANCE_PALETTE.get(sig_class, '#CCCCCC'), symbol=SIGNIFICANCE_MARKERS.get(sig_class, 'circle'), size=7, opacity=0.65, line=dict(width=0.5, color='#343A40')),name=sig_class, legendgroup="significance", legendgrouptitle_text="Variant Significance",customdata=group_df[['Position', 'Name', 'StandardType', 'SignificanceClass', 'PhenotypeListClean']],hovertemplate=("<b>Position:</b> %{customdata[0]}<br><b>Name:</b> %{customdata[1]}<br><b>Type:</b> %{customdata[2]}<br><b>Significance:</b> %{customdata[3]}<br><b>Phenotype:</b> %{customdata[4]}<extra></extra>")))
    
    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local)
    fig.update_layout(shapes=domain_shapes)
    if domains_data_local:
        add_domain_legend_traces(fig, domains_data_local)
        fig.add_annotation(
            text="Protein Position (aa)", showarrow=False, xref="paper", yref="paper",
            x=0.5, y=(DOMAIN_TRACK_PAPER_HEIGHT * 0.1), 
            xanchor='center', yanchor='top',
            font=dict(size=PLOT_FONT_SIZE_AXIS_TITLE-1, family=PLOT_FONT_FAMILY)
        )
    fig.add_annotation(**data_source_annotation_density)
    return fig

def create_origin_pie_chart(df_gene_variants_filtered, gene_name):
    if df_gene_variants_filtered.empty or 'OriginSimple' not in df_gene_variants_filtered.columns: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Variant Origin (No Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    origin_counts = df_gene_variants_filtered['OriginSimple'].value_counts()
    if origin_counts.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Variant Origin (No Origin Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    fig = go.Figure(data=[go.Pie(labels=origin_counts.index, values=origin_counts.values, hole=.3,marker_colors=px.colors.qualitative.Pastel2)])
    fig.update_traces(textposition='inside', textinfo='percent+label')
    fig.update_layout(title_text=f'{gene_name} Variant Origin Distribution',height=GENE_EXPLORER_OTHER_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,legend_title_text='Origin', legend=dict(x=1.02, y=1, xanchor='left', yanchor='top'), margin=dict(r=170, t=50))
    return fig

def create_type_significance_stacked_bar(df_gene_variants_filtered, gene_name):
    if df_gene_variants_filtered.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Type vs Significance (No Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    counts = df_gene_variants_filtered.groupby(['StandardType', 'SignificanceClass']).size().reset_index(name='Count')
    if counts.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Type vs Significance (No Counts)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    all_types = sorted(list(VARIANT_TYPE_PALETTE.keys())); all_significances = list(SIGNIFICANCE_PALETTE.keys())
    pivot_df = counts.pivot_table(index='StandardType', columns='SignificanceClass', values='Count', fill_value=0)
    pivot_df = pivot_df.reindex(index=all_types, columns=all_significances, fill_value=0).dropna(how='all', axis=0)
    fig = go.Figure()
    for sig_class in all_significances:
        if sig_class in pivot_df.columns and pivot_df[sig_class].sum() > 0: fig.add_trace(go.Bar(name=sig_class,x=pivot_df.index,y=pivot_df[sig_class],marker_color=SIGNIFICANCE_PALETTE.get(sig_class)))
    fig.update_layout(barmode='stack',title_text=f'{gene_name}: Clinical Significance by Variant Type',xaxis_title='Variant Type',yaxis_title='Number of Variants',legend_title_text='Clinical Significance',height=GENE_EXPLORER_OTHER_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,xaxis={'categoryorder':'array', 'categoryarray': [t for t in all_types if t in pivot_df.index]}, legend=dict(x=1.02, y=1, xanchor='left', yanchor='top'), margin=dict(r=170, t=50))
    return fig

def create_phenotype_type_stacked_bar(df_gene_variants_filtered, gene_name, normalize=False, top_n=TOP_N_PHENOTYPES_GENE_EXPLORER):
    if df_gene_variants_filtered.empty or 'PhenotypeListClean' not in df_gene_variants_filtered.columns: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Phenotype vs Type (No Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    df_pheno_exploded = df_gene_variants_filtered.copy(); df_pheno_exploded['PhenotypeSingle'] = df_pheno_exploded['PhenotypeListClean'].str.split('; ')
    df_pheno_exploded = df_pheno_exploded.explode('PhenotypeSingle').dropna(subset=['PhenotypeSingle']); df_pheno_exploded = df_pheno_exploded[df_pheno_exploded['PhenotypeSingle'] != 'N/A']
    if df_pheno_exploded.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Phenotype vs Type (No Phenotype Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    top_phenos_list = df_pheno_exploded['PhenotypeSingle'].value_counts().nlargest(top_n).index.tolist()
    df_pheno_exploded['DisplayPhenotype'] = df_pheno_exploded['PhenotypeSingle'].apply(lambda x: x if x in top_phenos_list else "Other Phenotypes")
    phenotype_order = top_phenos_list[:];
    if "Other Phenotypes" in df_pheno_exploded['DisplayPhenotype'].unique() and "Other Phenotypes" not in phenotype_order: phenotype_order.append("Other Phenotypes")
    counts = df_pheno_exploded.groupby(['DisplayPhenotype', 'StandardType']).size().reset_index(name='Count')
    if counts.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} Phenotype vs Type (No Counts)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
        return fig
    pivot_df = counts.pivot_table(index='DisplayPhenotype', columns='StandardType', values='Count', fill_value=0)
    if normalize: pivot_df = pivot_df.apply(lambda x: (x / x.sum() * 100) if x.sum() > 0 else x, axis=1)
    all_types_ordered = list(VARIANT_TYPE_PALETTE.keys())
    pivot_df = pivot_df.reindex(index=phenotype_order, columns=all_types_ordered, fill_value=0).dropna(how='all', axis=0)
    fig = go.Figure()
    for var_type in all_types_ordered:
        if var_type in pivot_df.columns and pivot_df[var_type].sum() > 0: fig.add_trace(go.Bar(name=var_type,x=pivot_df.index,y=pivot_df[var_type],marker_color=VARIANT_TYPE_PALETTE.get(var_type)))
    yaxis_title = "Percentage of Variants (%)" if normalize else "Number of Variant Associations"; plot_title = f'{gene_name}: Relative Variant Types by Phenotype (Top {top_n})' if normalize else f'{gene_name}: Variant Types by Phenotype (Top {top_n})'
    fig.update_layout(barmode='stack',title_text=plot_title,xaxis_title='Phenotype',yaxis_title=yaxis_title,legend_title_text='Variant Type',height=GENE_EXPLORER_OTHER_PLOT_HEIGHT + 100,font_family=PLOT_FONT_FAMILY,xaxis={'categoryorder':'array', 'categoryarray': [p for p in phenotype_order if p in pivot_df.index]},xaxis_tickangle=-30, legend=dict(x=1.02, y=1, xanchor='left', yanchor='top'), margin=dict(r=170, b=150, t=50)) 
    return fig

# HPO Plotting Functions
def plot_hpo_phenotype_categories_stacked_bar(df_hpo_phenotypes, gene_name):
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} HPO Phenotype Categories (No Data)", height=HPO_PLOT_HEIGHT)
        return fig
    df_copy = df_hpo_phenotypes.copy()
    df_copy['Category'] = df_copy['name'].apply(classify_hpo_phenotype_category)
    category_counts = df_copy.groupby('Category').size().reset_index(name='Count')
    category_order = category_counts.sort_values('Count', ascending=False)['Category'].tolist()
    fig = px.bar(category_counts, x='Category', y='Count', color='Category', title=f"HPO Phenotype Categories for {gene_name}",labels={'Count': 'Number of HPO Phenotypes', 'Category': 'Phenotype Category'},height=HPO_PLOT_HEIGHT,category_orders={"Category": category_order}, color_discrete_sequence=px.colors.qualitative.Pastel)
    fig.update_layout(font_family=PLOT_FONT_FAMILY, showlegend=False, xaxis_tickangle=-30, margin=dict(b=120, t=50)) 
    return fig

def plot_hpo_phenotype_frequency_bar(df_hpo_phenotypes, gene_name, top_n=15):
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} HPO Phenotype Frequency (No Data)", height=HPO_PLOT_HEIGHT + 100)
        return fig
    phenotype_counts_series = df_hpo_phenotypes['name'].value_counts()
    if phenotype_counts_series.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} HPO Phenotype Frequency (No Phenotypes Counted)", height=HPO_PLOT_HEIGHT + 100)
        return fig
    
    phenotype_counts_df = phenotype_counts_series.nlargest(top_n).reset_index()
    phenotype_counts_df.columns = ['Phenotype', 'Frequency']
    
    fig = px.bar(phenotype_counts_df.sort_values('Frequency', ascending=True), 
                 x='Frequency', y='Phenotype', orientation='h',
                 title=f"Top {top_n} Most Frequent HPO Phenotypes for {gene_name}",
                 labels={'Frequency': 'Frequency', 'Phenotype': 'HPO Phenotype Term'},
                 height=HPO_PLOT_HEIGHT + 50 + (min(top_n, len(phenotype_counts_df)) * 20), 
                 color_discrete_sequence=['#89CFF0']) 
    fig.update_layout(font_family=PLOT_FONT_FAMILY, yaxis={'categoryorder':'total ascending', 'automargin':True}, xaxis={'automargin':True}, plot_bgcolor='white', margin=dict(l=100, r=30, t=50)) 
    return fig

def create_hpo_top_phenotypes_table_df(df_hpo_phenotypes, top_n=15):
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns or 'id' not in df_hpo_phenotypes.columns:
        return pd.DataFrame(columns=['Phenotype Name', 'HPO ID', 'Frequency'])
    
    phenotype_counts = df_hpo_phenotypes['name'].value_counts().nlargest(top_n)
    if phenotype_counts.empty:
        return pd.DataFrame(columns=['Phenotype Name', 'HPO ID', 'Frequency'])
        
    top_phenotypes_list = []
    unique_pheno_names_with_ids = df_hpo_phenotypes[['name', 'id']].drop_duplicates(subset=['name'])
    
    for pheno_name, freq in phenotype_counts.items():
        hpo_id_series = unique_pheno_names_with_ids[unique_pheno_names_with_ids['name'] == pheno_name]['id']
        hpo_id = hpo_id_series.iloc[0] if not hpo_id_series.empty else "N/A"
        top_phenotypes_list.append({'Phenotype Name': pheno_name, 'HPO ID': hpo_id, 'Frequency': freq})
    
    return pd.DataFrame(top_phenotypes_list)

def create_hpo_wordcloud_image(df_hpo_phenotypes, gene_name):
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns: return None
    phenotype_names = df_hpo_phenotypes["name"].tolist()
    if not phenotype_names: return None
    text = " ".join(phenotype_names)
    try:
        wordcloud = WordCloud(width=800, height=400, background_color="white", colormap="viridis", font_path=None).generate(text) 
        img_path = f"{gene_name}_hpo_wordcloud.png" 
        wordcloud.to_file(img_path)
        return img_path
    except Exception as e: st.error(f"Error generating wordcloud for {gene_name}: {e}"); return None

def plot_hpo_gene_phenotype_network_plotly(df_hpo_phenotypes, df_hpo_diseases, gene_name, top_n_pheno=10, top_n_dis=5):
    if df_hpo_phenotypes.empty and df_hpo_diseases.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} HPO Network (No Data)", height=COMPLEX_PLOT_HEIGHT)
        return fig
        
    G = nx.Graph(); G.add_node(gene_name, type="gene", size=30, color='#FF6B6B', label=gene_name) 
    
    pheno_counts = df_hpo_phenotypes['name'].value_counts() if not df_hpo_phenotypes.empty else pd.Series()
    pheno_df_top = df_hpo_phenotypes[df_hpo_phenotypes['name'].isin(pheno_counts.nlargest(top_n_pheno).index)].drop_duplicates(subset=['id']) if not pheno_counts.empty else pd.DataFrame()

    for _, row in pheno_df_top.iterrows():
        pheno_label = row['name'][:30] + '...' if len(row['name']) > 30 else row['name']
        G.add_node(row['id'], type="hpo_phenotype", size=15, color='#4ECDC4', label=pheno_label); G.add_edge(gene_name, row['id']) 
    
    dis_df_top = df_hpo_diseases.drop_duplicates(subset=['id']).head(top_n_dis) if not df_hpo_diseases.empty else pd.DataFrame()

    for _, row in dis_df_top.iterrows():
        dis_label = row['name'][:30] + '...' if len(row['name']) > 30 else row['name']
        G.add_node(row['id'], type="hpo_disease", size=20, color='#FFE66D', label=dis_label); G.add_edge(gene_name, row['id']) 

    if not G.nodes() or len(G.nodes()) == 1: 
        fig = go.Figure()
        fig.update_layout(title_text=f"{gene_name} HPO Network (Not Enough Nodes)", height=COMPLEX_PLOT_HEIGHT)
        return fig
    
    pos = nx.spring_layout(G, k=0.9, iterations=70, seed=42); edge_x = []; edge_y = []
    for edge in G.edges(): x0, y0 = pos[edge[0]]; x1, y1 = pos[edge[1]]; edge_x.extend([x0, x1, None]); edge_y.extend([y0, y1, None])
    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.7, color='#B0B0B0'), hoverinfo='none', mode='lines')
    
    node_x = []; node_y = []; node_text_hover = []; node_text_display = []; node_size = []; node_color = []
    for node in G.nodes():
        x, y = pos[node]; node_x.append(x); node_y.append(y)
        node_type_display = G.nodes[node]['type'].replace('hpo_phenotype', 'Phenotype').replace('hpo_disease', 'Disease').capitalize()
        node_text_hover.append(f"<b>{G.nodes[node]['label']}</b><br>Type: {node_type_display}<br>ID: {node}")
        
        is_gene = G.nodes[node]['type'] == 'gene'
        label_text = G.nodes[node]['label']
        if not is_gene and len(label_text) > 15 : label_text = label_text[:12] + '...'
        node_text_display.append(label_text)

        node_size.append(G.nodes[node]['size']); node_color.append(G.nodes[node]['color'])
    
    node_trace = go.Scatter(x=node_x, y=node_y,mode='markers+text', hoverinfo='text',text= node_text_display, hovertext=node_text_hover, marker=dict(showscale=False, size=node_size, sizemode='diameter', color=node_color, line_width=1.5, line_color='black'),textfont=dict(size=9, color='black', family=PLOT_FONT_FAMILY), textposition="bottom center")
    
    fig = go.Figure(data=[edge_trace, node_trace], 
                    layout=go.Layout(
                        title=dict(text=f'Network of {gene_name} with Top HPO Phenotypes & Diseases', font=dict(size=PLOT_FONT_SIZE_MAIN_TITLE-1, family=PLOT_FONT_FAMILY)),
                        showlegend=False, hovermode='closest',height=COMPLEX_PLOT_HEIGHT,
                        margin=dict(b=20, l=5, r=5, t=50),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        plot_bgcolor='white'
                        ))
    return fig

# ─── Comparison Plotting Functions (CHD Family Tab) ─ CORRECTED FOR EMPTY DATA ─────────────────────────────
def plot_significance_distribution_per_gene(all_variants_df):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Significance Distribution", height=COMPARISON_PLOT_HEIGHT)
        return fig
    counts = all_variants_df.groupby(['Gene', 'SignificanceClass']).size().reset_index(name='Count')
    all_sig_classes = list(SIGNIFICANCE_PALETTE.keys())
    pivot_df = counts.pivot(index='Gene', columns='SignificanceClass', values='Count').fillna(0)
    for sig_class in all_sig_classes:
        if sig_class not in pivot_df.columns: pivot_df[sig_class] = 0
    pivot_df = pivot_df[all_sig_classes].reindex(ALL_GENES).fillna(0)
    fig = go.Figure()
    for sig_class in all_sig_classes:
        if sig_class in pivot_df.columns and pivot_df[sig_class].sum() > 0:
            fig.add_trace(go.Bar(name=sig_class,x=pivot_df.index,y=pivot_df[sig_class],marker_color=SIGNIFICANCE_PALETTE.get(sig_class)))
    fig.update_layout(barmode='stack',title_text='Clinical Significance Distribution per CHD Gene',xaxis_title='Gene',yaxis_title='Number of Variants',legend_title_text='Clinical Significance',height=COMPARISON_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY, legend=dict(x=1.02, y=1, xanchor='left', yanchor='top'), margin=dict(r=170, t=50))
    return fig

def plot_pathogenic_variant_types_heatmap(all_variants_df):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Pathogenic Variant Types Heatmap", height=COMPARISON_PLOT_HEIGHT)
        return fig
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]
    if patho_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No Pathogenic/Likely Pathogenic variants for Heatmap", height=COMPARISON_PLOT_HEIGHT)
        return fig
    counts = patho_df.groupby(['Gene', 'StandardType']).size().reset_index(name='Count')
    pivot_df = counts.pivot(index='Gene', columns='StandardType', values='Count').fillna(0).reindex(index=ALL_GENES).fillna(0)
    sorted_variant_types = sorted(pivot_df.columns.tolist())
    pivot_df = pivot_df[sorted_variant_types]
    fig = go.Figure(data=go.Heatmap(z=pivot_df.values,x=pivot_df.columns,y=pivot_df.index,colorscale='Blues',text=pivot_df.values,texttemplate="%{text}",hoverongaps=False))
    fig.update_layout(title_text='Heatmap of Pathogenic/Likely Pathogenic Variant Types per CHD Gene',xaxis_title='Variant Type',yaxis_title='Gene',height=COMPARISON_PLOT_HEIGHT + 100,font_family=PLOT_FONT_FAMILY,xaxis_tickangle=-45, yaxis_automargin=True, xaxis_automargin=True, margin=dict(r=50, b=120, l=70, t=50))
    return fig

def plot_pathogenic_variant_types_per_gene_stacked(all_variants_df):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Pathogenic Variant Types Distribution", height=COMPARISON_PLOT_HEIGHT)
        return fig
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]
    if patho_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No Pathogenic/Likely Pathogenic variants for Type Distribution", height=COMPARISON_PLOT_HEIGHT)
        return fig
    counts = patho_df.groupby(['Gene', 'StandardType']).size().reset_index(name='Count')
    all_variant_types = sorted(list(VARIANT_TYPE_PALETTE.keys()))
    pivot_df = counts.pivot(index='Gene', columns='StandardType', values='Count').fillna(0)
    for vt in all_variant_types:
        if vt not in pivot_df.columns: pivot_df[vt] = 0
    pivot_df = pivot_df.reindex(columns=all_variant_types, fill_value=0).reindex(index=ALL_GENES).fillna(0)
    fig = go.Figure()
    for var_type in all_variant_types:
        if var_type in pivot_df.columns and pivot_df[var_type].sum() > 0: 
             fig.add_trace(go.Bar(name=var_type,x=pivot_df.index,y=pivot_df[var_type],marker_color=VARIANT_TYPE_PALETTE.get(var_type, "#CCCCCC")))
    fig.update_layout(barmode='stack',title_text='Distribution of Pathogenic/Likely Pathogenic Variant Types per CHD Gene',xaxis_title='Gene',yaxis_title='Number of Pathogenic Variants',legend_title_text='Variant Type',height=COMPARISON_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY, legend=dict(x=1.02, y=1, xanchor='left', yanchor='top'), margin=dict(r=170, t=50))
    return fig

def plot_overall_pathogenic_phenotype_distribution(all_variants_df, top_n=DEFAULT_TOP_N_PHENOTYPES):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Phenotype Distribution", height=COMPARISON_PLOT_HEIGHT)
        return fig
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]
    if patho_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No Pathogenic/Likely Pathogenic variants for Phenotype Pie Chart", height=COMPARISON_PLOT_HEIGHT)
        return fig
    phenotypes_series = patho_df['PhenotypeListClean'].str.split('; ').explode().str.strip()
    phenotypes_series = phenotypes_series[phenotypes_series != 'N/A'].dropna()
    if phenotypes_series.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No valid phenotypes found for Pathogenic variants", height=COMPARISON_PLOT_HEIGHT)
        return fig
    phenotype_counts = Counter(phenotypes_series); top_phenotypes_data = phenotype_counts.most_common(top_n)
    if not top_phenotypes_data: 
        fig = go.Figure()
        fig.update_layout(title_text=f"No phenotypes to display for Top {top_n}", height=COMPARISON_PLOT_HEIGHT)
        return fig
    labels = [item[0] for item in top_phenotypes_data]; values = [item[1] for item in top_phenotypes_data]
    if len(phenotype_counts) > top_n: other_count = sum(count for pheno, count in phenotype_counts.items() if pheno not in labels); labels.append(f'Other ({len(phenotype_counts) - top_n} types)'); values.append(other_count)
    
    colors = px.colors.qualitative.Alphabet[:len(labels)] if len(labels) > len(PHENOTYPE_COLORS) else PHENOTYPE_COLORS[:len(labels)]

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3, marker_colors=colors, sort=False)])
    fig.update_traces(textposition='inside', textinfo='percent+label', pull=[0.05 if i==0 and len(labels)>1 else 0 for i in range(len(labels))])
    fig.update_layout(title_text=f'Overall Distribution of Top {top_n} Pathogenic ClinVar Phenotypes (All CHD Genes)',height=COMPARISON_PLOT_HEIGHT + 150,font_family=PLOT_FONT_FAMILY,legend_title_text='Phenotype',showlegend=True,
                      legend=dict(orientation="v", yanchor="top", y=0.95, xanchor="left", x=1.02, traceorder="normal"),  
                      margin=dict(r=200, t=80)) 
    return fig

def plot_gene_phenotype_heatmap(all_variants_df, top_n_phenotypes=DEFAULT_TOP_N_PHENOTYPES):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Gene-Phenotype Heatmap", height=COMPLEX_PLOT_HEIGHT)
        return fig
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)].copy()
    if patho_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No Pathogenic variants for Gene-Phenotype Heatmap", height=COMPLEX_PLOT_HEIGHT)
        return fig
    patho_df['PhenotypeSingle'] = patho_df['PhenotypeListClean'].str.split('; ')
    exploded_pheno_df = patho_df.explode('PhenotypeSingle'); exploded_pheno_df = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'] != 'N/A'].dropna(subset=['PhenotypeSingle'])
    if exploded_pheno_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No valid phenotypes for Gene-Phenotype Heatmap", height=COMPLEX_PLOT_HEIGHT)
        return fig
    top_phenos_list = exploded_pheno_df['PhenotypeSingle'].value_counts().nlargest(top_n_phenotypes).index.tolist()
    if not top_phenos_list: 
        fig = go.Figure()
        fig.update_layout(title_text=f"No phenotypes found for Top {top_n_phenotypes}", height=COMPLEX_PLOT_HEIGHT)
        return fig
    filtered_exploded_df = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'].isin(top_phenos_list)]
    counts = filtered_exploded_df.groupby(['Gene', 'PhenotypeSingle']).size().reset_index(name='Count')
    pivot_df = counts.pivot(index='Gene', columns='PhenotypeSingle', values='Count').fillna(0); pivot_df = pivot_df.reindex(index=ALL_GENES, columns=top_phenos_list).fillna(0)
    fig = go.Figure(data=go.Heatmap(z=pivot_df.values,x=pivot_df.columns,y=pivot_df.index,colorscale='YlOrRd',text=pivot_df.values,texttemplate="%{text}",hoverongaps=False,hovertemplate="Gene: %{y}<br>Phenotype: %{x}<br>Pathogenic Variants: %{z}<extra></extra>"))
    fig.update_layout(title_text=f'Heatmap of Pathogenic Variants by Gene and Top {top_n_phenotypes} ClinVar Phenotypes',xaxis_title='Phenotype',yaxis_title='Gene',height=COMPLEX_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,xaxis_tickangle=-45, yaxis_automargin=True, xaxis_automargin=True, margin=dict(r=50, b=150, l=70, t=80)) 
    return fig

def plot_pathogenic_bubble_chart_gene_phenotype(all_variants_df, top_n_phenotypes=DEFAULT_TOP_N_PHENOTYPES):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
        return fig
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)].copy()
    if patho_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No Pathogenic variants for Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
        return fig
    patho_df['PhenotypeSingle'] = patho_df['PhenotypeListClean'].str.split('; ')
    exploded_pheno_df = patho_df.explode('PhenotypeSingle'); exploded_pheno_df = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'] != 'N/A'].dropna(subset=['PhenotypeSingle'])
    if exploded_pheno_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No valid phenotypes for Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
        return fig
    top_phenos_list = exploded_pheno_df['PhenotypeSingle'].value_counts().nlargest(top_n_phenotypes).index.tolist()
    if not top_phenos_list: 
        fig = go.Figure()
        fig.update_layout(title_text=f"No phenotypes found for Top {top_n_phenotypes}", height=COMPLEX_PLOT_HEIGHT)
        return fig
    bubble_data = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'].isin(top_phenos_list)]
    counts_df = bubble_data.groupby(['Gene', 'PhenotypeSingle']).size().reset_index(name='PathogenicVariantCount')
    if counts_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"No data for selected Top {top_n_phenotypes} phenotypes in Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
        return fig
    fig = px.scatter(counts_df,x='Gene',y='PhenotypeSingle',size='PathogenicVariantCount',color='Gene',hover_name='PhenotypeSingle',size_max=60,title=f'Bubble Chart of Pathogenic Variants by Gene and Top {top_n_phenotypes} ClinVar Phenotypes',labels={'Gene': 'Gene', 'PhenotypeSingle': 'Phenotype', 'PathogenicVariantCount': 'Number of Pathogenic Variants'},color_discrete_sequence=GENE_COLORS)
    fig.update_layout(height=COMPLEX_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,yaxis={'categoryorder':'total ascending', 'automargin':True}, xaxis_automargin=True, legend=dict(x=1.02, y=1, xanchor='left', yanchor='top'), margin=dict(r=170, b=100, l=100, t=80)) 
    fig.update_traces(hovertemplate="Gene: %{x}<br>Phenotype: %{y}<br>Count: %{marker.size}<extra></extra>")
    return fig

def plot_sankey_variant_type_to_phenotype(all_variants_df, top_n_phenotypes=DEFAULT_TOP_N_PHENOTYPES):
    if all_variants_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No data for Sankey Diagram", height=COMPLEX_PLOT_HEIGHT)
        return fig
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)].copy()
    if patho_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No Pathogenic variants for Sankey Diagram", height=COMPLEX_PLOT_HEIGHT)
        return fig
    patho_df['PhenotypeSingle'] = patho_df['PhenotypeListClean'].str.split('; ')
    exploded_df = patho_df.explode('PhenotypeSingle'); exploded_df = exploded_df[exploded_df['PhenotypeSingle'] != 'N/A'].dropna(subset=['PhenotypeSingle', 'StandardType'])
    if exploded_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text="No valid phenotype/type data for Sankey", height=COMPLEX_PLOT_HEIGHT)
        return fig
    top_phenos_list = exploded_df['PhenotypeSingle'].value_counts().nlargest(top_n_phenotypes).index.tolist()
    if not top_phenos_list: 
        fig = go.Figure()
        fig.update_layout(title_text=f"No phenotypes found for Top {top_n_phenotypes}", height=COMPLEX_PLOT_HEIGHT)
        return fig
    sankey_df = exploded_df[exploded_df['PhenotypeSingle'].isin(top_phenos_list)]
    if sankey_df.empty: 
        fig = go.Figure()
        fig.update_layout(title_text=f"No data for Top {top_n_phenotypes} phenotypes in Sankey", height=COMPLEX_PLOT_HEIGHT)
        return fig
    
    all_variant_types = sorted(list(VARIANT_TYPE_PALETTE.keys()))
    all_phenotypes_sankey = sorted(sankey_df['PhenotypeSingle'].unique())
    labels = all_variant_types + all_phenotypes_sankey
    label_indices = {label: i for i, label in enumerate(labels)}
    
    source = []
    target = []
    value = []
    link_colors = [] 
    
    link_counts = sankey_df.groupby(['StandardType', 'PhenotypeSingle']).size().reset_index(name='Count')
    type_palette = VARIANT_TYPE_PALETTE
    
    link_alpha = 0.6 

    for i, row in link_counts.iterrows():
        if row['StandardType'] in label_indices and row['PhenotypeSingle'] in label_indices:
            src_idx = label_indices[row['StandardType']]
            tgt_idx = label_indices[row['PhenotypeSingle']]
            source.append(src_idx)
            target.append(tgt_idx)
            value.append(row['Count'])
            
            base_color_hex = type_palette.get(row['StandardType'], "#CCCCCC")
            link_colors.append(hex_to_rgba(base_color_hex, alpha=link_alpha)) 
            
    node_colors_final = [type_palette.get(vt, "#CCCCCC") for vt in all_variant_types]
    pheno_sankey_colors = px.colors.qualitative.Pastel[:len(all_phenotypes_sankey)]
    if len(all_phenotypes_sankey) > len(px.colors.qualitative.Pastel): 
        pheno_sankey_colors = px.colors.qualitative.Pastel * (len(all_phenotypes_sankey) // len(px.colors.qualitative.Pastel) +1)
    node_colors_final.extend(pheno_sankey_colors[:len(all_phenotypes_sankey)])
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,thickness=20,line=dict(color="black", width=0.5),
            label=labels,color=node_colors_final,
            hovertemplate='Node: %{label}<extra></extra>'
            ),
        link=dict(
            source=source,target=target,value=value,
            color=link_colors, 
            hovertemplate='From %{source.label} to %{target.label}: %{value} variants<extra></extra>'
            ))])
    fig.update_layout(title=dict(text=f'Sankey Diagram: Pathogenic Variant Types to Top {top_n_phenotypes} ClinVar Phenotypes',font=dict(size=PLOT_FONT_SIZE_MAIN_TITLE-1, family=PLOT_FONT_FAMILY)),height=COMPLEX_PLOT_HEIGHT + 100, font_family=PLOT_FONT_FAMILY, margin=dict(t=80, b=50, l=50, r=50)) 
    return fig

# ────────────────────────── STREAMLIT APP LAYOUT & LOGIC ──────────────────────────
# --- Sidebar Setup ---
if os.path.exists(LOGO_PATH):
    try:
        logo_image_sidebar = Image.open(LOGO_PATH)
        st.sidebar.image(logo_image_sidebar, use_container_width=True) 
    except Exception as e:
        st.sidebar.warning(f"Could not load logo: {e}")
else:
    st.sidebar.caption("Logo not found.")

st.sidebar.header("Gene Selection")
selected_gene = st.sidebar.radio(
    "Select Gene:",
    ALL_GENES,
    index=ALL_GENES.index("CHD1") if "CHD1" in ALL_GENES else 0,
    key="gene_selector_main"
)

st.sidebar.header(f"Filters for {selected_gene} Positional Plots")
df_gene_variants, gene_domains, gene_length, clinvar_load_errors = load_gene_data(selected_gene)
if clinvar_load_errors:
    for error_msg in clinvar_load_errors:
        st.sidebar.warning(error_msg)
min_pos_slider = 0
max_pos_slider = gene_length
if max_pos_slider <= min_pos_slider:
    max_pos_slider = min_pos_slider + 10
default_slider_val = [min_pos_slider, max_pos_slider]
slider_session_key = f'pos_slider_val_{selected_gene}'
current_slider_val = st.session_state.get(slider_session_key, default_slider_val)
if not (isinstance(current_slider_val, (list, tuple)) and len(current_slider_val) == 2 and
        current_slider_val[0] <= current_slider_val[1] and
        current_slider_val[0] >= min_pos_slider and current_slider_val[1] <= max_pos_slider):
    current_slider_val = list(default_slider_val)
pos_range_val = st.sidebar.slider(
    f"Position Range ({selected_gene}):",
    min_pos_slider, max_pos_slider,
    tuple(current_slider_val), 1,
    key=f"pos_slider_widget_{selected_gene}"
)
st.session_state[slider_session_key] = list(pos_range_val)
sig_options = sorted(list(SIGNIFICANCE_PALETTE.keys()), key=lambda x: list(SIGNIFICANCE_PALETTE.keys()).index(x))
default_sig = st.session_state.get(f'sel_sig_{selected_gene}', sig_options)
unique_sig_in_data = df_gene_variants['SignificanceClass'].dropna().unique() if 'SignificanceClass' in df_gene_variants else []
actual_sig_options = [s for s in sig_options if s in unique_sig_in_data]
default_sig = [s for s in default_sig if s in actual_sig_options] if default_sig else actual_sig_options

selected_significance_filter = st.sidebar.multiselect(
    f"Significance (Positional Plots):",
    actual_sig_options, default_sig,
    key=f"sig_multi_widget_{selected_gene}"
)
st.session_state[f'sel_sig_{selected_gene}'] = selected_significance_filter

type_options = sorted(list(VARIANT_TYPE_PALETTE.keys()))
default_type = st.session_state.get(f'sel_type_{selected_gene}', type_options)
unique_types_in_data = df_gene_variants['StandardType'].dropna().unique() if 'StandardType' in df_gene_variants else []
actual_type_options = [t for t in type_options if t in unique_types_in_data]
default_type = [t for t in default_type if t in actual_type_options] if default_type else actual_type_options

selected_types_filter = st.sidebar.multiselect(
    f"Variant Type (Positional Plots):",
    actual_type_options, default_type,
    key=f"type_multi_widget_{selected_gene}"
)
st.session_state[f'sel_type_{selected_gene}'] = selected_types_filter
st.sidebar.markdown("---")
st.sidebar.markdown("""<div class="sidebar-footer">
<b>Developed by:</b><br>
Zi-Hao Wang<br>
Instutute of Biomedical Sciences(IBS), Fudan University
</div>""", unsafe_allow_html=True)
st.sidebar.markdown("""<div class="sidebar-footer" style="margin-top:10px;">
<b>Data Sources:</b><br>
<a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank">NCBI ClinVar</a><br>
<a href="https://www.uniprot.org/" target="_blank">UniProt</a><br>
<a href="https://hpo.jax.org/" target="_blank">HPO</a>
</div>""", unsafe_allow_html=True)
contact_email_sidebar = "soap@fastemail.io" 
current_date_sidebar = datetime.now().strftime("%Y.%m.%d")
version_info_sidebar = f"v3.8.9" 
st.sidebar.markdown(f"""<div class="sidebar-footer" style="margin-top:10px;">
<b>Contact:</b> <a href="mailto:{contact_email_sidebar}">{contact_email_sidebar}</a><br>
<b>Last Updated:</b> {current_date_sidebar}<br>
<b>Version:</b> {version_info_sidebar}<br>
<i>Continuously Evolving...</i>
</div>""", unsafe_allow_html=True)


st.title("CHD Gene Family Variant Database")

df_hpo_phenotypes, df_hpo_diseases, hpo_load_error = load_hpo_data(selected_gene)
if hpo_load_error and "file not found" not in hpo_load_error.lower(): 
    st.warning(hpo_load_error)

main_tabs_list = ["Gene Explorer (ClinVar)", "CHD Family Comparison (ClinVar)", "HPO Phenotype Explorer"]
main_tab1, main_tab2, main_tab3 = st.tabs(main_tabs_list)

with main_tab1:
    st.header(f"{selected_gene} ClinVar Variant Explorer")
    if df_gene_variants.empty and not gene_domains :
        st.warning(f"No ClinVar variant or domain data could be loaded for {selected_gene}.")
    else:
        df_positions_filtered = pd.DataFrame()
        if 'Position' in df_gene_variants.columns and not df_gene_variants.empty:
            df_positions_filtered_temp = df_gene_variants.dropna(subset=['Position']).copy()
            if not df_positions_filtered_temp.empty:
                df_positions_filtered_temp['Position'] = pd.to_numeric(df_positions_filtered_temp['Position'], errors='coerce')
                df_positions_filtered_temp = df_positions_filtered_temp.dropna(subset=['Position'])
                if not df_positions_filtered_temp.empty: 
                    df_positions_filtered_temp['Position'] = df_positions_filtered_temp['Position'].astype(int)
                    
                    df_positions_filtered = df_positions_filtered_temp[
                        (df_positions_filtered_temp['Position'] >= pos_range_val[0]) &
                        (df_positions_filtered_temp['Position'] <= pos_range_val[1])
                    ].copy() 

                    if selected_significance_filter and not df_positions_filtered.empty:
                        df_positions_filtered = df_positions_filtered[df_positions_filtered['SignificanceClass'].isin(selected_significance_filter)]
                    elif not selected_significance_filter: 
                         df_positions_filtered = pd.DataFrame(columns=df_positions_filtered.columns)


                    if selected_types_filter and not df_positions_filtered.empty:
                        df_positions_filtered = df_positions_filtered[df_positions_filtered['StandardType'].isin(selected_types_filter)]
                    elif not selected_types_filter: 
                         df_positions_filtered = pd.DataFrame(columns=df_positions_filtered.columns)
        
        num_variants_for_plot_tab_title = len(df_positions_filtered) if not df_positions_filtered.empty else 0
        
        ge_tab_pos, ge_tab_type_sig, ge_tab_pheno, ge_tab_origin = st.tabs([
            f"Positional Plots & Search ({num_variants_for_plot_tab_title:,} variants in plots)",
            "Type & Significance",
            "ClinVar Phenotypes",
            "Variant Origin"
        ])

        with ge_tab_pos:
            st.subheader(f"Search Variant by Protein Position in {selected_gene}")
            searchable_df = pd.DataFrame()
            if 'Position' in df_gene_variants.columns and not df_gene_variants.empty:
                 df_gene_variants_copy = df_gene_variants.copy()
                 df_gene_variants_copy['Position'] = pd.to_numeric(df_gene_variants_copy['Position'], errors='coerce')
                 searchable_df = df_gene_variants_copy.dropna(subset=['Position']).copy()
                 if not searchable_df.empty:
                    searchable_df['Position'] = searchable_df['Position'].astype(int)

            search_col1_pos, search_col2_pos = st.columns([3, 1])
            with search_col1_pos:
                search_pos = search_col1_pos.number_input(
                    "Enter exact protein position (aa):", min_value=1,
                    max_value=gene_length if gene_length > 0 else 10000,
                    step=1, value=None, placeholder="e.g., 123",
                    key=f"search_pos_input_{selected_gene}", label_visibility="collapsed"
                )
            with search_col2_pos:
                search_button = st.button(f"🔍 Search Position", key=f"search_pos_button_{selected_gene}", use_container_width=True)

            if search_button and search_pos is not None:
                st.markdown("---")
                st.subheader(f"Search Results for Position {search_pos}")
                results_df_search = pd.DataFrame()
                if not searchable_df.empty:
                    results_df_search = searchable_df[searchable_df['Position'] == search_pos]
                    if not results_df_search.empty:
                        st.success(f"Found {len(results_df_search)} variant(s) at position {search_pos}:")
                        for index, row in results_df_search.iterrows():
                            with st.expander(f"**{row.get('Name', 'N/A')}** (ClinVar ID: {row.get('VariationID', 'N/A')})"):
                                st.markdown(f"""
                                *   **Type:** {row.get('StandardType', 'N/A')}
                                *   **Clinical Significance:** {row.get('SignificanceClass', 'N/A')}
                                *   **Reported Origin:** {row.get('OriginSimple', 'N/A')}
                                *   **Phenotypes (ClinVar):** {row.get('PhenotypeListClean', 'N/A')}
                                """)
                                if pd.notna(row.get('VariationID')):
                                     st.link_button("View on ClinVar", f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{row.get('VariationID')}/")
                    else:
                        st.info(f"No ClinVar variants found exactly at protein position {search_pos} for {selected_gene}.")
                else:
                     st.warning(f"No positional data available in the loaded ClinVar file for {selected_gene} to perform search.")
                st.markdown("---")
            elif search_button and search_pos is None:
                 st.warning("Please enter a protein position to search.")

            st.subheader(f"Positional Variant Plots for {selected_gene}")
            st.markdown(f"_Displaying variants matching sidebar filters (Position Range, Significance, Type). Total matching: {num_variants_for_plot_tab_title:,}._")

            fig_binned = create_binned_stacked_variant_plot(df_positions_filtered, gene_domains, gene_length, selected_gene)
            st.plotly_chart(fig_binned, use_container_width=True)
            
            if gene_domains:
                domain_details_md = "<b>Protein Domain Details:</b><br>"
                for i, domain in enumerate(gene_domains):
                    color = DOMAIN_COLORS[i % len(DOMAIN_COLORS)]
                    domain_details_md += f"<span style='color:{color}; font-weight:bold;'>■ {domain.get('name', f'Domain {i+1}')}</span>: {domain.get('start','N/A')} - {domain.get('end','N/A')} aa<br>"
                st.caption(domain_details_md, unsafe_allow_html=True)
            
            st.markdown("---") 

            if not df_positions_filtered.empty:
                # Checkbox placed directly under the binned plot and domain details
                show_details = st.checkbox(
                    "📊 Show Detailed Variant Plots (Lollipop, Waterfall, Individual Points)", 
                    value=False, 
                    key=f"show_details_{selected_gene}_main", # Ensure key is unique if used elsewhere
                    help="Display additional detailed plots for the filtered variants. These plots may take a moment to render."
                )

                if show_details:
                    st.markdown("---") 
                    if search_button and search_pos is not None and 'results_df_search' in locals() and not results_df_search.empty: 
                        if pos_range_val[0] <= search_pos <= pos_range_val[1]:
                            st.markdown(f"_*Note: Variants found at searched position {search_pos} are highlighted or visible in the plots below if they match the selected filters._")
                        else:
                            st.markdown(f"_*Note: Searched position {search_pos} is outside the current plot range [{pos_range_val[0]}, {pos_range_val[1]}]. Adjust the 'Position Range' filter to include it._")
                    
                    st.plotly_chart(create_interactive_lollipop(df_positions_filtered, gene_domains, gene_length, selected_gene), use_container_width=True)
                    st.markdown("---")
                    st.plotly_chart(create_interactive_waterfall(df_positions_filtered, gene_domains, gene_length, selected_gene), use_container_width=True)
                    st.markdown("---")
                    st.plotly_chart(create_interactive_density(df_positions_filtered, gene_domains, gene_length, selected_gene), use_container_width=True)
            
            elif df_gene_variants.empty and ('Position' not in df_gene_variants or df_gene_variants['Position'].isnull().all()):
                 pass 
            else: 
                 if not df_positions_filtered.empty: 
                    st.info(f"No variants match current sidebar filters for {selected_gene} to display detailed plots. Adjust filters or note that all variants might be filtered out.")


        with ge_tab_type_sig:
            st.subheader(f"Variant Type vs. Clinical Significance for {selected_gene} (All ClinVar Variants for this Gene)")
            st.plotly_chart(create_type_significance_stacked_bar(df_gene_variants, selected_gene), use_container_width=True)
        with ge_tab_pheno:
            st.subheader(f"ClinVar Phenotype Analysis for {selected_gene} (Top {TOP_N_PHENOTYPES_GENE_EXPLORER} Phenotypes, All ClinVar Variants for this Gene)")
            st.plotly_chart(create_phenotype_type_stacked_bar(df_gene_variants, selected_gene, normalize=False, top_n=TOP_N_PHENOTYPES_GENE_EXPLORER), use_container_width=True)
            st.markdown("---")
            st.plotly_chart(create_phenotype_type_stacked_bar(df_gene_variants, selected_gene, normalize=True, top_n=TOP_N_PHENOTYPES_GENE_EXPLORER), use_container_width=True)
        with ge_tab_origin:
            st.subheader(f"Variant Origin Distribution for {selected_gene} (All ClinVar Variants for this Gene)")
            st.plotly_chart(create_origin_pie_chart(df_gene_variants, selected_gene), use_container_width=True)

        st.markdown("---");
        st.markdown(f"Gene Length Used for Plotting {selected_gene}: **{gene_length} aa**. Domains: **{', '.join([d['name'] for d in gene_domains]) if gene_domains else 'N/A'}**.")
        st.markdown(f"Data source: ClinVar (2025年5月6日 release).") 
        st.markdown("---"); st.subheader(f"Explore {selected_gene} on External Resources")
        link_defs = [{"name": "PubMed", "url": f"https://pubmed.ncbi.nlm.nih.gov/?term={selected_gene}", "icon": "🔬"},{"name": "OMIM", "url": f"https://www.omim.org/search?index=entry&search={selected_gene}", "icon": "🧬"},{"name": "DECIPHER", "url": f"https://www.deciphergenomics.org/gene/{selected_gene}/overview/protein-genomic", "icon": "💡"},{"name": "gnomAD", "url": f"https://gnomad.broadinstitute.org/search/{selected_gene}", "icon": "📊"},{"name": "Protein Atlas", "url": f"https://www.proteinatlas.org/search/{selected_gene}", "icon": "🖼️"}]
        cols = st.columns(len(link_defs));
        for i, link in enumerate(link_defs): cols[i].link_button(f"{link['icon']} {link['name']}", link['url'], use_container_width=True)

with main_tab2:
    st.header("CHD Gene Family Comparison (ClinVar Data)")
    all_chd_variants_df, comparison_summary_df, all_genes_load_errors = load_all_chd_data_for_comparison()
    st.subheader("Basic Summary Statistics")
    if "UniProt ID" in comparison_summary_df.columns and "RefSeq NM" in comparison_summary_df.columns:
        summary_df_for_display = comparison_summary_df.copy()
        summary_df_for_display["UniProt Link"] = summary_df_for_display["UniProt ID"].apply(lambda x: f"https://www.uniprot.org/uniprotkb/{x}/entry" if x != "N/A" else "")
        summary_df_for_display["RefSeq Link"] = summary_df_for_display["RefSeq NM"].apply(lambda x: f"https://www.ncbi.nlm.nih.gov/nuccore/{x}" if x != "N/A" else "")
        st.dataframe(summary_df_for_display, use_container_width=True, hide_index=True,column_config={"UniProt ID": st.column_config.TextColumn("UniProt ID"),"UniProt Link": st.column_config.LinkColumn("UniProt Link", display_text="🔗"),"RefSeq NM": st.column_config.TextColumn("RefSeq NM"),"RefSeq Link": st.column_config.LinkColumn("RefSeq Link", display_text="🔗")},column_order=["Gene", "UniProt ID", "UniProt Link", "RefSeq NM", "RefSeq Link", "Protein Length (aa)", "Total Variants", "Pathogenic/Likely Path. Variants", "Data Loading Status"])
    else: st.dataframe(comparison_summary_df, use_container_width=True, hide_index=True)
    if all_genes_load_errors:
        st.subheader("Data Loading Issues:");
        for gene, errors in all_genes_load_errors.items():
            with st.expander(f"Details for {gene}"):
                for e in errors: st.warning(f"- {e}")
    st.markdown("---"); st.subheader("Comparative Visualizations")
    if not all_chd_variants_df.empty:
        top_n_pheno_slider_val = st.slider(label="Select Top N ClinVar Phenotypes for complex plots:", min_value=5, max_value=30, value=DEFAULT_TOP_N_PHENOTYPES, step=1, key="top_n_clinvar_pheno_slider")
        st.markdown("##### Pathogenic/Likely Pathogenic Variants per Gene")
        fig_patho_counts = px.bar(comparison_summary_df.sort_values("Pathogenic/Likely Path. Variants", ascending=False), x="Gene", y="Pathogenic/Likely Path. Variants",title="Pathogenic/Likely Pathogenic Variants per CHD Gene",labels={"Pathogenic/Likely Path. Variants": "Number of Patho/LP Variants"},height=COMPARISON_PLOT_HEIGHT,color_discrete_sequence=["#007AFF"])
        fig_patho_counts.update_layout(font_family=PLOT_FONT_FAMILY, margin=dict(r=170, t=50), legend=dict(x=1.02, y=1, xanchor='left', yanchor='top')); st.plotly_chart(fig_patho_counts, use_container_width=True); st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("##### Clinical Significance Distribution")
        fig_sig_dist = plot_significance_distribution_per_gene(all_chd_variants_df); st.plotly_chart(fig_sig_dist, use_container_width=True); st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("##### Pathogenic Variant Types Heatmap (Gene vs Type)")
        fig_patho_types_heatmap = plot_pathogenic_variant_types_heatmap(all_chd_variants_df); st.plotly_chart(fig_patho_types_heatmap, use_container_width=True); st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(f"##### Overall Pathogenic ClinVar Phenotype Distribution (Top {top_n_pheno_slider_val})")
        fig_overall_pheno_pie = plot_overall_pathogenic_phenotype_distribution(all_chd_variants_df, top_n=top_n_pheno_slider_val); st.plotly_chart(fig_overall_pheno_pie, use_container_width=True); st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("##### Pathogenic Variant Type Distribution per Gene (Stacked)")
        fig_patho_types_stacked = plot_pathogenic_variant_types_per_gene_stacked(all_chd_variants_df); st.plotly_chart(fig_patho_types_stacked, use_container_width=True); st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(f"##### Heatmap of Pathogenic Variants (Gene vs Top {top_n_pheno_slider_val} ClinVar Phenotypes)")
        fig_gene_pheno_heatmap = plot_gene_phenotype_heatmap(all_chd_variants_df, top_n_phenotypes=top_n_pheno_slider_val)
        st.plotly_chart(fig_gene_pheno_heatmap, use_container_width=True)
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(f"##### Bubble Chart of Pathogenic Variants (Gene vs Top {top_n_pheno_slider_val} ClinVar Phenotypes)")
        fig_bubble_gene_pheno = plot_pathogenic_bubble_chart_gene_phenotype(all_chd_variants_df, top_n_phenotypes=top_n_pheno_slider_val)
        st.plotly_chart(fig_bubble_gene_pheno, use_container_width=True)
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown(f"##### Sankey: Pathogenic Variant Types to Top {top_n_pheno_slider_val} ClinVar Phenotypes")
        fig_sankey_type_pheno = plot_sankey_variant_type_to_phenotype(all_chd_variants_df, top_n_phenotypes=top_n_pheno_slider_val)
        st.plotly_chart(fig_sankey_type_pheno, use_container_width=True)
        st.markdown("<br>", unsafe_allow_html=True)
    else: st.warning("No variant data available across all CHD genes to generate comparative plots.")
    st.markdown("---"); st.subheader("Further Comparison Ideas:")
    st.markdown("""*   **Comparative Domain Architecture Plot.**\n*   **Phenotype Similarity/Clustering.**\n*   **Variant Hotspot Comparison.**\n*   **Gene Expression Correlation (External Data).**\n*   **Temporal Analysis** (e.g., variant submission dates, if available).\n*   **Structural Variant Impact.**\n*   **Functional Impact Scores Comparison** (e.g., CADD, SIFT, PolyPhen from external annotations).\n*   **Conservation Analysis.**""")

with main_tab3:
    st.header(f"{selected_gene} HPO Phenotype Explorer")
    if df_hpo_phenotypes.empty and df_hpo_diseases.empty and hpo_load_error:
         st.warning(hpo_load_error) 
    elif df_hpo_phenotypes.empty and df_hpo_diseases.empty:
        st.warning(f"No HPO phenotype or disease data could be loaded or found for {selected_gene}. Please ensure '{selected_gene}_annotations.json' exists in '{HPO_FILE_DIR}'.")
    else:
        st.info(f"Displaying HPO (Human Phenotype Ontology) annotations for {selected_gene}.")
        st.markdown("---")
        
        hpo_tab1, hpo_tab2, hpo_tab3, hpo_tab4, hpo_tab5 = st.tabs([
            "Phenotype Overview", 
            "Disease Associations", 
            "Phenotype Categories",
            "Phenotype Wordcloud",
            "Gene-Phenotype Network"
        ])
        
        with hpo_tab1:
            st.subheader("Phenotype Overview")
            
            if not df_hpo_phenotypes.empty:
                pheno_count = len(df_hpo_phenotypes['id'].unique()) 
                st.metric(label=f"Unique HPO Phenotypes for {selected_gene}", value=f"{pheno_count:,}")
            else:
                st.metric(label=f"Unique HPO Phenotypes for {selected_gene}", value="0")
            
            st.markdown("##### Top HPO Phenotypes")
            hpo_top_n_slider_overview = st.slider(f"Select Top N HPO Phenotypes to display details:", 5, 50, 15, 1, key="hpo_top_n_slider_overview_tab1")
            
            if not df_hpo_phenotypes.empty:
                top_phenos_df = create_hpo_top_phenotypes_table_df(df_hpo_phenotypes, top_n=hpo_top_n_slider_overview)
                st.dataframe(top_phenos_df, use_container_width=True, hide_index=True, height=(len(top_phenos_df)+1)*35+3) 
            else:
                st.info("No HPO phenotype data to display in table.")


        with hpo_tab2:
            st.subheader("Associated Diseases (from HPO annotations)")
            if not df_hpo_diseases.empty:
                disease_count = len(df_hpo_diseases['id'].unique())
                st.metric(label=f"Unique HPO Disease Associations for {selected_gene}", value=f"{disease_count:,}")
                
                st.markdown("##### Disease Details Table")
                display_disease_df = df_hpo_diseases[['name', 'id', 'mondoId']].copy().drop_duplicates()
                display_disease_df.rename(columns={'name':'Disease Name', 'id':'HPO Disease ID', 'mondoId': 'Mondo ID'}, inplace=True)
                if 'Mondo ID' in display_disease_df.columns:
                    display_disease_df['Mondo Link'] = display_disease_df['Mondo ID'].apply(
                        lambda x: f"https://monarchinitiative.org/disease/{x}" if pd.notna(x) and x else ""
                    )
                    st.dataframe(display_disease_df, use_container_width=True, hide_index=True,
                                 column_config={
                                     "Mondo ID": st.column_config.TextColumn("Mondo ID"),
                                     "Mondo Link": st.column_config.LinkColumn("Link", display_text="🔗 Monarch")
                                 },
                                 column_order=["Disease Name", "HPO Disease ID", "Mondo ID", "Mondo Link"]
                                 )
                else:
                    st.dataframe(display_disease_df, use_container_width=True, hide_index=True)
            else: 
                st.metric(label=f"Unique HPO Disease Associations for {selected_gene}", value="0")
                st.info("No specific disease associations found in HPO data for this gene.")

        with hpo_tab3:
            st.subheader("Phenotype Categories (based on HPO terms)")
            st.plotly_chart(plot_hpo_phenotype_categories_stacked_bar(df_hpo_phenotypes, selected_gene), use_container_width=True)
        
        with hpo_tab4:
            st.subheader("Phenotype Word Cloud")
            if not df_hpo_phenotypes.empty:
                img_path = create_hpo_wordcloud_image(df_hpo_phenotypes, selected_gene)
                if img_path and os.path.exists(img_path):
                    st.image(img_path, caption=f"Word Cloud for {selected_gene} HPO Phenotypes", use_container_width=True)
                elif img_path is None and not df_hpo_phenotypes.empty : st.warning("Could not generate word cloud.")
                else: st.info("No HPO phenotype names available to generate a word cloud.") 
            else: st.info("No HPO phenotype data to generate a word cloud.")
        
        with hpo_tab5:
            st.subheader("Gene - HPO Phenotype - Disease Network (Simplified)")
            net_col1, net_col2 = st.columns(2)
            with net_col1:
                net_top_pheno = st.slider("Top N Phenotypes for Network:", 3, 20, 7, key="net_top_pheno")
            with net_col2:
                net_top_dis = st.slider("Top N Diseases for Network:", 1, 10, 3, key="net_top_dis")
            st.plotly_chart(plot_hpo_gene_phenotype_network_plotly(df_hpo_phenotypes, df_hpo_diseases, selected_gene, top_n_pheno=net_top_pheno, top_n_dis=net_top_dis), use_container_width=True)
