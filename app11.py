#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CHD1-9 ClinVar & HPO Explorer: Interactive Positional Analysis (Streamlit Version)
Modifications:
- Updated UniProt IDs in GENE_IDENTIFIERS based on provided image.
- Moved Contact and Version information to the bottom of the sidebar.
- ... (all previous modifications) ...
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

# ────────────────────────── PAGE CONFIGURATION (MUST BE FIRST STREAMLIT COMMAND) ─────────────────
st.set_page_config(layout="wide", page_title="CHD Gene Family Variant & HPO Explorer")

# --- Custom CSS ---
st.markdown("""
<style>
/* General Styles */
body {
    font-family: 'Segoe UI', 'Roboto', 'Open Sans', 'Lato', sans-serif;
    font-size: 16px; 
    color: #333; 
    background-color: #f9f9f9; 
}

/* Headings */
h1, .stTitle { 
    color: #2c3e50; 
    font-weight: 600;
}
h2, .stHeader { 
    color: #34495e; 
    font-weight: 500;
    margin-top: 2em;
    margin-bottom: 0.8em;
    border-bottom: 1px solid #e0e0e0; 
    padding-bottom: 0.3em;
}
h3, .stSubheader { 
    color: #34495e;
    font-weight: 500;
    margin-top: 1.5em;
    margin-bottom: 0.6em;
}

/* Sidebar specific styling */
.stSidebar > div:first-child { 
    background-color: #f0f2f6; 
    padding: 20px; 
    border-right: 1px solid #ddd; 
}
.stSidebar .stRadio > label { 
    font-size: 1.05em; 
}
.stSidebar .stSlider > label { 
     font-size: 1em;
     font-weight: 500;
     color: #2c3e50;
     margin-bottom: 0.5rem; 
}
.stSidebar .stMultiSelect > label { 
    font-size: 1em;
    font-weight: 500;
    color: #2c3e50;
    margin-bottom: 0.5rem;
}
.sidebar-footer { /* Custom class for sidebar footer */
    font-size: 0.8em;
    color: #555;
    padding-top: 15px;
    border-top: 1px solid #ddd;
    margin-top: 20px;
}
.sidebar-footer a {
    color: #007bff;
}


/* Tab styling */
div[data-baseweb="tab-list"] button[data-baseweb="tab"] {
    font-size: 1.1em; 
    padding: 12px 20px; 
    font-weight: 500; 
    color: #555; 
    border-bottom: 3px solid transparent; 
    transition: all 0.25s ease-in-out;
}
div[data-baseweb="tab-list"] button[data-baseweb="tab"][aria-selected="true"] {
    color: #007bff; 
    border-bottom: 3px solid #007bff; 
    font-weight: 600;
}

/* Multiselect tags */
.stMultiSelect [data-baseweb="tag"] {
    background-color: #007bff; 
    color: white; 
    border-radius: 15px; 
    padding: 4px 10px; 
    font-size: 0.9em;
}

/* Link buttons styling (targets st.link_button) */
.stButton > button[kind="secondaryLink"], .stButton > button[kind="secondary"] {
    border: 1px solid #007bff !important;
    color: #007bff !important;
    background-color: transparent !important;
    border-radius: 5px !important;
    padding: 0.4em 0.8em !important; 
    font-weight: 500 !important;
    transition: background-color 0.2s ease-in-out, color 0.2s ease-in-out;
}
.stButton > button[kind="secondaryLink"]:hover, .stButton > button[kind="secondary"]:hover {
    background-color: #007bff !important;
    color: white !important;
    border-color: #0056b3 !important;
}
.stButton > button[kind="secondaryLink"]:focus, .stButton > button[kind="secondary"]:focus {
    box-shadow: 0 0 0 0.2rem rgba(0,123,255,.5) !important; 
}


/* General Info/Warning boxes */
.stAlert { 
    border-radius: 6px;
    border-left-width: 5px;
    padding: 12px;
    font-size: 0.95em;
    margin-top: 1em;
    margin-bottom: 1em;
}

/* Markdown links */
a {
    color: #0066cc; 
    text-decoration: none;
}
a:hover {
    text-decoration: underline;
    color: #004C99;
}
</style>
""", unsafe_allow_html=True)


# ────────────────────────── CONFIGURATION ──────────────────────────
ALL_GENES = [f"CHD{i}" for i in range(1, 10)]
DEFAULT_GENE_LENGTH_FALLBACK = 1000

# UPDATED GENE_IDENTIFIERS based on your image
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
PLOT_FONT_FAMILY = "Arial, sans-serif"; PLOT_FONT_SIZE_MAIN_TITLE = 18; PLOT_FONT_SIZE_AXIS_TITLE = 12; PLOT_FONT_SIZE_AXIS_TICKS = 10; PLOT_FONT_SIZE_LEGEND_TITLE = 11; PLOT_FONT_SIZE_LEGEND_ITEMS = 10
DOMAIN_TRACK_PAPER_HEIGHT = 0.09; MAIN_PLOT_Y_DOMAIN_BOTTOM = DOMAIN_TRACK_PAPER_HEIGHT + 0.02; FIG_HEIGHT = 600; GENE_EXPLORER_OTHER_PLOT_HEIGHT = 450; COMPARISON_PLOT_HEIGHT = 500; COMPLEX_PLOT_HEIGHT = 700 
PATHOGENIC_SIGNIFICANCES = ["Pathogenic", "Likely pathogenic"]; DEFAULT_TOP_N_PHENOTYPES = 10 
TOP_N_PHENOTYPES_GENE_EXPLORER = 7 
HPO_PLOT_HEIGHT = 500
HPO_FILE_DIR = "chd_phenotype_data" 

# ────────────────────────── DATA PREPROCESSING HELPER FUNCTIONS ──────────────────────────
# (All helper functions as before)
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

# ────────────────────────── GENE-SPECIFIC DATA LOADING & PREPROCESSING ──────────────────────────
# (load_gene_data and load_all_chd_data_for_comparison as before - GENE_IDENTIFIERS will be used in load_all_chd_data_for_comparison)
@st.cache_data
def load_gene_data(gene_name):
    base_path = os.path.dirname(os.path.abspath(__file__))
    clinvar_file_path = os.path.join(base_path, f"{gene_name}_clinvar.txt"); domains_file_path = os.path.join(base_path, f"{gene_name}_domains.json")
    df_clinvar_data = pd.DataFrame(); domains_data_list = []; error_messages = []
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

# ─────────────────── HPO DATA LOADING & PREPROCESSING ───────────────────
# (load_hpo_data and classify_hpo_phenotype_category as before)
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
# (All ClinVar plotting functions for Gene Explorer tab as before)
def add_domain_legend_traces(fig, domains_data_local):
    for i, domain in enumerate(domains_data_local):
        color = DOMAIN_COLORS[i % len(DOMAIN_COLORS)]; name = domain.get("name", f"Domain {i+1}")
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',marker=dict(color=color, size=10, symbol='square'),name=name,legendgroup="domains", legendgrouptitle_text="Protein Domains"))
def create_domain_shapes_paper(domains_data_local, gene_length_local):
    shapes = [dict(type="rect", xref="x", yref="paper",x0=0, y0=0, x1=gene_length_local, y1=DOMAIN_TRACK_PAPER_HEIGHT,fillcolor="#F0F2F6", line_width=0.5, line_color="#D1D5DB", layer="below", opacity=1)]
    if not domains_data_local: return shapes
    for i, domain in enumerate(domains_data_local):
        color = DOMAIN_COLORS[i % len(DOMAIN_COLORS)]; start, end, name = domain["start"], domain["end"], domain.get("name", f"Domain {i+1}")
        if end - start > 0: shapes.append(dict(type="rect", xref="x", yref="paper",x0=start, y0=DOMAIN_TRACK_PAPER_HEIGHT * 0.1,x1=end,   y1=DOMAIN_TRACK_PAPER_HEIGHT * 0.9,fillcolor=color, line_color='#374151', line_width=1, layer="below", opacity=0.9,name=name))
    return shapes
def _setup_single_plot_with_domain_space(fig_title, y_axis_title, y_axis_range, gene_length_local):
    fig = make_subplots(rows=1, cols=1, x_title="Protein Position (aa)")
    fig.update_layout(title_text=fig_title, title_font_size=PLOT_FONT_SIZE_MAIN_TITLE,font_family=PLOT_FONT_FAMILY, plot_bgcolor='white',xaxis=dict(range=[-gene_length_local * 0.02, gene_length_local * 1.02],gridcolor='rgba(200,200,200,0.3)', zeroline=False,tickfont_size=PLOT_FONT_SIZE_AXIS_TICKS, title_font_size=PLOT_FONT_SIZE_AXIS_TITLE,title_standoff=15),yaxis=dict(title_text=y_axis_title,domain=[MAIN_PLOT_Y_DOMAIN_BOTTOM, 1.0],range=y_axis_range,gridcolor='rgba(200,200,200,0.3)', zeroline=False,tickfont_size=PLOT_FONT_SIZE_AXIS_TICKS, title_font_size=PLOT_FONT_SIZE_AXIS_TITLE, title_standoff=10),height=FIG_HEIGHT,margin=dict(t=50, b=80, l=80, r=220),legend=dict(yanchor="top", y=1, xanchor="left", x=1.01,bgcolor='rgba(255,255,255,0.95)', bordercolor="#ADB5BD", borderwidth=1,font_size=PLOT_FONT_SIZE_LEGEND_ITEMS, title_font_size=PLOT_FONT_SIZE_LEGEND_TITLE, tracegroupgap=15))
    return fig
def create_interactive_lollipop(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name): 
    if df_positions_filtered.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Lollipop (No Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
    significance_hierarchy = list(SIGNIFICANCE_MARKERS.keys()); df_copy = df_positions_filtered.copy(); df_copy.loc[:, 'SignificanceOrdered'] = pd.Categorical(df_copy['SignificanceClass'], categories=significance_hierarchy, ordered=True)
    lollipop_agg_data = df_copy.groupby('Position', as_index=False).agg(TotalCount=('VariationID', 'nunique'), DominantSignificance=('SignificanceOrdered', 'min'), VariantTypes=('StandardType', lambda x: ', '.join(sorted(list(x.unique())))), VariantNames=('Name', lambda x: '<br>'.join(sorted(list(x.unique()))[:5]) + (f'<br>...and {len(x.unique())-5} more' if len(x.unique()) > 5 else '')), Phenotypes=('PhenotypeListClean', lambda x: '<br>'.join(sorted(list(x.unique()))[:3]) + (f'<br>...and {len(x.unique())-3} more' if len(x.unique()) > 3 else '')))
    lollipop_agg_data.loc[:, 'DominantSignificance'] = lollipop_agg_data['DominantSignificance'].astype(str)
    if lollipop_agg_data.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Lollipop (No Aggregated Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
    max_total_count = lollipop_agg_data['TotalCount'].max() if not lollipop_agg_data.empty else 1; y_min = -max_total_count * 0.05 if max_total_count > 1 else -0.1; y_max = max_total_count * 1.15 if max_total_count > 0 else 1.15
    fig = _setup_single_plot_with_domain_space(f"{current_gene_name} Variant Lollipop Plot", "Unique Variants per Position", [y_min, y_max], gene_length_local); fig.update_layout(yaxis_zeroline=True, yaxis_zerolinecolor='rgba(200,200,200,0.5)')
    for _, row in lollipop_agg_data.iterrows(): fig.add_trace(go.Scatter(x=[row['Position'], row['Position']], y=[0, row['TotalCount']], mode='lines', line=dict(color='#BCCCDC', width=1), hoverinfo='none', showlegend=False))
    min_px_size, max_px_size = 5, 18; scaled_plotly_sizes = pd.Series(min_px_size, index=lollipop_agg_data.index)
    if max_total_count > 0: log_counts = np.log1p(lollipop_agg_data['TotalCount']); log_max_count = np.log1p(max_total_count); scaled_plotly_sizes = min_px_size + (max_px_size - min_px_size) * (log_counts / log_max_count) if log_max_count > 0 else pd.Series(min_px_size, index=lollipop_agg_data.index)
    for sig_class in significance_hierarchy:
        subset = lollipop_agg_data[lollipop_agg_data['DominantSignificance'] == sig_class]
        if not subset.empty: fig.add_trace(go.Scatter(x=subset['Position'], y=subset['TotalCount'], mode='markers', marker=dict(color=SIGNIFICANCE_PALETTE.get(sig_class, '#CCCCCC'), size=scaled_plotly_sizes.loc[subset.index], symbol=SIGNIFICANCE_MARKERS.get(sig_class, 'circle'), line=dict(width=0.7, color='#1F2937')), name=sig_class, legendgroup="significance", legendgrouptitle_text="Variant Significance", customdata=subset[['Position', 'TotalCount', 'DominantSignificance', 'VariantTypes', 'VariantNames', 'Phenotypes']], hovertemplate=("<b>Position:</b> %{customdata[0]}<br><b>Variant Count:</b> %{customdata[1]}<br><b>Dominant Significance:</b> %{customdata[2]}<br><b>Variant Types:</b> %{customdata[3]}<br><b>Names:</b><br>%{customdata[4]}<br><b>Phenotypes:</b><br>%{customdata[5]}<extra></extra>")))
    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local); fig.update_layout(shapes=domain_shapes);
    if domains_data_local: add_domain_legend_traces(fig, domains_data_local)
    return fig
def create_interactive_waterfall(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name): 
    if df_positions_filtered.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Waterfall (No Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
    agg_for_stacking = df_positions_filtered.groupby(["Position", "SignificanceClass"]).size().unstack(fill_value=0)
    sig_order_plot = [s for s in reversed(list(SIGNIFICANCE_PALETTE.keys())) if s in agg_for_stacking.columns] + [s for s in agg_for_stacking.columns if s not in SIGNIFICANCE_PALETTE.keys()]
    agg_for_stacking = agg_for_stacking.reindex(columns=sig_order_plot, fill_value=0); max_y_val_overall = agg_for_stacking.sum(axis=1).max() if not agg_for_stacking.empty else 0
    fig = _setup_single_plot_with_domain_space(f"{current_gene_name} Variant Waterfall Plot", "Cumulative Variant Count", [0, max_y_val_overall * 1.1 if max_y_val_overall > 0 else 1], gene_length_local); fig.update_layout(barmode='stack')
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
    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local); fig.update_layout(shapes=domain_shapes)
    if domains_data_local: add_domain_legend_traces(fig, domains_data_local)
    return fig
def create_interactive_density(df_positions_filtered, domains_data_local, gene_length_local, current_gene_name):
    if df_positions_filtered.empty: return go.Figure().update_layout(title_text=f"{current_gene_name} Density (No Data)", height=FIG_HEIGHT, xaxis_title="Protein Position (aa)", xaxis_range=[-gene_length_local * 0.02, gene_length_local * 1.02])
    fig = _setup_single_plot_with_domain_space(f"{current_gene_name} Variant Density Plot", "", [-0.5, 0.5], gene_length_local); fig.update_layout(yaxis_showgrid=False, yaxis_showticklabels=False, yaxis_zeroline=False)
    df_copy = df_positions_filtered.copy(); df_copy['y_jitter'] = np.random.uniform(-0.4, 0.4, size=len(df_copy))
    significance_plot_order = [s for s in SIGNIFICANCE_MARKERS.keys()]
    for sig_class in significance_plot_order:
        group_df = df_copy[df_copy['SignificanceClass'] == sig_class]
        if group_df.empty: continue
        fig.add_trace(go.Scatter(x=group_df['Position'], y=group_df['y_jitter'], mode='markers',marker=dict(color=SIGNIFICANCE_PALETTE.get(sig_class, '#CCCCCC'), symbol=SIGNIFICANCE_MARKERS.get(sig_class, 'circle'), size=7, opacity=0.65, line=dict(width=0.5, color='#343A40')),name=sig_class, legendgroup="significance", legendgrouptitle_text="Variant Significance",customdata=group_df[['Position', 'Name', 'StandardType', 'SignificanceClass', 'PhenotypeListClean']],hovertemplate=("<b>Position:</b> %{customdata[0]}<br><b>Name:</b> %{customdata[1]}<br><b>Type:</b> %{customdata[2]}<br><b>Significance:</b> %{customdata[3]}<br><b>Phenotype:</b> %{customdata[4]}<extra></extra>")))
    domain_shapes = create_domain_shapes_paper(domains_data_local, gene_length_local); fig.update_layout(shapes=domain_shapes)
    if domains_data_local: add_domain_legend_traces(fig, domains_data_local)
    return fig
def create_origin_pie_chart(df_gene_variants_filtered, gene_name):
    if df_gene_variants_filtered.empty or 'OriginSimple' not in df_gene_variants_filtered.columns: return go.Figure().update_layout(title_text=f"{gene_name} Variant Origin (No Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    origin_counts = df_gene_variants_filtered['OriginSimple'].value_counts()
    if origin_counts.empty: return go.Figure().update_layout(title_text=f"{gene_name} Variant Origin (No Origin Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    fig = go.Figure(data=[go.Pie(labels=origin_counts.index, values=origin_counts.values, hole=.3,marker_colors=px.colors.qualitative.Pastel2)])
    fig.update_traces(textposition='inside', textinfo='percent+label')
    fig.update_layout(title_text=f'{gene_name} Variant Origin Distribution',height=GENE_EXPLORER_OTHER_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,legend_title_text='Origin')
    return fig
def create_type_significance_stacked_bar(df_gene_variants_filtered, gene_name):
    if df_gene_variants_filtered.empty: return go.Figure().update_layout(title_text=f"{gene_name} Type vs Significance (No Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    counts = df_gene_variants_filtered.groupby(['StandardType', 'SignificanceClass']).size().reset_index(name='Count')
    if counts.empty: return go.Figure().update_layout(title_text=f"{gene_name} Type vs Significance (No Counts)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    all_types = sorted(list(VARIANT_TYPE_PALETTE.keys())); all_significances = list(SIGNIFICANCE_PALETTE.keys())
    pivot_df = counts.pivot_table(index='StandardType', columns='SignificanceClass', values='Count', fill_value=0)
    pivot_df = pivot_df.reindex(index=all_types, columns=all_significances, fill_value=0).dropna(how='all', axis=0)
    fig = go.Figure()
    for sig_class in all_significances:
        if sig_class in pivot_df.columns and pivot_df[sig_class].sum() > 0: fig.add_trace(go.Bar(name=sig_class,x=pivot_df.index,y=pivot_df[sig_class],marker_color=SIGNIFICANCE_PALETTE.get(sig_class)))
    fig.update_layout(barmode='stack',title_text=f'{gene_name}: Clinical Significance by Variant Type',xaxis_title='Variant Type',yaxis_title='Number of Variants',legend_title_text='Clinical Significance',height=GENE_EXPLORER_OTHER_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,xaxis={'categoryorder':'array', 'categoryarray': [t for t in all_types if t in pivot_df.index]})
    return fig
def create_phenotype_type_stacked_bar(df_gene_variants_filtered, gene_name, normalize=False, top_n=TOP_N_PHENOTYPES_GENE_EXPLORER):
    if df_gene_variants_filtered.empty or 'PhenotypeListClean' not in df_gene_variants_filtered.columns: return go.Figure().update_layout(title_text=f"{gene_name} Phenotype vs Type (No Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    df_pheno_exploded = df_gene_variants_filtered.copy(); df_pheno_exploded['PhenotypeSingle'] = df_pheno_exploded['PhenotypeListClean'].str.split('; ')
    df_pheno_exploded = df_pheno_exploded.explode('PhenotypeSingle').dropna(subset=['PhenotypeSingle']); df_pheno_exploded = df_pheno_exploded[df_pheno_exploded['PhenotypeSingle'] != 'N/A']
    if df_pheno_exploded.empty: return go.Figure().update_layout(title_text=f"{gene_name} Phenotype vs Type (No Phenotype Data)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    top_phenos_list = df_pheno_exploded['PhenotypeSingle'].value_counts().nlargest(top_n).index.tolist()
    df_pheno_exploded['DisplayPhenotype'] = df_pheno_exploded['PhenotypeSingle'].apply(lambda x: x if x in top_phenos_list else "Other Phenotypes")
    phenotype_order = top_phenos_list; 
    if "Other Phenotypes" in df_pheno_exploded['DisplayPhenotype'].unique() and "Other Phenotypes" not in phenotype_order: phenotype_order.append("Other Phenotypes")
    counts = df_pheno_exploded.groupby(['DisplayPhenotype', 'StandardType']).size().reset_index(name='Count')
    if counts.empty: return go.Figure().update_layout(title_text=f"{gene_name} Phenotype vs Type (No Counts)", height=GENE_EXPLORER_OTHER_PLOT_HEIGHT)
    pivot_df = counts.pivot_table(index='DisplayPhenotype', columns='StandardType', values='Count', fill_value=0)
    if normalize: pivot_df = pivot_df.apply(lambda x: (x / x.sum() * 100) if x.sum() > 0 else x, axis=1)
    all_types = sorted(list(VARIANT_TYPE_PALETTE.keys())); pivot_df = pivot_df.reindex(index=phenotype_order, columns=all_types, fill_value=0).dropna(how='all', axis=0)
    fig = go.Figure()
    for var_type in all_types:
        if var_type in pivot_df.columns and pivot_df[var_type].sum() > 0: fig.add_trace(go.Bar(name=var_type,x=pivot_df.index,y=pivot_df[var_type],marker_color=VARIANT_TYPE_PALETTE.get(var_type)))
    yaxis_title = "Percentage of Variants (%)" if normalize else "Number of Variant Associations"; plot_title = f'{gene_name}: Relative Variant Types by Phenotype (Top {top_n})' if normalize else f'{gene_name}: Variant Types by Phenotype (Top {top_n})'
    fig.update_layout(barmode='stack',title_text=plot_title,xaxis_title='Phenotype',yaxis_title=yaxis_title,legend_title_text='Variant Type',height=GENE_EXPLORER_OTHER_PLOT_HEIGHT + 100,font_family=PLOT_FONT_FAMILY,xaxis={'categoryorder':'array', 'categoryarray': [p for p in phenotype_order if p in pivot_df.index]},xaxis_tickangle=-30)
    return fig

# ─────────────────── PLOTLY FIGURE CREATION FUNCTIONS (HPO Phenotype Explorer Tab - NEW) ───────────────────
# (All HPO plotting functions as before)
def plot_hpo_phenotype_counts_bar(df_hpo_phenotypes, gene_name):
    if df_hpo_phenotypes.empty: return go.Figure().update_layout(title_text=f"{gene_name} HPO Phenotypes (No Data)", height=HPO_PLOT_HEIGHT)
    phenotype_count = len(df_hpo_phenotypes)
    fig = go.Figure(data=[go.Bar(x=[gene_name], y=[phenotype_count], text=[phenotype_count], textposition='auto', marker_color='skyblue',name='HPO Phenotype Count')])
    fig.update_layout(title_text=f"Total Number of HPO Phenotypes for {gene_name}",xaxis_title="Gene",yaxis_title="Number of HPO Phenotypes",height=HPO_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,showlegend=False)
    return fig
def plot_hpo_disease_associations_bar(df_hpo_diseases, gene_name):
    if df_hpo_diseases.empty: return go.Figure().update_layout(title_text=f"{gene_name} HPO Disease Associations (No Data)", height=HPO_PLOT_HEIGHT)
    disease_count = len(df_hpo_diseases)
    fig = go.Figure(data=[go.Bar(x=[gene_name], y=[disease_count], text=[disease_count], textposition='auto', marker_color='lightcoral')])
    fig.update_layout(title_text=f"Number of HPO Disease Associations for {gene_name}",xaxis_title="Gene",yaxis_title="Number of Associated Diseases (HPO)",height=HPO_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY)
    return fig
def plot_hpo_phenotype_categories_stacked_bar(df_hpo_phenotypes, gene_name):
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns: return go.Figure().update_layout(title_text=f"{gene_name} HPO Phenotype Categories (No Data)", height=HPO_PLOT_HEIGHT)
    df_hpo_phenotypes['Category'] = df_hpo_phenotypes['name'].apply(classify_hpo_phenotype_category)
    category_counts = df_hpo_phenotypes.groupby('Category').size().reset_index(name='Count')
    category_order = category_counts.sort_values('Count', ascending=False)['Category'].tolist()
    fig = px.bar(category_counts, x='Category', y='Count', color='Category', title=f"HPO Phenotype Categories for {gene_name}",labels={'Count': 'Number of HPO Phenotypes', 'Category': 'Phenotype Category'},height=HPO_PLOT_HEIGHT,category_orders={"Category": category_order})
    fig.update_layout(font_family=PLOT_FONT_FAMILY, showlegend=False)
    return fig
def plot_hpo_phenotype_frequency_bar(df_hpo_phenotypes, gene_name, top_n=15):
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns: return go.Figure().update_layout(title_text=f"{gene_name} HPO Phenotype Frequency (No Data)", height=HPO_PLOT_HEIGHT + 100)
    phenotype_counts = df_hpo_phenotypes['name'].value_counts().nlargest(top_n).reset_index(); phenotype_counts.columns = ['Phenotype', 'Frequency']
    fig = px.bar(phenotype_counts.sort_values('Frequency', ascending=True), x='Frequency', y='Phenotype', orientation='h',title=f"Top {top_n} HPO Phenotypes for {gene_name}",labels={'Frequency': 'Frequency', 'Phenotype': 'HPO Phenotype Term'},height=HPO_PLOT_HEIGHT + 50 + (top_n * 15),color_discrete_sequence=['skyblue'])
    fig.update_layout(font_family=PLOT_FONT_FAMILY, yaxis={'categoryorder':'total ascending'})
    return fig
def create_hpo_wordcloud_image(df_hpo_phenotypes, gene_name): 
    if df_hpo_phenotypes.empty or 'name' not in df_hpo_phenotypes.columns: return None
    phenotype_names = df_hpo_phenotypes["name"].tolist()
    if not phenotype_names: return None
    text = " ".join(phenotype_names)
    try:
        wordcloud = WordCloud(width=800, height=400, background_color="white", colormap="viridis").generate(text)
        img_path = f"{gene_name}_hpo_wordcloud.png" 
        wordcloud.to_file(img_path)
        return img_path
    except Exception as e: st.error(f"Error generating wordcloud for {gene_name}: {e}"); return None
def plot_hpo_gene_phenotype_network_plotly(df_hpo_phenotypes, df_hpo_diseases, gene_name, top_n_pheno=10, top_n_dis=5):
    if df_hpo_phenotypes.empty and df_hpo_diseases.empty: return go.Figure().update_layout(title_text=f"{gene_name} HPO Network (No Data)", height=COMPLEX_PLOT_HEIGHT)
    G = nx.Graph(); G.add_node(gene_name, type="gene", size=30, color='red', label=gene_name)
    pheno_df_top = df_hpo_phenotypes.head(top_n_pheno) if not df_hpo_phenotypes.empty else pd.DataFrame()
    for _, row in pheno_df_top.iterrows():
        pheno_label = row['name'][:30] + '...' if len(row['name']) > 30 else row['name']
        G.add_node(row['id'], type="hpo_phenotype", size=15, color='lightblue', label=pheno_label); G.add_edge(gene_name, row['id'])
    dis_df_top = df_hpo_diseases.head(top_n_dis) if not df_hpo_diseases.empty else pd.DataFrame()
    for _, row in dis_df_top.iterrows():
        dis_label = row['name'][:30] + '...' if len(row['name']) > 30 else row['name']
        G.add_node(row['id'], type="hpo_disease", size=20, color='lightgreen', label=dis_label); G.add_edge(gene_name, row['id'])
    if not G.nodes(): return go.Figure().update_layout(title_text=f"{gene_name} HPO Network (No Nodes)", height=COMPLEX_PLOT_HEIGHT)
    pos = nx.spring_layout(G, k=0.8, iterations=50); edge_x = []; edge_y = []
    for edge in G.edges(): x0, y0 = pos[edge[0]]; x1, y1 = pos[edge[1]]; edge_x.extend([x0, x1, None]); edge_y.extend([y0, y1, None])
    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')
    node_x = []; node_y = []; node_text_hover = []; node_text_display = []; node_size = []; node_color = []
    for node in G.nodes():
        x, y = pos[node]; node_x.append(x); node_y.append(y)
        node_text_hover.append(f"{G.nodes[node]['label']} ({G.nodes[node]['type']})")
        node_text_display.append(G.nodes[node]['label'] if G.nodes[node]['type'] == 'gene' else G.nodes[node]['label'][:10] + '...' if len(G.nodes[node]['label']) > 10 else G.nodes[node]['label'])
        node_size.append(G.nodes[node]['size']); node_color.append(G.nodes[node]['color'])
    node_trace = go.Scatter(x=node_x, y=node_y,mode='markers+text', hoverinfo='text',text= node_text_display, hovertext=node_text_hover, marker=dict(showscale=False, size=node_size, sizemode='diameter', color=node_color, line_width=1, line_color='black'),textfont=dict(size=8, color='black'), textposition="bottom center")
    fig = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(title=dict(text=f'Network of {gene_name} with Top HPO Phenotypes & Diseases',font=dict(size=16, family=PLOT_FONT_FAMILY)),showlegend=False, hovermode='closest',height=COMPLEX_PLOT_HEIGHT,margin=dict(b=20, l=5, r=5, t=40),xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
    return fig

# ─────────────────── PLOTLY FIGURE CREATION FUNCTIONS (Comparison Tab - Unchanged) ───────────────────
# (All comparison tab plotting functions as before)
def plot_significance_distribution_per_gene(all_variants_df): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Significance Distribution", height=COMPARISON_PLOT_HEIGHT)
    counts = all_variants_df.groupby(['Gene', 'SignificanceClass']).size().reset_index(name='Count')
    all_sig_classes = list(SIGNIFICANCE_PALETTE.keys())
    pivot_df = counts.pivot(index='Gene', columns='SignificanceClass', values='Count').fillna(0)
    for sig_class in all_sig_classes:
        if sig_class not in pivot_df.columns: pivot_df[sig_class] = 0
    pivot_df = pivot_df[all_sig_classes].reindex(ALL_GENES).fillna(0)
    fig = go.Figure()
    for sig_class in all_sig_classes:
        if sig_class in pivot_df.columns: fig.add_trace(go.Bar(name=sig_class,x=pivot_df.index,y=pivot_df[sig_class],marker_color=SIGNIFICANCE_PALETTE.get(sig_class)))
    fig.update_layout(barmode='stack',title_text='Clinical Significance Distribution per CHD Gene',xaxis_title='Gene',yaxis_title='Number of Variants',legend_title_text='Clinical Significance',height=COMPARISON_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY)
    return fig
def plot_pathogenic_variant_types_heatmap(all_variants_df): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Pathogenic Variant Types Heatmap", height=COMPARISON_PLOT_HEIGHT)
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]
    if patho_df.empty: return go.Figure().update_layout(title_text="No Pathogenic/Likely Pathogenic variants for Heatmap", height=COMPARISON_PLOT_HEIGHT)
    counts = patho_df.groupby(['Gene', 'StandardType']).size().reset_index(name='Count')
    pivot_df = counts.pivot(index='Gene', columns='StandardType', values='Count').fillna(0).reindex(index=ALL_GENES).fillna(0)
    sorted_variant_types = sorted(pivot_df.columns.tolist())
    pivot_df = pivot_df[sorted_variant_types]
    fig = go.Figure(data=go.Heatmap(z=pivot_df.values,x=pivot_df.columns,y=pivot_df.index,colorscale='Blues',text=pivot_df.values,texttemplate="%{text}",hoverongaps=False))
    fig.update_layout(title_text='Heatmap of Pathogenic/Likely Pathogenic Variant Types per CHD Gene',xaxis_title='Variant Type',yaxis_title='Gene',height=COMPARISON_PLOT_HEIGHT + 100,font_family=PLOT_FONT_FAMILY,xaxis_tickangle=-45)
    return fig
def plot_pathogenic_variant_types_per_gene_stacked(all_variants_df): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Pathogenic Variant Types Distribution", height=COMPARISON_PLOT_HEIGHT)
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]
    if patho_df.empty: return go.Figure().update_layout(title_text="No Pathogenic/Likely Pathogenic variants for Type Distribution", height=COMPARISON_PLOT_HEIGHT)
    counts = patho_df.groupby(['Gene', 'StandardType']).size().reset_index(name='Count')
    all_variant_types = sorted(patho_df['StandardType'].dropna().unique())
    pivot_df = counts.pivot(index='Gene', columns='StandardType', values='Count').fillna(0)
    for vt in all_variant_types:
        if vt not in pivot_df.columns: pivot_df[vt] = 0
    pivot_df = pivot_df[all_variant_types].reindex(ALL_GENES).fillna(0)
    fig = go.Figure()
    for i, var_type in enumerate(all_variant_types):
        fig.add_trace(go.Bar(name=var_type,x=pivot_df.index,y=pivot_df[var_type],marker_color=VARIANT_TYPE_PALETTE.get(var_type, "#CCCCCC"))) 
    fig.update_layout(barmode='stack',title_text='Distribution of Pathogenic/Likely Pathogenic Variant Types per CHD Gene',xaxis_title='Gene',yaxis_title='Number of Pathogenic Variants',legend_title_text='Variant Type',height=COMPARISON_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY)
    return fig
def plot_overall_pathogenic_phenotype_distribution(all_variants_df, top_n=DEFAULT_TOP_N_PHENOTYPES): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Phenotype Distribution", height=COMPARISON_PLOT_HEIGHT)
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)]
    if patho_df.empty: return go.Figure().update_layout(title_text="No Pathogenic/Likely Pathogenic variants for Phenotype Pie Chart", height=COMPARISON_PLOT_HEIGHT)
    phenotypes_series = patho_df['PhenotypeListClean'].str.split('; ').explode().str.strip()
    phenotypes_series = phenotypes_series[phenotypes_series != 'N/A'].dropna()
    if phenotypes_series.empty: return go.Figure().update_layout(title_text="No valid phenotypes found for Pathogenic variants", height=COMPARISON_PLOT_HEIGHT)
    phenotype_counts = Counter(phenotypes_series); top_phenotypes_data = phenotype_counts.most_common(top_n)
    if not top_phenotypes_data: return go.Figure().update_layout(title_text=f"No phenotypes to display for Top {top_n}", height=COMPARISON_PLOT_HEIGHT)
    labels = [item[0] for item in top_phenotypes_data]; values = [item[1] for item in top_phenotypes_data]
    if len(phenotype_counts) > top_n: other_count = sum(count for pheno, count in phenotype_counts.items() if pheno not in labels); labels.append(f'Other ( {len(phenotype_counts) - top_n} types)'); values.append(other_count)
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3, marker_colors=PHENOTYPE_COLORS[:len(labels)])])
    fig.update_traces(textposition='inside', textinfo='percent+label', pull=[0.05 if i==0 else 0 for i in range(len(labels))])
    fig.update_layout(title_text=f'Overall Distribution of Top {top_n} Pathogenic Phenotypes (All CHD Genes)',height=COMPARISON_PLOT_HEIGHT + 150,font_family=PLOT_FONT_FAMILY,legend_title_text='Phenotype',showlegend=True,legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    return fig
def plot_gene_phenotype_heatmap(all_variants_df, top_n_phenotypes=DEFAULT_TOP_N_PHENOTYPES): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Gene-Phenotype Heatmap", height=COMPLEX_PLOT_HEIGHT)
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)].copy()
    if patho_df.empty: return go.Figure().update_layout(title_text="No Pathogenic variants for Gene-Phenotype Heatmap", height=COMPLEX_PLOT_HEIGHT)
    patho_df['PhenotypeSingle'] = patho_df['PhenotypeListClean'].str.split('; ')
    exploded_pheno_df = patho_df.explode('PhenotypeSingle'); exploded_pheno_df = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'] != 'N/A'].dropna(subset=['PhenotypeSingle'])
    if exploded_pheno_df.empty: return go.Figure().update_layout(title_text="No valid phenotypes for Gene-Phenotype Heatmap", height=COMPLEX_PLOT_HEIGHT)
    top_phenos_list = exploded_pheno_df['PhenotypeSingle'].value_counts().nlargest(top_n_phenotypes).index.tolist()
    if not top_phenos_list: return go.Figure().update_layout(title_text=f"No phenotypes found for Top {top_n_phenotypes}", height=COMPLEX_PLOT_HEIGHT)
    filtered_exploded_df = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'].isin(top_phenos_list)]
    counts = filtered_exploded_df.groupby(['Gene', 'PhenotypeSingle']).size().reset_index(name='Count')
    pivot_df = counts.pivot(index='Gene', columns='PhenotypeSingle', values='Count').fillna(0); pivot_df = pivot_df.reindex(index=ALL_GENES, columns=top_phenos_list).fillna(0)
    fig = go.Figure(data=go.Heatmap(z=pivot_df.values,x=pivot_df.columns,y=pivot_df.index,colorscale='YlOrRd',text=pivot_df.values,texttemplate="%{text}",hoverongaps=False,hovertemplate="Gene: %{y}<br>Phenotype: %{x}<br>Pathogenic Variants: %{z}<extra></extra>"))
    fig.update_layout(title_text=f'Heatmap of Pathogenic Variants by Gene and Top {top_n_phenotypes} Phenotypes',xaxis_title='Phenotype',yaxis_title='Gene',height=COMPLEX_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,xaxis_tickangle=-45)
    return fig
def plot_pathogenic_bubble_chart_gene_phenotype(all_variants_df, top_n_phenotypes=DEFAULT_TOP_N_PHENOTYPES): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)].copy()
    if patho_df.empty: return go.Figure().update_layout(title_text="No Pathogenic variants for Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
    patho_df['PhenotypeSingle'] = patho_df['PhenotypeListClean'].str.split('; ')
    exploded_pheno_df = patho_df.explode('PhenotypeSingle'); exploded_pheno_df = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'] != 'N/A'].dropna(subset=['PhenotypeSingle'])
    if exploded_pheno_df.empty: return go.Figure().update_layout(title_text="No valid phenotypes for Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
    top_phenos_list = exploded_pheno_df['PhenotypeSingle'].value_counts().nlargest(top_n_phenotypes).index.tolist()
    if not top_phenos_list: return go.Figure().update_layout(title_text=f"No phenotypes found for Top {top_n_phenotypes}", height=COMPLEX_PLOT_HEIGHT)
    bubble_data = exploded_pheno_df[exploded_pheno_df['PhenotypeSingle'].isin(top_phenos_list)]
    counts_df = bubble_data.groupby(['Gene', 'PhenotypeSingle']).size().reset_index(name='PathogenicVariantCount')
    if counts_df.empty: return go.Figure().update_layout(title_text=f"No data for selected Top {top_n_phenotypes} phenotypes in Bubble Chart", height=COMPLEX_PLOT_HEIGHT)
    fig = px.scatter(counts_df,x='Gene',y='PhenotypeSingle',size='PathogenicVariantCount',color='Gene',hover_name='PhenotypeSingle',size_max=60,title=f'Bubble Chart of Pathogenic Variants by Gene and Top {top_n_phenotypes} Phenotypes',labels={'Gene': 'Gene', 'PhenotypeSingle': 'Phenotype', 'PathogenicVariantCount': 'Number of Pathogenic Variants'},color_discrete_sequence=GENE_COLORS)
    fig.update_layout(height=COMPLEX_PLOT_HEIGHT,font_family=PLOT_FONT_FAMILY,yaxis={'categoryorder':'total ascending'})
    fig.update_traces(hovertemplate="Gene: %{x}<br>Phenotype: %{y}<br>Count: %{marker.size}<extra></extra>")
    return fig
def plot_sankey_variant_type_to_phenotype(all_variants_df, top_n_phenotypes=DEFAULT_TOP_N_PHENOTYPES): 
    if all_variants_df.empty: return go.Figure().update_layout(title_text="No data for Sankey Diagram", height=COMPLEX_PLOT_HEIGHT)
    patho_df = all_variants_df[all_variants_df['SignificanceClass'].isin(PATHOGENIC_SIGNIFICANCES)].copy()
    if patho_df.empty: return go.Figure().update_layout(title_text="No Pathogenic variants for Sankey Diagram", height=COMPLEX_PLOT_HEIGHT)
    patho_df['PhenotypeSingle'] = patho_df['PhenotypeListClean'].str.split('; ')
    exploded_df = patho_df.explode('PhenotypeSingle'); exploded_df = exploded_df[exploded_df['PhenotypeSingle'] != 'N/A'].dropna(subset=['PhenotypeSingle', 'StandardType'])
    if exploded_df.empty: return go.Figure().update_layout(title_text="No valid phenotype/type data for Sankey", height=COMPLEX_PLOT_HEIGHT)
    top_phenos_list = exploded_df['PhenotypeSingle'].value_counts().nlargest(top_n_phenotypes).index.tolist()
    if not top_phenos_list: return go.Figure().update_layout(title_text=f"No phenotypes found for Top {top_n_phenotypes}", height=COMPLEX_PLOT_HEIGHT)
    sankey_df = exploded_df[exploded_df['PhenotypeSingle'].isin(top_phenos_list)]
    if sankey_df.empty: return go.Figure().update_layout(title_text=f"No data for Top {top_n_phenotypes} phenotypes in Sankey", height=COMPLEX_PLOT_HEIGHT)
    all_variant_types = sorted(sankey_df['StandardType'].unique()); all_phenotypes_sankey = sorted(sankey_df['PhenotypeSingle'].unique()) 
    labels = all_variant_types + all_phenotypes_sankey; label_indices = {label: i for i, label in enumerate(labels)}
    source = []; target = []; value = []; link_colors = []
    link_counts = sankey_df.groupby(['StandardType', 'PhenotypeSingle']).size().reset_index(name='Count')
    type_palette = VARIANT_TYPE_PALETTE 
    for i, row in link_counts.iterrows():
        src_idx = label_indices[row['StandardType']]; tgt_idx = label_indices[row['PhenotypeSingle']]
        source.append(src_idx); target.append(tgt_idx); value.append(row['Count'])
        link_colors.append(type_palette.get(row['StandardType'], "#CCCCCC")) 
    node_colors = [type_palette.get(vt, "#CCCCCC") for vt in all_variant_types] + PHENOTYPE_COLORS[:len(all_phenotypes_sankey)]; node_colors.extend(PHENOTYPE_COLORS * ( (len(labels) - len(node_colors)) // len(PHENOTYPE_COLORS) + 1)); node_colors = node_colors[:len(labels)]
    fig = go.Figure(data=[go.Sankey(node=dict(pad=15,thickness=20,line=dict(color="black", width=0.5),label=labels,color=node_colors),link=dict(source=source,target=target,value=value,color=link_colors,hovertemplate='From %{source.label} to %{target.label}: %{value} variants<extra></extra>'))])
    fig.update_layout(title=dict(text=f'Sankey Diagram: Pathogenic Variant Types to Top {top_n_phenotypes} Phenotypes',font=dict(size=16, family=PLOT_FONT_FAMILY)),height=COMPLEX_PLOT_HEIGHT + 100) 
    return fig


# ────────────────────────── STREAMLIT APP LAYOUT & LOGIC ──────────────────────────
# ... (之前的 import 和 CONFIGURATION 部分保持不变) ...
# ... (所有 DATA PREPROCESSING 和 PLOTTING 函数保持不变) ...
# --- Sidebar Setup ---
st.sidebar.header("Gene Selection")
selected_gene = st.sidebar.radio(
    "Select Gene:",
    ALL_GENES,
    index=ALL_GENES.index("CHD7") if "CHD7" in ALL_GENES else 0,
    key="gene_selector_main"
)

# Filters directly under Gene Selection for the selected gene
st.sidebar.header(f"Filters for {selected_gene} Positional Plots")

# Load ClinVar data needed for filters and gene_length
# Note: df_gene_variants is used for filter options and other plots later
df_gene_variants, gene_domains, gene_length, clinvar_load_errors = load_gene_data(selected_gene)
if clinvar_load_errors:
    for error_msg in clinvar_load_errors:
        st.sidebar.warning(error_msg) # Display ClinVar load errors in sidebar

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

sig_options = sorted(df_gene_variants['SignificanceClass'].dropna().unique(), key=lambda x: list(SIGNIFICANCE_PALETTE.keys()).index(x) if x in SIGNIFICANCE_PALETTE else float('inf')) if not df_gene_variants.empty else []
default_sig = st.session_state.get(f'sel_sig_{selected_gene}', sig_options)
default_sig = [s for s in default_sig if s in sig_options]
selected_significance_filter = st.sidebar.multiselect(
    f"Significance (Positional Plots):",
    sig_options, default_sig,
    key=f"sig_multi_widget_{selected_gene}"
)
st.session_state[f'sel_sig_{selected_gene}'] = selected_significance_filter

type_options = sorted(df_gene_variants['StandardType'].dropna().unique()) if not df_gene_variants.empty else []
default_type = st.session_state.get(f'sel_type_{selected_gene}', type_options)
default_type = [t for t in default_type if t in type_options]
selected_types_filter = st.sidebar.multiselect(
    f"Variant Type (Positional Plots):",
    type_options, default_type,
    key=f"type_multi_widget_{selected_gene}"
)
st.session_state[f'sel_type_{selected_gene}'] = selected_types_filter

# --- Sidebar Footer Information (Consolidated) ---
st.sidebar.markdown("---") # Separator before footer info
# ... (Sidebar footer code remains the same) ...
st.sidebar.markdown("""<div class="sidebar-footer">
<b>Developed by:</b><br>
Zihao Wang<br>
<i>Lab Feng</i><br>
IBS, Fudan University
</div>""", unsafe_allow_html=True)
st.sidebar.markdown("""<div class="sidebar-footer" style="margin-top:10px;">
<b>Data Sources:</b><br>
<a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank">NCBI ClinVar</a><br>
<a href="https://www.uniprot.org/" target="_blank">UniProt</a><br>
<a href="https://hpo.jax.org/" target="_blank">HPO</a>
</div>""", unsafe_allow_html=True)
contact_email_sidebar = "soap@fastemail.io"
current_date_sidebar = datetime.now().strftime("%Y.%m.%d")
version_info_sidebar = f"Wangzihao_CHD_Explorer_{current_date_sidebar}_v2.9" # Incremented version for search feature
st.sidebar.markdown(f"""<div class="sidebar-footer" style="margin-top:10px;">
<b>Contact:</b> <a href="mailto:{contact_email_sidebar}">{contact_email_sidebar}</a><br>
<b>Last Updated:</b> {current_date_sidebar}<br>
<b>Version:</b> {version_info_sidebar}<br>
<i>Continuously Evolving...</i>
</div>""", unsafe_allow_html=True)


# Load HPO data (can be done after sidebar setup if not needed for sidebar itself)
df_hpo_phenotypes, df_hpo_diseases, hpo_load_error = load_hpo_data(selected_gene)
if hpo_load_error and "file not found" not in hpo_load_error.lower():
    st.warning(hpo_load_error) # Show HPO load error in main area if critical

# Define main tabs
main_tabs_list = ["Gene Explorer (ClinVar)", "CHD Family Comparison (ClinVar)", "HPO Phenotype Explorer"]
main_tab1, main_tab2, main_tab3 = st.tabs(main_tabs_list)

# --- Main Tab 1: Gene Explorer (ClinVar) ---
with main_tab1:
    st.title(f"{selected_gene} ClinVar Variant Explorer")
    if df_gene_variants.empty and not gene_domains :
        st.warning(f"No ClinVar variant or domain data could be loaded for {selected_gene}.")
    else:
        # Prepare df_positions_filtered for positional plots using sidebar filters
        df_positions_filtered = df_gene_variants.dropna(subset=['Position']) if 'Position' in df_gene_variants.columns else pd.DataFrame()
        if not df_positions_filtered.empty:
            # Make sure Position is numeric for filtering
            df_positions_filtered['Position'] = pd.to_numeric(df_positions_filtered['Position'], errors='coerce')
            df_positions_filtered = df_positions_filtered.dropna(subset=['Position'])
            df_positions_filtered['Position'] = df_positions_filtered['Position'].astype(int)

            df_positions_filtered = df_positions_filtered[
                (df_positions_filtered['Position'] >= pos_range_val[0]) &
                (df_positions_filtered['Position'] <= pos_range_val[1])
            ]
            if selected_significance_filter:
                df_positions_filtered = df_positions_filtered[df_positions_filtered['SignificanceClass'].isin(selected_significance_filter)]
            else: # If no significance selected for positional, show no positional variants
                df_positions_filtered = pd.DataFrame(columns=df_positions_filtered.columns)

            if not df_positions_filtered.empty: # If still data after significance filter
                if selected_types_filter:
                    df_positions_filtered = df_positions_filtered[df_positions_filtered['StandardType'].isin(selected_types_filter)]
                else: # If no type selected for positional, show no positional variants
                    df_positions_filtered = pd.DataFrame(columns=df_positions_filtered.columns)

        # Define tabs for Gene Explorer (ClinVar)
        ge_tab_pos, ge_tab_type_sig, ge_tab_pheno, ge_tab_origin = st.tabs([
            f"Positional Plots & Search ({len(df_positions_filtered):,} variants in plots)", # Updated tab title
            "Type & Significance",
            "ClinVar Phenotypes",
            "Variant Origin"
        ])

        with ge_tab_pos:
            st.subheader(f"Search Variant by Protein Position in {selected_gene}")

            # Make sure the original df has numeric Position
            if 'Position' in df_gene_variants.columns:
                 df_gene_variants['Position'] = pd.to_numeric(df_gene_variants['Position'], errors='coerce')
                 # Keep original df intact, but filter out NaNs for searching
                 searchable_df = df_gene_variants.dropna(subset=['Position']).copy()
                 searchable_df['Position'] = searchable_df['Position'].astype(int)
            else:
                searchable_df = pd.DataFrame() # No position data to search

            search_col1, search_col2 = st.columns([3, 1])
            with search_col1:
                search_pos = search_col1.number_input(
                    "Enter exact protein position (aa):",
                    min_value=1,
                    max_value=gene_length if gene_length > 0 else 10000, # Use gene_length as max
                    step=1,
                    value=None, # Default to None, no value pre-filled
                    placeholder="e.g., 123",
                    key=f"search_pos_input_{selected_gene}",
                    label_visibility="collapsed" # Hide label, use placeholder
                )
            with search_col2:
                search_button = st.button(f"🔍 Search Position", key=f"search_pos_button_{selected_gene}", use_container_width=True)

            # --- Display Search Results ---
            if search_button and search_pos is not None:
                st.markdown("---") # Separator
                st.subheader(f"Search Results for Position {search_pos}")
                if not searchable_df.empty:
                    # Filter the *original* full gene df (before positional plot filters)
                    results_df = searchable_df[searchable_df['Position'] == search_pos]

                    if not results_df.empty:
                        st.success(f"Found {len(results_df)} variant(s) at position {search_pos}:")
                        for index, row in results_df.iterrows():
                            # Nicer display using expander and markdown
                            with st.expander(f"**{row.get('Name', 'N/A')}** (ClinVar ID: {row.get('VariationID', 'N/A')})"):
                                st.markdown(f"""
                                *   **Type:** {row.get('StandardType', 'N/A')}
                                *   **Clinical Significance:** {row.get('SignificanceClass', 'N/A')}
                                *   **Reported Origin:** {row.get('OriginSimple', 'N/A')}
                                *   **Phenotypes (ClinVar):** {row.get('PhenotypeListClean', 'N/A')}
                                """)
                                # Add ClinVar link if VariationID exists
                                if pd.notna(row.get('VariationID')):
                                     st.link_button("View on ClinVar", f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{row.get('VariationID')}/")
                    else:
                        st.info(f"No ClinVar variants found exactly at protein position {search_pos} for {selected_gene}.")
                else:
                     st.warning(f"No positional data available in the loaded ClinVar file for {selected_gene} to perform search.")
                st.markdown("---") # Separator after results
            elif search_button and search_pos is None:
                 st.warning("Please enter a protein position to search.")

            # --- Display Positional Plots ---
            st.subheader(f"Positional Variant Plots for {selected_gene}")
            st.markdown(f"_Displaying variants matching sidebar filters (Position Range, Significance, Type). Total matching: {len(df_positions_filtered)}._")

            if df_positions_filtered.empty:
                st.info(f"No variants with protein positions match current sidebar filters for {selected_gene}.")
            else:
                # Check if search results exist and if the searched position is within the current plot range
                if search_button and search_pos is not None and not results_df.empty:
                     if pos_range_val[0] <= search_pos <= pos_range_val[1]:
                           st.markdown(f"_*Note: Variants found at searched position {search_pos} are highlighted or visible in the plots below if they match the selected filters._")
                     else:
                           st.markdown(f"_*Note: Searched position {search_pos} is outside the current plot range [{pos_range_val[0]}, {pos_range_val[1]}]._")

                st.plotly_chart(create_interactive_lollipop(df_positions_filtered, gene_domains, gene_length, selected_gene), use_container_width=True)
                st.plotly_chart(create_interactive_waterfall(df_positions_filtered, gene_domains, gene_length, selected_gene), use_container_width=True)
                st.plotly_chart(create_interactive_density(df_positions_filtered, gene_domains, gene_length, selected_gene), use_container_width=True)

        # --- Other Gene Explorer Tabs (Type/Sig, Pheno, Origin) ---
        with ge_tab_type_sig:
            # ... (content remains the same) ...
            st.subheader(f"Variant Type vs. Clinical Significance for {selected_gene} (All ClinVar Variants for this Gene)")
            st.plotly_chart(create_type_significance_stacked_bar(df_gene_variants, selected_gene), use_container_width=True)

        with ge_tab_pheno:
            # ... (content remains the same) ...
            st.subheader(f"ClinVar Phenotype Analysis for {selected_gene} (Top {TOP_N_PHENOTYPES_GENE_EXPLORER} Phenotypes, All ClinVar Variants for this Gene)")
            st.plotly_chart(create_phenotype_type_stacked_bar(df_gene_variants, selected_gene, normalize=False, top_n=TOP_N_PHENOTYPES_GENE_EXPLORER), use_container_width=True)
            st.markdown("---")
            st.plotly_chart(create_phenotype_type_stacked_bar(df_gene_variants, selected_gene, normalize=True, top_n=TOP_N_PHENOTYPES_GENE_EXPLORER), use_container_width=True)

        with ge_tab_origin:
            # ... (content remains the same) ...
            st.subheader(f"Variant Origin Distribution for {selected_gene} (All ClinVar Variants for this Gene)")
            st.plotly_chart(create_origin_pie_chart(df_gene_variants, selected_gene), use_container_width=True)

        # --- Gene Summary & External Links ---
        st.markdown("---");
        st.markdown(f"Gene Length Used for Plotting {selected_gene}: **{gene_length} aa**. Domains: **{', '.join([d['name'] for d in gene_domains]) if gene_domains else 'N/A'}**.")
        st.markdown("---"); st.subheader(f"Explore {selected_gene} on External Resources")
        # ... (external links remain the same) ...
        link_defs = [{"name": "PubMed", "url": f"https://pubmed.ncbi.nlm.nih.gov/?term={selected_gene}", "icon": "🔬"},{"name": "OMIM", "url": f"https://www.omim.org/search?index=entry&search={selected_gene}", "icon": "🧬"},{"name": "DECIPHER", "url": f"https://www.deciphergenomics.org/gene/{selected_gene}/overview/protein-genomic", "icon": "💡"},{"name": "gnomAD", "url": f"https://gnomad.broadinstitute.org/search/{selected_gene}", "icon": "📊"},{"name": "Protein Atlas", "url": f"https://www.proteinatlas.org/search/{selected_gene}", "icon": "🖼️"}]
        cols = st.columns(len(link_defs));
        for i, link in enumerate(link_defs): cols[i].link_button(f"{link['icon']} {link['name']}", link['url'], use_container_width=True)


# --- Main Tab 2: CHD Family Comparison (ClinVar) ---
with main_tab2:
    # ... (Content of main_tab2 remains unchanged) ...
    st.title("CHD Gene Family Comparison (ClinVar Data)")
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
        fig_patho_counts = px.bar(comparison_summary_df.sort_values("Pathogenic/Likely Path. Variants", ascending=False), x="Gene", y="Pathogenic/Likely Path. Variants",title="Pathogenic/Likely Pathogenic Variants per CHD Gene",labels={"Pathogenic/Likely Path. Variants": "Number of Patho/LP Variants"},height=COMPARISON_PLOT_HEIGHT,color_discrete_sequence=["#007bff"])
        fig_patho_counts.update_layout(font_family=PLOT_FONT_FAMILY); st.plotly_chart(fig_patho_counts, use_container_width=True); st.markdown("<br>", unsafe_allow_html=True)
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
    st.markdown("""*   **Comparative Domain Architecture Plot.**\n*   **Phenotype Similarity/Clustering.**\n*   **Variant Hotspot Comparison.**\n*   **Gene Expression Correlation (External Data).**\n*   **Temporal Analysis**.\n*   **Structural Variant Impact.**\n*   **Functional Impact Scores Comparison.**\n*   **Conservation Analysis.**""")


# --- Main Tab 3: HPO Phenotype Explorer ---
with main_tab3:
    # ... (Content of main_tab3 remains unchanged) ...
    st.title(f"{selected_gene} HPO Phenotype Explorer")
    if df_hpo_phenotypes.empty and df_hpo_diseases.empty and hpo_load_error:
         st.warning(hpo_load_error)
    elif df_hpo_phenotypes.empty and df_hpo_diseases.empty:
        st.warning(f"No HPO phenotype or disease data could be loaded or found for {selected_gene}. Please ensure '{selected_gene}_annotations.json' exists in '{HPO_FILE_DIR}'.")
    else:
        st.info(f"Displaying HPO (Human Phenotype Ontology) annotations for {selected_gene}.")
        st.markdown("---")
        hpo_tab1, hpo_tab2, hpo_tab3, hpo_tab4, hpo_tab5 = st.tabs(["Phenotype Overview", "Disease Associations", "Phenotype Categories","Phenotype Wordcloud","Gene-Phenotype Network"])
        with hpo_tab1:
            st.subheader("Phenotype Counts and Frequency")
            col1, col2 = st.columns(2)
            with col1: st.plotly_chart(plot_hpo_phenotype_counts_bar(df_hpo_phenotypes, selected_gene), use_container_width=True)
            hpo_top_n_slider = st.slider(f"Select Top N HPO Phenotypes to display:", 5, 50, 15, 1, key="hpo_top_n_slider")
            st.plotly_chart(plot_hpo_phenotype_frequency_bar(df_hpo_phenotypes, selected_gene, top_n=hpo_top_n_slider), use_container_width=True)
        with hpo_tab2:
            st.subheader("Associated Diseases (from HPO annotations)")
            st.plotly_chart(plot_hpo_disease_associations_bar(df_hpo_diseases, selected_gene), use_container_width=True)
            if not df_hpo_diseases.empty:
                st.markdown("##### Disease Details Table")
                st.dataframe(df_hpo_diseases[['name', 'id', 'mondoId']].rename(columns={'name':'Disease Name', 'id':'HPO Disease ID'}), use_container_width=True, hide_index=True)
            else: st.info("No specific disease associations found in HPO data for this gene.")
        with hpo_tab3:
            st.subheader("Phenotype Categories (based on HPO terms)")
            st.plotly_chart(plot_hpo_phenotype_categories_stacked_bar(df_hpo_phenotypes, selected_gene), use_container_width=True)
        with hpo_tab4:
            st.subheader("Phenotype Word Cloud")
            if not df_hpo_phenotypes.empty:
                img_path = create_hpo_wordcloud_image(df_hpo_phenotypes, selected_gene)
                if img_path and os.path.exists(img_path):
                    st.image(img_path, caption=f"Word Cloud for {selected_gene} HPO Phenotypes", use_column_width=True)
                elif img_path is None and not df_hpo_phenotypes.empty : st.warning("Could not generate word cloud.")
                else: st.info("No HPO phenotype names available to generate a word cloud.")
            else: st.info("No HPO phenotype data to generate a word cloud.")
        with hpo_tab5:
            st.subheader("Gene - HPO Phenotype - Disease Network (Simplified)")
            net_top_pheno = st.slider("Top N Phenotypes for Network:", 3, 15, 7, key="net_top_pheno")
            net_top_dis = st.slider("Top N Diseases for Network:", 1, 10, 3, key="net_top_dis")
            st.plotly_chart(plot_hpo_gene_phenotype_network_plotly(df_hpo_phenotypes, df_hpo_diseases, selected_gene, top_n_pheno=net_top_pheno, top_n_dis=net_top_dis), use_container_width=True)
