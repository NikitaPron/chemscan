import streamlit as st
import pandas as pd
import numpy as np
import base64
from io import BytesIO
import math

st.set_page_config(
    page_title="ChemScan - Drug Screening Platform",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

LIPINSKI_LIMITS = {
    'MW': 500.0,
    'logP (pred.)': 5.0,
    'HBD': 5,
    'HBA': 10,
    'TPSA': 140.0,
}

def parse_bool(value):
    if pd.isna(value):
        return None
    return str(value).lower() in ['true', '1', 'yes', '+']

@st.cache_data
def load_data():
    try:
        df = pd.read_csv("table_drugs.csv", decimal=',')
        # Updated column names
        numeric_columns = ['pKa (basic)', 'pKa (acidic)', 'pKa (exp.)_1', 'pKa (exp.)_2', 
                          'logP (pred.)', 'logP (exp.)', 'DyeLeS score', 'First approval (year)', 'Oral',
                          'MW', 'TPSA', 'HBD', 'HBA', 'lipinski_violations']
        
        for col in numeric_columns:
            if col in df.columns:
                df[col] = df[col].astype(str).str.replace(',', '.').str.replace(' ', '')
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Renaming Fluorescent column for convenience
        if 'Fluorescent' in df.columns:
            df['is_fluorescent'] = df['Fluorescent'].apply(parse_bool)

        if 'lipinski_pass' in df.columns:
            df['lipinski_pass'] = df['lipinski_pass'].apply(parse_bool)
        
        return df
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return pd.DataFrame()

def apply_optional_range_filter(df, column, enabled, value_min, value_max):
    if not enabled or column not in df.columns:
        return df
    mask = df[column].isna() | df[column].between(value_min, value_max)
    return df[mask]

def get_download_link(df, filename="chemscan_results.csv"):
    csv = df.to_csv(index=False, encoding='utf-8-sig')
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">📥 Download results (CSV)</a>'
    return href

def create_histogram(df, column, title, bins=20):
    import altair as alt
    data = df[[column]].dropna().rename(columns={column: 'value'})
    if len(data) == 0:
        return None
    
    chart = alt.Chart(data).mark_bar().encode(
        x=alt.X('value:Q', bin=alt.Bin(maxbins=bins), title=title),
        y=alt.Y('count()', title='Count'),
        tooltip=[
            alt.Tooltip('value:Q', title='Value', format='.2f'),
            alt.Tooltip('count()', title='Count'),
        ],
    ).properties(
        height=300,
        title=title,
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=14,
    ).configure_title(
        fontSize=16,
        anchor='middle',
    )
    
    return chart

def create_histogram_from_values(values, title, bins=20):
    import altair as alt
    data = pd.DataFrame({'value': values})
    if len(data) == 0:
        return None
    return alt.Chart(data).mark_bar().encode(
        x=alt.X('value:Q', bin=alt.Bin(maxbins=bins), title=title),
        y=alt.Y('count()', title='Count'),
        tooltip=[
            alt.Tooltip('value:Q', title='Value', format='.2f'),
            alt.Tooltip('count()', title='Count'),
        ],
    ).properties(
        height=300,
        title=title,
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=14,
    ).configure_title(
        fontSize=16,
        anchor='middle',
    )

def render_distribution(df, column, title):
    if column not in df.columns or not df[column].notna().any():
        st.info(f"No data for {title}")
        return
    st.write(f"**{title}**")
    hist_chart = create_histogram(df, column, title)
    if hist_chart:
        st.altair_chart(hist_chart, use_container_width=True)
    series = df[column].dropna()
    st.caption(f"Mean: {series.mean():.2f}, Median: {series.median():.2f}, n={len(series)}")

def render_combined_pka_exp(df, columns, title):
    values = []
    for col in columns:
        if col in df.columns:
            values.extend(df[col].dropna().tolist())
    if not values:
        st.info(f"No data for {title}")
        return
    st.write(f"**{title}**")
    hist_chart = create_histogram_from_values(values, title)
    if hist_chart:
        st.altair_chart(hist_chart, use_container_width=True)
    series = pd.Series(values)
    st.caption(f"Mean: {series.mean():.2f}, Median: {series.median():.2f}, n={len(series)}")

def get_molecule_image(smiles, size=(400, 400)):
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        
        if pd.isna(smiles) or smiles == '' or str(smiles) == 'nan':
            return None
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        
        img = Draw.MolToImage(mol, size=size, kekulize=True)
        
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        return f'<img src="data:image/png;base64,{img_str}" style="width: 100%; max-width: 400px; height: auto; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">'
    except Exception as e:
        return None

@st.cache_data
def extract_atc_levels(code):
    if pd.isna(code):
        return None, None, None
    code = str(code).strip().upper().replace(' ', '')
    l1 = code[:1] if len(code) >= 1 and code[0].isalpha() else None
    l2 = code[:3] if len(code) >= 3 and code[1:3].isdigit() else None
    l3 = code[:4] if len(code) >= 4 and code[3].isalpha() else None
    return l1, l2, l3

def main():
    st.title("🔬 ChemScan")
    st.subheader("Drug Screening Platform for Subcellular Accumulation and Targeted Delivery")
    
    with st.expander("ℹ️ About the Platform"):
        st.markdown("""
        **ChemScan** is a computational platform for screening approved drugs based on their suitability for subcellular accumulation.
        
        ### Key Screening Criteria:
        - **Basic pKa**: 8-12 (recommended for platelet accumulation)
        - **Acid pKa**: > 8 (to avoid hindering cell membrane penetration)
        - **LogP**: 1-6 (facilitates penetration)
        
        ### Additional Features:
        - ATC classification filtering (3 levels)
        - Fluorescence property assessment (DyeLeS score)
        - Molecular structure visualization
        - Results export
        """)
    
    df = load_data()
    
    if df.empty:
        st.warning("⚠️ Data not loaded. Please check that the file table_drugs.csv exists.")
        return
    
    atc_code_col = None
    for col in ['ATC code(s)', 'ATC Code', 'atc_code']:
        if col in df.columns:
            atc_code_col = col
            break
    
    if atc_code_col and df[atc_code_col].notna().any():
        df['_atc_l1'] = df[atc_code_col].apply(lambda x: extract_atc_levels(x)[0])
        df['_atc_l2'] = df[atc_code_col].apply(lambda x: extract_atc_levels(x)[1])
        df['_atc_l3'] = df[atc_code_col].apply(lambda x: extract_atc_levels(x)[2])
    
    st.sidebar.header("🔍 Filters")
    
    st.sidebar.subheader("📊 pKa and LogP Value Type")
    value_type = st.sidebar.radio(
        "Select value type:",
        ["Predicted", "Experimental"],
        help="Predicted values obtained via Chemaxon, experimental values from literature"
    )
    
    is_experimental = value_type == "Experimental"
    if is_experimental:
        logp_col = 'logP (exp.)'
        exp_pka_cols = [col for col in ['pKa (exp.)_1', 'pKa (exp.)_2'] if col in df.columns]
    else:
        bpka_col = 'pKa (basic)'
        apka_col = 'pKa (acidic)'
        logp_col = 'logP (pred.)'
        exp_pka_cols = []

    display_atc_col = atc_code_col if atc_code_col else 'ATC therapeutic group'
    
    st.sidebar.subheader("⚗️ pKa Parameters")
    
    enable_bpka_filter = False
    enable_apka_filter = False
    enable_exp_pka_filter = False
    basic_pka_min, basic_pka_max = 0.0, 15.0
    acid_pka_min, acid_pka_max = 0.0, 35.0
    exp_pka_min, exp_pka_max = 0.0, 15.0

    if is_experimental:
        exp_pka_values = pd.concat(
            [df[col] for col in exp_pka_cols],
            ignore_index=True
        ).dropna()
        if len(exp_pka_values) > 0:
            enable_exp_pka_filter = st.sidebar.checkbox("Enable pKa (exp.) filter", value=True)
            min_exp_pka = float(exp_pka_values.min())
            max_exp_pka = float(exp_pka_values.max())
            if enable_exp_pka_filter:
                exp_pka_min, exp_pka_max = st.sidebar.slider(
                    "pKa (exp.) range",
                    min_value=min_exp_pka,
                    max_value=max_exp_pka,
                    value=(8.0, 12.0),
                    help="Keeps drugs where any experimental pKa value falls in this range"
                )
            else:
                exp_pka_min, exp_pka_max = min_exp_pka, max_exp_pka
        else:
            st.sidebar.warning("No experimental pKa data available")
    else:
        enable_bpka_filter = st.sidebar.checkbox("Enable Basic pKa filter", value=True)
        if bpka_col in df.columns and not df[bpka_col].isna().all():
            min_bpka = float(df[bpka_col].min())
            max_bpka = float(df[bpka_col].max())
            
            if enable_bpka_filter:
                basic_pka_min, basic_pka_max = st.sidebar.slider(
                    "Basic pKa range",
                    min_value=min_bpka,
                    max_value=max_bpka,
                    value=(8.0, 12.0),
                    help="Recommended range: 8-12"
                )
            else:
                basic_pka_min, basic_pka_max = min_bpka, max_bpka
        else:
            st.sidebar.warning("No Basic pKa data available")
            basic_pka_min, basic_pka_max = 0.0, 15.0
            enable_bpka_filter = False
            
        enable_apka_filter = st.sidebar.checkbox("Enable Acid pKa filter", value=True)
        if apka_col in df.columns and not df[apka_col].isna().all():
            min_apka = float(df[apka_col].min())
            max_apka = float(df[apka_col].max())
            
            if enable_apka_filter:
                acid_pka_min, acid_pka_max = st.sidebar.slider(
                    "Acid pKa range",
                    min_value=min_apka,
                    max_value=max_apka,
                    value=(8.0, 35.0),
                    help="Recommended range: > 8"
                )
            else:
                acid_pka_min, acid_pka_max = min_apka, max_apka
        else:
            st.sidebar.warning("No Acid pKa data available")
            acid_pka_min, acid_pka_max = 0.0, 35.0
            enable_apka_filter = False
    
    st.sidebar.subheader("🧪 LogP")
    enable_logp_filter = st.sidebar.checkbox("Enable LogP filter", value=True)
    if logp_col in df.columns and not df[logp_col].isna().all():
        min_logp = float(df[logp_col].min())
        max_logp = float(df[logp_col].max())
        
        if enable_logp_filter:
            logp_min, logp_max = st.sidebar.slider(
                "LogP range",
                min_value=min_logp,
                max_value=max_logp,
                value=(1.0, 6.0),
                help="Recommended range: 1-6"
            )
        else:
            logp_min, logp_max = min_logp, max_logp
    else:
        st.sidebar.warning("No LogP data available")
        logp_min, logp_max = -5.0, 10.0
        enable_logp_filter = False
    
    st.sidebar.subheader("💡 Fluorescence")
    fluorescence_filter = st.sidebar.selectbox(
        "Fluorescence filter",
        ["All", "Fluorescent only", "Non-fluorescent only"],
        help="Filter by DyeLeS score (>0.5 = fluorescent)"
    )
    
    st.sidebar.subheader("📋 ATC Classification")
    
    selected_l1, selected_l2, selected_l3 = [], [], []
    
    if atc_code_col and df[atc_code_col].notna().any():
        l1_options = sorted([x for x in df['_atc_l1'].dropna().unique() if pd.notna(x)])
        selected_l1 = st.sidebar.multiselect(
            "🔸 Level 1: Anatomical group",
            options=l1_options,
            default=[],
            help="Example: A — Alimentary tract, C — Cardiovascular system, N — Nervous system"
        )
        
        if selected_l1:
            l2_candidates = df[df['_atc_l1'].isin(selected_l1)]['_atc_l2'].dropna().unique()
            l2_options = sorted([x for x in l2_candidates if pd.notna(x)])
            selected_l2 = st.sidebar.multiselect(
                "🔹 Level 2: Therapeutic subgroup",
                options=l2_options,
                default=[],
                help="Example: A02 — Drugs for acid-related disorders"
            )
        else:
            selected_l2 = []
        
        if selected_l2:
            l3_candidates = df[df['_atc_l2'].isin(selected_l2)]['_atc_l3'].dropna().unique()
            l3_options = sorted([x for x in l3_candidates if pd.notna(x)])
            selected_l3 = st.sidebar.multiselect(
                "▫️ Level 3: Pharmacological subgroup",
                options=l3_options,
                default=[],
                help="Example: A02BC — Proton pump inhibitors"
            )
        else:
            selected_l3 = []
            if selected_l1:
                st.sidebar.info("↪️ Select a therapeutic subgroup above")
    else:
        st.sidebar.warning("⚠️ ATC code column not found. ATC filtering disabled.")
    
    st.sidebar.subheader("💊 Administration")
    oral_only = st.sidebar.checkbox("Oral drugs only", value=False)

    lipinski_col = 'logP (pred.)'
    lipinski_filters = {
        'MW': {'enabled': False, 'min': 0.0, 'max': 5000.0},
        lipinski_col: {'enabled': False, 'min': -10.0, 'max': 10.0},
        'HBD': {'enabled': False, 'min': 0, 'max': 50},
        'HBA': {'enabled': False, 'min': 0, 'max': 50},
        'TPSA': {'enabled': False, 'min': 0.0, 'max': 2000.0},
    }

    with st.sidebar.expander("📐 Lipinski-like (optional)", expanded=False):
        st.caption(
            "Rule of Five: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10, TPSA ≤ 140 Å². "
            "Sliders default to these limits when enabled."
        )

        lipinski_specs = [
            ('MW', 'MW range (Da)', 'Molecular weight', 1.0),
            (lipinski_col, 'LogP range (predicted)', 'Lipinski uses predicted LogP', 0.1),
            ('HBD', 'H-bond donors', 'NumHDonors (RDKit)', 1.0),
            ('HBA', 'H-bond acceptors', 'NumHAcceptors (RDKit)', 1.0),
            ('TPSA', 'TPSA range (Å²)', 'Topological polar surface area', 1.0),
        ]

        for col, label, help_text, step in lipinski_specs:
            if col not in df.columns or df[col].isna().all():
                st.warning(f"No {col} data available")
                continue

            col_min = float(df[col].min())
            col_max = float(df[col].max())
            lipinski_max = float(LIPINSKI_LIMITS[col])
            default_max = min(lipinski_max, col_max)
            default_min = col_min

            lipinski_filters[col]['enabled'] = st.checkbox(
                f"Enable {label.split('(')[0].strip()} filter",
                value=False,
                key=f"enable_lipinski_{col}",
            )

            if lipinski_filters[col]['enabled']:
                if col in ('HBD', 'HBA'):
                    slider_min = int(col_min)
                    slider_max = int(col_max)
                    default_min_i = int(default_min)
                    default_max_i = int(default_max)
                    selected = st.slider(
                        label,
                        min_value=slider_min,
                        max_value=slider_max,
                        value=(default_min_i, default_max_i),
                        step=int(step),
                        help=help_text,
                        key=f"lipinski_slider_{col}",
                    )
                    lipinski_filters[col]['min'] = float(selected[0])
                    lipinski_filters[col]['max'] = float(selected[1])
                else:
                    selected = st.slider(
                        label,
                        min_value=col_min,
                        max_value=col_max,
                        value=(default_min, default_max),
                        step=step,
                        help=help_text,
                        key=f"lipinski_slider_{col}",
                    )
                    lipinski_filters[col]['min'] = selected[0]
                    lipinski_filters[col]['max'] = selected[1]
            else:
                lipinski_filters[col]['min'] = col_min
                lipinski_filters[col]['max'] = col_max
    
    filtered_df = df.copy()
    
    if is_experimental and enable_exp_pka_filter and exp_pka_cols:
        has_exp_pka = pd.Series(False, index=filtered_df.index)
        in_exp_pka_range = pd.Series(False, index=filtered_df.index)
        for col in exp_pka_cols:
            has_exp_pka |= filtered_df[col].notna()
            in_exp_pka_range |= filtered_df[col].between(exp_pka_min, exp_pka_max)
        filtered_df = filtered_df[~has_exp_pka | in_exp_pka_range]
    else:
        if enable_bpka_filter and bpka_col in filtered_df.columns:
            mask_bpka = filtered_df[bpka_col].isna() | filtered_df[bpka_col].between(basic_pka_min, basic_pka_max)
            filtered_df = filtered_df[mask_bpka]
        
        if enable_apka_filter and apka_col in filtered_df.columns:
            mask_apka = filtered_df[apka_col].isna() | filtered_df[apka_col].between(acid_pka_min, acid_pka_max)
            filtered_df = filtered_df[mask_apka]
    
    if enable_logp_filter and logp_col in filtered_df.columns:
        mask_logp = filtered_df[logp_col].isna() | filtered_df[logp_col].between(logp_min, logp_max)
        filtered_df = filtered_df[mask_logp]
    
    if 'is_fluorescent' in filtered_df.columns and fluorescence_filter != "All":
        if fluorescence_filter == "Fluorescent only":
            filtered_df = filtered_df[filtered_df['is_fluorescent'] == True]
        elif fluorescence_filter == "Non-fluorescent only":
            filtered_df = filtered_df[filtered_df['is_fluorescent'] == False]
    
    if atc_code_col and atc_code_col in filtered_df.columns:
        mask_atc = pd.Series([True] * len(filtered_df), index=filtered_df.index)
        
        if selected_l3:
            prefixes = tuple(str(p).strip() for p in selected_l3 if pd.notna(p))
            if prefixes:
                mask_atc = filtered_df[atc_code_col].astype(str).str.strip().str.upper().str.startswith(prefixes)
        
        elif selected_l2:
            prefixes = tuple(str(p).strip() for p in selected_l2 if pd.notna(p))
            if prefixes:
                mask_atc = filtered_df[atc_code_col].astype(str).str.strip().str.upper().str.startswith(prefixes)
        
        elif selected_l1:
            prefixes = tuple(str(p).strip() for p in selected_l1 if pd.notna(p))
            if prefixes:
                mask_atc = filtered_df[atc_code_col].astype(str).str.strip().str.upper().str.startswith(prefixes)
        
        filtered_df = filtered_df[mask_atc]
    
    if 'Oral' in filtered_df.columns and oral_only:
        filtered_df = filtered_df[filtered_df['Oral'] == 1.0]

    for col, settings in lipinski_filters.items():
        filtered_df = apply_optional_range_filter(
            filtered_df,
            col,
            settings['enabled'],
            settings['min'],
            settings['max'],
        )
    
    exp_signal_cols = [
        col for col in ['pKa_comment']
        if col in filtered_df.columns
    ]
    if exp_signal_cols:
        filtered_df['_has_exp_data'] = False
        for col in exp_signal_cols:
            filtered_df['_has_exp_data'] |= filtered_df[col].notna()
    else:
        filtered_df['_has_exp_data'] = False

    sort_cols = ['_has_exp_data']
    if is_experimental and exp_pka_cols:
        filtered_df['_has_pka'] = False
        for col in exp_pka_cols:
            filtered_df['_has_pka'] |= filtered_df[col].notna()
        sort_cols.append('_has_pka')
    elif (not is_experimental) and ('pKa (basic)' in filtered_df.columns):
        filtered_df['_has_bpka'] = filtered_df['pKa (basic)'].notna()
        sort_cols.append('_has_bpka')

    filtered_df = filtered_df.sort_values(by=sort_cols, ascending=[False] * len(sort_cols), kind='mergesort')
    filtered_df = filtered_df.drop(columns=[c for c in ['_has_exp_data', '_has_pka', '_has_bpka'] if c in filtered_df.columns])
    
    st.header("📊 Search Results")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total in database", len(df))
    with col2:
        st.metric("Found", len(filtered_df))
    with col3:
        percentage = (len(filtered_df) / len(df) * 100) if len(df) > 0 else 0
        if 0 < percentage < 0.1:
            percentage_display = 0.1
        else:
            percentage_display = round(percentage, 1)
        st.metric("Percentage", f"{percentage_display:.1f}%")
    
    tab1, tab2, tab3 = st.tabs(["💊 Molecules", "📈 Statistics", "📋 Table"])
    
    with tab1:
        st.subheader("Molecular Cards")
        
        if len(filtered_df) > 0:
            display_count = min(20, len(filtered_df))
            
            for idx in range(display_count):
                row = filtered_df.iloc[idx]
                
                # Using container with fixed structure
                with st.container():
                    # Creating two columns with 1:2 ratio, responsive design
                    col1, col2 = st.columns([1, 2], gap="medium")
                    
                    with col1:
                        if 'SMILES' in row and pd.notna(row['SMILES']):
                            mol_img = get_molecule_image(row['SMILES'])
                            if mol_img:
                                st.markdown(mol_img, unsafe_allow_html=True)
                            else:
                                st.info("Structure not available")
                        else:
                            st.info("SMILES not available")
                    
                    with col2:
                        # Adding top margin for better alignment
                        st.markdown(f"**{row.get('Name', 'N/A')}**")
                        atc_display = row.get(atc_code_col, 'N/A') if atc_code_col else row.get('ATC therapeutic group', 'N/A')
                        st.caption(f"ATC: {atc_display}")
                        
                        param_col1, param_col2 = st.columns(2)
                        
                        with param_col1:
                            st.caption("Predicted (Chemaxon)")
                            if pd.notna(row.get('pKa (basic)')):
                                st.write(f"**Basic pKa:** {float(row.get('pKa (basic)')):.2f}")
                            else:
                                st.markdown('<span style="color: red;">**Basic pKa:** —</span>', unsafe_allow_html=True)
                            if pd.notna(row.get('pKa (acidic)')):
                                st.write(f"**Acid pKa:** {float(row.get('pKa (acidic)')):.2f}")
                            else:
                                st.markdown('<span style="color: red;">**Acid pKa:** —</span>', unsafe_allow_html=True)
                            if pd.notna(row.get('logP (pred.)')):
                                st.write(f"**LogP:** {float(row.get('logP (pred.)')):.2f}")
                            else:
                                st.markdown('<span style="color: red;">**LogP:** —</span>', unsafe_allow_html=True)
                        
                        with param_col2:
                            st.caption("Experimental (literature)")
                            if pd.notna(row.get('pKa (exp.)_1')):
                                st.write(f"**pKa (exp.) 1:** {float(row.get('pKa (exp.)_1')):.2f}")
                            else:
                                st.markdown('<span style="color: red;">**pKa (exp.) 1:** —</span>', unsafe_allow_html=True)
                            if pd.notna(row.get('pKa (exp.)_2')):
                                st.write(f"**pKa (exp.) 2:** {float(row.get('pKa (exp.)_2')):.2f}")
                            else:
                                st.markdown('<span style="color: red;">**pKa (exp.) 2:** —</span>', unsafe_allow_html=True)
                            if pd.notna(row.get('pKa_comment')):
                                st.write(f"**pKa comment:** {row['pKa_comment']}")
                            if pd.notna(row.get('logP (exp.)')):
                                st.write(f"**LogP (exp.):** {float(row.get('logP (exp.)')):.2f}")
                            else:
                                st.markdown('<span style="color: red;">**LogP (exp.):** —</span>', unsafe_allow_html=True)
                            if 'is_fluorescent' in row:
                                st.write(f"**Fluorescent:** {'✅ Yes' if row['is_fluorescent'] else '❌ No'}")
                            if 'Oral' in row:
                                if pd.notna(row.get('Oral')):
                                    st.write(f"**Oral:** {'✅ Yes' if row['Oral'] == 1.0 else '❌ No'}")
                                else:
                                    st.markdown('<span style="color: red;">**Oral:** —</span>', unsafe_allow_html=True)
                        
                        lipinski_parts = []
                        if pd.notna(row.get('MW')):
                            lipinski_parts.append(f"MW: {float(row['MW']):.1f}")
                        if pd.notna(row.get('TPSA')):
                            lipinski_parts.append(f"TPSA: {float(row['TPSA']):.1f}")
                        if pd.notna(row.get('HBD')):
                            lipinski_parts.append(f"HBD: {int(row['HBD'])}")
                        if pd.notna(row.get('HBA')):
                            lipinski_parts.append(f"HBA: {int(row['HBA'])}")
                        if lipinski_parts:
                            st.caption(" · ".join(lipinski_parts))
                        if 'lipinski_violations' in row and pd.notna(row.get('lipinski_violations')):
                            violations = int(row['lipinski_violations'])
                            status = "✅ Lipinski-like" if violations <= 1 else f"⚠️ {violations} Ro5 violations"
                            st.caption(status)

                        if 'First approval (year)' in row and pd.notna(row['First approval (year)']):
                            st.caption(f"📅 Approved: {int(row['First approval (year)'])}")
                        
                        if 'SMILES' in row and pd.notna(row['SMILES']):
                            with st.expander("🔬 SMILES and Identifiers"):
                                st.code(row['SMILES'])
                                if 'ChEMBL ID' in row and pd.notna(row['ChEMBL ID']):
                                    st.write(f"**ChEMBL ID:** {row['ChEMBL ID']}")
                                if 'DrugBank ID' in row and pd.notna(row['DrugBank ID']):
                                    st.write(f"**DrugBank ID:** {row['DrugBank ID']}")
                    
                    st.divider()
            
            if len(filtered_df) > 10:
                st.info(f"Showing first 20 of {len(filtered_df)} molecules. Use the table to view all.")
            
            st.markdown(get_download_link(filtered_df), unsafe_allow_html=True)
        else:
            st.warning("⚠️ No molecules found with the selected filters. Try expanding the ranges.")
    
    with tab2:
        st.subheader("📈 Statistics")
        
        if len(filtered_df) > 0:
            col_basic, col_acid, col_exp = st.columns(3)
            with col_basic:
                render_distribution(filtered_df, 'pKa (basic)', 'Basic pKa (predicted)')
            with col_acid:
                render_distribution(filtered_df, 'pKa (acidic)', 'Acid pKa (predicted)')
            with col_exp:
                render_combined_pka_exp(
                    filtered_df,
                    ['pKa (exp.)_1', 'pKa (exp.)_2'],
                    'pKa (experimental)',
                )

            col_logp_pred, col_logp_exp = st.columns(2)
            with col_logp_pred:
                render_distribution(filtered_df, 'logP (pred.)', 'LogP (predicted)')
            with col_logp_exp:
                render_distribution(filtered_df, 'logP (exp.)', 'LogP (experimental)')

            with st.expander("📐 Lipinski descriptors", expanded=False):
                col_mw, col_logp_ro5, col_tpsa = st.columns(3)
                with col_mw:
                    render_distribution(filtered_df, 'MW', 'Molecular Weight (Da)')
                with col_logp_ro5:
                    render_distribution(filtered_df, 'logP (pred.)', 'LogP (predicted, Lipinski)')
                with col_tpsa:
                    render_distribution(filtered_df, 'TPSA', 'TPSA (Å²)')
                col_hbd, col_hba, col_pass = st.columns(3)
                with col_hbd:
                    render_distribution(filtered_df, 'HBD', 'H-bond donors')
                with col_hba:
                    render_distribution(filtered_df, 'HBA', 'H-bond acceptors')
                with col_pass:
                    if 'lipinski_pass' in filtered_df.columns:
                        st.write("**Lipinski-like (≤1 violation)**")
                        pass_counts = filtered_df['lipinski_pass'].value_counts()
                        pass_df = pd.DataFrame({
                            'Status': ['Pass', 'Fail'],
                            'Count': [pass_counts.get(True, 0), pass_counts.get(False, 0)],
                        })
                        st.bar_chart(pass_df.set_index('Status'), use_container_width=True)

            col_fluo, col_atc = st.columns(2)
            with col_fluo:
                if 'is_fluorescent' in filtered_df.columns:
                    st.write("**Fluorescence Distribution**")
                    fluo_counts = filtered_df['is_fluorescent'].value_counts()
                    fluo_df = pd.DataFrame({
                        'Type': ['Fluorescent', 'Non-fluorescent'],
                        'Count': [fluo_counts.get(True, 0), fluo_counts.get(False, 0)],
                    })
                    st.bar_chart(fluo_df.set_index('Type'), use_container_width=True)
            with col_atc:
                if display_atc_col in filtered_df.columns:
                    st.write("**Top 10 ATC Groups**")
                    atc_counts = filtered_df[display_atc_col].value_counts().head(10)
                    st.bar_chart(atc_counts, use_container_width=True)
            
            st.write("**📊 Summary Statistics**")
            
            if is_experimental:
                numeric_cols = exp_pka_cols + [logp_col, 'DyeLeS score', 'MW', 'TPSA', 'HBD', 'HBA', 'lipinski_violations', 'First approval (year)']
            else:
                numeric_cols = [bpka_col, apka_col, logp_col, 'DyeLeS score', 'MW', 'TPSA', 'HBD', 'HBA', 'lipinski_violations', 'First approval (year)']
            available_numeric_cols = [col for col in numeric_cols if col in filtered_df.columns and not filtered_df[col].isna().all()]
            
            if available_numeric_cols:
                stats_df = filtered_df[available_numeric_cols].describe()
                st.dataframe(
                    stats_df.round(2),
                    use_container_width=True,
                    height=300
                )
        else:
            st.warning("⚠️ No data available for statistics")
    
    with tab3:
        st.subheader("📋 Full Table")
        
        if len(filtered_df) > 0:
            show_smiles = st.checkbox("Show SMILES", value=False)
            show_ids = st.checkbox("Show identifiers", value=False)
            if is_experimental:
                show_pred = st.checkbox("Show predicted values", value=False)
            else:
                show_pred = False
                show_exp = st.checkbox("Show experimental values", value=True)
            show_atc_levels = st.checkbox("Show ATC levels", value=False)
            
            lipinski_table_cols = ['MW', 'HBD', 'HBA', 'TPSA', 'lipinski_violations', 'lipinski_pass']
            if is_experimental:
                table_columns = ['Name', display_atc_col] + exp_pka_cols + ['pKa_comment', logp_col] + lipinski_table_cols + ['DyeLeS score', 'is_fluorescent', 'Oral', 'First approval (year)']
                if show_pred:
                    pred_columns = ['pKa (basic)', 'pKa (acidic)', 'logP (pred.)']
                    table_columns.extend(col for col in pred_columns if col in filtered_df.columns)
            else:
                table_columns = ['Name', display_atc_col, bpka_col, apka_col, logp_col] + lipinski_table_cols + ['DyeLeS score', 'is_fluorescent', 'Oral', 'First approval (year)']
                if show_exp:
                    exp_detail_columns = ['pKa (exp.)_1', 'pKa (exp.)_2', 'pKa_comment', 'logP (exp.)']
                    table_columns.extend(
                        col for col in exp_detail_columns
                        if col in filtered_df.columns and col not in table_columns
                    )
            
            if show_atc_levels and '_atc_l1' in filtered_df.columns:
                table_columns.extend(['_atc_l1', '_atc_l2', '_atc_l3'])
            
            if show_smiles and 'SMILES' in filtered_df.columns:
                table_columns.append('SMILES')
            
            if show_ids:
                id_columns = ['ChEMBL ID', 'DrugBank ID', 'ligand_IPTM']
                table_columns.extend([col for col in id_columns if col in filtered_df.columns])
            
            table_columns = list(dict.fromkeys(
                col for col in table_columns if col in filtered_df.columns
            ))
            
            st.dataframe(
                filtered_df[table_columns],
                use_container_width=True,
                height=500
            )
            
            st.markdown(get_download_link(filtered_df[table_columns]), unsafe_allow_html=True)
        else:
            st.warning("⚠️ No data to display")
    
    st.divider()
    st.markdown("""
    ### 📚 About the Project
    
    The ChemScan platform is designed to support decision-making when selecting drugs not only for platelet delivery but also for a wide range of targeted therapy strategies. The system allows modeling compound accumulation in various acidic subcellular compartments (such as lysosomes, endosomes, or mitochondria) of different cell types using physicochemical parameters (pKa, logP) and ion trapping mechanisms.
    
    **Contact:** nick-pronn@yandex.ru | **Version:** 1.0
    """)

if __name__ == "__main__":
    main()