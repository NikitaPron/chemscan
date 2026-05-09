import streamlit as st
import pandas as pd
import numpy as np
import base64
from io import BytesIO
import math

st.set_page_config(
    page_title="ChemScan - Платформа подбора лекарств",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@st.cache_data
def load_data():
    try:
        df = pd.read_csv("table_drugs.csv", decimal=',')
        numeric_columns = ['acd_most_apka', 'acd_most_bpka', 'pka exp', 'pka exp 2', 
                          'LogP', 'LogP_exp', 'dyeles_score', 'first_approval', 'oral']
        
        for col in numeric_columns:
            if col in df.columns:
                df[col] = df[col].astype(str).str.replace(',', '.').str.replace(' ', '')
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        return df
    except Exception as e:
        st.error(f"Ошибка загрузки данных: {e}")
        return pd.DataFrame()

def get_download_link(df, filename="chemscan_results.csv"):
    csv = df.to_csv(index=False, encoding='utf-8-sig')
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">📥 Скачать результаты (CSV)</a>'
    return href

def create_histogram(df, column, title, bins=20):
    import altair as alt
    data = df[[column]].dropna()
    if len(data) == 0:
        return None
    
    chart = alt.Chart(data).mark_bar().encode(
        x=alt.X(f'{column}:Q', bin=alt.Bin(maxbins=bins), title=title),
        y=alt.Y('count()', title='Количество'),
        tooltip=[alt.Tooltip(f'{column}:Q', title='Значение'),
                alt.Tooltip('count()', title='Количество')]
    ).properties(
        width=400,
        height=300,
        title=title
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=14
    ).configure_title(
        fontSize=16,
        anchor='middle'
    )
    
    return chart

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
        
        return f'<img src="data:image/png;base64,{img_str}" style="max-width: 400px; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">'
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
    st.subheader("Платформа подбора лекарственных средств для субклеточного накопления и таргетной доставки")
    
    with st.expander("ℹ️ О платформе"):
        st.markdown("""
        **ChemScan** - это вычислительная платформа для скрининга одобренных лекарственных средств 
        по критериям пригодности к субклеточному накоплению.
        
        ### Ключевые критерии отбора:
        - **Basic pKa**: 8-12 (рекомендуемый диапазон для тромбоцитарного накопления)
        - **Acid pKa**: > 8 (чтобы не препятствовать проникновению через клеточную мембрану)
        - **LogP**: 1-6 (способствует проникновению)
        
        ### Дополнительные возможности:
        - Фильтрация по АТХ-классификации (3 уровня)
        - Оценка флуоресцентных свойств (DyeLeS score)
        - Визуализация молекулярных структур
        - Экспорт результатов
        """)
    
    df = load_data()
    
    if df.empty:
        st.warning("⚠️ Данные не загружены. Проверьте наличие файла table_drugs.csv")
        return
    
    atc_code_col = None
    for col in ['atc_code', 'ATC Code', 'atc_therapeutic_group']:
        if col in df.columns:
            atc_code_col = col
            break
    
    if atc_code_col and df[atc_code_col].notna().any():
        df['_atc_l1'] = df[atc_code_col].apply(lambda x: extract_atc_levels(x)[0])
        df['_atc_l2'] = df[atc_code_col].apply(lambda x: extract_atc_levels(x)[1])
        df['_atc_l3'] = df[atc_code_col].apply(lambda x: extract_atc_levels(x)[2])
    
    st.sidebar.header("🔍 Фильтры")
    
    st.sidebar.subheader("📊 Тип значений pKa и LogP")
    value_type = st.sidebar.radio(
        "Выберите тип значений:",
        ["Predicted", "Experimental"],
        help="Расчетные значения получены через ACD/Labs, экспериментальные - из литературы"
    )
    
    if value_type == "Predicted":
        bpka_col = 'acd_most_bpka'
        apka_col = 'acd_most_apka'
        logp_col = 'LogP'
    else:
        bpka_col = 'pka exp'
        apka_col = 'pka exp 2'
        logp_col = 'LogP_exp'
    
    st.sidebar.subheader("⚗️ pKa параметры")
    
    enable_bpka_filter = st.sidebar.checkbox("Включить фильтр Basic pKa", value=True)
    if bpka_col in df.columns and not df[bpka_col].isna().all():
        min_bpka = float(df[bpka_col].min())
        max_bpka = float(df[bpka_col].max())
        
        if enable_bpka_filter:
            basic_pka_min, basic_pka_max = st.sidebar.slider(
                "Basic pKa диапазон",
                min_value=min_bpka,
                max_value=max_bpka,
                value=(8.0, 12.0),
                help="Рекомендуемый диапазон: 8-12"
            )
        else:
            basic_pka_min, basic_pka_max = min_bpka, max_bpka
    else:
        st.sidebar.warning("Нет данных по Basic pKa")
        basic_pka_min, basic_pka_max = 0.0, 15.0
        enable_bpka_filter = False
        
    enable_apka_filter = st.sidebar.checkbox("Включить фильтр Acid pKa", value=True)
    if apka_col in df.columns and not df[apka_col].isna().all():
        min_apka = float(df[apka_col].min())
        max_apka = float(df[apka_col].max())
        
        if enable_apka_filter:
            acid_pka_min, acid_pka_max = st.sidebar.slider(
                "Acid pKa диапазон",
                min_value=min_apka,
                max_value=max_apka,
                value=(8.0, 15.0),
                help="Рекомендуемый диапазон: > 8"
            )
        else:
            acid_pka_min, acid_pka_max = min_apka, max_apka
    else:
        st.sidebar.warning("Нет данных по Acid pKa")
        acid_pka_min, acid_pka_max = 0.0, 20.0
        enable_apka_filter = False
    
    st.sidebar.subheader("🧪 LogP")
    enable_logp_filter = st.sidebar.checkbox("Включить фильтр LogP", value=True)
    if logp_col in df.columns and not df[logp_col].isna().all():
        min_logp = float(df[logp_col].min())
        max_logp = float(df[logp_col].max())
        
        if enable_logp_filter:
            logp_min, logp_max = st.sidebar.slider(
                "LogP диапазон",
                min_value=min_logp,
                max_value=max_logp,
                value=(1.0, 6.0),
                help="Рекомендуемый диапазон: 1-6"
            )
        else:
            logp_min, logp_max = min_logp, max_logp
    else:
        st.sidebar.warning("Нет данных по LogP")
        logp_min, logp_max = -5.0, 10.0
        enable_logp_filter = False
    
    st.sidebar.subheader("💡 Флуоресценция")
    fluorescence_filter = st.sidebar.selectbox(
        "Фильтр по флуоресценции",
        ["Все", "Только флуоресцентные", "Только не флуоресцентные"],
        help="Фильтрация по DyeLeS score (>0.5 = флуоресцентный)"
    )
    
    st.sidebar.subheader("📋 АТХ-классификация")
    
    selected_l1, selected_l2, selected_l3 = [], [], []
    
    if atc_code_col and df[atc_code_col].notna().any():
        l1_options = sorted([x for x in df['_atc_l1'].dropna().unique() if pd.notna(x)])
        selected_l1 = st.sidebar.multiselect(
            "🔸 Уровень 1: Анатомическая группа",
            options=l1_options,
            default=[],
            help="Пример: A — ЖКТ, C — Сердечно-сосудистая, N — Нервная система"
        )
        
        if selected_l1:
            l2_candidates = df[df['_atc_l1'].isin(selected_l1)]['_atc_l2'].dropna().unique()
            l2_options = sorted([x for x in l2_candidates if pd.notna(x)])
            selected_l2 = st.sidebar.multiselect(
                "🔹 Уровень 2: Терапевтическая подгруппа",
                options=l2_options,
                default=[],
                help="Пример: A02 — Препараты для лечения кислотозависимых заболеваний"
            )
        else:
            selected_l2 = []
        
        if selected_l2:
            l3_candidates = df[df['_atc_l2'].isin(selected_l2)]['_atc_l3'].dropna().unique()
            l3_options = sorted([x for x in l3_candidates if pd.notna(x)])
            selected_l3 = st.sidebar.multiselect(
                "▫️ Уровень 3: Фармакологическая подгруппа",
                options=l3_options,
                default=[],
                help="Пример: A02BC — Ингибиторы протонной помпы"
            )
        else:
            selected_l3 = []
            if selected_l1:
                st.sidebar.info("↪️ Выберите терапевтическую подгруппу выше")
    else:
        st.sidebar.warning("⚠️ Колонка с АТХ-кодами не найдена. Фильтрация по АТХ отключена.")
    
    st.sidebar.subheader("💊 Применение")
    oral_only = st.sidebar.checkbox("Только пероральные препараты", value=False)
    
    filtered_df = df.copy()
    
    if enable_bpka_filter and bpka_col in filtered_df.columns:
       if value_type == "Experimental" and 'pka exp' in filtered_df.columns and 'pka exp 2' in filtered_df.columns:
           mask_bpka = (
               (filtered_df['pka exp'].between(basic_pka_min, basic_pka_max)) |    
               (filtered_df['pka exp 2'].between(basic_pka_min, basic_pka_max))    
           )
           filtered_df = filtered_df[mask_bpka]
       else:
           mask_bpka = filtered_df[bpka_col].isna() | filtered_df[bpka_col].between(basic_pka_min, basic_pka_max)
           filtered_df = filtered_df[mask_bpka]
    
    if enable_apka_filter and apka_col in filtered_df.columns:
        mask_apka = filtered_df[apka_col].isna() | filtered_df[apka_col].between(acid_pka_min, acid_pka_max)
        filtered_df = filtered_df[mask_apka]
    
    if enable_logp_filter and logp_col in filtered_df.columns:
        mask_logp = filtered_df[logp_col].isna() | filtered_df[logp_col].between(logp_min, logp_max)
        filtered_df = filtered_df[mask_logp]
    
    if 'is_fluorescent' in filtered_df.columns and fluorescence_filter != "Все":
        if fluorescence_filter == "Только флуоресцентные":
            filtered_df = filtered_df[filtered_df['is_fluorescent'] == True]
        elif fluorescence_filter == "Только не флуоресцентные":
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
    
    if 'oral' in filtered_df.columns and oral_only:
        filtered_df = filtered_df[filtered_df['oral'] == 1.0]
    
    if bpka_col in filtered_df.columns:
        filtered_df['_has_bpka'] = filtered_df[bpka_col].notna()
        filtered_df = filtered_df.sort_values(by=['_has_bpka'], ascending=False)
        filtered_df = filtered_df.drop(columns=['_has_bpka'])
    
    st.header("📊 Результаты поиска")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Всего в базе", len(df))
    with col2:
        st.metric("Найдено", len(filtered_df))
    with col3:
        percentage = (len(filtered_df) / len(df) * 100) if len(df) > 0 else 0
        if 0 < percentage < 0.1:
            percentage_display = 0.1
        else:
            percentage_display = round(percentage, 1)
        st.metric("Доля", f"{percentage_display:.1f}%")
    
    tab1, tab2, tab3 = st.tabs(["💊 Молекулы", "📈 Статистика", "📋 Таблица"])
    
    with tab1:
        st.subheader("Карточки молекул")
        
        if len(filtered_df) > 0:
            display_count = min(20, len(filtered_df))
            
            for idx in range(display_count):
                row = filtered_df.iloc[idx]
                
                with st.container():
                    col1, col2 = st.columns([1, 2])
                    
                    with col1:
                        if 'SMILES' in row and pd.notna(row['SMILES']):
                            mol_img = get_molecule_image(row['SMILES'])
                            if mol_img:
                                st.markdown(mol_img, unsafe_allow_html=True)
                            else:
                                st.info("Структура недоступна")
                        else:
                            st.info("SMILES не доступен")
                    
                    with col2:
                        st.markdown(f"**{row.get('Name', 'N/A')}**")
                        atc_display = row.get(atc_code_col, 'N/A') if atc_code_col else row.get('atc_therapeutic_group', 'N/A')
                        st.caption(f"ATX: {atc_display}")
                        
                        if pd.notna(row.get(bpka_col)):
                            st.write(f"**Basic pKa:** {float(row.get(bpka_col)):.2f}")
                        else:
                            st.markdown('<span style="color: red;">**Basic pKa:** Unknown</span>', unsafe_allow_html=True)
                        if pd.notna(row.get(apka_col)):
                            st.write(f"**Acid pKa:** {float(row.get(apka_col)):.2f}")
                        else:
                            st.markdown('<span style="color: red;">**Acid pKa:** Unknown</span>', unsafe_allow_html=True)
                        if pd.notna(row.get(logp_col)):
                            st.write(f"**LogP:** {float(row.get(logp_col)):.2f}")
                        
                        if 'pka exp' in row and pd.notna(row['pka_comment']):
                            st.write(f"**pKa (exp):** {row['pka_comment']}")
                        if 'LogP_exp' in row and pd.notna(row['LogP_exp']):
                            st.write(f"**LogP (exp):** {float(row['LogP_exp']):.2f}")
                        
                        if 'is_fluorescent' in row:
                            st.write(f"**Флуоресцентный:** {'✅ Да' if row['is_fluorescent'] else '❌ Нет'}")
                        
                        if 'SMILES' in row and pd.notna(row['SMILES']):
                            with st.expander("SMILES"):
                                st.code(row['SMILES'])
                    
                    st.divider()
            
            if len(filtered_df) > 10:
                st.info(f"Показано первые 20 из {len(filtered_df)} молекул. Используйте таблицу для просмотра всех.")
            
            st.markdown(get_download_link(filtered_df), unsafe_allow_html=True)
        else:
            st.warning("⚠️ По заданным фильтрам ничего не найдено. Попробуйте расширить диапазоны.")
    
    with tab2:
        st.subheader("📈 Статистика")
        
        if len(filtered_df) > 0:
            col1, col2 = st.columns(2)
            
            with col1:
                if bpka_col in filtered_df.columns and not filtered_df[bpka_col].isna().all():
                    st.write(f"**Распределение Basic pKa**")
                    hist_chart = create_histogram(filtered_df, bpka_col, "Basic pKa")
                    if hist_chart:
                        st.altair_chart(hist_chart, use_container_width=True)
                    st.caption(f"Среднее: {filtered_df[bpka_col].mean():.2f}, Медиана: {filtered_df[bpka_col].median():.2f}")
            
            with col2:
                if apka_col in filtered_df.columns and not filtered_df[apka_col].isna().all():
                    st.write(f"**Распределение Acid pKa**")
                    hist_chart = create_histogram(filtered_df, apka_col, "Acid pKa")
                    if hist_chart:
                        st.altair_chart(hist_chart, use_container_width=True)
                    st.caption(f"Среднее: {filtered_df[apka_col].mean():.2f}, Медиана: {filtered_df[apka_col].median():.2f}")
            
            if logp_col in filtered_df.columns and not filtered_df[logp_col].isna().all():
                st.write("**Распределение LogP**")
                hist_chart = create_histogram(filtered_df, logp_col, "LogP")
                if hist_chart:
                    st.altair_chart(hist_chart, use_container_width=True)
                st.caption(f"Среднее: {filtered_df[logp_col].mean():.2f}, Медиана: {filtered_df[logp_col].median():.2f}")
            
            display_atc_col = atc_code_col if atc_code_col else 'atc_therapeutic_group'
            if display_atc_col in filtered_df.columns:
                st.write("**Распределение по АТХ группам**")
                atc_counts = filtered_df[display_atc_col].value_counts().head(10)
                st.bar_chart(atc_counts, use_container_width=True)
            
            st.write("**📊 Сводная статистика**")
            
            numeric_cols = [bpka_col, apka_col, logp_col, 'dyeles_score']
            available_numeric_cols = [col for col in numeric_cols if col in filtered_df.columns and not filtered_df[col].isna().all()]
            
            if available_numeric_cols:
                stats_df = filtered_df[available_numeric_cols].describe()
                st.dataframe(
                    stats_df.round(2),
                    use_container_width=True,
                    height=300
                )
        else:
            st.warning("⚠️ Нет данных для отображения статистики")
    
    with tab3:
        st.subheader("📋 Полная таблица")
        
        if len(filtered_df) > 0:
            show_smiles = st.checkbox("Показать SMILES", value=False)
            show_ids = st.checkbox("Показать идентификаторы", value=False)
            show_exp = st.checkbox("Показать экспериментальные значения", value=True)
            show_atc_levels = st.checkbox("Показать уровни АТХ", value=False)
            
            table_columns = ['Name', display_atc_col, bpka_col, apka_col, logp_col, 'dyeles_score', 'is_fluorescent', 'oral', 'first_approval']
            
            if show_exp:
                if 'pka exp' in filtered_df.columns:
                    table_columns.append('pka exp')
                if 'pka exp 2' in filtered_df.columns:
                    table_columns.append('pka exp 2')
                if 'LogP_exp' in filtered_df.columns:
                    table_columns.append('LogP_exp')
            
            if show_atc_levels and '_atc_l1' in filtered_df.columns:
                table_columns.extend(['_atc_l1', '_atc_l2', '_atc_l3'])
            
            if show_smiles and 'SMILES' in filtered_df.columns:
                table_columns.append('SMILES')
            
            if show_ids:
                id_columns = ['molecule_chembl_id', 'DrugBank ID', 'CAS Number', 'PubChem Compound ID']
                table_columns.extend([col for col in id_columns if col in filtered_df.columns])
            
            table_columns = [col for col in table_columns if col in filtered_df.columns]
            
            st.dataframe(
                filtered_df[table_columns],
                use_container_width=True,
                height=500
            )
            
            st.markdown(get_download_link(filtered_df[table_columns]), unsafe_allow_html=True)
        else:
            st.warning("⚠️ Нет данных для отображения")
    
    st.divider()
    st.markdown("""
    ### 📚 О проекте
    
    Платформа ChemScan разработана для поддержки принятия решений при выборе лекарственных препаратов не только для тромбоцитарной доставки, но и для широкого спектра стратегий таргетной терапии. Система позволяет моделировать накопление соединений в различных кислых субклеточных компартментах (таких как лизосомы, эндосомы или митохондрии) разных типов клеток, используя физико-химические параметры (pKa, logP) и механизмы ионного захвата.
    
    **Контакты:** nick-pronn@yandex.ru | **Версия:** 1.1
    """)

if __name__ == "__main__":
    main()