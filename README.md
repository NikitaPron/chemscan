# 🔬 ChemScan

> A computational platform for drug screening based on subcellular accumulation and targeted delivery criteria

## 🎯 About the Project

**ChemScan** is a computational platform for selecting approved drugs based on their physicochemical properties suitable for accumulation in acidic subcellular compartments (lysosomes, endosomes, platelet granules) through ion trapping.

### 🔑 Key Features

- ✅ **Basic pKa** filtering (recommended range: 8–12)
- ✅ **Acid pKa** and **LogP** filtering
- ✅ Support for calculated and experimental values
- ✅ Cascade filtering by **ATC classification** (3 levels)
- ✅ Fluorescence property assessment via **DyeLeS score**
- ✅ **Lipinski's Rule of Five** descriptors (MW, LogP, HBD, HBA, TPSA) with optional filters
- ✅ Molecular structure visualization (RDKit)
- ✅ CSV export of results
- ✅ Interactive histograms and statistics
- ✅ Drug lookup by name, ChEMBL ID, or DrugBank ID

## 📦 Dataset (`table_drugs.csv`)

The ChemScan database contains **2,063 unique approved small-molecule drugs** assembled from public **ChEMBL** and **DrugBank** snapshots (retrieved February 2026).

### Inclusion criteria

- Regulatory approval in at least one major jurisdiction (`max_phase = 4` in ChEMBL; approved annotation in DrugBank)
- Valid canonical **SMILES** string (entries without SMILES or with unparseable structures were excluded)
- Deduplication by canonical SMILES

### Descriptor calculation

| Field | Source / method |
|-------|-----------------|
| Basic and acidic pKa (predicted) | ChemAxon Marvin |
| logP, MW, TPSA, HBD, HBA (predicted) | RDKit |
| pKa, logP (experimental) | ChEMBL / DrugBank annotations where available |
| DyeLeS score | Pre-computed fluorescence-likelihood metric |
| `ABCB1_Pgp_CHEMBL` | Curated ChEMBL bioactivity records against human ABCB1 (CHEMBL4302) |
| ATC codes, oral flag, first-approval year | Repository metadata |

Experimental coverage in the current snapshot: predicted basic pKa and logP for all retained compounds; experimental pKa for ~21% and experimental logP for ~52% of entries.

The main screening table is `table_drugs.csv` in the project root. The Streamlit app loads this file at startup. Where available, `ABCB1_Pgp_CHEMBL` entries are curated experimental bioactivity annotations imported from ChEMBL (not model-predicted values).

### Recent platform updates (post peer review)

- Optional **Lipinski-like filters** in the sidebar (MW, predicted logP, HBD, HBA, TPSA; off by default)
- **ABCB1_Pgp_CHEMBL** annotations shown verbatim on each drug card
- **Drug lookup** by name, ChEMBL ID, or DrugBank ID (single search field)
- Lipinski descriptor panels in the **Statistics** tab and additional columns in the **Table** export

## 🚀 Quick Start

### Local Launch

```bash
# 1. Clone the repository
git clone https://github.com/NikitaPron/chemscan.git
cd chemscan

# 2. Create a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: .\venv\Scripts\Activate.ps1

# 3. Install dependencies
pip install -r requirements.txt

# 4. Prepare data
# Place the table_drugs.csv file in the project root

# 5. Run the application
streamlit run app.py