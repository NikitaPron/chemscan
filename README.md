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
- ✅ Molecular structure visualization (RDKit)
- ✅ CSV export of results
- ✅ Interactive histograms and statistics

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