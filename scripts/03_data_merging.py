from pathlib import Path
import pandas as pd

# Base path (ruta del script) — evita depender del CWD
BASE = Path(__file__).resolve().parent.parent
FILE1 = BASE / "data/processed/clinical_clean.tsv"
FILE2 = BASE / "data/processed/TCGA-PAAD_merged.maf"
OUTPUT = BASE / "data/processed/merged_clinical_TCGA-PAAD.tsv"

# Leer tablas (lanza FileNotFoundError si faltan)
df1 = pd.read_csv(FILE1, sep="\t")
df2 = pd.read_csv(FILE2, sep="\t")

# Tomar los primeros X caracteres de la columna necesaria en df2
X = 12
df2['Tumor_Sample_Barcode'] = df2['Tumor_Sample_Barcode'].str[:X]

# Unir los DataFrames (ajusta left_on/right_on según columnas reales)
merged = df1.merge(df2, left_on='cases.submitter_id', right_on='Tumor_Sample_Barcode', how='inner')

# Guardar resultado
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(OUTPUT, sep="\t", index=False)
print(f"Merge completado: {merged.shape[0]} filas -> {OUTPUT}")
