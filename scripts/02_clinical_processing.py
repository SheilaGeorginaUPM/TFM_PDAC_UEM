#!/usr/bin/env python3
"""
Script 02: Procesamiento de datos clínicos TCGA-PAAD
Extrae submitter_id y estadio patológico, quedándose con
el estadio más avanzado por paciente.
"""

import pandas as pd

# Orden de estadios de menor a mayor
stage_order = {
    'Stage 0': 0, 'Stage 0a': 1, 'Stage 0is': 2,
    'Stage I': 3, 'Stage IA': 4, 'Stage IB': 5,
    'Stage II': 6, 'Stage IIA': 7, 'Stage IIB': 8,
    'Stage III': 9, 'Stage IV': 10
}

# Cargar solo las columnas que necesitamos
df = pd.read_csv(
    'data/raw/clinical.tsv',
    sep='\t',
    usecols=['cases.submitter_id', 'diagnoses.ajcc_pathologic_stage']
)

# Eliminar filas sin estadio
df = df.dropna(subset=['diagnoses.ajcc_pathologic_stage'])
df = df[df['diagnoses.ajcc_pathologic_stage'] != "'--"]

# Añadir columna numérica de orden
df['stage_rank'] = df['diagnoses.ajcc_pathologic_stage'].map(stage_order)

# Quedarse con el estadio más avanzado por paciente
df = df.sort_values('stage_rank', ascending=False)
df = df.drop_duplicates(subset='cases.submitter_id', keep='first')

# Eliminar columna auxiliar y guardar
df = df.drop(columns='stage_rank')
df.to_csv('data/processed/clinical_clean.tsv', sep='\t', index=False)

print(f"Pacientes únicos con estadio: {len(df)}")
print(df['diagnoses.ajcc_pathologic_stage'].value_counts())
