#!/bin/bash
# Script 01: Descarga y preparación de datos TCGA-PAAD
# TFM - Firmas mutacionales en PDAC
# Autor: Sheila Georgina Jimenez

# ── 0. DESCARGA DESDE GDC DATA PORTAL ─────────────────────────
# Requisito: descargar el GDC Data Transfer Tool desde
# https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
#
# Pasos manuales previos:
# 1. Acceder a https://portal.gdc.cancer.gov/repository
# 2. Filtros: Project=TCGA-PAAD, Data Category=Simple Nucleotide Variation,
#    Data Type=Masked Somatic Mutation, Workflow Type=Aliquot Ensemble
# 3. Añadir todos los archivos al carrito y descargar el manifest (gdc_manifest.txt)
# 4. Ejecutar la descarga con el cliente GDC:
#
# ./gdc-client download -m gdc_manifest.txt -d data/raw/
#
# Alternativamente, descargar el .tar.gz directamente desde el carrito
# y descomprimir en data/raw/:
#
# tar -xzf gdc_download_*.tar.gz -C data/raw/

# ── 1. FUSIÓN DE MAF INDIVIDUALES ──────────────────────────────
# Extraer cabecera del primer archivo MAF
primera=$(ls -d tcga_paad/*/ | head -1)
archivo=$(ls "$primera"*.maf.gz)
gzcat "$archivo" | grep -v "^#" | sed -n '1p' > data/processed/cabecera.tsv

# Concatenar cuerpo de los 180 MAF (sin cabecera)
for muestra in tcga_paad/*/*.maf.gz; do
    gzcat "$muestra" | grep -v "^#" | sed -n '2,$p'
done >> data/processed/cuerpo.tsv

# Unir cabecera + cuerpo en un único MAF fusionado
cat data/processed/cabecera.tsv data/processed/cuerpo.tsv > data/processed/TCGA-PAAD_merged.maf

# ── 2. VERIFICACIÓN DE SANIDAD ──────────────────────────────────
# Contar mutaciones y muestras únicas
wc -l data/processed/TCGA-PAAD_merged.maf
cut -f16 data/processed/TCGA-PAAD_merged.maf | sed -n '2,$p' | sort -u | wc -l

# Top 10 genes más mutados
cut -f1 data/processed/TCGA-PAAD_merged.maf | sed -n '2,$p' | sort | uniq -c | sort -rn | head -10
