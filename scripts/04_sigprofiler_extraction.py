from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matGen.SigProfilerMatrixGeneratorFunc(
    project="TFM_PDAC",      # nombre del proyecto
    reference_genome="GRCh38",       # genoma de referencia
    path_to_input_files="data/raw/maf_clean/",     # ruta a la carpeta con los vcf
    exome=True,        # True si es exoma, False si es genoma completo
    bed_file=None,     # None porque no hay un archivo BED específico
    chrom_based=False,  # False por defecto
    plot=True,         # True para generar gráficos
)
