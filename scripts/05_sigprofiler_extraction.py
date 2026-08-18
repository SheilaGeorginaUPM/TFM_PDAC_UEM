from SigProfilerExtractor import sigpro as sig

if __name__ == '__main__':
	sig.sigProfilerExtractor(
   	 input_type="matrix",      # formato del input
    	output="results/sigprofiler_extraction/",          # carpeta donde guardar los resultados
    	input_data="data/processed/matrices/SBS/TFM_PDAC.SBS96.exome",      # ruta al archivo de input
    	reference_genome="GRCh38", # genoma de referencia
    	minimum_signatures=1, # mínimo de firmas a probar
    	maximum_signatures=10, # máximo de firmas a probar
    	nmf_replicates=100   # número de repeticiones para estabilidad
)
