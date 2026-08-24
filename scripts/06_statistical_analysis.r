# Instalar paquetes necesarios
install.packages(c("tidyverse","ggplot2","epitools","ggpubr"))

# Cargar librerías
library(tidyverse)

# Crear carpeta para figuras
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Cargar archivos
clinical <- read_tsv("data/processed/clinical_clean.tsv")
activities <- read_tsv("results/sigprofiler_extraction/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt")

# Trim del barcode
activities <- activities %>%
  mutate(Samples = substr(Samples, 1, 12))

# Merge
merged <- clinical %>%
  inner_join(activities, by = c("cases.submitter_id" = "Samples"))

nrow(merged)

# Reorganizar datos en formato largo para ggplot
merged_long <- merged %>%
  pivot_longer(cols = c(SBS1, SBS5, SBS10a, SBS10b, SBS15, SBS107),
               names_to = "Signature",
               values_to = "Activity")

# Boxplot de actividad por estadio (con outliers)
ggplot(merged_long, aes(x = diagnoses.ajcc_pathologic_stage, 
                         y = Activity, 
                         fill = Signature)) +
  geom_boxplot() +
  facet_wrap(~Signature, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Signature activity by pathological stage",
       x = "Stage", y = "Number of mutations")

ggsave("results/figures/boxplot_by_stage_with_outliers.png", width = 12, height = 8, dpi = 300)

# Ver los pacientes con más mutaciones totales
merged %>%
  mutate(total = SBS1 + SBS5 + SBS10a + SBS10b + SBS15 + SBS107) %>%
  arrange(desc(total)) %>%
  select(cases.submitter_id, diagnoses.ajcc_pathologic_stage, total, SBS10a, SBS10b) %>%
  head(10)

# Identificar y excluir hipermutados (>500 mutaciones totales)
merged <- merged %>%
  mutate(total_mutations = SBS1 + SBS5 + SBS10a + SBS10b + SBS15 + SBS107)

# Separar hipermutados
hypermutated <- merged %>% filter(total_mutations > 500)
merged_filtered <- merged %>% filter(total_mutations <= 500)

cat("Pacientes excluidos (hipermutados):", nrow(hypermutated), "\n")
cat("Pacientes en análisis principal:", nrow(merged_filtered), "\n")

# Reorganizar datos en formato largo para ggplot (sin outliers)
merged_without_outliers <- merged_filtered %>%
  pivot_longer(cols = c(SBS1, SBS5, SBS10a, SBS10b, SBS15, SBS107),
               names_to = "Signature",
               values_to = "Activity")

# Boxplot de actividad por estadio (sin outliers)
ggplot(merged_without_outliers, aes(x = diagnoses.ajcc_pathologic_stage, 
                                    y = Activity, 
                                    fill = Signature)) +
  geom_boxplot() +
  facet_wrap(~Signature, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Signature activity by pathological stage (hypermutated excluded)",
       x = "Stage", y = "Number of mutations")

ggsave("results/figures/boxplot_by_stage_filtered.png", width = 12, height = 8, dpi = 300)

# Clasificar en tempranos vs tardíos
merged_filtered <- merged_filtered %>%
  mutate(stage_group = case_when(
    diagnoses.ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage IIA") ~ "Early",
    diagnoses.ajcc_pathologic_stage %in% c("Stage IIB", "Stage III", "Stage IV") ~ "Late",
    TRUE ~ NA_character_
  ))

table(merged_filtered$stage_group)

# Test de Wilcoxon para cada firma
signatures <- c("SBS1", "SBS5", "SBS10a", "SBS10b", "SBS15", "SBS107")

results <- data.frame(Signature = signatures,
                      p_value = NA,
                      median_early = NA,
                      median_late = NA)

for (i in seq_along(signatures)) {
  sig <- signatures[i]
  early <- merged_filtered %>% filter(stage_group == "Early") %>% pull(!!sym(sig))
  late <- merged_filtered %>% filter(stage_group == "Late") %>% pull(!!sym(sig))
  
  results$p_value[i] <- wilcox.test(early, late)$p.value
  results$median_early[i] <- median(early)
  results$median_late[i] <- median(late)
}

results$p_adjusted <- p.adjust(results$p_value, method = "bonferroni")
print(results)

# Para cada firma, crear tabla 2x2 y test de Fisher
results_fisher <- data.frame(
  Signature = signatures,
  p_value = NA,
  prop_early = NA,
  prop_late = NA
)

for (i in seq_along(signatures)) {
  sig <- signatures[i]
  
  early_active <- merged_filtered %>% 
    filter(stage_group == "Early") %>% 
    summarise(active = sum(!!sym(sig) > 0), total = n())
  
  late_active <- merged_filtered %>% 
    filter(stage_group == "Late") %>% 
    summarise(active = sum(!!sym(sig) > 0), total = n())
  
  tabla <- matrix(c(early_active$active, 
                    early_active$total - early_active$active,
                    late_active$active,
                    late_active$total - late_active$active),
                  nrow = 2)
  
  results_fisher$p_value[i] <- fisher.test(tabla)$p.value
  results_fisher$prop_early[i] <- early_active$active / early_active$total
  results_fisher$prop_late[i] <- late_active$active / late_active$total
}

results_fisher$p_adjusted <- p.adjust(results_fisher$p_value, method = "bonferroni")
print(results_fisher)


# Preparar datos para el gráfico de proporciones
results_fisher_long <- results_fisher %>%
  pivot_longer(cols = c(prop_early, prop_late),
               names_to = "stage",
               values_to = "proportion") %>%
  mutate(group = ifelse(stage == "prop_early", "Early", "Late"))

# Gráfico de barras
ggplot(results_fisher_long, aes(x = Signature, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(title = "Proportion of Active Signatures by Stage",
       x = "Signature", y = "Proportion")

ggsave("results/figures/proportion_active_early_vs_late.png", width = 10, height = 6, dpi = 300)

# Heatmap de actividad por paciente ordenado por estadio
merged_filtered_ordered <- merged_filtered %>%
  arrange(diagnoses.ajcc_pathologic_stage) %>%
  mutate(cases.submitter_id = factor(cases.submitter_id, levels = cases.submitter_id))

merged_heatmap <- merged_filtered_ordered %>%
  pivot_longer(cols = c(SBS1, SBS5, SBS10a, SBS10b, SBS15, SBS107),
               names_to = "Signature",
               values_to = "Activity")

ggplot(merged_heatmap, aes(x = cases.submitter_id, y = Signature, fill = log1p(Activity))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Mutational signature activity per patient (ordered by stage)",
       x = "Patients", y = "Signature", fill = "log(Activity+1)")

ggsave("results/figures/heatmap_signatures_by_patient.png", width = 14, height = 5, dpi = 300)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)

# Preparar matriz
mat <- merged_filtered_ordered %>%
  select(SBS1, SBS5, SBS10a, SBS10b, SBS15, SBS107) %>%
  as.matrix() %>%
  t()

colnames(mat) <- merged_filtered_ordered$cases.submitter_id

# Aplicar log1p
mat_log <- log1p(mat)

# Anotación de estadio
stage_colors <- c(
  "Stage I" = "#fee8c8",
  "Stage IA" = "#fdbb84", 
  "Stage IB" = "#e34a33",
  "Stage II" = "#ffffb2",
  "Stage IIA" = "#fecc5c",
  "Stage IIB" = "#fd8d3c",
  "Stage III" = "#f03b20",
  "Stage IV" = "#bd0026"
)

ha <- HeatmapAnnotation(
  Stage = merged_filtered_ordered$diagnoses.ajcc_pathologic_stage,
  col = list(Stage = stage_colors)
)

# Heatmap
ht <- Heatmap(mat_log,
        name = "log(Activity+1)",
        top_annotation = ha,
        show_column_names = FALSE,
        row_names_side = "left",
        col = colorRamp2(c(0, 2, 4), c("white", "salmon", "darkred")),
        column_title = "Patients ordered by stage",
        row_title = "Signature")

png("results/figures/heatmap_signatures_annotated.png", width = 1400, height = 600, res = 150)
draw(ht)
dev.off()

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

gse57495 <- getGEO(filename = "data/raw/GSE57495_series_matrix.txt.gz", 
                   getGPL = FALSE)

# Ver metadatos clínicos
pdata <- pData(gse57495)
colnames(pdata)

table(pdata$`vital.status:ch1`)
table(pdata$`Stage:ch1`)
head(pdata$`overall survival (month):ch1`)

install.packages("pROC")
library(pROC)

# Preparar datos de supervivencia
survival_data <- pdata %>%
  select(
    os_months = `overall survival (month):ch1`,
    vital_status = `vital.status:ch1`,
    stage = `Stage:ch1`
  ) %>%
  mutate(
    os_months = as.numeric(os_months),
    deceased = ifelse(vital_status == "DEAD", 1, 0),
    stage_numeric = case_when(
      stage == "1" ~ 1,
      stage == "1B" ~ 2,
      stage == "2A" ~ 3,
      stage == "2B" ~ 4,
      TRUE ~ NA_real_
    )
  )

# Curva ROC: estadio como predictor de mortalidad
roc_curve <- roc(survival_data$deceased, 
                 survival_data$stage_numeric,
                 na.rm = TRUE)

# Ver AUC
auc(roc_curve)

# Graficar
png("results/figures/roc_curve_stage_survival.png", width = 800, height = 800, res = 150)
plot(roc_curve, 
     main = "ROC curve: stage as predictor of mortality (GSE57495)",
     print.auc = TRUE,
     col = "darkred",
     lwd = 2)
dev.off()

# Tabla resumen de firmas
tabla_resumen <- data.frame(
  Signature = c("SBS1", "SBS5", "SBS10a", "SBS10b", "SBS15", "SBS107"),
  Aetiology = c("Spontaneous deamination (ageing)",
                "Unknown (ageing-related)",
                "POLE exonuclease domain mutations",
                "POLE exonuclease domain mutations", 
                "Defective DNA mismatch repair",
                "Unknown"),
  Present_in_PDAC = c("Yes", "Yes", "Yes (rare)", "Yes (rare)", "Yes", "Yes (rare)"),
  Prop_Early = round(results_fisher$prop_early, 3),
  Prop_Late = round(results_fisher$prop_late, 3),
  Fisher_p_adjusted = round(results_fisher$p_adjusted, 3),
  Stage_association = c("None", "None", "Late only", "Late > Early", "Early > Late", "Early > Late")
)

print(tabla_resumen)
s
# Guardar como CSV
write.csv(tabla_resumen, "results/signature_summary_table.csv", row.names = FALSE)
