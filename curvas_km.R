# ==========================================================
#  1️⃣Librerias
# ==========================================================
library(tidyr);
library(flextable);
library(officer);
library(openxlsx);
library(dplyr);
library(survival);# analisis de supervivencia
library(survminer); # visualización de las curvas de supervivencia
library(viridis); # visualización
library(ggsci) # visaulización (paquete de colores)
library(survMisc)
# ==========================================================
# 🔹 3️⃣ Ajustar el modelo Kaplan–Meier
# ==========================================================
# Calculate survival estimates stratified by Gene Fusion
fit_km <- survfit(surv_EFS ~ gen_fusion, data = my_data)
# View the Life Table (time, n.risk, survival probability)
summary(fit_km)
# Test de Tarone-Ware (librería survMisc)
fit_ten <- ten(fit_km)
comp(fit_ten)
# ==========================================================
# 🔹 4️⃣ Test log-rank (comparación entre curvas)
# ==========================================================
# Run the Log-Rank test comparing groups defined by 'gen_fusion'
# Note: Ensure 'surv_obj' is the Surv() object created earlier (e.g., surv_OS)
log_rank <- survdiff(surv_EFS ~ aneuplo, data = my_data)
# Print the standard output (shows Observed vs Expected events)
log_rank
# --- Manual P-value Extraction ---
# Extract Chi-square and calculate P-value manually for precise reporting
chisq_stat <- log_rank$chisq
df <- length(log_rank$n) - 1 # Degrees of freedom = (Groups - 1)
p_val <- 1 - pchisq(chisq_stat, df)
# Output the result cleanly
cat("Log-Rank Test Result:\n")
cat("Chi-Square:", round(chisq_stat, 2), "\n")
cat("P-value:", round(p_val, 4), "\n")
# ==========================================================
# 🔹 5️⃣ Calcular la mediana de supervivencia por grupo
# ==========================================================
# Extraer la mediana de supervivencia y los intervalos de confianza (IC 95%)
median_surv <- surv_median(fit_km)
# Interpretación rápida para el usuario:
# - 'median': Tiempo en el que el 50% de la población ha tenido el evento.
# - Si es NA, significa que más del 50% sigue vivo al final del estudio.
median_surv
# ==========================================================
# 🔹 6️⃣ Graficar la curva Kaplan–Meier
# ==========================================================
# 1. Convertir la columna de texto a "factor"
my_data$gen_fusion <- as.factor(my_data$gen_fusion)
# 2. Ahora sí, pregúntale los niveles
levels(my_data$gen_fusion)
# cambio de nombre de las etiquetas
mis_etiquetas <- c("BCR::ABL1", 
                   "ETV6::RUNX1", 
                   "KMT2A::AFF1", 
                   "NEGATIVO", 
                   "TCF3::PBX1")
# Generar el gráfico
plot_km <- ggsurvplot(
  fit_km,
  data = my_data,
  # --- 1. Estilo de las líneas y censura ---
  size = 2,              # Grosor de la línea
  linetype = "solid",     # Diferentes tipos de línea por grupo: opciones 'dashed','dotted','dotdash'
  censor.shape = 124,      # Forma de la censura (barra vertical '|')
  censor.size = 3,         # Tamaño de la marca de censura
  # --- 2. Intervalos y Estadísticos ---
  conf.int = FALSE,        # Ocultar intervalo de confianza (sombreado)
  surv.median.line = "none", # "hv" dibuja líneas vert/horiz en la mediana (Opcional, tenías NULL)
  # --- 3. P-Valor ---
  pval = TRUE,             # Mostrar p-valor del Log-Rank
  pval.size = 4,
  pval.coord = c(0, 0.1),  # Coordenadas (x,y) dentro del gráfico. Ajustar según los datos.
  # --- 4. Colores y Tema ---
  palette = c("#B22222", "#4682B4", "#2E8B57", "#D2691E", "#696969"), # Paleta Viridis, 5 número de categorías, opciones: D, A, C, E
  ggtheme = theme_bw(base_size = 14) +             # Tema base blanco y negro
    theme(
      panel.grid = element_blank(),                # Quitar rejilla de fondo
      axis.line = element_line(size = 0.4),        # Líneas de los ejes finas
      legend.position = "right",                   # Leyenda a la derecha
      # ---- AQUÍ MODIFICAS LA LEYENDA ----
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 13),
      legend.key.size = unit(1.2, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),# subtítulo centrado
      axis.title = element_text(face = "bold")     # Títulos de ejes en negrita
    ),
  # --- 5. Tabla de Riesgo (Risk Table) ---
  risk.table = TRUE,
  risk.table.fontsize = 3.5,
  risk.table.height = 0.25,      # Altura de la tabla respecto al gráfico
  risk.table.title = "Pacientes en riesgo",
  risk.table.y.text = FALSE,     # Ocultar nombres de grupos en el eje Y de la tabla (más limpio)
  # --- 6. Etiquetas y Textos ---
  title = "Curvas de supervivencia por Gen de fusión",
  subtitle = 'Supervivencia Libre de Eventos: Cohorte 2017 - 2024',
  xlab = "Tiempo (meses)",
  ylab = "Probabilidad de supervivencia",
  legend.labs = mis_etiquetas, # Reemplaza "gen_fusion=..." por tus textos limpios
  legend.title = "Gen de fusión",
)
# Imprimir el gráfico
print(plot_km)
# ==========================================================
# 7. EXPORTAR GRÁFICO (Alta Resolución .JPG)
# ==========================================================

# 1. Definir el nombre del archivo y sus propiedades
jpeg(filename = "Curva_SLE_GF.jpg", # Nombre del archivo
     width = 12,      # Ancho
     height = 10,     # Alto
     units = "in",    # Unidades (pulgadas)
     res = 300)       # Resolución: 300 dpi (Estándar para publicación/impresión)

# 2. Imprimir el objeto del gráfico
# Es IMPORTANTE usar print() explícitamente para que 'survminer'
# ensamble la curva + la tabla de riesgos correctamente.
print(plot_km)

# 3. Cerrar el dispositivo (Esto es lo que guarda el archivo en tu carpeta)
dev.off()


# ==========================================================
#  SUPERVIVENCIA A TIEMPOS FIJOS (12, 24)
# ==========================================================
# --------------------------------------------------
# Resumen KM a tiempos fijos
# --------------------------------------------------
km_summary <- summary(fit_km, times = c(12, 24), extend = TRUE)
km_table_wide2 <- tibble(
  Grupo        = km_summary$strata,
  Tiempo_meses = km_summary$time,
  Supervivencia = round(km_summary$surv * 100, 1),
  IC_inf       = round(km_summary$lower * 100, 1),
  IC_sup       = round(km_summary$upper * 100, 1)
) %>%
  mutate(
    Grupo = gsub("gen_fusion=", "", Grupo),
    IC = paste0(IC_inf, "–", IC_sup)
  ) %>%
  select(Grupo, Tiempo_meses, Supervivencia, IC) %>%
  pivot_wider(
    names_from  = Tiempo_meses,
    values_from = c(Supervivencia, IC),
    names_glue  = "{.value}_{Tiempo_meses}"
  )
# --------------------------------------------------
# Tabla para Word
# --------------------------------------------------
km_ft2 <- flextable(km_table_wide2) %>%
  set_header_labels(
    Grupo = "Genes de fusión",
    Supervivencia_12 = "12 meses (%)",
    IC_12 = "IC 95% (12 m)",
    Supervivencia_24 = "24 meses (%)",
    IC_24 = "IC 95% (24 m)"
  ) %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  set_caption("Supervivencia libre de eventos a 12 y 24 meses según Genes de Fusión")
# --------------------------------------------------
# Exportación a Word
# --------------------------------------------------
doc <- read_docx() %>%
  body_add_par(
    "Supervivencia estimada por Kaplan–Meier en tiempos fijos",
    style = "heading 1"
  ) %>%
  body_add_par(
    "La tabla muestra la probabilidad de supervivencia acumulada a 12 y 24 meses,
    estratificada por grupo, junto con sus intervalos de confianza al 95%.",
    style = "Normal"
  ) %>%
  body_add_flextable(km_ft2)
print(doc, target = "SLE_KM_tiempos_GF.docx")
# ==========================================================
#  Tabla de Tarone-Ware test
# ==========================================================
salida <- capture.output(comp(fit_ten))
# Identificar dónde empieza la segunda tabla
fin_tabla1 <- grep("\\$tft", salida)[1] - 1
tabla1_txt <- salida[1:fin_tabla1]
# convertir a data frame
tabla_tests <- read.table(
  text = tabla1_txt,
  header = TRUE,
  row.names = 1
)
tabla_tests <- tibble::rownames_to_column(tabla_tests, "Test")
tabla_tests
# extrear Tarone-Ware
tabla_tarone <- dplyr::filter(tabla_tests, Test == "sqrtN") |>
  dplyr::mutate(Test = "Tarone-Ware")
tabla_tarone
# Flextable
flextable(tabla_tarone) |>
  set_header_labels(
    Test   = "Test",
    chiSq  = "Chi-cuadrado",
    df     = "Grados de libertad",
    pChisq = "Valor p"
  ) |>
  colformat_num(j = c("chiSq", "pChisq"), digits = 4) |>
  theme_booktabs() |>
  autofit()
# Exportar a word
doc <- read_docx() |>
  body_add_par("Comparación de curvas de supervivencia", style = "heading 1") |>
  body_add_par("Prueba de Tarone-Ware", style = "heading 2") |>
  body_add_flextable(ft_tarone)
print(doc, target = "tarone_SLE_GF.docx")
# ==========================================================
#  Tabla de Log-rank test
# ==========================================================
# Extract statistics
chisq_stat <- log_rank$chisq
n_groups   <- length(log_rank$n)
df         <- n_groups - 1
p_val      <- 1 - pchisq(chisq_stat, df)
#--------------------------------------------------
# CREATE SUMMARY TABLE (publication-friendly)
#--------------------------------------------------
logrank_table <- data.frame(
  variable = "gen_fusion",
  n_grupos = n_groups,
  chi2 = round(chisq_stat, 2),
  gl = df,
  p_valor = ifelse(p_val < 0.001, "< 0.001", round(p_val, 4))
)
#--------------------------------------------------
# CONVERT TO FLEXTABLE
#--------------------------------------------------
logrank_ft <- flextable(logrank_table) %>%
  set_header_labels(
    variable = "Variable",
    n_grupos = "N° de grupos",
    chi2 = "χ²",
    gl = "gl",
    p_valor = "p-valor"
  ) %>%
  autofit() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  set_caption("Comparación de curvas de supervivencia mediante el test de log-rank")
#--------------------------------------------------
# EXPORT TO WORD
#--------------------------------------------------
doc <- read_docx() %>%
  body_add_par("Resultados del test de log-rank", style = "heading 2") %>%
  body_add_flextable(logrank_ft)
print(doc, target = "LogRank_SLE_GF.docx")





mis_etiquetas <- c("BCR::ABL1", 
                   "ETV6::RUNX1", 
                   "KMT2A::AFF1", 
                   "NEGATIVO", 
                   "TCF3::PBX1")
mis_etiquetas <- c("46 chr", 
                   "46 chr + Alt", 
                   "Hiperdiploidía", 
                   "Hipodiploidía", 
                   "No crecimiento")


