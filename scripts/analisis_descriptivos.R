# ==========================================================
# 📌 ANÁLISIS DESCRIPTIVO
# Autor: S.Caro
# Fecha: 13/01/2026
# Objetivo: análisis descriptivo: tablas y gráficos
# ==========================================================

# 1️⃣ Librerías
# 2️⃣ Carga de datos
# 3️⃣ Limpieza y preparación
# 4️⃣ Grafico:histograma
# 5️⃣ Modelos estadísticos

# ===================================================================
# 1️⃣ Librerías
# ===================================================================
# 📊 Visualización de datos
library(forcats)
library(ggplot2)
library(ggdist)
library(ggrain)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
# 📋 Tablas y presentación de resultados
library(gt)
# 📁 Manipulación y lectura de datos
library(dplyr)
library(readr)
# 📤 Importación y Exportación de datos
library(openxlsx)

# ===================================================================
# 2️⃣ Carga de datos
# ===================================================================
my_data <- read.xlsx('clean_data.xlsx')

# ===================================================================
# 3️⃣ Limpieza y preparación
# ===================================================================
# Correción del formato de fechas
vars_fechas <- c("fecha_nacimiento","fecha_dx","fecha_recaida",
                 "fecha_TPH","fecha_ultima","fecha_galenos","diferimiento_i","diferimiento_f")
my_data <- my_data %>%
  mutate(across(all_of(vars_fechas),~ as.Date(.x, origin = "1899-12-30")))
# Consistencia de las fechas
checks_fechas <- my_data %>%
  mutate(
    regla_nac_dx_ult = # la fecha de dx es menor a la fecha ultima
      fecha_nacimiento <= fecha_dx &
      fecha_dx <= fecha_ultima,
    regla_recaida_dx = # la fecha de dx es menor a la fecha de recaida
      is.na(fecha_recaida) | fecha_recaida >= fecha_dx,
    regla_TPH_dx = # la fecha de dx es menor a la fecha de TPH
      is.na(fecha_TPH) | fecha_TPH >= fecha_dx
  )
# Resumen del analisis
tabla_checks <- checks_fechas %>%
  reframe(
    Regla = c(
      "Nacimiento ≤ Dx ≤ Última","Recaída ≥ Dx (o NA)","TPH ≥ Dx (o NA)"
    ),
    Errores = c(
      sum(!regla_nac_dx_ult, na.rm = TRUE),
      sum(!regla_recaida_dx, na.rm = TRUE),
      sum(!regla_TPH_dx, na.rm = TRUE)
    )
  )
tabla_checks
# Revisar valores
checks_fechas %>%
  filter(!regla_recaida_dx | !regla_TPH_dx) %>%
  select(fecha_nacimiento,fecha_dx,fecha_recaida,fecha_TPH,fecha_ultima)

# Para analizar los datos dentro de las variables
my_data %>%
  select(all_of(vars_demo)) %>%
  str()
# Para obtener un resumen estadístico de las variables
my_data %>%
  select(all_of(vars_demo)) %>%
  summary()

# ------------------------------------------------------------------
# GRÁFICO HISTOGRAMA
# ------------------------------------------------------------------
n_total  <- nrow(my_data)
n_valid  <- sum(!is.na(my_data$talla_dx))
ggplot(my_data, aes(x = talla_dx)) + # creación del objeto ggplot
# ----------------------------------------------------------
# Capa 1: geom_histogram() construye el histograma de la variable
# ----------------------------------------------------------
  geom_histogram(
    binwidth = 5, # ancho de las barras
    fill = "#4DB6AC", # color de relleno de las barras
    color = "white", # color del borde de las barras
    alpha = 0.8, # transparencia del relleno
    na.rm = TRUE # eliminar valores NA antes de graficar
  ) +
# ----------------------------------------------------------
# Capa 2: geom_vline() agrega una línea vertical al gráfico
# ----------------------------------------------------------
  geom_vline(
    aes(
      xintercept = median(talla_dx, na.rm = TRUE),
      linetype = "Mediana"),
    color = "#8a1708",linewidth = 1.5
  ) +
  scale_linetype_manual(
    name = "Línea",values = c("Mediana" = "dashed")
  ) +
# ----------------------------------------------------------
# Etiquetas del gráfico: labs() define títulos y nombres de ejes
# ----------------------------------------------------------
  labs(
    title = "Distribución del AMO33 al diagnóstico",
    subtitle = paste0("Total: ", n_total," | Válidos: ", n_valid),
    x = "Células leucémicas (%)",
    y = "Frecuencia (n)"
  ) +
# ----------------------------------------------------------
# Tema base del gráfico: theme_classic() aplica un estilo limpio
# ----------------------------------------------------------
  guides(
  linetype = guide_legend(nrow = 1)
  ) +
  theme_classic(base_size = 16) +
# ----------------------------------------------------------
# Ajustes finos del tema: theme() modifica elementos específicos del gráfico
# ----------------------------------------------------------
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), # centra el titulo y titulo en negrita
    plot.subtitle = element_text(size = 14, hjust = 0.5), # subtitulo centrado
    axis.title = element_text(size = 12, face = "bold"), # titulos y ejes en negrita
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )


# ------------------------------------------------------------------
# GRÁFICO BARRAS SIMPLE
# ------------------------------------------------------------------
n_total  <- nrow(my_data)
n_valid  <- sum(!is.na(my_data$aneuplo))
# ------------------------------------------------------------------
df_prop <- my_data |>
  dplyr::filter(!is.na(aneuplo)) |>
  dplyr::count(aneuplo) |>
  dplyr::mutate(prop = n / sum(n))
df_prop <- df_prop |>
  dplyr::mutate(
    aneuplo = forcats::fct_recode(
      aneuplo,
      "Hiperdiploidía" = "HIPERDIPLOIDIA",
      "Hipodiploidía" = "HIPODIPLOIDIA",
      "No Crecimiento" = "NO CRECIMIENTO",
      "46 chr" = "CARIOTIPO NORMAL",
      "46 chr + Alt" = "CARIOTIPO NORMAL + ALTERACIONES",
    ))
df_prop <- df_prop |>
  dplyr::mutate(
    aneuplo = factor(
      aneuplo,
      levels = aneuplo[order(-prop)]
    ))
# ------------------------------------------------------------------
n_total  <- nrow(my_data)
n_valid  <- sum(!is.na(my_data$gen_fusion))
df_prop <- my_data |>
  dplyr::filter(!is.na(gen_fusion)) |>
  dplyr::count(gen_fusion) |>
  dplyr::mutate(prop = n / sum(n))
df_prop <- df_prop |>
  dplyr::mutate(
    gen_fusion = forcats::fct_recode(
      gen_fusion,
      "BCR::ABL1"    = "BCR::ABL1",
      "ETV6::RUNX1"  = "ETV6::RUNX1",
      "KMT2A::AFF1"  = "KMT2A::AFF1",
      "TCF3::PBX1"   = "TCF3::PBX1",
      "Negativo"     = "NEGATIVO"
    ))
df_prop <- df_prop |>
  dplyr::mutate(
    gen_fusion = factor(
      gen_fusion,
      levels = gen_fusion[order(-prop)]
    ))
# ------------------------------------------------------------------
ggplot(df_prop, aes(x = gen_fusion, y = prop)) +
  geom_col(
    fill = "#4DB6AC",
    color = "white",
    width = 0.7
  ) +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 1)),
    vjust = -0.3,
    size = 5,
    fontface = "bold"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.1))  # aumento leve para que no corte el texto
  ) +
  labs(
    title = "Genes de Fusión",
    subtitle = paste0("Total: ", n_total," | Válidos: ", n_valid),
    x = "Categorías",
    y = "Frecuencia (%)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

# ------------------------------------------------------------------
# GRÁFICO BOXPLOT
# ------------------------------------------------------------------
ggplot(my_data, aes(x = l_con, y = gen_fusion)) +
  # Patrón de circulos
  geom_quasirandom( size = 3, color = "#428d99", alpha = 0.4, shape = 16) +
  # Figura de violin transparente
  geom_violin(fill = "#4DB6AC", alpha = 0.3, color = NA) +
  # Figura de caja
  geom_boxplot(fill = "#4DB6AC",color = "#2C3E50",linewidth = 1.1,width = 0.2,outlier.shape = NA) +
  # Punto de la media
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 4, color = "#8a1708") +
  # Agregar el p-value
  #ggpubr::stat_compare_means(method = "wilcox.test",
                             #label.y = max(my_data$edad_anios) + 1) +
  # Renombrar los valores
  scale_x_discrete(
    labels = c("F" = "Femenino", 
               "M" = "Masculino")
  ) +
  # Etiquetas
  labs(
    title = "Distribución de leucocitos según Genes de fusión",
    x = "Genes de fusión",
    y = "Leucocitos (años)"
  ) +
  # Tema de los titulos y ejes
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ------------------------------------------------------------------
# 6️⃣ Gráfico Scatterplot
# -------------------------------------------------
ggplot(my_data, aes(x = DHL_dx, y = blastos_dx)) +
  geom_point(alpha = 0.7, size = 3, color = "#1f77b4") +   # puntos profesionales
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) + # línea de tendencia
  labs(
    title = "Relación entre Variable X y Variable Y",
    subtitle = "Pacientes pediátricos con B-ALL",
    x = "Variable X (unidades)",
    y = "Variable Y (unidades)"
  ) +
  # Títulos
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    # Etiquetas de ejes
    axis.title = element_text(face = "bold", size = 16),
    # Texto de ejes
    axis.text = element_text(size = 14, color = "black"),
    # Líneas del tema
    axis.line = element_line(linewidth = 0.8, color = "black"),
    # Márgenes
    plot.margin = margin(10, 10, 10, 10)
    )
# ------------------------------------------------------------------
# 7️⃣ Gráfico de Barras apiladas
# ------------------------------------------------------------------
df_prop <- my_data %>%
  group_by(grupo_OMS, infiltracion_extramedular_dx) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

ggplot(df_prop, aes(
  x = grupo_OMS,
  y = prop,
  fill = infiltracion_extramedular_dx
)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#4C72B0", "#55A868")) +
  labs(
    title = "Proporción de infiltración extramedular según fusión genética",
    x = "Fusión genética",
    y = "Proporción (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# ------------------------------------------------------------------
# 8️⃣ Gráfico de Heatmaps

# ------------------------------------------------------------------
# 9️⃣ Gráfico rainplot
# ------------------------------------------------------------------
library(ggplot2)
library(ggdist)
p1 <- ggplot(my_data, aes(x = sexo, y = edad_anios, fill = sexo)) +
  # Half violin (densidad)
  stat_halfeye(adjust = 0.6, width = 0.6, justification = -0.3, .width = 0) +
  # Boxplot
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  # Puntos individuales (jitter)
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.3) +
  
  theme_minimal() +
  labs(
    title = "Raincloud plot de edad según sexo",
    x = "Sexo",
    y = "Edad (años)"
  )

p1
# ------------------------------------------------------------------