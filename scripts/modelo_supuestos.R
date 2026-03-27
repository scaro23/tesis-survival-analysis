# ==========================================================
#  1️⃣ Librerias
# ==========================================================
library(dplyr)
library(flextable)
library(car)
# ==========================================================
cox_multi1 <- coxph(
  surv_EFS ~ 
    gen_fusion + aneuplo + edad_anios_exacto + SNC_dx +
    TESTICULO_dx + log2_leuco + plaq_100 + log2_urico + s_down +
    sexo,
  data = my_data)
# ==========================================================
# 1️⃣ Supuesto de riesgos proporcionales
# ==========================================================
ph_test <- cox.zph(cox_multi1)
print(ph_test)
# -------------------------------
# Gráficos de residuos de Schoenfeld
# -------------------------------
plot(
  ph_test,
  var = "aneuplo",
  col = "steelblue",
  lwd = 2,
  xlab = "Tiempo de seguimiento (años)",
  ylab = "",
  yaxt = "n",
  main = "Evaluación del supuesto PH: Grupo de Ploidía")
axis(2)
mtext("Residuos de Schoenfeld escalados", side = 2, line = 3)
abline(h = 0, lty = 2, col = "red")
legend(
  "topright",
  legend = c(
    "Residuos de Schoenfeld (suavizado)",
    "Referencia: riesgo proporcional"),
  col = c("steelblue", "red"),
  lty = c(1, 2),
  lwd = c(2, 1),
  bty = "n")
# ==========================================================
# 2️⃣ RESIDUOS DE MARTINGALE (evaluacion de la linealidad)
# ==========================================================
library(splines)
cox_spline_leuco <- coxph(
  surv_OS ~ 
    gen_fusion + aneuplo + edad_anios_exacto + SNC_dx +
    TESTICULO_dx + ns(log2_leuco, df = 3) +
    plaq_100 + log2_urico + s_down + sexo,
  data = my_data)
anova(cox_multi1, cox_spline_leuco)
ph_test_2 <- cox.zph(cox_spline_leuco)
print(ph_test_2)
print(ph_test)

mf <- model.frame(cox_multi1)
martingale_resid <- residuals(cox_multi1, type = "martingale")
plot(
  mf$edad_anios_exacto,
  martingale_resid,
  pch = 19,
  col = "gray40",
  xlab = "Edad (años)",
  ylab = "",
  yaxt = "n",
  main = "Evaluación de la forma funcional: Edad"
)
# Eje Y manual
axis(2, las = 1)
mtext("Residuos de Martingale", side = 2, line = 3)
# Línea de referencia
abline(h = 0, lty = 2, col = "black")
# Suavizado
lines(
  lowess(mf$edad_anios_exacto, martingale_resid),
  col = "steelblue",
  lwd = 2
)
# Leyenda
legend(
  "topright",
  legend = c(
    "Residuos de Martingale",
    "Suavizado (LOESS)",
    "Referencia (0)"
  ),
  col = c("gray40", "steelblue", "black"),
  pch = c(19, NA, NA),
  lty = c(NA, 1, 2),
  lwd = c(NA, 2, 1),
  bty = "n",
  cex = 0.9
)
# ==========================================================
# regla de Events Per Variable (EPV)
# ==========================================================
# Número de eventos (1 = evento, 0 = censura)
n_eventos <- sum(my_data$status_EFS == 1, na.rm = TRUE)
n_eventos
# Contar parámetros del modelo
fit <- cox_multi1
# Número de coeficientes estimados
n_parametros <- length(coef(fit))
n_parametros
# Calcular EPV
EPV <- n_eventos / n_parametros
EPV
# Interpretación
if (EPV >= 10) {
  print("EPV adecuado (bajo riesgo de sobreajuste)")
} else if (EPV >= 5) {
  print("EPV intermedio (interpretar con cautela)")
} else {
  print("EPV bajo (alto riesgo de sobreajuste, modelo inestable)")
}
# ==========================================================
# Multicolinealidad → VIF
# ==========================================================
vif_values <- vif(cox_multi1)
print(vif_values)
# ==========================================================
# ANALISIS OPCIONALES
# ==========================================================
# dfbeta
dfb <- residuals(cox_multi1, type = "dfbeta")
matplot(dfb, type = "p")
abline(h = c(-0.2, 0.2), lty = 2)
# Residuos de devianza
dev_res <- residuals(cox_multi1, type = "deviance")
plot(dev_res, ylab = "Residuos de devianza")
abline(h = 0, lty = 2)


# ==========================================================
# CREAR UNA TABLA (1)
# ==========================================================
vif_raw <- vif(cox_multi1)

if (is.matrix(vif_raw)) {
  vif_values <- vif_raw[, "GVIF^(1/(2*Df))"]
} else {
  vif_values <- vif_raw
}
tabla_vif <- tibble(
  Variable = names(vif_values),
  VIF      = round(as.numeric(vif_values), 2)
) %>%
  mutate(
    Interpretación = case_when(
      VIF < 2  ~ "Sin colinealidad relevante",
      VIF < 5  ~ "Colinealidad moderada",
      VIF < 10 ~ "Colinealidad alta",
      TRUE     ~ "Problema serio de colinealidad"
    )
  )
ft_vif <- flextable(tabla_vif) %>%
  autofit() %>%
  theme_booktabs() %>%
  set_caption("Evaluación de multicolinealidad entre variables del modelo (VIF ajustado)")

doc <- read_docx() %>%
  body_add_par("Evaluación del modelo multivariado de Cox", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_vif)

print(doc, target = "supuestos_sle.docx")
# ==========================================================
# CREAR UNA TABLA (3)
# ==========================================================
# Convertir a data.frame
ph_df <- as.data.frame(ph_test$table) %>%
  tibble::rownames_to_column("Variable") %>%
  rename(
    `Chi-cuadrado` = chisq,
    `gl` = df,
    `p-valor` = p
  )
# Redondear a data.frame
ph_df_round <- ph_df
ph_df_round[, c("Chi-cuadrado", "p-valor")] <-
  round(ph_df_round[, c("Chi-cuadrado", "p-valor")], 2)
# Crear la tabla con flextable
ft_ph <- flextable(ph_df_round) %>%
  set_header_labels(
    Variable = "Variable",
    `Chi-cuadrado` = "χ²",
    gl = "gl",
    `p-valor` = "p-valor"
  ) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "Variable", align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    caption = "Evaluación del supuesto de riesgos proporcionales mediante residuos de Schoenfeld"
  )
# exportar a word
doc <- read_docx() %>%
  body_add_par(
    "Evaluación del supuesto de riesgos proporcionales mediante residuos de Schoenfeld",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft_ph)
print(doc, target = "hr_sle.docx")


vars_cont <- my_data %>%
  select(edad_anios_exacto, log2_leuco, plaq_100, log2_urico)

cor_matrix <- cor(vars_cont, use = "complete.obs")
print(cor_matrix)




