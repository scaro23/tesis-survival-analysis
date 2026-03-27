# ===================================================================
# LIBRERIAS
# ====================================================================
# Análisis de supervivencia
library(survival)
library(cmprsk)
# Manipulación y lectura de datos
library(dplyr)
library(readr)
library(openxlsx)
library(lubridate)
library(purrr)
# Modelado y outputs ordenados
library(broom)
# ===================================================================
# OVERALL SURVIVAL (OS)
# ===================================================================
my_data <- my_data %>%
  mutate(
    # ===============================
    # 🔹 Status de supervivencia OS
    # ===============================
    status_OS = case_when(
      Estado_vital == "FALLECIDO" ~ 1,   # evento
      Estado_vital %in% c("VIVO", "EXTRANJERO") ~ 0),
    # ===============================
    # 🔹 Tiempo de supervivencia
    # ===============================
    fecha_evento_OS = if_else(
      status_OS == 1,
      fecha_ultima,    # fallecido → fecha de muerte
      fecha_galenos    # censurado → última fecha de contacto
    ),
    time_OS_months = as.numeric(fecha_evento_OS - fecha_dx) / 30.44)
# ===============================
# 🔹 Objeto de supervivencia (Overall survival)
# ===============================
surv_OS <- with(my_data, Surv(time_OS_months, status_OS))
cox_model1 <- coxph(surv_OS ~ gen_fusion,data = my_data) # ejemplo
cox.zph(cox_model1)
# ===================================================================
# EVENT FREE SURVIVAL (EFS)
# ===================================================================
my_data <- my_data %>%
  mutate(
    fecha_fallo_induccion = if_else(
      respuesta_induc == "Refractaria / persistente",
      fecha_dx + 33,
      as.Date(NA)
    ),
    fecha_recaida_evento = if_else(
      recaida == "SI" & !is.na(fecha_recaida),
      fecha_recaida,
      as.Date(NA)
    ),
    fecha_muerte_evento = if_else(
      Estado_vital == "FALLECIDO",
      fecha_ultima,
      as.Date(NA)
    ),
    fecha_evento_real = pmin(
      fecha_fallo_induccion,
      fecha_recaida_evento,
      fecha_muerte_evento,
      na.rm = TRUE
    ),
    fecha_evento_real = if_else(
      is.infinite(fecha_evento_real),
      as.Date(NA),
      fecha_evento_real
    ),
    # 🔹 Fecha real de censura
    fecha_censura_real = case_when(
      abandono == 1 ~ diferimiento_i,
      TRUE ~ fecha_galenos
    ),
    # 🔹 Aplicar censura correctamente
    status_EFS = if_else(
      is.na(fecha_evento_real) | fecha_evento_real > fecha_censura_real,
      0L, 1L
    ),
    fecha_evento_EFS = if_else(
      status_EFS == 1,
      fecha_evento_real,
      fecha_censura_real
    ),
    tipo_evento_EFS = case_when(
      status_EFS == 0 ~ "Censura",
      !is.na(fecha_fallo_induccion) & fecha_evento_real == fecha_fallo_induccion ~ "Fallo de inducción",
      !is.na(fecha_recaida_evento) & fecha_evento_real == fecha_recaida_evento ~ "Recaída",
      !is.na(fecha_muerte_evento) & fecha_evento_real == fecha_muerte_evento ~ "Muerte"
    ),
    time_EFS_months = as.numeric(fecha_evento_EFS - fecha_dx) / 30.44)
# ===============================
# 🔹 Objeto de supervivencia (Event free survival)
# ===============================
surv_EFS <- with(my_data, Surv(time_EFS_months, status_EFS))
cox_model2 <- coxph(surv_EFS ~ gen_fusion,data = my_data) # ejemplo
cox.zph(cox_model2)

# ==========================================================
# CUMULATIVE COMPETING RISK (CCR)
# ==========================================================
my_data_ccr <- my_data %>%
  # 🔹 Excluir fallos de inducción (solo pacientes en remisión)
  filter(respuesta_induc != "Refractaria / persistente") %>%
  mutate(
    # --------------------------------------------
    # 1 = Recaída (evento de interés)
    # 2 = Muerte sin recaída (evento competitivo)
    # 0 = Censura
    # --------------------------------------------
    evento_CCR = case_when(
      recaida == "SI" ~ 1,
      Estado_vital == "FALLECIDO" & recaida != "SI" ~ 2,
      TRUE ~ 0
    ),
    # 🔹 Fecha real del evento
    fecha_evento_CCR = case_when(
      evento_CCR == 1 ~ fecha_recaida,
      evento_CCR == 2 ~ fecha_ultima,
      TRUE ~ fecha_censura_real
    ),
    # 🔹 Tiempo en días
    tiempo_CCR_dias = as.numeric(fecha_evento_CCR - fecha_dx)
  ) %>%
  filter(!is.na(tiempo_CCR_dias) & tiempo_CCR_dias >= 0)