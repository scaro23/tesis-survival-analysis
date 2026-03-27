# ==========================================================
#  1️⃣ Load libraries
# ==========================================================
library(tidyr);
library(flextable);
library(officer);
library(openxlsx);
library(dplyr);
library(survival);
library(survminer);
library(viridis);
library(ggsci);
library(survMisc)
# ==========================================================
# 2️⃣ Kaplan–Meier
# ==========================================================
fit_km <- survfit(surv_EFS ~ gen_fusion, data = my_data)
summary(fit_km)
fit_ten <- ten(fit_km)
comp(fit_ten)
# ==========================================================
# 3️⃣ Log-rank tests
# ==========================================================
log_rank <- survdiff(surv_EFS ~ aneuplo, data = my_data)
log_rank
# --- Manual P-value Extraction ---
chisq_stat <- log_rank$chisq
df <- length(log_rank$n) - 1 # Degrees of freedom = (Groups - 1)
p_val <- 1 - pchisq(chisq_stat, df)
cat("Log-Rank Test Result:\n")
cat("Chi-Square:", round(chisq_stat, 2), "\n")
cat("P-value:", round(p_val, 4), "\n")
# ==========================================================
# 4️⃣ Calculate the la median survival by group
# ==========================================================
median_surv <- surv_median(fit_km)
median_surv
# ==========================================================
# 5️⃣ Kaplan–Meier curves
# ==========================================================
my_data$gen_fusion <- as.factor(my_data$gen_fusion)
levels(my_data$gen_fusion)
# change of the name of the levels
mis_etiquetas <- c("BCR::ABL1", 
                   "ETV6::RUNX1", 
                   "KMT2A::AFF1", 
                   "NEGATIVO", 
                   "TCF3::PBX1")
# --- 5.1 Generate the graph ---
plot_km <- ggsurvplot(
  fit_km,
  data = my_data,
  # --- 5.2 Line styles ---
  size = 2,              
  linetype = "solid",     
  censor.shape = 124,      
  censor.size = 3,         
  # --- 5.3 Intervals and statistics ---
  conf.int = FALSE,        
  surv.median.line = "none",
  # --- 5.4 P-value ---
  pval = TRUE,             
  pval.size = 4,
  pval.coord = c(0, 0.1), 
  # --- 5.5 Colors and theme ---
  palette = c("#B22222", "#4682B4", "#2E8B57", "#D2691E", "#696969"),
  ggtheme = theme_bw(base_size = 14) +             
    theme(
      panel.grid = element_blank(),                
      axis.line = element_line(size = 0.4),        
      legend.position = "right",                   
      # ---- 5.6 Modify the legends of the graph ----
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 13),
      legend.key.size = unit(1.2, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold")
    ),
  # --- 5.7 Risk Table ---
  risk.table = TRUE,
  risk.table.fontsize = 3.5,
  risk.table.height = 0.25,
  risk.table.title = "Pacientes en riesgo",
  risk.table.y.text = FALSE,
  # --- 5.8 Labels in the graph ---
  title = "Curvas de supervivencia por Gen de fusión",
  subtitle = 'Supervivencia Libre de Eventos: Cohorte YYYY - YYYY',
  xlab = "Tiempo (meses)",
  ylab = "Probabilidad de supervivencia",
  legend.labs = mis_etiquetas,
  legend.title = "Gen de fusión",
)
print(plot_km)
# ==========================================================
# 6️⃣ Export the graphic
# ==========================================================
jpeg(filename = "Curve_###_##.jpg",
     width = 12,      
     height = 10,
     units = "in",
     res = 300)
print(plot_km)
dev.off()
