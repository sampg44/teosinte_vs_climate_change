# =============================================================================
# shap.R
# Objetivo: calcular e interpretar SHAP values del modelo final de red neuronal
# =============================================================================



# --- 0. Paquetes -------------------------------------------------------------

library(kernelshap)
library(shapviz)
library(keras3)
library(ggplot2)
library(dplyr)

# --- 1. Rutas ----------------------------------------------------------------

RUTA_RES <- "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026"

# --- 2. Cargar modelo y objetos ---------------------------------------------

model_final <- load_model(
  file.path(RUTA_RES, "model_final.keras")
)

load(
  file.path(RUTA_RES, "modelo_final_objetos.RData")
)

cat("Modelo y objetos cargados.\n")
cat("Mejor configuración:", mejor_config, "\n")
cat("Variables usadas:", paste(vars_final, collapse = ", "), "\n")

# Verificación importante: que las variables coincidan con los parámetros guardados
cat("Número de variables en vars_final:", length(vars_final), "\n")
cat("Número de medias guardadas:", length(x_means_final), "\n")
cat("Número de sds guardadas:", length(x_sds_final), "\n")

if (length(vars_final) != length(x_means_final)) {
  stop("Inconsistencia: vars_final no coincide con x_means_final.")
}

if (length(vars_final) != length(x_sds_final)) {
  stop("Inconsistencia: vars_final no coincide con x_sds_final.")
}

# --- 3. Preparar X para SHAP -------------------------------------------------

X_shap_raw <- as.data.frame(scores[, vars_final])

# --- 4. Función de predicción ------------------------------------------------

pred_fun <- function(model, newdata) {
  
  X_new <- as.matrix(newdata)
  
  X_new_s <- scale(
    X_new,
    center = x_means_final,
    scale  = x_sds_final
  )
  
  pred_std <- as.vector(predict(model, X_new_s, verbose = 0))
  
  pred_pi <- pred_std * y_sd + y_mean
  
  return(pred_pi)
}

# --- 5. Calcular SHAP --------------------------------------------------------

set.seed(44)

# Si tarda mucho, puedes usar una muestra de fondo:
# bg_idx <- sample(seq_len(nrow(X_shap_raw)), size = min(80, nrow(X_shap_raw)))
# bg_data <- X_shap_raw[bg_idx, ]

# Si no tarda mucho, usa todos los datos como background:
bg_data <- X_shap_raw

shap_result <- kernelshap(
  object   = model_final,
  X        = X_shap_raw,
  bg_X     = bg_data,
  pred_fun = pred_fun,
  verbose  = TRUE
)

sv <- shapviz(shap_result)

save(
  shap_result,
  sv,
  X_shap_raw,
  file = file.path(RUTA_RES, "shap_result.RData")
)

cat("SHAP calculado y guardado.\n")

# --- 6. Importancia global ---------------------------------------------------

importancia_shap <- data.frame(
  variable = names(colMeans(abs(shap_result$S))),
  mean_abs_shap = as.numeric(colMeans(abs(shap_result$S)))
) |>
  arrange(desc(mean_abs_shap))

print(importancia_shap)

write.csv(
  importancia_shap,
  file.path(RUTA_RES, "shap_importancia_global.csv"),
  row.names = FALSE
)

# --- 7. Gráficas SHAP --------------------------------------------------------

p_bar <- sv_importance(sv, kind = "bar") +
  labs(
    title = "Importancia global de variables",
    subtitle = "Media de valores absolutos de SHAP",
    x = "Importancia media absoluta",
    y = "Variable"
  )

p_bar


p_beeswarm <- sv_importance(sv, kind = "beeswarm") +
  labs(
    title = "Efecto de cada variable sobre π",
    subtitle = "Color = valor de la variable"
  )

p_beeswarm

# --- 8. Dependence plots para las variables más importantes ------------------

top_vars <- importancia_shap$variable[1:min(3, nrow(importancia_shap))]

for (v in top_vars) {
  
  p_dep <- sv_dependence(sv, v) +
    labs(
      title = paste("SHAP dependence:", v),
      x = v,
      y = "SHAP value"
    )
  
  ggsave(
    filename = file.path(RUTA_RES, paste0("shap_dependence_", v, ".png")),
    plot = p_dep,
    width = 6,
    height = 4,
    dpi = 300
  )
}

# --- 9. Guardar SHAP values como tabla --------------------------------------

shap_df <- as.data.frame(shap_result$S)

shap_df$Taxon   <- scores$Taxon
shap_df$pi      <- scores$pi
shap_df$cell_id <- scores$cell_id

write.csv(
  shap_df,
  file.path(RUTA_RES, "shap_values.csv"),
  row.names = FALSE
)

cat("SHAP values guardados.\n")
cat("Resumen de importancia global:\n")
importancia_shap$mean_abs_shap <- round(importancia_shap$mean_abs_shap, 5)
print(importancia_shap)
