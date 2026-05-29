# =============================================================================
# red_neuronal.R
# Objetivo: entrenar una red neuronal para predecir π (diversidad nucleotídica)
#           a partir de componentes principales de variables bioclimáticas,
#           comparando 3 configuraciones de inputs con k-fold CV (k=10).
#
# Configuraciones:
#   A: PC1-PC3 + es_parviglumis  (87.7% varianza bio)
#   B: PC1-PC4 + es_parviglumis  (91.9% varianza bio)
#   C: PC1-PC5 + es_parviglumis  (95.3% varianza bio)
#
# Samantha Melissa Pacheco Gómez
# =============================================================================


# --- 0. Instalación (solo la primera vez) ------------------------------------
# keras usa TensorFlow como backend. Instalar una vez:
#install.packages("keras3")
rm(list=ls())

# --- 1. Paquetes -------------------------------------------------------------

library(keras3)
library(dplyr)

# --- 2. Cargar datos ---------------------------------------------------------

scores <- read.csv(
  "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_modelo.csv"
)

cat("Dimensiones:", nrow(scores), "filas x", ncol(scores), "columnas\n")
cat("Distribución de subespecies:\n")
print(table(scores$Taxon))
cat("Resumen de π:\n")
print(summary(scores$pi))


# --- 3. Definir las 3 configuraciones ----------------------------------------

configs <- list(
  A = c("PC1_temperatura", "PC2_estacionalidad",
        "PC3_precipitacion", "es_parviglumis"),
  
  B = c("PC1_temperatura", "PC2_estacionalidad",
        "PC3_precipitacion", "PC4_rango_diurno",
        "es_parviglumis"),
  
  C = c("PC1_temperatura", "PC2_estacionalidad",
        "PC3_precipitacion", "PC4_rango_diurno",
        "PC5", "es_parviglumis")
)

# Variable respuesta
y_raw <- scores$pi



# --- 4. Estandarizar ---------------------------------------------------------
# IMPORTANTE: estandarizamos X e y POR SEPARADO.
# Los parámetros de estandarización se calculan en el set de entrenamiento
# de cada fold y se aplican al set de validación — nunca al revés.
# Aquí guardamos los parámetros globales para el modelo final.

# Variable respuesta sin estandarizar
y_raw <- scores$pi

# y_mean <- mean(y_raw)
# y_sd   <- sd(y_raw)
# y_std  <- (y_raw - y_mean) / y_sd   # π estandarizado




# --- 5. Función para construir la red ----------------------------------------


  # Arquitectura conservadora para n=204:
  #   - 2 capas ocultas con pocas neuronas (evita overfitting)
  #   - Dropout 0.2 como regularización
  #   - Activación ReLU en capas ocultas
  #   - Activación lineal en la salida (regresión)
  #   - Optimizador Adam (eficiente, robusto a learning rate)
  #   - Loss: MSE (equivalente a minimizar RMSE)
  

# 1er modelo tenía overfitting pipipi

# build_model <- function(n_inputs) {
#   
#   inputs <- keras_input(shape = n_inputs)
#   
#   outputs <- inputs |>
#     layer_dense(units = 32, activation = "relu") |>
#     layer_dropout(rate = 0.2) |>
#     layer_dense(units = 16, activation = "relu") |>
#     layer_dropout(rate = 0.2) |>
#     layer_dense(units = 1, activation = "linear")
#   
#   model <- keras_model(inputs = inputs, outputs = outputs)
#   
#   model |> compile(
#     optimizer = optimizer_adam(learning_rate = 0.001),
#     loss      = "mse",
#     metrics   = list("mae")
#   )
#   
#   return(model)
# }


# intento 2 tenía underfitting

# build_model <- function(n_inputs) {
#   
#   inputs <- keras_input(shape = n_inputs)
#   
#   outputs <- inputs |>
#     layer_dense(units = 16, activation = "relu",
#                 kernel_regularizer = regularizer_l2(0.01)) |>
#     layer_dropout(rate = 0.3) |>
#     layer_dense(units = 8, activation = "relu",
#                 kernel_regularizer = regularizer_l2(0.01)) |>
#     layer_dropout(rate = 0.3) |>
#     layer_dense(units = 1, activation = "linear")
#   
#   model <- keras_model(inputs = inputs, outputs = outputs)
#   
#   model |> compile(
#     optimizer = optimizer_adam(learning_rate = 0.001),
#     loss      = "mse",
#     metrics   = list("mae")
#   )
#   
#   return(model)
# }

# intento 3 tenía un mejor r square pero aun con algo de overfitting
build_model <- function(n_inputs) {

  inputs <- keras_input(shape = n_inputs)

  outputs <- inputs |>
    layer_dense(units = 24, activation = "relu",
                kernel_regularizer = regularizer_l2(0.005)) |>
    layer_dropout(rate = 0.25) |>
    layer_dense(units = 12, activation = "relu",
                kernel_regularizer = regularizer_l2(0.005)) |>
    layer_dropout(rate = 0.25) |>
    layer_dense(units = 1, activation = "linear")

  model <- keras_model(inputs = inputs, outputs = outputs)

  model |> compile(
    optimizer = optimizer_adam(learning_rate = 0.0005),
    loss      = "mse",
    metrics   = list("mae")
  )

  return(model)
}

# perceptron simple no fue suficiente
# build_model <- function(n_inputs) {
#   inputs  <- keras_input(shape = n_inputs)
#   outputs <- inputs |>
#     layer_dense(units = 1, activation = "linear")  # sin capas ocultas
#   
#   model <- keras_model(inputs = inputs, outputs = outputs)
#   model |> compile(
#     optimizer = optimizer_adam(learning_rate = 0.001),
#     loss      = "mse", metrics = list("mae")
#   )
#   return(model)
# }

#modelo red sencillo no aprendió nada predice todo cerca de la media siempre
# build_model <- function(n_inputs) {
#   inputs  <- keras_input(shape = n_inputs)
#   outputs <- inputs |>
#     layer_dense(units = 8, activation = "relu",
#                 kernel_regularizer = regularizer_l2(0.005)) |>
#     layer_dense(units = 1, activation = "linear")
#   
#   model <- keras_model(inputs = inputs, outputs = outputs)
#   model |> compile(
#     optimizer = optimizer_adam(learning_rate = 0.001),
#     loss      = "mse", metrics = list("mae")
#   )
#   return(model)
# }

# # Prueba rápida con configuración B para probar que los codigos de cada modelo sirven

# vars_prueba <- c("PC1_temperatura", "PC2_estacionalidad",
#                  "PC3_precipitacion", "PC4_rango_diurno",
#                  "es_parviglumis")
# X_prueba <- scale(as.matrix(scores[, vars_prueba]))
# y_prueba <- scale(scores$pi)[, 1]
# 
# model_test <- build_model(5)
# summary(model_test)
# 
# history_test <- model_test |> fit(
#   x = X_prueba, y = y_prueba,
#   epochs = 10, batch_size = 32,
#   verbose = 1
# )

# --- 6. K-fold CV (k=10) para cada configuración ----------------------------

set.seed(44)
k       <- 10
n       <- nrow(scores)
fold_id <- sample(rep(1:k, length.out = n))   # asignar folds aleatoriamente

# Tabla de resultados
resultados <- data.frame(
  config     = character(),
  fold       = integer(),
  rmse_train = numeric(),
  mae_train  = numeric(),
  r2_train   = numeric(),
  acc_train  = numeric(),
  rmse_val   = numeric(),
  mae_val    = numeric(),
  r2_val     = numeric(),
  acc_val    = numeric(),
  stringsAsFactors = FALSE
)

cat("=== K-FOLD CV (k=10) ===\n\n")

for (cfg_name in names(configs)) {
  
  vars    <- configs[[cfg_name]]
  X_raw   <- as.matrix(scores[, vars])
  cat("--- Configuración", cfg_name,
      "(", paste(vars, collapse=", "), ") ---\n")
  
  fold_rmse <- numeric(k)
  
  for (fold in 1:k) {
    
    # Separar train / validación
    idx_val   <- which(fold_id == fold)
    idx_train <- which(fold_id != fold)
    
    X_train <- X_raw[idx_train, , drop = FALSE]
    X_val   <- X_raw[idx_val,   , drop = FALSE]
    
    # y en escala original
    y_train_raw <- y_raw[idx_train]
    y_val_raw   <- y_raw[idx_val]
    
    # Estandarizar y usando SOLO train
    y_mean_fold <- mean(y_train_raw)
    y_sd_fold   <- sd(y_train_raw)
    
    y_train <- (y_train_raw - y_mean_fold) / y_sd_fold
    y_val   <- (y_val_raw   - y_mean_fold) / y_sd_fold
    
    # Estandarizar X con parámetros del train (evitar data leakage)
    x_means <- colMeans(X_train)
    x_sds   <- apply(X_train, 2, sd)
    x_sds[x_sds == 0] <- 1   # evitar división por cero (es_parviglumis)
    
    X_train_s <- scale(X_train, center = x_means, scale = x_sds)
    X_val_s   <- scale(X_val,   center = x_means, scale = x_sds)
    
    # Construir y entrenar modelo
    model <- build_model(n_inputs = length(vars))
    
    history <- model |> fit(
      x          = X_train_s,
      y          = y_train,
      epochs     = 200,
      batch_size = 32,
      validation_data = list(X_val_s, y_val),
      callbacks  = list(
        callback_early_stopping(
          monitor   = "val_loss",
          patience  = 20,          # detener si no mejora en 20 épocas
          restore_best_weights = TRUE
        )
      ),
      verbose = 1   #0 es silencioso; cambiar a 1 para ver el entrenamiento
    )
    
    # Predicciones en train y validación
    y_pred_train_std <- as.vector(predict(model, X_train_s))
    y_pred_val_std   <- as.vector(predict(model, X_val_s))
    
    # Regresar a escala original de π usando los parámetros del fold
    y_pred_train_pi <- y_pred_train_std * y_sd_fold + y_mean_fold
    y_pred_val_pi   <- y_pred_val_std   * y_sd_fold + y_mean_fold
    
    y_train_pi <- y_train_raw
    y_val_pi   <- y_val_raw
    
    # Métricas TRAIN
    rmse_train <- sqrt(mean((y_train_pi - y_pred_train_pi)^2))
    mae_train  <- mean(abs(y_train_pi - y_pred_train_pi))
    r2_train   <- 1 - sum((y_train_pi - y_pred_train_pi)^2) /
      sum((y_train_pi - mean(y_train_pi))^2)
    
    # Métricas VALIDACIÓN / TEST FOLD
    rmse_val <- sqrt(mean((y_val_pi - y_pred_val_pi)^2))
    mae_val  <- mean(abs(y_val_pi - y_pred_val_pi))
    r2_val   <- 1 - sum((y_val_pi - y_pred_val_pi)^2) /
      sum((y_val_pi - mean(y_val_pi))^2)
    
    # Desempeño porcentual aproximado, "accuracy aproximado"
    # Debido a que el modelo es de regresión, el desempeño principal se evaluó con
    #RMSE, MAE y R². Además, se calculó una métrica porcentual de desempeño basada
    # en el RMSE relativo al rango observado de π, reportada como 
    # “accuracy aproximada” solo para facilitar la interpretación
    
    acc_train <- 100 * (1 - rmse_train / diff(range(y_raw)))
    acc_val   <- 100 * (1 - rmse_val   / diff(range(y_raw)))
    
    
    
    fold_rmse[fold] <- rmse_val
    
    resultados <- rbind(resultados, data.frame(
      config     = cfg_name,
      fold       = fold,
      rmse_train = round(rmse_train, 5),
      mae_train  = round(mae_train, 5),
      r2_train   = round(r2_train, 4),
      acc_train  = round(acc_train, 2),
      rmse_val   = round(rmse_val, 5),
      mae_val    = round(mae_val, 5),
      r2_val     = round(r2_val, 4),
      acc_val    = round(acc_val, 2)
    ))
    
    cat("  Fold", fold,
        "| RMSE train:", round(rmse_train, 5),
        "| RMSE val:", round(rmse_val, 5),
        "| R² train:", round(r2_train, 4),
        "| R² val:", round(r2_val, 4),
        "| acc train:", round(acc_train, 2), "%",
        "| acc val:", round(acc_val, 2), "%\n")
  }
  
  cat("  → RMSE val medio:", round(mean(fold_rmse), 5),
      "± SD:", round(sd(fold_rmse), 5), "\n\n")
}


# --- 7. Comparar configuraciones ---------------------------------------------

cat("=== RESUMEN COMPARATIVO ===\n")
  resumen <- resultados |>
    group_by(config) |>
    summarise(
      RMSE_train_medio = round(mean(rmse_train), 5),
      RMSE_val_medio   = round(mean(rmse_val), 5),
      RMSE_val_sd      = round(sd(rmse_val), 5),
      MAE_train_medio  = round(mean(mae_train), 5),
      MAE_val_medio    = round(mean(mae_val), 5),
      R2_train_medio   = round(mean(r2_train), 4),
      R2_val_medio     = round(mean(r2_val), 4),
      acc_train_media  = round(mean(acc_train), 2),
      acc_val_media    = round(mean(acc_val), 2),
      .groups = "drop"
    ) |>
    arrange(RMSE_val_medio)

print(resumen)

mejor_config <- resumen$config[1]
cat("\nMejor configuración:", mejor_config,
    "(menor RMSE en validación)\n")





# --- 8. Modelo final con la mejor configuración ------------------------------
# Entrena sobre TODOS los datos (sin CV) para el modelo definitivo.



cat("\n=== ENTRENANDO MODELO FINAL (config", mejor_config, ") ===\n")

vars_final  <- configs[[mejor_config]]
X_final_raw <- as.matrix(scores[, vars_final])


# Estandarizar X con todos los datos
x_means_final <- colMeans(X_final_raw)
x_sds_final   <- apply(X_final_raw, 2, sd)
x_sds_final[x_sds_final == 0] <- 1
X_final_s     <- scale(X_final_raw,
                       center = x_means_final,
                       scale  = x_sds_final)

model_final <- build_model(n_inputs = length(vars_final))

# Estandarización final de y usando todos los datos
# Solo para entrenar el modelo definitivo después del CV
y_mean <- mean(y_raw)
y_sd   <- sd(y_raw)
y_std  <- (y_raw - y_mean) / y_sd

history_final <- model_final |> fit(
  x          = X_final_s,
  y          = y_std,
  epochs     = 300,
  batch_size = 32,
  validation_split = 0.1,   # 10% interno para early stopping
  callbacks  = list(
    callback_early_stopping(
      monitor  = "val_loss",
      patience = 30,
      restore_best_weights = TRUE
    )
  ),
  verbose = 1
)

# --------------diagnósticos

# Verificar que los datos están bien
cat("Resumen de X para el modelo:\n")

X_check <- as.matrix(scores[, vars_final])
X_check_s <- scale(X_check)

cat("Medias después de scale:\n")
print(round(colMeans(X_check_s), 4))
cat("SDs después de scale:\n")
print(round(apply(X_check_s, 2, sd), 4))
cat("Proporción parviglumis:", mean(scores$es_parviglumis), "\n")

# Ver si y también está bien estandarizado
y_s <- scale(scores$pi)[,1]
cat("Media y:", round(mean(y_s), 4), 
    "| SD y:", round(sd(y_s), 4), "\n")







# --- 9. Evaluación del modelo final ------------------------------------------


y_pred_final_std <- as.vector(predict(model_final, X_final_s))
y_pred_final_pi  <- y_pred_final_std * y_sd + y_mean


rmse_final_train <- sqrt(mean((y_raw - y_pred_final_pi)^2))
mae_final_train  <- mean(abs(y_raw - y_pred_final_pi))
r2_final_train   <- 1 - sum((y_raw - y_pred_final_pi)^2) /
  sum((y_raw - mean(y_raw))^2)

cat("\nMétricas del modelo final entrenado con todos los datos:\n")
cat("  RMSE entrenamiento completo:", round(rmse_final_train, 5), "\n")
cat("  MAE entrenamiento completo: ", round(mae_final_train, 5), "\n")
cat("  R² entrenamiento completo:  ", round(r2_final_train, 4), "\n")
cat("  Nota: para generalización usar las métricas del CV.\n")

# Curva de entrenamiento
plot(history_final,
     main = paste("Curva de entrenamiento — config", mejor_config))


# --- 10. Diagnósticos visuales -----------------------------------------------

# Predicho vs. real
plot(y_raw, y_pred_final_pi,
     pch  = 19, cex = 0.7,
     col  = ifelse(scores$Taxon == "Zea mays ssp. parviglumis",
                   "mediumorchid3", "aquamarine3"),
     xlab = "π observado",
     ylab = "π predicho",
     main = paste("predicho vs. real"))
abline(0, 1, col = "hotpink2", lty = 2, lwd = 1.5)
legend("topleft",
       legend = c("parviglumis", "mexicana"),
       col    = c("mediumorchid3", "aquamarine3"),
       pch    = 19, bty = "n", cex = 0.85)

# Residuos vs. predicho
residuos <- y_raw - y_pred_final_pi
plot(y_pred_final_pi, residuos,
     pch  = 19, cex = 0.7, col = "aquamarine3",
     xlab = "π predicho",
     ylab = "Residuo (observado − predicho)",
     main = "residuos vs. predicho")
abline(h = 0, col = "mediumorchid3", lty = 2, lwd = 1.5)



# --- 11. Guardar resultados --------------------------------------------------
#resumen comparativo
write.csv(resumen,
          "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/cv_resumen_configuraciones.csv",
          row.names = FALSE)

# Tabla de métricas por fold
write.csv(resultados,
          "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/cv_resultados.csv",
          row.names = FALSE)

# Predicciones del modelo final con metadatos
pred_final <- scores[, c("cell_id", "lon_pixel", "lat_pixel",
                         "n_ind", "Taxon", "pi")]
pred_final$pi_predicho <- round(y_pred_final_pi, 6)
pred_final$residuo     <- round(residuos, 6)

write.csv(pred_final,
          "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/predicciones_finales.csv",
          row.names = FALSE)

# Parámetros de estandarización (necesarios para predecir en datos nuevos)
params_std <- data.frame(
  variable = c("pi_mean", "pi_sd", vars_final,
               paste0(vars_final, "_sd")),
  valor    = c(y_mean, y_sd, x_means_final, x_sds_final)
)
write.csv(params_std,
          "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/params_estandarizacion.csv",
          row.names = FALSE)










# Cuantifica el sesgo en los extremos
pred_final <- read.csv(
  "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/predicciones_finales.csv"
)

# Comparar error en pi bajo vs pi típico
pred_final$grupo <- ifelse(pred_final$pi < 0.19, 
                           "pi bajo (<0.19)",
                           "pi típico (≥0.19)")

pred_final %>%
  group_by(grupo) %>%
  summarise(
    n    = n(),
    RMSE = round(sqrt(mean(residuo^2)), 5),
    sesgo_medio = round(mean(residuo), 5)
  )


print(resumen)



# --- Guardar modelo final y objetos necesarios para SHAP ---------------------

# guardar modelo keras
save_model(
  model_final,
  "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/model_final.keras"
)

# guardar objetos de R necesarios para interpretar y reutilizar el modelo
save(
  scores,
  configs,
  mejor_config,
  vars_final,
  X_final_raw,
  X_final_s,
  x_means_final,
  x_sds_final,
  y_mean,
  y_sd,
  y_std,
  y_raw,
  y_pred_final_pi,
  pred_final,
  resultados,
  resumen,
  history_final,
  file = "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/modelo_final_objetos.RData"
)

cat("Modelo final guardado como .keras\n")
cat("Objetos del modelo guardados como modelo_final_objetos.RData\n")



# =============================================================================