# 
# PCA sobre las 19 variables bioclimáticas (WorldClim) de las
# poblaciones de teosinte
#
# 27 ago 2026
#
# tomo la tabla que produce 1_extraer_valores_ambientales.R 
#
# procesito
#   1. carga la tabla de valores crudos
#   2. estandariza cada variable (z-score) ANTES del PCA
#   3. corre el PCA
#   4. guarda varianza explicada, loadings y scores
#
# no aplica la fórmula de "diferenciación ambiental" de BayeScEnv (|valor - referencia| / sd). 
# esa se aplicará después, sobre el componente elegido, no sobre las 19 variables
#


# ============
# 1. rutas
# ==========

ruta_valores_crudos  <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/valores_ambientales_crudos.csv"
ruta_salida_loadings <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/pca_19_bio/loadings_pca_19_bio.csv"
ruta_salida_varianza <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/pca_19_bio/varianza_explicada_pca_19_bio.csv"
ruta_salida_scores   <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/pca_19_bio/scores_pca_19_bio.csv"

nombres_variables <- paste0("bio", 1:19)  #debe coincidir con el script 1


# =================
# 2. cargar datos
# ===============

tabla <- read.csv(ruta_valores_crudos, stringsAsFactors = FALSE)

faltantes <- apply(is.na(tabla[, nombres_variables]), 1, any)
if (any(faltantes)) {
  message(sum(faltantes), " población(es) con valores faltantes en al menos ",
          "una variable -- el PCA las va a excluir por completo.")
  print(tabla[faltantes, c("pop_index", "Accession")])
}

completos <- tabla[!faltantes, ]


# ============================================================
# 3. estandarizar (z-score) y correr PCA
# ============================================================

X <- as.matrix(completos[, nombres_variables])
pca <- prcomp(X, center = TRUE, scale. = TRUE)

scores <- pca$x  # poblaciones x componentes -- ya estandarizados por construcción


# ============================================================
# 4. guardar resultados
# ============================================================

n_comp <- length(nombres_variables)
nombres_pc <- paste0("PC", 1:n_comp)

varianza_explicada <- (pca$sdev^2) / sum(pca$sdev^2) * 100
varianza <- data.frame(
  componente = nombres_pc,
  varianza_explicada_pct = varianza_explicada,
  varianza_acumulada_pct = cumsum(varianza_explicada)
)
write.csv(varianza, ruta_salida_varianza, row.names = FALSE)

# loadings: cuánto pesa cada variable original en cada componente
loadings <- as.data.frame(pca$rotation)
loadings <- cbind(variable = rownames(loadings), loadings)
write.csv(loadings, ruta_salida_loadings, row.names = FALSE)

# scores: el valor de cada población en cada componente, EN EL
# MISMO ORDEN que orden_pop.csv (porque nunca se reordenó la tabla)
tabla_scores <- cbind(completos[, c("pop_index", "Accession")], as.data.frame(scores))
colnames(tabla_scores) <- c("pop_index", "Accession", nombres_pc)
write.csv(tabla_scores, ruta_salida_scores, row.names = FALSE)

message("Varianza explicada guardada en: ", ruta_salida_varianza)
message("Loadings guardados en: ", ruta_salida_loadings)
message("Scores por población guardados en: ", ruta_salida_scores)

message("\nVarianza explicada por los primeros 5 componentes:")
print(head(varianza, 5))

message("\nLoadings de PC1 y PC2 (para empezar a interpretar qué representan):")
print(loadings[order(-abs(loadings$PC1)), c("variable", "PC1", "PC2")])


# ============
# 5. chequeo 
# ===========

message("Media de cada componente (debe ser ~0 por construcción del PCA):")
print(round(colMeans(scores)[1:3], 6))
message("Desvío estándar de cada componente (PC1 debe ser el mayor, decreciendo):")
print(round(apply(scores, 2, sd)[1:5], 4))

