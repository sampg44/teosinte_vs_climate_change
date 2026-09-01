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
print(loadings[order(-abs(loadings$PC1)), c("PC1", "PC2")])


# ============
# 5. chequeo 
# ===========

message("Media de cada componente (debe ser ~0 por construcción del PCA):")
print(round(colMeans(scores)[1:3], 6))
message("Desvío estándar de cada componente (PC1 debe ser el mayor, decreciendo):")
print(round(apply(scores, 2, sd)[1:5], 4))


# 6. gráficas
# ===============

# si falta alguno: install.packages(c("ggplot2", "tidyr"))
library(ggplot2)
library(tidyr)   # solo para pasar los loadings de formato ancho a largo

ruta_varianza <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/pca_19_bio/varianza_explicada_pca_19_bio.csv"
ruta_loadings <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/pca_19_bio/loadings_pca_19_bio.csv"


varianza <- read.csv(ruta_varianza, stringsAsFactors = FALSE)

varianza$componente <- factor(varianza$componente, levels = varianza$componente)

# las barras (% individual, máximo ~46 en tus datos) y la línea
# (% acumulado, siempre llega a 100) comparten el mismo eje si no se
# hace nada -- eso aplasta las barras contra abajo. "factor_escala"
# reescala la línea para que ocupe visualmente el mismo rango que las
# barras; sec_axis() se encarga de que el eje derecho siga mostrando
# el % acumulado real (0-100), no el valor reescalado.
max_individual <- max(varianza$varianza_explicada_pct)
factor_escala <- max_individual / 100


color_barras <- "lightgreen"
color_linea  <- "mediumpurple1"

scree_plot <- ggplot(varianza, aes(x = componente)) +

  geom_col(aes(y = varianza_explicada_pct, fill = "% individual"), alpha = 0.85) +
  geom_line(aes(y = varianza_acumulada_pct * factor_escala, color = "% acumulado", group = 1),
            linewidth = 1) +
  geom_point(aes(y = varianza_acumulada_pct * factor_escala, color = "% acumulado"), size = 2) +
  scale_y_continuous(
    name = "% individual",
    sec.axis = sec_axis(~ . / factor_escala, name = "% acumulado")
  ) +
  scale_fill_manual(name = "", values = c("% individual" = color_barras)) +
  scale_color_manual(name = "", values = c("% acumulado" = color_linea)) +
  labs(title = "Varianza explicada por componente", x = "Componente") +
  theme_minimal() +
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5))

print(scree_plot)  


#heatmap de loadings (variables x componentes)


n_componentes_mostrar <- 8  


loadings <- read.csv(ruta_loadings, stringsAsFactors = FALSE)

nombres_pc_mostrar <- paste0("PC", 1:n_componentes_mostrar)

# pivot_longer: pasa de "una columna por componente" (formato ancho,
# como está guardado el CSV) a "una fila por combinación
# variable-componente" (formato largo), que es lo que necesita
# ggplot para poder colorear un heatmap por el valor de loading
loadings_largo <- pivot_longer(
  loadings[, c("variable", nombres_pc_mostrar)],
  cols = all_of(nombres_pc_mostrar),
  names_to = "componente", values_to = "loading"
)


loadings_largo$variable <- factor(loadings_largo$variable, levels = paste0("bio", 1:19))
loadings_largo$componente <- factor(loadings_largo$componente, levels = nombres_pc_mostrar)

heatmap_loadings <- ggplot(loadings_largo, aes(x = componente, y = variable, fill = loading)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(loading, 2)), size = 3) +
  # gradiente divergente centrado en 0: azul = negativo, rojo = positivo,
  # blanco = cerca de cero -- así se ve de un vistazo qué variables
  # cargan fuerte (de cualquier signo) en cada componente
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Loadings del PCA ambiental", x = "", y = "", fill = "Loading") +
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))

print(heatmap_loadings)

