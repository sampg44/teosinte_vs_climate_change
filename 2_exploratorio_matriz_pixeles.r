
# =============================================================================
# exploratorio_matriz_pixeles.R
# Objetivo: entender la distribución de π y las 19 variables bioclimáticas,
#           detectar correlaciones entre predictores, y decidir qué entra
#           a la red neuronal.
#
# Decisiones ya tomadas:
#   - Solo Zea mays ssp. parviglumis y mexicana
#   - Resolución 30s
#   - Umbral mínimo 10 individuos por píxel  →  ~204 observaciones (pixeles)
#
# Samantha Melissa Pacheco Gómez
# =============================================================================

rm(list=ls())

# --- 0. Paquetes -------------------------------------------------------------
# install.packages(c("corrplot", "ggplot2", "car", "terra", "dplyr"))

library(dplyr)
library(ggplot2)
library(corrplot)
library(car)      
library(terra)    


# --- 1. Cargar y filtrar la matriz -------------------------------------------

mat_raw <- read.csv(
  "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_pi_bio.csv"
)

# Aplicar filtros definitivos
mat <- mat_raw[mat_raw$n_ind >= 10, ]

cat("Observaciones después de filtrar (n_ind >= 10):", nrow(mat), "\n")
cat("Columnas disponibles:", paste(colnames(mat), collapse = ", "), "\n\n")

BIO_COLS <- paste0("bio", 1:19)


# =============================================================================
# BLOQUE 1: distribución de π
# =============================================================================

cat("--- Resumen de π ---\n")
print(summary(mat$pi))
cat("CV de π (desv/media):", round(sd(mat$pi) / mean(mat$pi), 3), "\n\n")

# Histograma de π
hist(mat$pi,
     col    = "steelblue",
     border = "white",
     main   = "Distribución de π por píxel",
     xlab   = "π (diversidad nucleotídica)",
     ylab   = "Frecuencia")

# ¿Está sesgada? Si la cola es larga, log(π) puede ayudar al modelo.
# no ayudóporque esa transformación es para arreglar colas hacia la derecha y la tenía a la izq
# Compara el histograma con y sin log.
hist(log(mat$pi),
     #breaks = 30,
     col    = "steelblue3",
     border = "white",
     main   = "Distribución de log(π)",
     xlab   = "log(π)",
     ylab   = "Frecuencia")


# ¿Parviglumis y mexicana tienen π distintos?
# 
# # Cargar metadata original
# Data <- read.csv(
#   "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv",
#   header = TRUE, sep = ",", stringsAsFactors = FALSE
# )
# Data$Latitude  <- as.numeric(as.character(Data$Latitude))
# Data$Longitude <- as.numeric(as.character(Data$Longitude))
# 
# # Cargar matriz ya filtrada (n_ind >= 10, parv+mex)
# mat <- read.csv(
#   "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_pi_bio.csv"
# )
# mat <- mat[mat$n_ind >= 10, ]
# 
# # Para cada pixel, obtener el taxon mayoritario de sus individuos
# library(terra)
# ref <- rast("/home/sam/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif")[[1]]
# 
# Data$cell_id <- cellFromXY(ref, cbind(Data$Longitude, Data$Latitude))
# 
# taxon_por_pixel <- aggregate(
#   Taxon ~ cell_id,
#   data  = Data[!is.na(Data$cell_id), ],
#   FUN   = function(x) names(sort(table(x), decreasing = TRUE))[1]
# )
# 
# mat <- merge(mat, taxon_por_pixel, by = "cell_id", all.x = TRUE)
# 
# # Ahora sí:
# tapply(mat$pi, mat$Taxon, summary)
# 
# boxplot(pi ~ Taxon, data = mat,
#         col  = c("steelblue", "tomato"),
#         main = "π por subespecie",
#         ylab = "π", xlab = "")
# 
# 
# # Conteo de individuos por taxon en cada pixel
# library(dplyr)
# 
# taxon_conteo <- Data %>%
#   filter(!is.na(cell_id), cell_id %in% mat$cell_id) %>%
#   count(cell_id, Taxon) %>%
#   tidyr::pivot_wider(names_from = Taxon, 
#                      values_from = n, 
#                      values_fill = 0) %>%
#   rename(n_mexicana   = `Zea mays ssp. mexicana`,
#          n_parviglumis = `Zea mays ssp. parviglumis`)
# 
# # Agregar a mat
# mat <- merge(mat, taxon_conteo, by = "cell_id", all.x = TRUE)
# 
# # ¿Cuántos píxeles son "puros" vs mixtos?
# mat$tipo_pixel <- case_when(
#   mat$n_mexicana == 0                    ~ "solo parviglumis",
#   mat$n_parviglumis == 0                 ~ "solo mexicana",
#   mat$n_mexicana > mat$n_parviglumis     ~ "mayoría mexicana",
#   mat$n_parviglumis > mat$n_mexicana     ~ "mayoría parviglumis",
#   TRUE                                   ~ "empate"
# )
# 
# table(mat$tipo_pixel)


# agregar taxon como variable

# =============================================================================
# BLOQUE 2: ¿n_ind confunde la estimación de π?
# =============================================================================

# Si los píxeles con más individuos tienen sistemáticamente más o menos π,
# el tamaño de muestra está sesgando tu variable respuesta.

cat("--- Correlación n_ind ~ π ---\n")
cor_test <- cor.test(mat$n_ind, mat$pi, method = "spearman")
print(cor_test)
cat("ρ de Spearman:", round(cor_test$estimate, 3),
    " | p-valor:", round(cor_test$p.value, 4), "\n\n")

# Visualización
plot(mat$n_ind, mat$pi,
     pch  = 19, cex = 0.7, col = "steelblue",
     xlab = "n individuos por píxel",
     ylab = "π",
     main = "¿Más individuos = más π?")
abline(lm(pi ~ n_ind, data = mat), col = "tomato", lwd = 2)

# Interpretación:
# |ρ| < 0.2  → n_ind no confunde, umbral ≥10 es suficiente.
# |ρ| 0.2-0.4 → efecto moderado, menciónalo como limitación.
# |ρ| > 0.4  → hay sesgo real; considera subir el umbral o corregir π.

# esta ok el umbral, no confunde
# =============================================================================
# BLOQUE 3: distribución de cada variable bioclimática
# =============================================================================

# Resumen rápido de las 19 bios
cat("--- Resumen de variables bioclimáticas ---\n")
print(summary(mat[, BIO_COLS]))

# Histogramas en una sola figura (4 filas × 5 columnas)
par(mfrow = c(4, 5), mar = c(3, 2, 2, 1))
for (bio in BIO_COLS) {
  hist(mat[[bio]],
       breaks = 15,
       col    = "slategray3",
       border = "white",
       main   = bio,
       xlab   = "",
       ylab   = "")
}
par(mfrow = c(1, 1))  # restablecer layout

# Lo que buscas: ¿alguna variable tiene distribución muy rara (bimodal,
# casi constante, con outliers extremos)? Esas merecen atención especial.

# la bimodalidad en las variables de temperatura no es un problema que haya que corregir 
# es exactamente la diferencia climática entre subespecies que el modelo necesita aprender.
# =============================================================================
# BLOQUE 4: correlación entre las 19 variables bioclimáticas
# =============================================================================

# Matriz de correlaciones de Pearson entre bios
cor_bio <- cor(mat[, BIO_COLS], use = "complete.obs")

# Mapa de calor de correlaciones
corrplot(cor_bio,
         method  = "color",
         type    = "upper",
         tl.cex  = 0.75,
         tl.col  = "black",
         addCoef.col = "black",
         number.cex  = 0.45,
         col     = colorRampPalette(c("steelblue", "white", "tomato"))(200),
         title   = "Correlación entre variables bioclimáticas",
         mar     = c(0, 0, 1.5, 0))

# Las 19 variables de WorldClim suelen agruparse en ~3-4 bloques:
#   • Temperatura media y derivadas (bio1, bio2, bio3, bio5, bio6, bio7, bio10, bio11)
#   • Estacionalidad de temperatura (bio4, bio7)
#   • Precipitación total (bio12, bio13, bio14, bio16, bio17, bio18, bio19)
#   • Estacionalidad de precipitación (bio15)
# Si ves correlaciones > 0.9 entre pares, hay redundancia severa.


# =============================================================================
# BLOQUE 5: PCA para después reducir dimesion de las variables
# =============================================================================

# PCA sobre las 19 bios
pca_result <- prcomp(mat[, paste0("bio", 1:19)], scale. = TRUE)

# Varianza explicada acumulada — tabla y scree plot juntos
var_exp <- summary(pca_result)$importance

# Tabla limpia
tabla_pca <- data.frame(
  PC            = paste0("PC", 1:10),
  SD            = round(var_exp[1, 1:10], 3),
  Var_propia    = round(var_exp[1, 1:10]^2, 3),
  Prop_varianza = round(var_exp[2, 1:10] * 100, 1),
  Acumulada     = round(var_exp[3, 1:10] * 100, 1)
)
print(tabla_pca)

# Scree plot con varianza acumulada superpuesta
par(mar = c(4, 4, 3, 4))
vp <- var_exp[2, 1:10] * 100
vc <- var_exp[3, 1:10] * 100

bp <- barplot(vp,
              names.arg = paste0("PC", 1:10),
              col  = "hotpink1",
              border = "white",
              ylim = c(0, 100),
              ylab = "% varianza explicada",
              main = "PCA — varianza por componente y acumulada",
              las  = 1)

# Línea de varianza acumulada
lines(bp, vc, col = "aquamarine3", lwd = 2, type = "b", pch = 19)

# Líneas de referencia
abline(h = 80, lty = 2, col = "mediumorchid3")
abline(h = 90, lty = 2, col = "mediumpurple4")

axis(side = 4, las = 1)
mtext("% acumulado", side = 4, line = 3)

# Leyenda reubicada manualmente
legend(x = bp[8], y = 78,
       legend = c("% por componente", "% acumulado", "80%", "90%"),
       col    = c("hotpink1", "aquamarine3", "mediumorchid3", "mediumpurple4"),
       lty    = c(NA, 1, 2, 2),
       pch    = c(15, 19, NA, NA),
       bty    = "n",
       cex    = 0.85,
       pt.cex = 1.1,
       y.intersp = 0.9)

# Contribución de cada bio a PC1, PC2, PC3
loadings_df <- as.data.frame(round(pca_result$rotation[, 1:8], 3))
loadings_df$bio <- rownames(loadings_df)
loadings_df <- loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ]
print(loadings_df)

# =============================================================================
# BLOQUE 6: correlación de cada bio con π (relaciones bivariadas)
# =============================================================================

cor_con_pi <- sapply(BIO_COLS, function(bio) {
  cor(mat[[bio]], mat$pi, use = "complete.obs", method = "spearman")
})

cor_pi_df <- data.frame(
  variable    = names(cor_con_pi),
  rho_con_pi  = round(cor_con_pi, 3)
) |> arrange(desc(abs(rho_con_pi)))

cat("\n--- Correlación Spearman de cada bio con π (ordenada por |ρ|) ---\n")
print(cor_pi_df)

# Gráfico de barras de correlaciones con π
barplot(cor_pi_df$rho_con_pi,
        names.arg = cor_pi_df$variable,
        las       = 2,
        col       = ifelse(cor_pi_df$rho_con_pi > 0, "#9DE0E6", "mediumorchid4"),
        main      = "Correlación Spearman: bio ~ π",
        ylab      = "ρ de Spearman",
        cex.names = 0.75)
abline(h = 0, lwd = 1)


# =============================================================================
# BLOQUE 7: scatterplots de las 3-4 bios más correlacionadas con π
# =============================================================================

# Tomar las 4 variables con mayor |ρ| con π
top4 <- cor_pi_df$variable[1:4]

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (bio in top4) {
  plot(mat[[bio]], mat$pi,
       pch  = 19, cex = 0.7, col = "#9DE0E6",
       xlab = bio,
       ylab = "π",
       main = paste("π ~", bio))
  lines(lowess(mat[[bio]], mat$pi), col = "mediumorchid4", lwd = 2)
}
par(mfrow = c(1, 1))

# La línea roja (lowess) muestra la tendencia sin asumir linealidad.
# Si es curva → la relación no es lineal, la red neuronal es especialmente
# apropiada para capturarla.


# =============================================================================
# BLOQUE 8: mapa de π sobre Mesoamérica
# =============================================================================

library(terra)
library(scales)  # para alpha()

# Cargar raster solo para el extent
ref <- rast("/home/sam/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif")[[1]]

# Definir extent ajustado a tus datos con un poco de margen
ext_zoom <- ext(-108, -92, 14, 23)
ref_zoom  <- crop(ref, ext_zoom)

# Paleta de colores para π
pal      <- hcl.colors(100, "YlOrRd", rev = TRUE)
pi_range <- range(mat$pi, na.rm = TRUE)
col_idx  <- findInterval(mat$pi,
                         seq(pi_range[1], pi_range[2], length.out = 101),
                         all.inside = TRUE)
pt_cols  <- pal[col_idx]

# Graficar
plot(ref_zoom,
     col     = gray.colors(50, start = 0.85, end = 0.95),
     legend  = FALSE,
     main    = "Diversidad nucleotídica (π) por píxel en Mesoamérica",
     xlab    = "Longitud (°O)",
     ylab    = "Latitud (°N)",
     mar     = c(3, 3, 2, 6))

# Puntos coloreados por π
points(mat$lon_pixel, mat$lat_pixel,
       pch = 19, cex = 0.9, col = pt_cols)

# Leyenda manual
legend_image <- as.raster(matrix(rev(pal), ncol = 1))

rasterImage(legend_image,
            xleft = -93.3, ybottom = 14.7,
            xright = -92.9, ytop = 19.3)

text(x = -92.75,
     y = c(14.7, 17.0, 19.3),
     labels = round(c(pi_range[1], mean(pi_range), pi_range[2]), 3),
     cex = 0.7, adj = 0)

text(x = -93.1, y = 19.7, labels = "π", cex = 0.8)


# =============================================================================
# BLOQUE 9: preparar matriz completa para modelo con ia
# =============================================================================
# Recuperar Taxon en mat 

Data <- read.csv(
  "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv",
  header = TRUE, sep = ",", stringsAsFactors = FALSE
)
Data$Latitude  <- as.numeric(as.character(Data$Latitude))
Data$Longitude <- as.numeric(as.character(Data$Longitude))

library(terra)
ref      <- rast("/home/sam/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif")[[1]]
Data$cell_id <- cellFromXY(ref, cbind(Data$Longitude, Data$Latitude))

taxon_por_pixel <- aggregate(
  Taxon ~ cell_id,
  data = Data[!is.na(Data$cell_id), ],
  FUN  = function(x) names(sort(table(x), decreasing = TRUE))[1]
)

mat <- merge(mat, taxon_por_pixel, by = "cell_id", all.x = TRUE)

# Verificar que quedó
cat("¿Taxon en mat?", "Taxon" %in% colnames(mat), "\n")
table(mat$Taxon)

# Preparar matriz completa con PC1-PC6 (más que suficiente para las 3 opciones)
scores <- as.data.frame(pca_result$x[, 1:6])
colnames(scores) <- c("PC1_temperatura",
                      "PC2_estacionalidad", 
                      "PC3_precipitacion",
                      "PC4_rango_diurno",
                      "PC5",
                      "PC6")

scores$es_parviglumis <- as.integer(
  mat$Taxon == "Zea mays ssp. parviglumis"
)
scores$pi     <- mat$pi
scores$log_pi <- log(mat$pi)      # guardamos ambas versiones del target

# Metadata para interpretar resultados después (no entra al modelo)
scores$cell_id    <- mat$cell_id
scores$lon_pixel  <- mat$lon_pixel
scores$lat_pixel  <- mat$lat_pixel
scores$n_ind      <- mat$n_ind
scores$Taxon      <- mat$Taxon



write.csv(scores,
          "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_modelo.csv",
          row.names = FALSE)

save(
  scores,
  mat,
  pca_result,
  file = "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_modelo_pca.RData"
)

cat("Matriz modelo guardada en CSV y RData.\n")



# =============================================================================

