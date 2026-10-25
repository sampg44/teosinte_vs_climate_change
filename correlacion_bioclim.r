# ==============================================================================
# correlacion_bioclim.r
# 
# Correlación entre las 19 variables bioclimáticas, con los valores de los
# puntos de muestreo (un valor por población).
# ==============================================================================

library(ggplot2)

#  1. rutas y así 
ruta_bioclim    <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/bioclim_por_pop.csv"
# es un csv porque el .tif de recortadas deja algunos pixeles conmás d euna pobalción y aquí se usa un csv que tiene un valor por pob aunque se repitan algunos pixeles
carpeta_salida  <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/bioclim/"
metodo          <- "pearson"   # "pearson" (relación lineal) o "spearman" (relación monótona, por rangos)
umbral_r        <- 0.7         # |r| >= umbral se considera muy relacionadas (regla práctica de Dormann et al. 2013)
dir.create(carpeta_salida, recursive = TRUE, showWarnings = FALSE)

morado <- "darkorchid3"; amarillo <- "gold"
tema <- theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = "gray20"))

# 2. leer archivo de variables (columnas bio 1 a 19, filas son poblaciones)
bio <- read.csv(ruta_bioclim, stringsAsFactors = FALSE)
cols_bio <- grep("^bio[0-9]+$", colnames(bio), value = TRUE)
cols_bio <- cols_bio[order(as.integer(sub("bio", "", cols_bio)))]   # bio1, bio2, ..., bio19
cat("Variables encontradas:", length(cols_bio), ", poblaciones (filas):", nrow(bio), "\n")
X <- bio[, cols_bio]
cat("Filas con algún NA:", sum(!complete.cases(X)), "\n")

#  3. correlación

# r(a,b) = covarianza(a,b) / (desv.est(a) * desv.est(b)), sobre las filas (poblaciones)

r <- cor(X, use = "pairwise.complete.obs", method = metodo)
write.csv(round(r, 4), file.path(carpeta_salida, "correlacion_bioclim_poblaciones.csv"))

#  4. visualización 

cl <- expand.grid(v1 = cols_bio, v2 = cols_bio, stringsAsFactors = FALSE)
cl$r <- as.vector(r[cbind(match(cl$v1, cols_bio), match(cl$v2, cols_bio))])
cl$v1 <- factor(cl$v1, levels = rev(cols_bio))   # filas: bio1 arriba
cl$v2 <- factor(cl$v2, levels = cols_bio)        # columnas: bio1 a la izquierda
cl$txt <- ifelse(abs(cl$r) > 0.7, "white", "gray1")
g <- ggplot(cl, aes(v2, v1, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", r), colour = txt), size = 3.2) +
  scale_colour_identity() +
  scale_fill_gradient2(low = amarillo, mid = "#f1f0ee", high = morado, midpoint = 0,
                       limits = c(-1, 1), name = "r") +
  scale_x_discrete(position = "top") + coord_equal() +
  labs(title = "Correlación entre variables bioclimáticas", x = NULL, y = NULL) + tema
ggsave(file.path(carpeta_salida, "correlacion_bioclim_poblaciones.png"), g,
       width = 11, height = 10, dpi = 150, bg = "white")

#  5. pares muy relacionadas y grupos

# Agrupamiento jerárquico con enlace completo: dentro de cada grupo, todos los
# pares tienen |r| >= umbral_r. 
# elegir cuál conservar de cada grupo es decisión ya mía

grupos <- cutree(hclust(as.dist(1 - abs(r)), method = "complete"), h = 1 - umbral_r)
tabla_grupos <- data.frame(variable = cols_bio, grupo = grupos[cols_bio])
tabla_grupos$n_en_grupo <- ave(tabla_grupos$grupo, tabla_grupos$grupo, FUN = length)
tabla_grupos <- tabla_grupos[order(tabla_grupos$grupo, match(tabla_grupos$variable, cols_bio)), ]
write.csv(tabla_grupos, file.path(carpeta_salida, "grupos_bioclim_correlacionadas.csv"), row.names = FALSE)
cat("\nGrupos de variables con |r| >=", umbral_r, "entre todos sus pares:\n")
for (g_id in unique(tabla_grupos$grupo)) {
  v <- tabla_grupos$variable[tabla_grupos$grupo == g_id]
  cat(sprintf("  grupo %d: %s\n", g_id, paste(v, collapse = ", ")))
}
cat("\noutputs guardados en:", carpeta_salida, "\n")
