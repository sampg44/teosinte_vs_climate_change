# ==============================================================================
# scatter_snp_bioclim.r
# 24 oct 2026

# este script ni lo tengo corrido pero si ya gasté agua pidiendo ayuda a la IA pues lo voy a guardar
# y también me da miedo que me regañen entonces mejor documentarlo
# cada punto es una población.
#   eje x esfrecuencia del alelo contado en esa población
#   eje y es valor de la variable bioclimática en esa población
# muestra los SNPs más significativos del GLM y, algunos que no lo son
# ==============================================================================
library(ggplot2)

# 1.- rutas y parámetros

ruta_matriz_frecuencia      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_frecuencia_poblacional_253.rds"
ruta_matriz_copias_posibles <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_copias_posibles_253.rds"
ruta_bioclim                <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/bioclim_por_pop.csv"
ruta_metadata               <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"

variable_bioclim <- "bio1"
ruta_glm_csv     <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio1.csv"
carpeta_salida   <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/24_oct/figuras_glm"

n_snps_hits       <- 6     # los N SNPs con menor p-valor
n_snps_contraste  <- 3     # N SNPs claramente no significativos (p > 0.5), elegidos al azar
semilla           <- 1
snps_manual       <- NULL  # c("S1_992727", "S3_123456") para graficar SNPs concretos
dir.create(carpeta_salida, recursive = TRUE, showWarnings = FALSE)

verde <- "#2f9e6e"; lila <- "#9b7fd8"; amarillo <- "#c98500"
tema <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = "#52514e"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "top")

# 2.- leer archivos
frec   <- readRDS(ruta_matriz_frecuencia)
copias <- readRDS(ruta_matriz_copias_posibles)
stopifnot(identical(dimnames(frec), dimnames(copias)))
bio    <- read.csv(ruta_bioclim, stringsAsFactors = FALSE)
meta   <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
stopifnot(variable_bioclim %in% colnames(bio))
glm    <- read.csv(ruta_glm_csv, stringsAsFactors = FALSE)

valor_bio <- bio[[variable_bioclim]][match(rownames(frec), bio$Accession)]
taxon     <- meta$Taxon[match(rownames(frec), meta$Accession)]
grupo <- ifelse(taxon == "Zea mays ssp. parviglumis", "parviglumis",
                ifelse(taxon == "Zea mays ssp. mexicana", "mexicana", "otros taxones"))
grupo <- factor(grupo, levels = c("parviglumis", "mexicana", "otros taxones"))

# 3.- elegir snps 
glm <- glm[is.finite(glm$p_valor) & glm$SNP %in% colnames(frec), ]
glm$p_bonf <- pmin(1, glm$p_valor * sum(is.finite(glm$p_valor)))
if (is.null(snps_manual)) {
  top <- head(glm$SNP[order(glm$p_valor)], n_snps_hits)
  set.seed(semilla)
  candidatos <- glm$SNP[glm$p_valor > 0.5]
  contraste <- sample(candidatos, min(n_snps_contraste, length(candidatos)))
  elegidos <- c(top, contraste)
} else {
  elegidos <- snps_manual
}
fila_glm <- glm[match(elegidos, glm$SNP), ]

# 4.- tabla una fila por población x SNP 
largo <- do.call(rbind, lapply(seq_along(elegidos), function(i) {
  s <- elegidos[i]
  data.frame(SNP = s, poblacion = rownames(frec), frecuencia = frec[, s],
             copias = copias[, s], bio = valor_bio, grupo = grupo,
             panel = sprintf("%s\ncoef = %.2f | p Bonferroni = %s", s, fila_glm$coef_bioclim[i],
                             formatC(fila_glm$p_bonf[i], format = "g", digits = 2)),
             stringsAsFactors = FALSE)
}))
largo$panel <- factor(largo$panel, levels = unique(largo$panel))
largo <- largo[!is.na(largo$frecuencia) & !is.na(largo$bio), ]

# 5.- figura
g <- ggplot(largo, aes(frecuencia, bio, colour = grupo, shape = grupo, size = copias)) +
  geom_point(alpha = 0.75) +
  facet_wrap(~ panel, ncol = 3) +
  scale_colour_manual(values = c(parviglumis = verde, mexicana = lila, `otros taxones` = amarillo), name = NULL) +
  scale_shape_manual(values = c(parviglumis = 16, mexicana = 17, `otros taxones` = 15), name = NULL) +
  scale_size_area(max_size = 3.5, name = "Copias\nposibles") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(title = paste0(variable_bioclim, " contra frecuencia alélica en cada población"),
       x = "Frecuencia del alelo contado en la población", y = variable_bioclim) + tema
alto <- 2.4 * ceiling(length(elegidos) / 3) + 1.2
ggsave(file.path(carpeta_salida, paste0("scatter_", variable_bioclim, "_snps.png")), g,
       width = 11, height = alto, dpi = 150, bg = "white")
cat("Figura guardada en:", carpeta_salida, "\n")
