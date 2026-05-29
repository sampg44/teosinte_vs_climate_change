# =============================================================================
# pi_por_pixel.R
# Objetivo: calcular diversidad nucleotídica (π) por píxel de WorldClim
#           y extraer las 19 variables bioclimáticas para construir la
#           matriz final individuos-pixeles vs variables+pi.
#
# Samantha Melissa Pacheco Gómez
# Requiere: SNPRelate, terra, dplyr
# =============================================================================

rm(list=ls())
# --- 0. Paquetes -------------------------------------------------------------

library(SNPRelate, lib.loc = "/home/sam/R/x86_64-pc-linux-gnu-library/4.5")
library(terra)      
library(dplyr)


# --- 1. Rutas ----------------------------------------------------------------

RUTA_META  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"
RUTA_FAM   <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.fam"
RUTA_BED   <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.bed"
RUTA_BIM   <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.bim"
RUTA_GDS   <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/teosinte.gds"
# TIF de Mesoamérica: un solo archivo multibanda con las 19 variables bio.
RUTA_MESO  <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif"


RUTA_OUT   <- "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_pi_bio.csv"


# --- 2. Parámetro de resolución ----------------------------------------------

RESOLUCION <- "30s"

# Mínimo de individuos por píxel para incluirlo en el análisis.
# Con 1 individuo π = heterocigosidad intraindividual (válido pero diferente
# interpretación). Subir el umbral si se quieren estimaciones más robustas.
MIN_IND_POR_PIXEL <- 10  


# --- 3. Metadata -------------------------------------------------------------

Data <- read.csv(RUTA_META, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Convertir coordenadas a numérico (en el csv están como character/factor)
Data$Latitude  <- as.numeric(as.character(Data$Latitude))
Data$Longitude <- as.numeric(as.character(Data$Longitude))

# Eliminar individuos sin coordenadas
n_antes <- nrow(Data)
Data <- Data[!is.na(Data$Latitude) & !is.na(Data$Longitude), ]
cat("Individuos con coordenadas:", nrow(Data), "de", n_antes, "\n")

# NOTA sobre taxones: los datos mezclan 7 taxones (Zea mays ssp., Z. diploperennis,
# Z. luxurians, Z. nicaraguensis, Z. perennis). Las poblaciones con prefijo "Z"
# (ZDNA, ZDJA, ZLOX…) tienen π mucho más bajo porque son especies distintas.
# filtrar solo parviglumis y mexicana:
Data <- Data[Data$Taxon %in% c("Zea mays ssp. parviglumis",
                                "Zea mays ssp. mexicana"), ]
cat("Taxones en el análisis:\n")
print(table(Data$Taxon))


# --- 4. Asignar individuos a píxeles -----------------------------------------

# El TIF de Mesoamérica es un solo archivo multibanda (19 bandas = bio1..bio19).
if (!file.exists(RUTA_MESO)) {
  stop("No se encontró el TIF de Mesoamérica en: ", RUTA_MESO,
       "\nDescárgalo del cluster antes de continuar (ver instrucciones al final).")
}

bio_stack  <- rast(RUTA_MESO)                          # SpatRaster con 19 bandas
names(bio_stack) <- paste0("bio", 1:nlyr(bio_stack))  # nombrar bio1 a bio19
cat("Bandas en el raster de Mesoamérica:", nlyr(bio_stack), "\n")

ref_raster <- bio_stack[[1]]   # primera banda como referencia de cuadrícula

# Coordenadas como matriz lon, lat (el orden importa en terra)
coords_mat <- cbind(lon = Data$Longitude, lat = Data$Latitude)

# cellFromXY: devuelve el número de celda del píxel donde cae cada individuo.
# Individuos fuera del extent del raster devuelven NA.
Data$cell_id <- cellFromXY(ref_raster, coords_mat)

n_fuera <- sum(is.na(Data$cell_id))
if (n_fuera > 0) {
  warning(n_fuera, " individuos fuera del extent del raster; se excluyen.")
  Data <- Data[!is.na(Data$cell_id), ]
}

cat("Píxeles únicos ocupados:", length(unique(Data$cell_id)), "\n")


# --- 5. Unir .fam con metadata -----------------------------------------------

fam <- read.table(RUTA_FAM, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

# IID del .fam == Sample_name en la metadata
meta_fam <- merge(fam, Data[, c("Sample_name", "cell_id", "Taxon",
                                "POB_CODE", "Latitude", "Longitude")],
                  by.x = "IID", by.y = "Sample_name",
                  all.x = TRUE)

cat("Individuos con match completo:", sum(!is.na(meta_fam$cell_id)), "\n")


# --- 6. Abrir archivo GDS ---------------------------------------------------
# uso archivo gds existente de mis análisis previos

if (!file.exists(RUTA_GDS)) {
  cat("Convirtiendo PLINK -> GDS (una sola vez)...\n")
  snpgdsBED2GDS(bed.fn = RUTA_BED,
                bim.fn = RUTA_BIM,
                fam.fn = RUTA_FAM,
                out.gdsfn = RUTA_GDS)
  cat("GDS creado en:", RUTA_GDS, "\n")
} else {
  cat("Usando GDS existente:", RUTA_GDS, "\n")
}

genofile <- snpgdsOpen(RUTA_GDS)


# --- 7. obtener el orden de sample.id en el GDS -----------------------------

samp_gds <- read.gdsn(index.gdsn(genofile, "sample.id"))

# cell_id para cada individuo en el orden del GDS
cell_gds <- meta_fam$cell_id[match(samp_gds, meta_fam$IID)]


# --- 8. calcular π por pixel -----------------------------------

calc_pi <- function(genofile, sample_ids) {
  # sample_ids: vector de IDs de individuos en ese pixel
  # Devuelve π (diversidad nucleotídica media sobre todos los SNPs)
  
  if (length(sample_ids) < 1) return(NA_real_)
  
  af <- snpgdsSNPRateFreq(genofile, sample.id = sample_ids,
                          with.id = FALSE)
  # af$AlleleFreq : frecuencia del alelo A1 (referencia) por SNP
  # af$MissingRate: fracción de genotipos faltantes por SNP
  
  p      <- af$AlleleFreq
  n_ind  <- length(sample_ids)
  
  # n_haplotipos con dato (diploides: 2 por individuo)
  nh <- 2 * n_ind * (1 - af$MissingRate)
  
  # Corrección por tamaño finito de muestra
  # (importante cuando hay pocos individuos por pixel)
  correction <- ifelse(nh > 1, nh / (nh - 1), NA_real_)
  
  pi_por_locus <- correction * 2 * p * (1 - p)
  
  mean(pi_por_locus, na.rm = TRUE)
}


# --- 9. Loop por pixel: calcular π ------------------------------------------

pixeles_unicos <- unique(cell_gds[!is.na(cell_gds)])

# filtrar pixeles con suficientes individuos
n_ind_por_pixel <- table(cell_gds)
pixeles_ok <- as.integer(names(n_ind_por_pixel[n_ind_por_pixel >= MIN_IND_POR_PIXEL]))

cat("Píxeles con >=", MIN_IND_POR_PIXEL, "individuo(s):", length(pixeles_ok), "\n")

# tabla de resultados
res_pi <- data.frame(
  cell_id = pixeles_ok,
  n_ind   = as.integer(n_ind_por_pixel[as.character(pixeles_ok)]),
  pi      = NA_real_
)

cat("Calculando π por pixel...\n")
t0 <- proc.time()

for (i in seq_len(nrow(res_pi))) {
  cid  <- res_pi$cell_id[i]
  idx  <- which(cell_gds == cid)
  ids  <- samp_gds[idx]
  res_pi$pi[i] <- calc_pi(genofile, ids)
  
  # progreso cada 10 píxeles
  if (i %% 10 == 0) {
    elapsed <- round((proc.time() - t0)[3], 1)
    cat("  Pixel", i, "de", nrow(res_pi), "| tiempo:", elapsed, "s\n")
  }
}

snpgdsClose(genofile)
cat("π calculado para", nrow(res_pi), "pixeles.\n")
print(summary(res_pi$pi))


# --- 10. extraer coordenadas del centro de cada pixel -----------------------

# xyFromCell devuelve las coordenadas lon/lat del centro del píxel
coords_pixel <- as.data.frame(xyFromCell(ref_raster, res_pi$cell_id))
colnames(coords_pixel) <- c("lon_pixel", "lat_pixel")
res_pi <- cbind(res_pi, coords_pixel)


# --- 11. Extraer las 19 variables bioclimáticas ----------------------------

# definimos el vector de nombres de las variables bioclimáticas
nombres_bio <- names(bio_stack)

# extraer valores en los centros de los píxeles de interés
coords_extraccion <- cbind(res_pi$lon_pixel, res_pi$lat_pixel)
bio_vals <- extract(bio_stack, coords_extraccion)
# bio_vals: data.frame con columnas bio1..bio19, una fila por pixel

# verificar NAs en las variables bioclimáticas
n_na_bio <- sum(is.na(bio_vals$bio1))
if (n_na_bio > 0) {
  warning(n_na_bio, " pixeles sin valor en bio1 (probablemente en mar o borde).")
}


# --- 12. Construir la matriz final ------------------------------------------

matriz_final <- cbind(res_pi, bio_vals)

# reordenar columnas para claridad
# filas = píxeles, columnas = cell_id, n_ind, lon, lat, π, bio1..bio19
matriz_final <- matriz_final[, c("cell_id", "n_ind", "lon_pixel", "lat_pixel",
                                 "pi", nombres_bio)]

cat("\nDimensiones de la matriz final:", dim(matriz_final), "\n")
cat("Primeras filas:\n")
print(head(matriz_final))

cat("\nResumen de π:\n")
print(summary(matriz_final$pi))

cat("\nNAs por columna:\n")
print(colSums(is.na(matriz_final)))


# --- 13. Guardar ------------------------------------------------------------


write.csv(matriz_final, RUTA_OUT, row.names = FALSE)

save(
  matriz_final,
  file = "/home/sam/Documents/sur_ecoevo_lab/results/mayo_2026/matriz_pi_bio.RData"
)

cat("\nMatriz guardada en CSV y RData.\n")

# =============================================================================
# exploración
# =============================================================================

# Ver distribución de individuos por pixel:
 hist(matriz_final$n_ind, main="Distribución de ind. por pixel",
      xlab = "Número de individuos por pixel",
      ylab= "Frecuencia",
      col = "#9DE0E6",
      border = "white")



# Ver distribución de pi:
 hist(matriz_final$pi, main="Distribución de π por pixel",
      xlab = "polimorfismo (pi)",
      ylab= "Frecuencia",
      col = "#9DE0E6",
      border = "white")

 
# Ver relación pi ~ bio1 (temperatura media anual):
 plot(matriz_final$bio1, matriz_final$pi,
      main= "Relación de π respecto a la Bio1",
      xlab="BIO1 (Temperatura media anual, °C×10)",
      ylab="π (diversidad nucleotídica o polimorfismo)", pch=19, cex=0.6, col="steelblue")
 abline(lm(pi ~ bio1, data=matriz_final), col="red")

# =============================================================================
# NOTAS PARA EL ANÁLISIS
# =============================================================================
# 1. RESOLUCIÓN: este script usa los TIFs de 30s (~1 km).
#
# 2. MIN_IND_POR_PIXEL: con 1 individuo, π = proporción de sitios heterocigotos
#    (diversidad intraindividual). Con ≥2 ya incluye diversidad entre individuos.
#
# 3. TAXONES MEZCLADOS: si incluiste todos los taxones, un pixel cerca de la
#    frontera de Z. diploperennis podría tener π muy bajo simplemente por ser
#    otra especie, no por el clima. 
#     por eso filtré solo parviglumis y mexicana
#
# 4. ESTANDARIZACIÓN: antes de correr los modelos, estandariza todas las
#    variables de entrada (bio1..bio19) con scale().
#
# 5. CORRELACIÓN ENTRE BIO VARS: las 19 variables de WorldClim están altamente
#    correlacionadas entre sí. Antes del modelado, considerar PCA o selección
#    de variables  para evitar multicolinealidad.
# =============================================================================


#
# Verificar que el archivo descargado tiene las 19 bandas esperadas en R:
  library(terra)
  r <- rast("~/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif")
  nlyr(r)   # debe imprimir 19
  names(r)  # muestra los nombres actuales de cada banda
# =============================================================================