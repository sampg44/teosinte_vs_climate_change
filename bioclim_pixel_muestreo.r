# ============================================================
# bioclim_pixel_muestreo.r
# 12 oct 2026
#
# quedarse con el pixel donde cae el punto de muestreo de cada población
# un pixel (~1 km2) por punto
#
# hacer un .tif con esos pixeles seleccionados y un csv para los modelos lineales, con la misma información
# pero por poblaciÓn
# ============================================================

library(terra)


# ========
# 1. rutas
# ========

ruta_tif_mesoamerica   <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif"
ruta_matriz_frecuencia <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_frecuencia_poblacional_253.rds"
ruta_metadata          <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"

ruta_salida_tif <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/meso_pixeles.tif"
ruta_salida_csv <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/bioclim_por_pop.csv" #porque


# ===========================
# 2. coordenadas de cada pop
# ==========================
# La lista de poblaciones sale de la matriz de frecuencias
# Las coordenadas salen de data_teosinte.csv, cada Accession tiene 1 solo par de coordenadas o se detiene.

matriz_frecuencia <- readRDS(ruta_matriz_frecuencia)
poblaciones <- rownames(matriz_frecuencia)

metadatos <- read.csv(ruta_metadata, stringsAsFactors = FALSE)

n_coordenadas_por_pob <- tapply(paste(metadatos$Latitude, metadatos$Longitude),
                                metadatos$Accession, function(x) length(unique(x)))
stopifnot(all(n_coordenadas_por_pob[poblaciones] == 1))

una_fila_por_pob <- metadatos[!duplicated(metadatos$Accession), ]
coordenadas <- una_fila_por_pob[match(poblaciones, una_fila_por_pob$Accession),
                                c("Accession", "Latitude", "Longitude")]
stopifnot(!anyNA(coordenadas$Accession), !anyNA(coordenadas$Latitude), !anyNA(coordenadas$Longitude))
cat("Poblaciones con coordenadas:", nrow(coordenadas), "\n")

# terra pide (x, y) = (LONGITUD, latitud) -- al reves de "lat, lon"
puntos <- as.matrix(data.frame(x = as.numeric(coordenadas$Longitude),
                               y = as.numeric(coordenadas$Latitude)))


# ================
# 3. leer el tif
# ===============

bio <- rast(ruta_tif_mesoamerica)

cat("\nnúmero de capas del tif: ", nlyr(bio), " con nombres:\n")
print(names(bio))
cat("Tienen una extensión de:", as.vector(ext(bio)), "(xmin, xmax, ymin, ymax)\n")

# los nombres deben terminar en el numero de la variable (ej. "..._bio_12" o "bio12")
numero_de_capa <- suppressWarnings(as.integer(sub("^.*bio_?([0-9]+)$", "\\1", names(bio))))
if (anyNA(numero_de_capa) || !setequal(numero_de_capa, 1:19)) {
  stop("No se puede identificar las capas a partir del nombre de la variable ")
}
bio <- bio[[order(numero_de_capa)]]       # ordenar por núm de variable
names(bio) <- paste0("bio", 1:19)


# ============================================================
# 4. en que pixel cae cada punto
# ============================================================

celda_de_cada_punto <- cellFromXY(bio, puntos)
fuera <- poblaciones[is.na(celda_de_cada_punto)]
if (length(fuera) > 0) {
  stop(length(fuera), " poblaciones caen fuera del tif: ", paste(fuera, collapse = ", "))
}

#utilizar el centro del pixel en vez de la coord original
centros_de_pixel <- xyFromCell(bio, celda_de_cada_punto)

# poblaciones que comparten pixel (reciben exactamente los mismos valores)
repetidas <- celda_de_cada_punto[duplicated(celda_de_cada_punto)]
if (length(repetidas) > 0) {
  cat("\n[AVISO] poblaciones que caen en el mismo pixel:\n")
  for (celda in unique(repetidas)) {
    cat("  ", paste(poblaciones[celda_de_cada_punto == celda], collapse = " + "), "\n")
  }
}


# =====================================
# 5. dejar valores solo en esos pixeles
# =====================================

# recortar el mapa a la caja que contiene todos los puntos (con un pixel de margen) para no cargar mas de lo necesario.
margen <- max(res(bio)) #tamaÑo del pixel en grados
caja <- ext(min(centros_de_pixel[, 1]) - margen, max(centros_de_pixel[, 1]) + margen,
            min(centros_de_pixel[, 2]) - margen, max(centros_de_pixel[, 2]) + margen)
bio_caja <- crop(bio, caja) #cortar el mapa a ese rectÁngulo imaginario


#hacer una plantilla vacÍa como si fuera papel calca y marcar con 1 los pixeles donde cae algun punto y lo demÁs es NA.
plantilla <- rast(bio_caja, nlyrs = 1)
values(plantilla) <- NA
plantilla[cellFromXY(bio_caja, centros_de_pixel)] <- 1

# todo lo que no este marcado en la plantilla se vuelve NA en el original (creo)
bio_pixeles <- mask(bio_caja, plantilla)


# ============================================================
# 6. revisar el resultado y guardar
# ============================================================

#checar que sÍ hayan 249 pixeles distintos
n_pixeles_con_dato <- global(!is.na(bio_pixeles[["bio1"]]), "sum")[1, 1]
cat("\npixeles con valor (bio1):", n_pixeles_con_dato,
    ", pixeles distintos donde caen los puntos:", length(unique(celda_de_cada_punto)), "\n")

#leer valores por pop (tabla de 253 filas en teorÍa)
valores <- extract(bio_pixeles, centros_de_pixel)
valores <- valores[, paste0("bio", 1:19), drop = FALSE]

#por si hay alguna poblaciÓn sin dato
sin_dato <- poblaciones[rowSums(is.na(valores)) > 0]
if (length(sin_dato) > 0) {
  cat("\n", length(sin_dato), "poblaciones quedaron con NA en alguna capa (pixel sin dato, ej. costa o mar):\n  ",
      paste(sin_dato, collapse = ", "), "\n")
}

#guardar el tif
dir.create(dirname(ruta_salida_tif), recursive = TRUE, showWarnings = FALSE)
writeRaster(bio_pixeles, ruta_salida_tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
cat("Guardado:", ruta_salida_tif, "\n")

#guardar la tabla
tabla_bioclim <- data.frame(coordenadas, celda_en_tif = celda_de_cada_punto, valores, row.names = NULL)
write.csv(tabla_bioclim, ruta_salida_csv, row.names = FALSE)
cat("Guardado:", ruta_salida_csv, "\n")

