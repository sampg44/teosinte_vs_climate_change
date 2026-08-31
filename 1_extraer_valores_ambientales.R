# ____________________________________________________________________________________________
# Extraer valores ambientales de WorldClim en las coordenadas de las 276 poblaciones de teosinte 
#
# 27 ago 2026
# ___________________________________________________________________________________________
# proceso
#   1. leer orden_pop.csv (pop_index -> Accession).
#   2. leer coordenadas (Latitude, Longitude, Altitude) por Accession
#      desde data_teosinte.csv 
#   3. abrir el raster .tif (19 bandas bioclimáticas) y extraer el
#      valor de cada banda en cada coordenada
#   4. guardar una tabla de referencia: una fila por población, en
#      el MISMO orden que orden_pop.csv, con los valores crudos
#      (sin transformar) de las 19 variables.
#
# Fuente de coordenadas: data_teosinte.csv

# (passport_data.xlsx le falta una Accession (CIM9479) que sí está en data_teosinte.csv)
# ============================================================

# si falta alguno: install.packages(c("terra", "readxl", "dplyr"))
library(terra)    # para leer el raster y extraer valores puntuales
library(readxl)   # para leer passport_data.xlsx (verificación opcional)
library(dplyr)


# _____________rutas____________________

ruta_orden_pop <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/orden_pop.csv"
ruta_metadata  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"        
ruta_passport  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/passport_data.xlsx"  #opcional solo por si acaso
ruta_tif       <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/mesoamerica_30s.tif" 
ruta_salida    <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/valores_ambientales_crudos.csv"


nombres_variables <- paste0("bio", 1:19)  # <-- CONFIRMAR

# Correcciones manuales de Accession, por si en el futuro aparece algún caso de nombre que no aparezca ni siquiera en data_teosinte.csv.
# formato: c("Accession_en_orden_pop" = "Accession_en_data_teosinte")
correcciones_manuales_accession <- c()


# =========================================
# 2. cargar orden canónico de poblaciones
# ===================================

orden <- read.csv(ruta_orden_pop, stringsAsFactors = FALSE)
orden$Accession_norm <- trimws(as.character(orden$Accession))

# ============================================================
# 3. Cargar coordenadas y unir preservando el orden de orden_pop.csv
# ============================================================

meta <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
meta$Accession_norm <- trimws(as.character(meta$Accession))

# chequeo: ¿alguna Accession tiene más de una coordenada distinta
# entre sus individuos? si la hay, hay que investigarlo a mano, NO
# promediar coordenadas ciegamente !!!
chequeo_var <- meta %>%
  group_by(Accession_norm) %>%
  summarise(n_lat = n_distinct(Latitude), n_lon = n_distinct(Longitude))
inconsistentes <- chequeo_var %>% filter(n_lat > 1 | n_lon > 1)
if (nrow(inconsistentes) > 0) {
  message(nrow(inconsistentes), " Accession(es) tienen más de una coordenada ",
          "distinta entre sus individuos en ", ruta_metadata,
          ". Revisar a mano, NO se está promediando automáticamente:")
  print(inconsistentes)
}

meta_pop <- meta %>%
  group_by(Accession_norm) %>%
  summarise(
    Latitude  = first(Latitude),
    Longitude = first(Longitude),
    Altitude  = first(Altitude)
  )

orden$Accession_busqueda <- orden$Accession_norm
if (length(correcciones_manuales_accession) > 0) {
  idx <- match(orden$Accession_busqueda, names(correcciones_manuales_accession))
  reemplazo <- correcciones_manuales_accession[idx]
  orden$Accession_busqueda[!is.na(reemplazo)] <- reemplazo[!is.na(reemplazo)]
}

tabla <- left_join(orden, meta_pop,
                    by = c("Accession_busqueda" = "Accession_norm"))

faltantes <- tabla[is.na(tabla$Latitude), ]
if (nrow(faltantes) > 0) {
  message(nrow(faltantes), " población(es) de orden_pop.csv NO tienen ",
          "coordenadas en ", ruta_metadata, ":")
  print(faltantes[, c("pop_index", "Accession")])
  message("Esto no debería pasar (data_teosinte.csv es la fuente de la que ",
          "salió orden_pop.csv) -- si aparece, revisar con cuidado.")
}


# ============================================================
# 3b. Verificación cruzada OPCIONAL contra passport_data.xlsx
# ============================================================

verificar_contra_passport <- function(tabla, ruta_passport) {
  if (!file.exists(ruta_passport)) {
    message("No se encontró ", ruta_passport, ", se omite la verificación cruzada.")
    return(invisible(NULL))
  }
  passport <- read_excel(ruta_passport)
  passport$Accession_norm <- trimws(as.character(passport$Accession))
  
  # poblaciones que no se pudieron comparar (su Accession no aparece
  # con el mismo nombre exacto en passport_data.xlsx -- no significa
  # que falten, solo que ese archivo usa otro string para esa Accession)
  
  sin_comparar <- anti_join(tabla, passport[, "Accession_norm", drop = FALSE],
                            by = "Accession_norm")
  if (nrow(sin_comparar) > 0) {
    message(nrow(sin_comparar), " población(es) no se pudieron comparar porque su ",
            "Accession no aparece con el mismo nombre en ", ruta_passport, ":")
    print(sin_comparar[, c("pop_index", "Accession")])
  }

  cruce <- inner_join(tabla, passport[, c("Accession_norm", "LATITUDE", "LONGITUDE")],
                       by = "Accession_norm")
  if (nrow(cruce) == 0) return(invisible(NULL))

  # distancia haversine, en km
  R_tierra <- 6371
  lat1 <- cruce$Latitude * pi / 180;  lon1 <- cruce$Longitude * pi / 180
  lat2 <- cruce$LATITUDE  * pi / 180; lon2 <- cruce$LONGITUDE  * pi / 180
  a <- sin((lat2 - lat1) / 2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
  dist_km <- 2 * R_tierra * asin(sqrt(a))

  message("\n--- Verificación cruzada contra ", ruta_passport, " ---")
  message("Poblaciones comparadas: ", nrow(cruce),
          ". Diferencia máxima encontrada: ", round(max(dist_km), 6), " km")
  sospechosas <- cruce[dist_km > 1, ]
  if (nrow(sospechosas) > 0) {
    message(nrow(sospechosas), " población(es) con más de 1 km de diferencia -- revisar:")
    print(sospechosas[, c("pop_index", "Accession")])
  } else {
    message("Sin diferencias relevantes (todo por debajo de 1 km).")
  }
}

verificar_contra_passport(tabla, ruta_passport)


# ============================================================
# 4. Extraer valores del raster en cada coordenada
# ============================================================

con_coords <- tabla[!is.na(tabla$Latitude) & !is.na(tabla$Longitude), ]

r <- rast(ruta_tif)
message("Raster abierto: ", ruta_tif)
message("  CRS: ", crs(r))  # imprime el WKT completo; alcanza para confirmar que es geográfico (grados) y no proyectado
message("  Bandas: ", nlyr(r))
message("  Nombres de banda guardados en el archivo: ", paste(names(r), collapse = ", "))

if (nlyr(r) != length(nombres_variables)) {
  stop("El raster tiene ", nlyr(r), " bandas pero nombres_variables tiene ",
       length(nombres_variables), " nombres. Revisar la config.")
}

puntos <- vect(con_coords, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
valores <- terra::extract(r, puntos)  # data.frame: columna ID + una columna por banda
valores$ID <- NULL
colnames(valores) <- nombres_variables

con_coords <- cbind(con_coords, valores)


con_na <- con_coords[apply(is.na(con_coords[, nombres_variables]), 1, any), ]
if (nrow(con_na) > 0) {
  message(nrow(con_na), " población(es) con NA en al menos una banda tras la extracción:")
  print(con_na[, c("pop_index", "Accession", "Latitude", "Longitude")])
}
valor_nodata <- NAflag(r)  # puede devolver NaN si no está seteado en el archivo
if (!is.na(valor_nodata)) {
  con_nodata_literal <- con_coords[apply(con_coords[, nombres_variables] == valor_nodata, 1, any, na.rm = TRUE), ]
  if (nrow(con_nodata_literal) > 0) {
    message(nrow(con_nodata_literal), " población(es) con el valor nodata (",
            valor_nodata, ") literal (no convertido a NA) en al menos una banda:")
    print(con_nodata_literal[, c("pop_index", "Accession", "Latitude", "Longitude")])
  }
}


# ==========================
# 5. guardar tabla de referencia
# ===============================

columnas <- c("pop_index", "Accession", "Latitude", "Longitude", "Altitude", nombres_variables)
write.csv(con_coords[, columnas], ruta_salida, row.names = FALSE)
message("\nTabla de referencia guardada en: ", ruta_salida)
message("Filas: ", nrow(con_coords), " (esperado: 276, salvo que haya excluido alguna)")


# ==========
# 6. Chequeo 

message("\n--- Chequeo ---")
message("Comparamos contra Altitude (dato de campo, independiente de WorldClim) para confirmar que la extracción tomó el pixel correcto en cada coordenada.")

validos <- con_coords[!is.na(con_coords$Altitude), ]
correlaciones <- sapply(nombres_variables, function(var) {
  if (sum(!is.na(validos[[var]])) > 2) cor(validos$Altitude, validos[[var]], use = "complete.obs") else NA
})
print(round(correlaciones, 3))

r_bio1 <- correlaciones["bio1"]
message("\nChequeo principal -- bio1 (temperatura media anual) vs Altitude: r = ", round(r_bio1, 3))
if (is.na(r_bio1) || r_bio1 > -0.5) {
  message("[AVISO] Se esperaba una correlación fuerte y negativa (< -0.7 aprox). ",
          "Esto sí ameritaría revisar coordenadas, CRS o la extracción.")
} else {
  message("Correlación fuerte y negativa, como se espera -- confirma que la",
          " extracción está tomando el pixel correcto por población.")
}
message("\nQue bio5/6/8/9/10/11 (también temperatura) salgan con correlación fuerte",
        " no es un problema -- es lo esperado, responden al mismo gradiente de",
        " elevación. Las de precipitación (bio12-19) no tienen una relación",
        " universal esperada con la altitud; que salgan moderadas o fuertes en",
        " cualquier dirección no es, por sí solo, señal de error.")


