# ============================================================
# generar_archivo_ambiental_bayescenv.r
# Generar el archivo ambiental de BayeScEnv a partir de un
# componente del PCA (por defecto PC1)

#
# fórmula (de Villemereuil & Gaggiotti 2015): g = |valor - referencia| / sd
#


# ========
# 1. rutas
# ========

# debe apuntar al orden_pop.csv de las 276 poblaciones completas, sin importar si vas a generar el archivo completo o uno
# de prueba (de ahí sale la media y sd que definen la escala de g).

# Así, el valor de g de una población es el mismo tanto en la corrida de prueba como en la corrida real (no se recalcula la escala cada
# vez con menos poblaciones). (si se quiere eso habría que hacer un script modificado)

ruta_orden_pop_referencia <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/orden_pop.csv"

# definir qué poblaciones y en qué orden se escriben en el archivo de salida. 

# para el archivo ambiental completo: el mismo que el de arriba (es como si solo hubera 1 varibale y no 2)
# para un archivo de prueba (subset): el orden_pop.csv NUEVO que genera convertir_a_bayescenv.r cuando
# corro con poblaciones_incluir puesto (ver seleccionar_poblaciones_prueba.R).
ruta_orden_pop_subset <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/orden_pop.csv"


ruta_scores <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/pca_19_bio/scores_pca_19_bio.csv"
componente_elegido <- "PC1"  # cambiar de componente si se necesita
ruta_salida <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/ambiental_bayescenv_PC1.txt"


# ============================================================
# 2. Calcular la referencia (media y sd) SIEMPRE sobre las 276
# ============================================================

orden_ref <- read.csv(ruta_orden_pop_referencia, stringsAsFactors = FALSE)
scores    <- read.csv(ruta_scores, stringsAsFactors = FALSE)

ref <- merge(orden_ref, scores[, c("Accession", componente_elegido)], by = "Accession")
if (nrow(ref) != nrow(orden_ref)) {
  warning(nrow(orden_ref) - nrow(ref), " población(es) del orden_pop.csv de referencia no se encontraron en scores_pca -- la media/sd se está calculando con menos de 276 poblaciones. Revisar antes de seguir.")
}

media_ref <- mean(ref[[componente_elegido]])
sd_ref    <- sd(ref[[componente_elegido]])
message("Referencia (", componente_elegido, ") calculada sobre ", nrow(ref),
        " poblaciones: media = ", round(media_ref, 4), ", sd = ", round(sd_ref, 4))


# ============================================================
# 3. tabla de salida (qué poblaciones y en qué orden)
# ============================================================

orden_salida <- read.csv(ruta_orden_pop_subset, stringsAsFactors = FALSE)
salida <- merge(orden_salida, scores[, c("Accession", componente_elegido)], by = "Accession")

# CRÍTICO: merge() no garantiza mantener el orden de pop_index --
# hay que reordenar a mano, porque el archivo de genotipos de
# BayeScEnv espera las poblaciones en ese orden exacto.

salida <- salida[order(salida$pop_index), ] # calcula la posición en la que debe ur cada fila para quedar ordenada por pop_index

faltantes <- salida[is.na(salida[[componente_elegido]]), ]
if (nrow(faltantes) > 0) {
  message(nrow(faltantes), " población(es) de ", ruta_orden_pop_subset,
          " no tienen valor de ", componente_elegido, " en scores_pca -- ",
          "revisar antes de escribir el archivo:")
  print(faltantes[, c("pop_index", "Accession")])
  stop("No se escribió el archivo -- resolver las poblaciones faltantes primero.")
}

if (nrow(salida) != nrow(orden_salida)) {
  stop("La tabla de salida tiene ", nrow(salida), " filas pero orden_salida tiene ", nrow(orden_salida), " -- algo no coincidió en el merge. Revisar.")
}


# ============================================================
# 4. calcular g y escribir el archivo en formato BayeScEnv
# ============================================================

g <- abs(salida[[componente_elegido]] - media_ref) / sd_ref
# [[]] siempre devuelve una sola columna como vector simple, me deja usar una variable para decir quiero la columna que se llama com o el VALOOR y si solo usara [], puede devolver varias columnas, (tipo data frame).

# una sola línea, valores separados por espacio, en el orden de
# orden_salida$pop_index (que es el mismo orden que debe tener el archivo de genotipos correspondiente)

writeLines(paste(g, collapse = " "), ruta_salida)

message("\nArchivo ambiental guardado en: ", ruta_salida)
message("Poblaciones escritas: ", nrow(salida), " | Componente usado: ", componente_elegido)
message("Rango de g: ", round(min(g), 3), " a ", round(max(g), 3),
        " (el paper señala que valores extremos normalmente no deberían pasar mucho de 2-3;",
        " si se ve algo mucho más grande, vale la pena revisar antes de correr el programa).")
