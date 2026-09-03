# ============================================================
# Seleccionar un subset de poblaciones para una corrida de prueba
# seleccionar_poblaciones_prueba.r
# 
# este script lo menciona convertir_a_--------bayescenv.r en la línea 36 aproximadamente
# ============================================================

ruta_orden_pop <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/orden_pop.csv"

n_poblaciones_prueba <- 20  # ajustar cuántas poblaciones para la prueba

ruta_salida_lista <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/bayescenv/prueba1_02.09.2026/pobs_prueba.csv"

# ==============================================
# Selección aleatoria (simple, no estratificada)

orden <- read.csv(ruta_orden_pop, stringsAsFactors = FALSE)

set.seed(44) # reproducible jaj
seleccion <- sort(sample(orden$Accession, n_poblaciones_prueba))

dir.create(dirname(ruta_salida_lista), recursive = TRUE, showWarnings = FALSE) #solo forza la creación de la carpeta si noe xiste
write.csv(data.frame(Accession = seleccion), ruta_salida_lista, row.names = FALSE)
message(length(seleccion), " poblaciones seleccionadas, guardadas en: ", ruta_salida_lista)

# ------------------------------------------------------------
# pegar en convertir_a_bayescenv.r, en la línea de "poblaciones_incluir <- NULL":
# ------------------------------------------------------------
cat("\nPegar esto en convertir_a_bayescenv.r (reemplaza la línea de poblaciones_incluir):\n\n")
cat("poblaciones_incluir <- c(", paste0('"', seleccion, '"', collapse = ", "), ")\n\n")

# IMPORTANTE, antes de volver a correr convertir_a_bayescenv.r !!!!

# Además de pegar el vector de arriba, hay que cambiar ESTAS DOS rutas " dentro de convertir_a_bayescenv.r a una carpeta distinta (por ejemplo ",'bayescenv/prueba/'), 
# o el script va a SOBRESCRIBIR el archivo de genotipos y el orden_pop.csv de las 276 poblaciones ya validadas:
# ruta_salida_geno  <- '.../bayescenv/prueba/subset_bayesformato.txt'"
# ruta_salida_orden <- '.../bayescenv/prueba/orden_pop_subset.csv'"
# El orden_pop_subset.csv que se genere ahí es el que hay que usar como'ruta_orden_pop_salida' en generar_archivo_ambiental_bayescenv.R",
# para que el archivo ambiental de prueba tenga las mismas poblaciones, en el mismo orden, que el archivo de genotipos de prueba.")
