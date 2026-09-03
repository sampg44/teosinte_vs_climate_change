# ============================================================
# Convertir genotipos PLINK (.bed/.bim/.fam) al formato de entrada de BayeScEnv
#
# formato de salida esperado:
#   [loci]=100
#   [populations]=16
#   [pop]=1
#   1 40 2 29 11
#   2 40 2 4 36
#   ...
# columnas por línea: locus  tamaño_muestra(copias de gen)  n_alelos  conteo_alelo1  conteo_alelo2
#
# Usa Accession (276 niveles) como población, NO POB_CODE (74),
# consistente con la corrección hecha en basic_analysis_teosinte.pdf
# 
# última versión: 2 sept 2026
# ============================================================

library(SNPRelate)
library(dplyr)

# --- 1. Rutas ---------
ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"
ruta_fam      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.fam"
ruta_bim      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.bim"
ruta_bed      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.bed"
ruta_gds      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/teosinte.gds"

#cambiar rutas cuando se hagan pruebas  (a su propia carpeta) o se puede sobreescribir 
ruta_salida_geno  <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/276_pob_bayesformato.txt"  # input para bayescenv
ruta_salida_orden <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/orden_pop.csv"             # mapping pop_index -> Accession
ruta_salida_diccionario <- file.path(dirname(ruta_salida_geno), "diccionario_alelos.csv")

#Deja NULL para correr con las 276 poblaciones completas.
poblaciones_incluir <- NULL  

# --- 2. Metadata y asignación de población (Accession, no POB_CODE) ------
Data <- read.csv2(ruta_metadata, header = TRUE, sep = ",", as.is = FALSE)

fam <- read.table(ruta_fam, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

meta_fam <- merge(fam, Data, by.x = "IID", by.y = "Sample_name", all.x = TRUE)
stopifnot(sum(is.na(meta_fam$Accession)) == 0)  # todos los individuos deben tener población asignada

# --- 3. Abrir (o crear) el GDS --------------------------------------------
if (!file.exists(ruta_gds)) {
  snpgdsBED2GDS(bed.fn = ruta_bed, bim.fn = ruta_bim, fam.fn = ruta_fam,
                out.gdsfn = ruta_gds)
}
genofile <- snpgdsOpen(ruta_gds)

samp   <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

poblacion <- meta_fam$Accession[match(samp, meta_fam$IID)]
stopifnot(sum(is.na(poblacion)) == 0)

# Orden fijo de poblaciones

poblaciones <- sort(unique(poblacion))

if (!is.null(poblaciones_incluir)) {
  faltantes <- setdiff(poblaciones_incluir, poblaciones)
  if (length(faltantes) > 0) warning("Estos Accession no existen en los datos: ",
                                     paste(faltantes, collapse = ", "))
  poblaciones <- poblaciones[poblaciones %in% poblaciones_incluir]
}
n_pops <- length(poblaciones)
if (is.null(poblaciones_incluir)) stopifnot(n_pops == 276)

write.csv(data.frame(pop_index = seq_along(poblaciones), Accession = poblaciones),
          ruta_salida_orden, row.names = FALSE)
message("Orden de poblaciones guardado en: ", ruta_salida_orden,
        "  usar mismo orden para el archivo ambiental")

# --- 4. Extraer matriz de genotipos (individuos x SNPs) -------------------
# Valores: 0/1/2 = copias del alelo fijado por SNPRelate para ese locus
# NA = dato faltante

geno <- snpgdsGetGeno(genofile, snp.id = snp_id, sample.id = samp,
                      with.id = FALSE, verbose = TRUE)

# dim(geno) debe ser: length(samp) x length(snp_id)

n_loci <- length(snp_id)  

# --- 5. Diccionario de identidad real de los alelos -----------------------
# "conteo_alelo1"/"conteo_alelo2" en el archivo de salida NO son A1/A2 del
# .bim directamente
# leo el campo snp.allele DEL MISMO GDS que usamos para contar para garantizar consistencia

bim <- read.table(ruta_bim, header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim")

alelos_gds <- read.gdsn(index.gdsn(genofile, "snp.allele"))  # formato "A/B", p.ej. "C/T"
alelos_split <- strsplit(alelos_gds, "/")

# snp_id NO son posiciones numéricas -- son los nombres de SNP (texto, como
# "S1_992727"), iguales a la columna SNP del bim. Hay que buscar la POSICIÓN
# de cada uno dentro de bim$SNP antes de poder indexar CHR/A1_bim/A2_bim.
posiciones_bim <- match(snp_id, bim$SNP)
stopifnot(sum(is.na(posiciones_bim)) == 0)  # todo snp_id debería encontrarse en el bim

diccionario_alelos <- data.frame(
  locus_index_bayescenv = seq_len(n_loci),        # el índice que aparece en el archivo de BayeScEnv
  snp_id_gds            = snp_id,
  CHR                    = bim$CHR[posiciones_bim],
  SNP                    = bim$SNP[posiciones_bim],
  A1_bim                 = bim$A1_bim[posiciones_bim],
  A2_bim                 = bim$A2_bim[posiciones_bim],
  alelo1_conteo          = vapply(alelos_split, `[`, character(1), 1),  # el que se cuenta en "conteo_alelo1"
  alelo2_conteo          = vapply(alelos_split, `[`, character(1), 2)   # el que se cuenta en "conteo_alelo2"
)
write.csv(diccionario_alelos, ruta_salida_diccionario, row.names = FALSE)
message("Diccionario de alelos (qué nucleótido es alelo1/alelo2 por locus) guardado en: ", ruta_salida_diccionario)


# Diagnóstico primero: ¿de dónde vienen los NA? Corre esto y revisa antes de
# confiar en el porcentaje de abajo.
cat("NAs en alelo1_conteo:", sum(is.na(diccionario_alelos$alelo1_conteo)), "\n")
cat("NAs en A1_bim:", sum(is.na(diccionario_alelos$A1_bim)), "\n")
# ejemplo de las primeras filas afectadas, para inspeccionar a mano:
print(head(diccionario_alelos[is.na(diccionario_alelos$alelo1_conteo) | is.na(diccionario_alelos$A1_bim), ]))

# Chequeo rápido de cordura: alelo1_conteo/alelo2_conteo deben coincidir con
# A1_bim/A2_bim en la mayoría de los loci -- si NO coinciden en muchos casos,
# eso confirma que SNPRelate sí hizo un swap respecto al bim, y hay que usar
# SIEMPRE diccionario_alelos.csv (nunca el bim original) para interpretar
# después cuáles nucleótidos están detrás de cada locus outlier.
# na.rm = TRUE evita que un solo NA vuelva todo el resultado NA (bug de la
# versión anterior); en cambio, reportamos aparte cuántos loci no se pudieron
# comparar.
n_na <- sum(is.na(diccionario_alelos$alelo1_conteo) | is.na(diccionario_alelos$A1_bim))
mismatch <- mean(diccionario_alelos$alelo1_conteo != diccionario_alelos$A1_bim, na.rm = TRUE)
message(n_na, " de ", n_loci, " loci no se pudieron comparar (NA en alelo1_conteo o A1_bim). ",
        "De los que sí: ", round(100 * mismatch, 1), "% tienen alelo1_conteo distinto de A1_bim -- ",
        "usa siempre diccionario_alelos.csv, no el bim, para interpretar resultados.")

# --- 6. escribir el archivo en formato BayeScEnv --
poblaciones_sin_datos <- character(0)  # aquí se van a acumular, si las hay

con <- file(ruta_salida_geno, "w")
writeLines(paste0("[loci]=", n_loci), con)
writeLines(paste0("[populations]=", n_pops), con)

for (p_index in seq_along(poblaciones)) {
  p   <- poblaciones[p_index]
  idx <- which(poblacion == p)
  genotipos_poblacion_actual <- geno[idx, , drop = FALSE]
  #genotipos_poblacion_actual es el recorte de la matriz de genotipos a su idx para quedarme con las filas (individuos)
  # de esa poblacion que está tomando el for en ese momento 
  
  n_ind_por_locus <- colSums(!is.na(genotipos_poblacion_actual))       # individuos genotipados en esta pob, por locus
  sample_size     <- 2 * n_ind_por_locus              # copias de gen = tamaño de muestra que pide el formato
  conteo_alelo1   <- colSums(genotipos_poblacion_actual, na.rm = TRUE)  # suma de dosis = copias del "alelo 1"
  conteo_alelo2   <- sample_size - conteo_alelo1
  
  # --- Chequeos automáticos
  stopifnot(all(conteo_alelo1 + conteo_alelo2 == sample_size))  #simple chequeo de que voy bien en el código jiji
  if (all(sample_size == 0)) poblaciones_sin_datos <- c(poblaciones_sin_datos, p)
  
  writeLines(paste0("[pop]=", p_index), con)
  bloque_poblacion <- paste(seq_len(n_loci), sample_size, 2L, conteo_alelo1, conteo_alelo2)
  writeLines(bloque_poblacion, con)
  # bloque_poblacion es el conjunto de líneas de texto (1 por locus)  que corresponden a cada población
}
close(con)
snpgdsClose(genofile)

message("Archivo de genotipos para BayeScEnv escrito en: ", ruta_salida_geno)

# --- 7. Resultado de los chequeos automáticos ------------------------------
if (length(poblaciones_sin_datos) > 0) {
  warning(length(poblaciones_sin_datos), " poblacion(es) sin NINGÚN locus genotipado (sample_size 0 en todos): ",
          paste(poblaciones_sin_datos, collapse = ", "))
} else {
  message("Chequeo OK: todas las poblaciones tienen al menos un locus con datos, ",
          "y en todos los loci alelo1 + alelo2 = sample_size.")
}

