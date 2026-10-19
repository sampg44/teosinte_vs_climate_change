# ============================================================
# matriz_frecuencias_poblacionales.r
# 1 oct 2026
# sam
#
# tomar matriz de genotipos (individuos x SNPs, con valores 0/1/2) y convertirla en una matriz
# de frecuencias alélicas (poblaciones x SNPs, con valores entre 0 y 1),
# promediando los individuos de cada población.
#
# ============================================================


# ========
# 1. rutas -- AJUSTAR a tus archivos reales
# ========

# del data original completo
# ruta_bfile    <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all"
# ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"
# ruta_salida   <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_frecuencia_poblacional.rds"
# ruta_salida_copias  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_copias_posibles.rds"
# ruta_salida_alelos  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/alelo_contado_por_snp.csv"




#del filtrado

ruta_bfile    <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/podado_ld"
ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"

# frecuencias (0 a 1).  copias posibles = 2 x individuos genotipados, por poblacion y SNP
ruta_salida         <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_frecuencia_poblacional_253.rds"
# copias posibles = 2 x individuos genotipados, por poblacion y SNP
ruta_salida_copias  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_copias_posibles_253.rds" #o sea de cada poblacion para el glm

# qué alelo es el que se esta contando en cada SNP
ruta_salida_alelos  <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/alelo_contado_por_snp_253.csv"


# ============================================================
# 2. obtener la matriz de genotipos (individuos x SNPs)
# ============================================================

ruta_raw <- paste0(ruta_bfile, ".raw")

if (!file.exists(ruta_raw)) {
  message("No existe el .raw todavia -- generandolo con plink2 (--export A)...")
  resultado <- system2("plink2", c("--bfile", ruta_bfile, "--export", "A", "--out", ruta_bfile))
  if (resultado != 0) stop("plink2 fallo generando el .raw -- revisa el .log en la misma carpeta del bfile.")
} else {
  message("Ya existe ", ruta_raw, ", no se vuelve a generar (borralo si quieres forzar que se rehaga).")
}

# leemos el .raw. Las primeras 6 columnas son datos del individuo (familia, ID, padre, madre, sexo, fenotipo). 
# Los SNPs empiezan en la columna 7

raw <- read.table(ruta_raw, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

genotipos <- as.matrix(raw[, -(1:6)])
rownames(genotipos) <- raw$IID

# el nombre de columna que pone plink2 es algo como "snp123_A" (le pega el alelo contado al final). 
# ese sufijo dice a cuál alelo corresponden los 0/1/2 (a cuál alelo corresponde la frecuencia).
# lo guardo en una tabla antes de quitarlo, para interpretar un resultado futuro.
alelo_contado <- sub("^.*_([A-Za-z]+)$", "\\1", colnames(genotipos))
nombre_snp    <- sub("_[A-Za-z]+$", "", colnames(genotipos))
write.csv(data.frame(SNP = nombre_snp, alelo_contado = alelo_contado),
          ruta_salida_alelos, row.names = FALSE)

# quitamos el sufijo de los nombres de columna para quedarnos solo con el
# nombre del SNP
colnames(genotipos) <- nombre_snp

cat("Matriz de genotipos: ", nrow(genotipos), " individuos x ", ncol(genotipos), " SNPs\n", sep = "")


# ============================================================
# 3. averiguar a que población pertenece cada individuo
# ============================================================
# usar la metadata para saber, para cada individuo (Sample_name), a qué poblacion pertenece (Accession)

metadatos <- read.csv(ruta_metadata, stringsAsFactors = FALSE)

poblacion_de_cada_individuo <- metadatos$Accession[match(rownames(genotipos), metadatos$Sample_name)]

n_sin_poblacion <- sum(is.na(poblacion_de_cada_individuo))
if (n_sin_poblacion > 0) {
  stop(n_sin_poblacion, " individuos de la matriz de genotipos no se encontraron en data_teosinte.csv -- ",
       "revisa que los IID del .fam coincidan exactamente con Sample_name.")
}

poblaciones <- sort(unique(poblacion_de_cada_individuo))
cat("Poblaciones encontradas:", length(poblaciones), "\n")


# ============================================================
# 4. construir la matriz de frecuencias (poblaciones x SNPs)
# ============================================================
# La idea, para 1 poblacion y 1 SNP:
#   - tomar solo los individuos de esa pop
#   - sumo sus valores 0/1/2 en ese SNP (copias del alelo que se estan contando)
#   - el numero de copias posibles es 2 por cada individuo genotipado
#     en ese SNP, osea 2 * (individuos sin NA en ese SNP)
#   - frecuencia = copias contadas / copias posibles
#
# cuento individuos sin NA y no atodos los individuos porque puede haber datos faltantes 
# y a esos no los cuento para no sesgar la frecuencia hacia 0.
#
# lo hago para cada pop, y se va llenando una matriz nueva: poblaciones en filas, mismos SNPs en columnas que
# la matriz de genotipos

n_pob <- length(poblaciones)
n_snps <- ncol(genotipos)

matriz_frecuencia <- matrix(NA_real_, nrow = n_pob, ncol = n_snps,
                            dimnames = list(poblaciones, colnames(genotipos)))

# luego hago otra matriz como de conteo d cuantas copias del gen se pudieron contar en
# cada población y SNP (el denominador de cada frecuencia) para la posteridad 
matriz_copias_posibles <- matrix(NA_real_, nrow = n_pob, ncol = n_snps,
                                 dimnames = list(poblaciones, colnames(genotipos)))

for (i in seq_along(poblaciones)) {
  pob_actual <- poblaciones[i]
  
  # filas (individuos) de esta población
  filas_de_esta_poblacion <- genotipos[poblacion_de_cada_individuo == pob_actual, , drop = FALSE]
  
  copias_contadas  <- colSums(filas_de_esta_poblacion, na.rm = TRUE)       # suma de 0/1/2 por SNP
  individuos_genotipados <- colSums(!is.na(filas_de_esta_poblacion))       # cuantos no son NA, por SNP
  copias_posibles  <- 2 * individuos_genotipados                           # diploide
  
  frecuencia <- copias_contadas / copias_posibles
  # si copias_posibles es 0 (nadie de esta poblacion tiene dato en ese
  # SNP), la division da NaN pero no quiere decir que la frecuencia es cero
  frecuencia[copias_posibles == 0] <- NA
  
  matriz_frecuencia[i, ] <- frecuencia
  matriz_copias_posibles[i, ] <- copias_posibles
}

cat("Matriz de frecuencias: ", nrow(matriz_frecuencia), " poblaciones x ",
    ncol(matriz_frecuencia), " SNPs\n", sep = "")

n_na <- sum(is.na(matriz_frecuencia))
if (n_na > 0) {
  cat( n_na, " celdas (de ", length(matriz_frecuencia), ") quedaron NA ",
      "poblaciones sin ningun individuo genotipado en ese SNP en particular.\n", sep = "")
}


# ===========
# 5. guardar
# ==========

saveRDS(matriz_frecuencia, ruta_salida)
saveRDS(matriz_copias_posibles, ruta_salida_copias)
cat("\nGuardada en:", ruta_salida, "\n")
cat("Copias posibles guardadas en:", ruta_salida_copias, "\n")
cat("Para volver a usarla despues: matriz_frecuencia <- readRDS('", ruta_salida, "')\n", sep = "")

