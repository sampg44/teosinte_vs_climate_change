# ============================================================
# hacer_matriz_freq_pops.r
# 1 oct 2026
#
# tomar matriz de genotipos (individuos x SNPs, con valores 0/1/2) y convertirla en una de
# de frecuencias alélicas (poblaciones x SNPs, con valores entre 0 y 1),
# promediando los individuos de cada poblacion.
#
# ============================================================


# ========
# 1. rutas
# ========

ruta_bfile    <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all"
ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"
ruta_salida   <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_frecuencia_poblacional.rds"


# ============================================================
# 2. obtener la matriz de genotipos (individuos x SNPs)
# ============================================================


ruta_raw <- paste0(ruta_bfile, ".raw")

if (!file.exists(ruta_raw)) {
  message("No existe el .raw todavia, generándolo con plink2")
  resultado <- system2("plink2", c("--bfile", ruta_bfile, "--export", "A", "--out", ruta_bfile))
  if (resultado != 0) stop("plink2 falló generando el .raw, revisar el .log en la misma carpeta del bfile.")
} else {
  message("Ya existe ", ruta_raw, ", no se vuelve a generar.")
}

# leemos el .raw. Las primeras 6 columnas son datos del individuo (familia, ID, padre, madre, sexo, fenotipo). 
# Los SNPs empiezan en la columna 7

raw <- read.table(ruta_raw, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

genotipos <- as.matrix(raw[, -(1:6)])   # quitamos esas primeras 6 columnas y ya es la mtariz de genotipos

rownames(genotipos) <- raw$IID

# el nombre de columna que pone plink2 es algo como "snp123_A" (le pega el alelo contado al final)
#lo quitamos para quedarnos solo con el nombre del SNP
colnames(genotipos) <- sub("_[A-Za-z]+$", "", colnames(genotipos))

cat("Matriz de genotipos: ", nrow(genotipos), " individuos x ", ncol(genotipos), " SNPs\n", sep = "")


# ============================================================
# 3. averiguar a qué poblacion pertenece cada individuo
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
# La idea, para 1 población y 1 SNP:
#   - tomamos solo los individuos de esa población
#   - sumamos sus valores 0/1/2 en ese SNP -> eso es el numero total
#     de copias del alelo que se estan contando, en esa poblacion
#   - el numero de copias posibls es 2 por cada individuo genotipado
#     en ese SNP (porque tengo puros diploides) osea 2 * (individuos sin NA en ese SNP)
#   - frecuencia = copias contadas / copias posibles
#
# individuos sin NA y no simplemente todos los individuos
# de la poblacion porque puede haber datos faltantes de genotipación, a esos no los
# contamos ni arriba ni abajo de la division, para no sesgar la frecuencia hacia 0.
#
# Hacemos esto para cada poblacion por separado, y vamos llenando una matriz nueva con
# poblaciones en filas, los mismos SNPs en columnas que la matriz de genotipos 
# (por eso se puede comparar celda a celda entre las dos matrices).

n_pob <- length(poblaciones)
n_snps <- ncol(genotipos)

matriz_frecuencia <- matrix(NA_real_, nrow = n_pob, ncol = n_snps,
                            dimnames = list(poblaciones, colnames(genotipos)))

for (i in seq_along(poblaciones)) {
  pob_actual <- poblaciones[i]
  
  # filas (individuos) de esa pob
  filas_de_esta_poblacion <- genotipos[poblacion_de_cada_individuo == pob_actual, , drop = FALSE]
  
  copias_contadas  <- colSums(filas_de_esta_poblacion, na.rm = TRUE)       # suma de 0/1/2 por SNP
  individuos_genotipados <- colSums(!is.na(filas_de_esta_poblacion))       # cuántos NO son NA, por SNP
  copias_posibles  <- 2 * individuos_genotipados                           # diploide
  frecuencia <- copias_contadas / copias_posibles
  # si copias_posibles es 0 (nadie de esta poblacion tiene dato en ese SNP), lo dejamos como NA explicito, porque
  # no es que sea 0 sino que la desconocemos
  frecuencia[copias_posibles == 0] <- NA
  
  matriz_frecuencia[i, ] <- frecuencia
}

cat("Matriz de frecuencias: ", nrow(matriz_frecuencia), " poblaciones x ",
    ncol(matriz_frecuencia), " SNPs\n", sep = "")

n_na <- sum(is.na(matriz_frecuencia))
if (n_na > 0) {
  cat("[AVISO] ", n_na, " celdas (de ", length(matriz_frecuencia), ") quedaron NA -- ",
      "poblaciones sin ningun individuo genotipado en ese SNP en particular.\n", sep = "")
}


# ============
# 5. guardar
# ===========

saveRDS(matriz_frecuencia, ruta_salida)
cat("\nGuardada en:", ruta_salida, "\n")
cat("Para volver a usarla despues: matriz_frecuencia <- readRDS('", ruta_salida, "')\n", sep = "")
