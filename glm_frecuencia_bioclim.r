# ============================================================
# glm_frecuencia_bioclim.r
# 10 oct 2026
#
# Regresion logística (GLM binomial, link=logit) por locus:
#   cbind(conteo_alelo1, conteo_alelo2) ~ variable_bioclim + PCs
#
# corre 1 variable bioclimática a la vez
#
# las PCs de correccion por estructura poblacional se calculan sobre matriz de frecuencias poblacionales (prcomp), no
# promediando las PCs individuales del PCA de plink2.
#
#no sé si sería mejor usar el promedio por población de las PCs individuales (de plink2 --pca)
#
# ============================================================


# ========
# 1. rutas y parámetros
# ========

ruta_matriz_frecuencia      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_frecuencia_poblacional_253.rds"
ruta_matriz_copias_posibles <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/matriz_copias_posibles_253.rds"


ruta_bioclim <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/bioclim_por_pop.csv"

variable_bioclim <- "bio1"   # cambio en cada bio (no lo automaticé porque quiero observar cada resultado 1 x 1, en mi mente me da paz)
n_pcs_correccion <- 10        # cuantas PCs de estructura poblacional incluir como covariable, no sé cuántas convendrían

# TRUE  = quasibinomial: igual que binomial pero estima cuanta variacion "extra" hay entre poblaciones (deriva, historia compartida) y ajusta los p-valores.
# FALSE = binomial puro (supone que SOLO el muestreo de alelos genera variacion).

usar_dispersion <- TRUE

ruta_salida <- paste0("/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_", variable_bioclim, ".csv")


# =============================
# 2. cargar matrices y bioclim
# =============================

matriz_frecuencia      <- readRDS(ruta_matriz_frecuencia)
matriz_copias_posibles <- readRDS(ruta_matriz_copias_posibles)

# las dos matrices deben tener exactamente las mismas poblaciones (filas) y SNPs (columnas), en el mismo orden
stopifnot(identical(dimnames(matriz_frecuencia), dimnames(matriz_copias_posibles)))
cat("Matriz leida:", nrow(matriz_frecuencia), "poblaciones x", ncol(matriz_frecuencia), "SNPs\n")

bioclim <- read.csv(ruta_bioclim, stringsAsFactors = FALSE)
stopifnot(variable_bioclim %in% colnames(bioclim))

# emparejar bioclim con el orden de poblaciones de la matriz (por nombre)
bioclim_ordenado <- bioclim[match(rownames(matriz_frecuencia), bioclim$Accession), ]
n_sin_bioclim <- sum(is.na(bioclim_ordenado[[variable_bioclim]]))
if (n_sin_bioclim > 0) {
  warning(n_sin_bioclim, " poblaciones de la matriz no encontraron bioclim -- se van a omitir de este GLM.")
}


# ============================================================
# 3. PCs de estructura poblacional (sobre la matriz de frecuencias)
# ============================================================
# prcomp no acepta NA -- si matriz_frecuencia tiene celdas NA (poblaciones
# sin individuos genotipados en algun locus), se usan solo los loci sin ningun NA para
# calcular las PCs de estructura. Esto es solo para las PCs de corrección;
# los GLM de cada locus (abajo) sí usan todas las poblaciones disponibles
# para ese locus en particular, no se descarta nada ahi.
#
# Ademas, con scale. = TRUE prcomp se detiene si algun SNP tiene EXACTAMENTE
# la misma frecuencia en todas las poblaciones (desviacion estandar 0, no se
# puede dividir entre 0). Esos SNPs no aportan nada a la estructura, asi que
# se quitan antes.

loci_completos <- colSums(is.na(matriz_frecuencia)) == 0
desv_est <- apply(matriz_frecuencia[, loci_completos, drop = FALSE], 2, sd)
loci_para_pca <- loci_completos
loci_para_pca[loci_completos] <- desv_est > 0
cat(sum(loci_completos), "de", ncol(matriz_frecuencia), "loci sin NA;",
    sum(loci_completos) - sum(loci_para_pca), "sin variacion entre poblaciones (se quitan);",
    sum(loci_para_pca), "usados para calcular las PCs de estructura\n")

pca_pob <- prcomp(matriz_frecuencia[, loci_para_pca], scale. = TRUE)
cat("Varianza explicada por las primeras", n_pcs_correccion, "PCs:",
    round(sum(summary(pca_pob)$importance[2, 1:n_pcs_correccion]) * 100, 1), "%\n\n")

pcs_estructura <- as.data.frame(pca_pob$x[, 1:n_pcs_correccion, drop = FALSE])
colnames(pcs_estructura) <- paste0("PC", 1:n_pcs_correccion)


# =========================
# 4. GLM binomial por locus
# =========================

familia_glm <- if (usar_dispersion) quasibinomial(link = "logit") else binomial(link = "logit")
cat("Familia del GLM:", familia_glm$family, "\n")

n_loci <- ncol(matriz_frecuencia)
resultados <- data.frame(
  SNP = colnames(matriz_frecuencia),
  n_pop_usadas = NA_integer_,
  coef_bioclim = NA_real_,
  error_estandar = NA_real_,
  p_valor = NA_real_,
  dispersion = NA_real_,     # 1 = solo variacion de muestreo; >1 = hay variación extra entre poblaciones
  convergio = FALSE
)

t0 <- Sys.time()
for (j in seq_len(n_loci)) {
  conteo_alelo1 <- round(matriz_frecuencia[, j] * matriz_copias_posibles[, j])
  conteo_alelo2 <- matriz_copias_posibles[, j] - conteo_alelo1
  
  df <- data.frame(
    exito = conteo_alelo1,
    fracaso = conteo_alelo2,
    bioclim = bioclim_ordenado[[variable_bioclim]],
    pcs_estructura
  )
  df <- df[complete.cases(df) & matriz_copias_posibles[, j] > 0, ]
  resultados$n_pop_usadas[j] <- nrow(df)
  
  if (nrow(df) < (3 + n_pcs_correccion)) next  # muy pocas poblaciones para ajustar el modelo
  
  formula_glm <- as.formula(paste("cbind(exito, fracaso) ~ bioclim +",
                                  paste(colnames(pcs_estructura), collapse = " + ")))
  modelo <- tryCatch(
    suppressWarnings(glm(formula_glm, data = df, family = familia_glm)),
    error = function(e) NULL
  )
  if (is.null(modelo)) next
  
  s <- summary(modelo)$coefficients
  if ("bioclim" %in% rownames(s)) {
    resultados$coef_bioclim[j]   <- s["bioclim", "Estimate"]
    resultados$error_estandar[j] <- s["bioclim", "Std. Error"]
    resultados$p_valor[j]        <- s["bioclim", 4]    # columna 4 = p-valor (se llama Pr(>|z|) en binomial y Pr(>|t|) en quasibinomial)
    resultados$dispersion[j]     <- summary(modelo)$dispersion
    resultados$convergio[j]      <- modelo$converged
  }
  
  if (j %% 5000 == 0) cat("  ", j, "de", n_loci, "loci procesados...\n")
}
t1 <- Sys.time()
cat("\nTiempo total:", round(as.numeric(difftime(t1, t0, units = "secs")), 1), "seg para", n_loci, "loci\n")


# =====================================
# 5. corrección por pruebas múltiples
# =====================================
# BH (Benjamini-Hochberg) = control de FDR
# Bonferroni es mas estricto (controla FWER, no FDR)

resultados$q_valor_BH <- p.adjust(resultados$p_valor, method = "BH")
resultados$p_bonferroni <- p.adjust(resultados$p_valor, method = "bonferroni")

n_fallidos <- sum(is.na(resultados$p_valor))
if (n_fallidos > 0) {
  cat("aviso !!!", n_fallidos, "loci no se pudieron ajustar (muy pocas poblaciones o no convergio) -- excluidos del conteo de significativos\n")
}

cat("\nSignificativos a q(BH)<0.05:", sum(resultados$q_valor_BH < 0.05, na.rm = TRUE), "de", n_loci, "\n")
cat("Significativos a p(Bonferroni)<0.05:", sum(resultados$p_bonferroni < 0.05, na.rm = TRUE), "de", n_loci, "\n")

write.csv(resultados, ruta_salida, row.names = FALSE)
cat("\nResultados guardados en:", ruta_salida, "\n")


# intento de visualización

res <- read.csv(ruta_salida)
p <- res$p_valor[!is.na(res$p_valor)]
mean(p < 0.05)                                                    # debería acercarse a 0.05
median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)   # lambda, debería acercarse a 1
hist(p, breaks = 20)                                              # debería verse casi plano

