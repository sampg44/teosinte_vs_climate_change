# =========================
# script_pca_poda.r
# 24 sep 2026
# sam
# después del ld prunning (en realidad no borró nada pero pues para tener la figura y no me regañen)
# =========================


library(adegenet)

reportar_tiempo <- function(etiqueta, t_referencia) {
  transcurrido <- as.numeric(difftime(Sys.time(), t_referenccia, units = "mins"))
  cat(sprintf("[tiempo] %s: %.2f min\n", etiqueta, transcurrido))
  Sys.time()
}
t0 <- Sys.time()
t_inicio_total <- t0


# ........
# 1. rutas
# ..........

# cluster
ruta_bfile <- "/mnt/data/sur/users/spacheco/data/teosinte/podado_ld" 
ruta_metadata <- "/mnt/data/sur/users/spacheco/data/teosinte/data_teosinte.csv"
ruta_carpeta_salida <- "/mnt/data/sur/users/spacheco/results/sep_2026/teo/24_sep/poda/"



dir.create(ruta_carpeta_salida, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(ruta_carpeta_salida)) {
  stop("No se pudo crear/acceder a la carpeta de salida: ", ruta_carpeta_salida,
       ", revisar mayúsculas/minúsculas y permisos antes de seguir.")
}


# .................................................................................
# 2. Convertir .bed/.bim/.fam a .raw con plink2 (tabs -> espacios, porque adegenet::read.PLINK necesita espacios)
# .................................................................................

ruta_raw_original <- paste0(ruta_bfile, ".raw")
ruta_raw_convertido <- paste0(ruta_bfile, "_convertido.raw")

if (!file.exists(ruta_raw_convertido)) {
  message("Generando .raw con plink2 --export A...")
  resultado <- system2("plink2", c("--bfile", ruta_bfile, "--export", "A", "--out", ruta_bfile))
  if (resultado != 0 || !file.exists(ruta_raw_original)) {
    stop("plink2 no corrió correctamente. No se puede continuar sin el .raw generado.")
  }
  message("Convirtiendo delimitador (tabs -> espacios)...")
  lineas <- readLines(ruta_raw_original)
  writeLines(gsub("\t", " ", lineas), ruta_raw_convertido)
} else {
  message(ruta_raw_convertido, " ya existe, no se regenera.")
}
t0 <- reportar_tiempo("conversión a .raw", t0)


# ................................................................
# 3. leer como genlight y asignar población (Accession)
# .......................................

gl <- read.PLINK(ruta_raw_convertido, quiet = TRUE)
cat("Individuos leídos:", nInd(gl), "| Loci:", nLoc(gl), "\n")

meta <- read.csv(ruta_metadata, stringsAsFactors = FALSE)
orden_gl <- data.frame(IID = indNames(gl))
cruce <- merge(orden_gl, meta[, c("Sample_name", "Accession")],
               by.x = "IID", by.y = "Sample_name", all.x = TRUE, sort = FALSE)
cruce <- cruce[match(indNames(gl), cruce$IID), ]

faltantes <- sum(is.na(cruce$Accession))
if (faltantes > 0) {
  message("[AVISO] ", faltantes, " individuo(s) del genlight no encontraron Accession -- revisar.")
}

pop(gl) <- cruce$Accession
cat("Poblaciones (Accession) distintas asignadas:", length(unique(pop(gl))), "\n")
t0 <- reportar_tiempo("leer genlight + asignar Accession", t0)

# guardar el genlight YA CON la población asignada -- esto es lo que
# dapc_estructura.R va a leer, para no repetir nada de lo anterior
saveRDS(gl, paste0(ruta_carpeta_salida, "gl.rds"))


# .......
# 4. PCA 
# .........

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
pca <- glPca(gl, nf = 10, parallel = (n_cores > 1), n.cores = n_cores)

varianza_pca <- pca$eig / sum(pca$eig) * 100
write.csv(data.frame(eje = paste0("PC", seq_along(varianza_pca)), varianza_pct = varianza_pca)[1:10, ],
          paste0(ruta_carpeta_salida, "pca_varianza.csv"), row.names = FALSE)
write.csv(data.frame(Accession = pop(gl), pca$scores),
          paste0(ruta_carpeta_salida, "pca_scores.csv"), row.names = FALSE)
cat("Varianza explicada, primeros 3 ejes:", round(varianza_pca[1:3], 2), "\n")
t0 <- reportar_tiempo("PCA (glPca)", t0)

library(ggplot2)
n_pops_distintas <- length(unique(pop(gl)))
df_pca <- data.frame(PC1 = pca$scores[, 1], PC2 = pca$scores[, 2], Accession = pop(gl))
p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Accession)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(x = paste0("PC1 (", round(varianza_pca[1], 1), "%)"),
       y = paste0("PC2 (", round(varianza_pca[2], 1), "%)"),
       title = "PCA -- coloreado por Accession") +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())
if (n_pops_distintas > 30) {
  p_pca <- p_pca + theme(legend.position = "none")
}
ggsave(paste0(ruta_carpeta_salida, "pca_plot.png"), p_pca, width = 9, height = 7, dpi = 150)

cat(sprintf("\n[tiempo] TOTAL de pca_estructura.R: %.2f min\n",
            as.numeric(difftime(Sys.time(), t_inicio_total, units = "mins"))))
cat("\nListo. Para el DAPC, corre dapc_estructura.R apuntando a esta misma ruta_carpeta_salida.\n")

