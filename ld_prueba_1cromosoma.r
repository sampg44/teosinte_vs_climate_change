# =========================
# ld_prueba_1cromosoma.r
# 30 sep 2026
# prueba si sirve esto para todo un chr
# ========================


library(snpStats)
library(LDheatmap)
library(ggplot2)


# ========
# 1.rutas (locales)
# ========

datasets <- list(
  r2_0.2 = list(
    ruta_bfile = "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/25_sep/r2_0.2/podado_ld",
    etiqueta   = "r2_0.2"
  ),
  data_inicial = list(
    ruta_bfile = "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all",
    etiqueta   = "data_inicial_T3604_33929"
  )
)

ruta_carpeta_salida <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/30_sep/ld_vis_prueba/"  # LOCAL, no de cluster

dir.create(ruta_carpeta_salida, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. función: procesar 1 cromosoma de 1 bfile,
#    sin límite de vecinos, matriz de r2 completa del cromosoma
# ============================================================

procesar_un_bfile <- function(ruta_bfile, etiqueta) {
  cat("\n.....\n")
  cat("Dataset:", etiqueta, "\n")
  cat("\n.....\n")
  
  t_lectura_0 <- Sys.time()
  datos_plink <- read.plink(paste0(ruta_bfile, ".bed"),
                            paste0(ruta_bfile, ".bim"),
                            paste0(ruta_bfile, ".fam"))
  t_lectura_1 <- Sys.time()
  cat("Tiempo leyendo el bfile:", round(as.numeric(difftime(t_lectura_1, t_lectura_0, units = "secs")), 1), "seg\n")
  
  genotipos_todos <- datos_plink$genotypes
  mapa_todos <- datos_plink$map
  cat("SNPs totales:", ncol(genotipos_todos), "| Individuos:", nrow(genotipos_todos), "\n")
  
  conteo_por_chr <- table(mapa_todos$chromosome)
  cat("SNPs por cromosoma:\n")
  print(conteo_por_chr)
  
  chr_mas_chico <- names(conteo_por_chr)[which.min(conteo_por_chr)]
  cat("\n-> Cromosoma elegido para la prueba (el de menos SNPs):", chr_mas_chico,
      "(", conteo_por_chr[chr_mas_chico], "SNPs )\n\n")
  
  idx_chr <- which(mapa_todos$chromosome == chr_mas_chico)
  idx_chr <- idx_chr[order(mapa_todos$position[idx_chr])]
  genotipos <- genotipos_todos[, idx_chr]
  mapa <- mapa_todos[idx_chr, ]
  n_snps_chr <- ncol(genotipos)
  
  # -------------------------------
  # 2a. r2 completo del cromosoma
  # -------------------------------
  cat("Calculando r2 completo (", n_snps_chr, "x", n_snps_chr, "pares =",
      format(n_snps_chr * (n_snps_chr - 1) / 2, big.mark = ","), "pares unicos ) ...\n")
  t0 <- Sys.time()
  r2_completo <- ld(genotipos, depth = n_snps_chr - 1, stats = "R.squared", symmetric = TRUE)
  r2_completo <- as.matrix(r2_completo)
  diag(r2_completo) <- 1
  t1 <- Sys.time()
  tiempo_r2 <- as.numeric(difftime(t1, t0, units = "secs"))
  cat("Tiempo calculando r2 completo del cromosoma:", round(tiempo_r2, 1), "seg\n")
  
  # ------------
  # 2b. heatmap 
  # ------------
  t0 <- Sys.time()
  png(paste0(ruta_carpeta_salida, "ld_heatmap_", etiqueta, "_chr", chr_mas_chico, ".png"),
      width = 5000, height = 5000, pointsize = 150)
  tryCatch({
    LDheatmap(r2_completo,
              genetic.distances = mapa$position,
              distances = "physical",
              LDmeasure = "r",
              title = paste0("LD heatmap ", etiqueta, " -- chr", chr_mas_chico, " completo (", n_snps_chr, " SNPs)"),
              color = colorRampPalette(c("blue", "magenta", "red"))(100))
  }, error = function(e) {
    cat("[ERROR en LDheatmap]:", conditionMessage(e), "\n")
  })
  dev.off()
  t1 <- Sys.time()
  cat("Tiempo dibujando el heatmap:", round(as.numeric(difftime(t1, t0, units = "secs")), 1), "seg\n")
  
  # --------------------------------------------------------
  # 2c. decay : todos los pares del cromosoma, agrupados en bins
  # --------------------------------------------------------
  pares <- which(upper.tri(r2_completo), arr.ind = TRUE)
  df_decay <- data.frame(
    distancia_bp = abs(mapa$position[pares[, 1]] - mapa$position[pares[, 2]]),
    r2 = r2_completo[pares]
  )
  
  ancho_bin <- 1e6
  df_decay$bin <- cut(df_decay$distancia_bp,
                      breaks = seq(0, max(df_decay$distancia_bp) + ancho_bin, by = ancho_bin))
  bins <- aggregate(r2 ~ bin, data = df_decay, FUN = mean)
  niveles <- levels(df_decay$bin)
  a <- as.numeric(sub("\\((.+),.*", "\\1", niveles))
  b <- as.numeric(sub(".*,(.+)\\]", "\\1", niveles))
  limites <- cbind(a, b)
  bins$distancia_media <- rowMeans(limites[match(bins$bin, niveles), , drop = FALSE])
  
  p_decay <- ggplot(bins, aes(x = distancia_media, y = r2)) + geom_point(color = "grey40", size = 1.5)
  if (nrow(bins) >= 10) {
    p_decay <- p_decay + geom_smooth(method = "loess", span = 0.3, color = "red", se = FALSE)
  } else {
    cat("[AVISO] solo", nrow(bins), "bins -- se omite la curva suavizada.\n")
  }
  p_decay <- p_decay +
    labs(title = paste0("LD decay ", etiqueta, " -- chr", chr_mas_chico, " completo"),
         x = "Distancia fisica (bp)", y = expression(R^2)) +
    theme_minimal(base_size = 12)
  ggsave(paste0(ruta_carpeta_salida, "ld_decay_", etiqueta, "_chr", chr_mas_chico, ".png"),
         p_decay, width = 8, height = 6, dpi = 150)
  
  cat("\n>>> RESUMEN", etiqueta, "chr", chr_mas_chico, ": ", n_snps_chr, "SNPs |",
      "tiempo r2:", round(tiempo_r2, 1), "seg <<<\n")
  
  invisible(list(etiqueta = etiqueta, chr = chr_mas_chico, n_snps = n_snps_chr, tiempo_r2_seg = tiempo_r2))
}


# ===============================
# 3. correr para los dos datasets
# ===============================

resultados <- lapply(datasets, function(d) procesar_un_bfile(d$ruta_bfile, d$etiqueta))

cat("\n\n============================================\n")
cat("antes de decidir si hacer los 10 cromosomas\n")
cat("============================================\n")
for (r in resultados) {
  cat(r$etiqueta, "| chr", r$chr, "|", r$n_snps, "SNPs |", round(r$tiempo_r2_seg, 1), "seg en calcular r2\n")
}


