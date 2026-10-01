# ============================================================
# ld_vis_por_cromosoma.r
# 30 sep 2026
# 
# version de ld_vis.r adaptada para correr sobre todo el genoma podado a r2=0.2 para compararlo con el data inicial
# gráficas por cromosoma
# ============================================================


library(snpStats)
library(LDheatmap)
library(ggplot2)


# ========
# 1. rutas
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

ruta_carpeta_salida <- "/home/sam/Documents/sur_ecoevo_lab/exp/sep_2026/teocintle/30_sep/ld_vis/"

dir.create(ruta_carpeta_salida, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. función: procesar un cromosoma completo (sin límite de vecinos)
# ============================================================

procesar_un_cromosoma <- function(genotipos_todos, mapa_todos, chr_id, etiqueta) {
  idx_chr <- which(mapa_todos$chromosome == chr_id)
  idx_chr <- idx_chr[order(mapa_todos$position[idx_chr])]
  genotipos <- genotipos_todos[, idx_chr]
  mapa <- mapa_todos[idx_chr, ]
  n_snps_chr <- ncol(genotipos)
  
  cat("--- ", etiqueta, " chr", chr_id, " (", n_snps_chr, " SNPs) ---\n", sep = "")
  
  if (n_snps_chr < 3) {
    cat("  muy pocos SNPs, se omite\n")
    return(invisible(list(etiqueta = etiqueta, chr = chr_id, n_snps = n_snps_chr,
                          tiempo_r2_seg = NA, tiempo_heatmap_seg = NA)))
  }
  
  t0 <- Sys.time()
  r2_completo <- ld(genotipos, depth = n_snps_chr - 1, stats = "R.squared", symmetric = TRUE)
  r2_completo <- as.matrix(r2_completo)
  diag(r2_completo) <- 1
  t1 <- Sys.time()
  tiempo_r2 <- as.numeric(difftime(t1, t0, units = "secs"))
  cat("  tiempo r2:", round(tiempo_r2, 1), "seg\n")
  
  # guardar la matriz, para poder rehacer plots despues (ej. zoom local, bins mas finos) sin volver a calcular nada
  saveRDS(list(r2 = r2_completo, position = mapa$position, snp_id = mapa$snp.name),
          paste0(ruta_carpeta_salida, "r2_matriz_", etiqueta, "_chr", chr_id, ".rds"))
  
  t0 <- Sys.time()
  png(paste0(ruta_carpeta_salida, "ld_heatmap_", etiqueta, "_chr", chr_id, ".png"),
      width = 5000, height = 5000, pointsize = 150)
  tryCatch({
    LDheatmap(r2_completo,
              genetic.distances = mapa$position,
              distances = "physical",
              LDmeasure = "r",
              title = paste0("LD heatmap ", etiqueta, " -- chr", chr_id, " completo (", n_snps_chr, " SNPs)"),
              color = colorRampPalette(c("blue", "magenta", "red"))(100))
  }, error = function(e) {
    cat("  [ERROR en LDheatmap]:", conditionMessage(e), "\n")
  })
  dev.off()
  t1 <- Sys.time()
  tiempo_heatmap <- as.numeric(difftime(t1, t0, units = "secs"))
  cat("  tiempo heatmap:", round(tiempo_heatmap, 1), "seg\n")
  
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
  }
  p_decay <- p_decay +
    labs(title = paste0("LD decay ", etiqueta, " -- chr", chr_id, " completo"),
         x = "Distancia fisica (bp)", y = expression(R^2)) +
    theme_minimal(base_size = 12)
  ggsave(paste0(ruta_carpeta_salida, "ld_decay_", etiqueta, "_chr", chr_id, ".png"),
         p_decay, width = 8, height = 6, dpi = 150)
  write.csv(df_decay, paste0(ruta_carpeta_salida, "ld_decay_pares_", etiqueta, "_chr", chr_id, ".csv"),
            row.names = FALSE)
  
  invisible(list(etiqueta = etiqueta, chr = chr_id, n_snps = n_snps_chr,
                 tiempo_r2_seg = tiempo_r2, tiempo_heatmap_seg = tiempo_heatmap))
}


# ============================================================
# 3. correr para los 10 cromosomas de cada dataset
# ============================================================

resumen <- list()

for (d in datasets) {
  cat("\n.....\n")
  cat("Dataset:", d$etiqueta, "\n")
  cat(".....\n")
  
  datos_plink <- read.plink(paste0(d$ruta_bfile, ".bed"),
                            paste0(d$ruta_bfile, ".bim"),
                            paste0(d$ruta_bfile, ".fam"))
  genotipos_todos <- datos_plink$genotypes
  mapa_todos <- datos_plink$map
  cromosomas <- sort(unique(mapa_todos$chromosome))
  cat("SNPs totales:", ncol(genotipos_todos), "| Individuos:", nrow(genotipos_todos),
      "| Cromosomas:", paste(cromosomas, collapse = ","), "\n\n")
  
  for (chr_id in cromosomas) {
    r <- procesar_un_cromosoma(genotipos_todos, mapa_todos, chr_id, d$etiqueta)
    resumen[[paste0(d$etiqueta, "_chr", chr_id)]] <- r
  }
}

df_resumen <- do.call(rbind, lapply(resumen, as.data.frame))
write.csv(df_resumen, paste0(ruta_carpeta_salida, "resumen_tiempos.csv"), row.names = FALSE)


print(df_resumen)
