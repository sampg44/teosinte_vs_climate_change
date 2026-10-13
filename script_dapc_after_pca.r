# .................................
# script_dapc_after_pca.r
# 24 sep 2026
#
# depende de que script_pca.r ya haya corrido sobre la misma carpeta de salida y haya dejado gl.rdsya armado (con Accession ya asignada).
#
# ..........................................


library(adegenet)

reportar_tiempo <- function(etiqueta, t_referencia) {
  transcurrido <- as.numeric(difftime(Sys.time(), t_referencia, units = "mins"))
  cat(sprintf("[tiempo] %s: %.2f min\n", etiqueta, transcurrido))
  Sys.time()
}
t0 <- Sys.time()
t_inicio_total <- t0

graficar_en_png <- function(ruta, ancho, alto, expr_grafico) {
  png(ruta, width = ancho, height = alto)
  on.exit(dev.off())
  expr_grafico()
}

buscar_K_con_salvaguarda <- function(gl, n_pca, k_max, etiqueta, carpeta_salida) {
  grupos <- find.clusters(gl, n.pca = n_pca, max.n.clust = k_max, choose.n.clust = FALSE)
  
  bic_df <- data.frame(K = seq_along(grupos$Kstat), BIC = as.numeric(grupos$Kstat))
  write.csv(bic_df, paste0(carpeta_salida, "bic_vs_K.csv"), row.names = FALSE)
  graficar_en_png(paste0(carpeta_salida, "bic_vs_K.png"), 800, 500, function() {
    plot(bic_df$K, bic_df$BIC, type = "b", xlab = "K", ylab = "BIC", main = "BIC vs K")
  })
  
  k_elegido <- length(unique(grupos$grp))
  if (k_elegido == k_max) {
    message("aviso: El K elegido automáticamente (", k_elegido, ") choca exactamente ",
            "contra el techo (max.n.clust=", k_max, "). El BIC probablemente no encontró un ",
            "mínimo real dentro del rango probado -- revisa bic_vs_K.png antes de confiar en este K.")
  }
  grupos
}


# ........................
# 1. rutas y parámetros
# ........................

# misma ruta_carpeta_salida que en script_pca.r (de ahí lee gl.rds)
ruta_carpeta_salida <- "/mnt/data/sur/users/spacheco/results/sep_2026/teo/24_sep/"

k_max <- 400

# NULL = correr xvalDapc para decidir el número de PCs (checkpoint 1, línea base).
# Un número = saltarse xvalDapc y usar ese número (checkpoints siguientes).
n_pcs_fijo <- 180

ruta_gl <- paste0(ruta_carpeta_salida, "gl.rds")
if (!file.exists(ruta_gl)) {
  stop("No existe ", ruta_gl, " , correr primero el script de pca sobre esta misma carpeta de salida.")
}
gl <- readRDS(ruta_gl)
cat("gl.rds cargado -- Individuos:", nInd(gl), "| Loci:", nLoc(gl),
    "| Poblaciones:", length(unique(pop(gl))), "\n")
t0 <- reportar_tiempo("cargar gl.rds", t0)


# ................................................................
# 2. Número de PCs para el DAPC: fijo, o decidido por xvalDapc
# ................................................................

if (is.null(n_pcs_fijo)) {
  cat("\n--- n_pcs_fijo es NULL -- corriendo xvalDapc para decidir ---\n")
  mat <- as.matrix(gl)
  mat[is.na(mat)] <- 0
  xval <- xvalDapc(mat, pop(gl), n.pca.max = xval_n_pca_max, n.rep = xval_n_rep, xval.plot = FALSE)
  n_pcs <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
  cat("Número de PCs sugerido por xvalDapc:", n_pcs, "\n")
  t0 <- reportar_tiempo(paste0("xvalDapc (hasta ", xval_n_pca_max, " PCs, ", xval_n_rep, " rep.)"), t0)
} else {
  n_pcs <- n_pcs_fijo
  cat("\n--- Usando n_pcs_fijo =", n_pcs, "(sin correr xvalDapc) ---\n")
}



#........
# 3. DAPC
# ........

grupos <- buscar_K_con_salvaguarda(gl, n_pcs, k_max, "dapc", ruta_carpeta_salida)
t0 <- reportar_tiempo(paste0("find.clusters (", n_pcs, " PCs, K=1 a ", k_max, ")"), t0)

if (length(unique(grupos$grp)) < 2) {
  stop("K encontrado fue 1, el DAPC no puede correr con un solo grupo. Revisar bic_vs_K.png.")
}

dapc_obj <- dapc(gl, pop = grupos$grp, n.pca = n_pcs, n.da = length(unique(grupos$grp)) - 1)
saveRDS(dapc_obj, paste0(ruta_carpeta_salida, "dapc_", n_pcs, "pcs.rds"))

tabla_asignaciones <- data.frame(
  IID = indNames(gl),
  Accession = pop(gl),
  grupo_dapc = grupos$grp
)
write.csv(tabla_asignaciones, paste0(ruta_carpeta_salida, "asignaciones_dapc.csv"), row.names = FALSE)

k_final <- length(unique(grupos$grp))
cat("\nDAPC con", n_pcs, "PCs -- K encontrado:", k_final, "\n")
cat("Proporción de reasignación correcta:", round(summary(dapc_obj)$assign.prop, 4), "\n")
t0 <- reportar_tiempo(paste0("dapc (", n_pcs, " PCs)"), t0)

paleta <- rainbow(k_final)
graficar_en_png(paste0(ruta_carpeta_salida, "dapc_scatter.png"), 900, 700, function() {
  scatter(dapc_obj, col = paleta, bg = "white", cstar = 0,
          legend = (k_final <= 30), posi.leg = "topright",
          scree.pca = TRUE, posi.pca = "bottomleft",
          main = paste0("DAPC (", n_pcs, " PCs, K=", k_final, ")"))
})
graficar_en_png(paste0(ruta_carpeta_salida, "compoplot.png"), 1200, 600, function() {
  compoplot(dapc_obj, col = paleta, legend = (k_final <= 30),
            main = paste0("Compoplot tipo ADMIXTURE (", n_pcs, " PCs, K=", k_final, ")"))
})

cat(sprintf("\n[tiempo] TOTAL de dapc_estructura.R: %.2f min\n",
            as.numeric(difftime(Sys.time(), t_inicio_total, units = "mins"))))
