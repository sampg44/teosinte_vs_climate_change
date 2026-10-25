# ==============================================================================
# visualizar_glm.r
# versión nueva (24 oct )
# ay llevo bastantes horas sin comer ayuda
# Visualiza los resultados de glm_frecuencia_bioclim.r
# ==============================================================================

library(ggplot2)


# 1. rutas y parámtero


archivos_glm <- c(
  bio1 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio1.csv",
  bio2 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio2.csv",
  bio3 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio3.csv",
  bio4 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio4.csv",
  bio5 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio5.csv",
  bio6 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio6.csv",
  bio7 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio7.csv",
  bio8 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio8.csv",
  bio9 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio9.csv",
  bio10 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio10.csv",
  bio11 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio11.csv",
  bio12 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio12.csv",
  bio13 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio13.csv",
  bio14 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio14.csv",
  bio15 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio15.csv",
  bio16 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio16.csv",
  bio17 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio17.csv",
  bio18 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio18.csv",
  bio19 = "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_bio19.csv"
)
carpeta_figuras <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/24_oct/figuras_glm"

alfa <- 0.05                  # nivel de significancia ANTES de corregir por pruebas múltiples
excluir_degenerados <- TRUE   # quita SNPs con ajuste inestable (ver sección 2)
hacer_qq     <- TRUE          # QQ-plot por variable
hacer_volcan <- FALSE          # volcán por variable
dir.create(carpeta_figuras, recursive = TRUE, showWarnings = FALSE)

# Paleta colores
verde_1 <- "#6fb890"; verde_2 <- "#b5dcc4"   # no significativos (cromosomas alternos)
verde   <- "#2f9e6e"                          # verde "fuerte" para puntos aislados
morado  <- "#6250d6"                          # supera Bonferroni
lila    <- "#9b7fd8"; amarillo <- "#c98500"
col_texto <- "#52514e"

tema <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = col_texto),
        legend.position = "top")


# 2. leer archivos

lambda_gc <- function(p) {      # factor de inflación genómica (Devlin & Roeder 1999)
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

preparar <- function(ruta, etiqueta) {
  d <- read.csv(ruta, stringsAsFactors = FALSE)
  m_pruebas <- sum(is.finite(d$p_valor))          # número de pruebas hechas (para Bonferroni)
  # nombre de SNP tipo "S2_5269870" (o "S2_5269870_-"): cromosoma y posición
  m <- regmatches(d$SNP, regexec("^S([0-9]+)_([0-9]+)", d$SNP))
  d$cromosoma <- as.integer(vapply(m, function(x) if (length(x)) x[2] else NA_character_, ""))
  d$posicion  <- as.numeric(vapply(m, function(x) if (length(x)) x[3] else NA_character_, ""))
  d$z <- d$coef_bioclim / d$error_estandar        # estadístico t (con signo)
  # convergio puede venir como TRUE/FALSE o como texto "True"/"False"
  d$convergio <- toupper(as.character(d$convergio)) %in% c("TRUE", "T")
  d$degenerado <- !d$convergio | d$error_estandar > 1 | d$dispersion > 50
  n_deg <- sum(d$degenerado, na.rm = TRUE)
  if (excluir_degenerados) d <- d[!d$degenerado, ]
  d <- d[is.finite(d$p_valor) & !is.na(d$cromosoma), ]
  d$variable <- etiqueta
  # Bonferroni: p ajustado = min(1, p * m)
  d$p_bonf <- pmin(1, d$p_valor * m_pruebas)
  d$hit <- d$p_bonf < alfa
  attr(d, "m_pruebas") <- m_pruebas
  cat(sprintf("%-8s pruebas: %d | usados: %d (excluidos por ajuste inestable: %d) | %%p<0.05: %.1f | lambda: %.2f | hits Bonferroni: %d\n",
              etiqueta, m_pruebas, nrow(d), if (excluir_degenerados) n_deg else 0L,
              100 * mean(d$p_valor < 0.05), lambda_gc(d$p_valor), sum(d$hit)))
  d
}
res <- Map(preparar, archivos_glm, names(archivos_glm))

# ------------------------
# 3. figuras por variable
# -------------------------
figuras_por_variable <- function(d, etiqueta) {
  m_pruebas <- attr(d, "m_pruebas")
  umbral_p  <- alfa / m_pruebas               # p mínimo para superar Bonferroni
  umbral_y  <- -log10(umbral_p)
  lam <- lambda_gc(d$p_valor)
  d$observado <- -log10(d$p_valor)
  d$grupo_hit <- factor(ifelse(d$hit, "si", "no"), levels = c("no", "si"))
  etiqueta_umbral <- "Bonferroni"
  
  
  #  QQ-plot 
  if (hacer_qq) {
    dq <- d[order(d$p_valor), ]
    m <- nrow(dq)
    dq$esperado <- -log10((seq_len(m) - 0.5) / m)
    lim <- max(dq$esperado, dq$observado)
    g1 <- ggplot(dq, aes(esperado, observado, colour = grupo_hit)) +
      geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed") +
      geom_point(size = 0.8, alpha = 0.7) +
      annotate("text", x = 0, y = lim, label = sprintf("lambda = %.2f", lam),
               hjust = 0, vjust = 1, size = 4.2, colour = "black") +
      scale_colour_manual(values = c(no = verde, si = morado), guide = "none") +
      coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
      labs(title = paste0(etiqueta, ": QQ-plot de p-valores"),
           x = expression("-log"[10]*"(p) esperado por azar"),
           y = expression("-log"[10]*"(p) observado")) + tema
    ggsave(file.path(carpeta_figuras, paste0(etiqueta, "_1_qq.png")), g1, width = 6.5, height = 6.8, dpi = 150, bg = "white")
  }
  
  #  Manhattan plot
  largo <- tapply(d$posicion, d$cromosoma, max)
  inicio <- c(0, cumsum(as.numeric(largo)))[seq_along(largo)]
  names(inicio) <- names(largo)
  d$x <- d$posicion + inicio[as.character(d$cromosoma)]
  centros <- inicio + as.numeric(largo) / 2
  d$grupo <- ifelse(d$hit, "hit", ifelse(d$cromosoma %% 2 == 0, "par", "impar"))
  d$grupo <- factor(d$grupo, levels = c("impar", "par", "hit"))
  g2 <- ggplot(d, aes(x, observado, colour = grupo)) +
    geom_point(size = 0.7, alpha = 0.8) +
    geom_hline(aes(yintercept = umbral_y, linetype = "u"), colour = "black") +
    scale_colour_manual(values = c(impar = verde_1, par = verde_2, hit = morado), guide = "none") +
    scale_linetype_manual(values = c(u = "dashed"), labels = etiqueta_umbral, name = NULL) +
    scale_x_continuous(breaks = centros, labels = names(largo), expand = c(0.01, 0)) +
    labs(title = paste0(etiqueta, ": Manhattan plot"),
         x = "Cromosoma", y = expression("-log"[10]*"(p)")) +
    tema + theme(panel.grid.major.x = element_blank())
  ggsave(file.path(carpeta_figuras, paste0(etiqueta, "_2_manhattan.png")), g2, width = 11, height = 4.8, dpi = 150, bg = "white")
  
  # Volcano plot 
  if (hacer_volcan) {
    g3 <- ggplot(d, aes(coef_bioclim, observado, colour = grupo_hit)) +
      geom_point(size = 0.7, alpha = 0.7) +
      geom_vline(xintercept = 0, colour = col_texto, linewidth = 0.3) +
      geom_hline(aes(yintercept = umbral_y, linetype = "u"), colour = "black") +
      scale_colour_manual(values = c(no = verde, si = morado), guide = "none") +
      scale_linetype_manual(values = c(u = "dashed"), labels = etiqueta_umbral, name = NULL) +
      labs(title = paste0(etiqueta, ": tamaño del efecto contra significancia"),
           x = "coef_bioclim (cambio en log-odds por +1 unidad de la variable)",
           y = expression("-log"[10]*"(p)")) + tema
    ggsave(file.path(carpeta_figuras, paste0(etiqueta, "_3_volcan.png")), g3, width = 7.5, height = 6, dpi = 150, bg = "white")
  }
}
invisible(Map(figuras_por_variable, res, names(res)))

# ------------------------------------------------------------------------------
# 4. figuras entre variables (>= 2 variables)
# ------------------------------------------------------------------------------
if (length(res) >= 2) {
  vars <- names(res)
  
  #  λ y proporción de p < 0.05 por variable 
  resumen <- do.call(rbind, lapply(res, function(d) data.frame(
    variable = d$variable[1], n_snps = nrow(d),
    prop_p05 = mean(d$p_valor < 0.05), lambda = lambda_gc(d$p_valor),
    hits_bonferroni = sum(d$hit))))
  write.csv(resumen, file.path(carpeta_figuras, "resumen_por_variable.csv"), row.names = FALSE)
  print(resumen, row.names = FALSE)
  g4 <- ggplot(resumen, aes(lambda, reorder(variable, lambda))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "black") +
    geom_point(colour = verde, size = 3.5) +
    geom_text(aes(label = sprintf("%.2f  (p<0.05: %.0f%%)", lambda, 100 * prop_p05)),
              hjust = -0.15, size = 3.4, colour = col_texto) +
    scale_x_continuous(limits = c(0, max(resumen$lambda) * 1.5)) +
    labs(title = "λ por variable bioclimática", x = "λ (lambda)", y = NULL) + tema
  ggsave(file.path(carpeta_figuras, "lambda_por_variable.png"), g4, width = 7, height = 1.5 + 0.35 * length(vars), dpi = 150, bg = "white")
  
  # utilidades para matrices (variables arriba y a la izquierda)
  ejes <- function(df) {
    df$v1 <- factor(df$v1, levels = rev(vars)); df$v2 <- factor(df$v2, levels = vars); df
  }
  tam <- function(k) c(4.5 + 0.4 * k, 4.2 + 0.4 * k)
  
  # Correlación entre los resultados del GLM de cada variable 
  snps_comunes <- Reduce(intersect, lapply(res, function(d) d$SNP))
  zmat <- sapply(res, function(d) d$z[match(snps_comunes, d$SNP)])
  colnames(zmat) <- vars
  cc <- cor(zmat, use = "pairwise.complete.obs")
  cl <- ejes(expand.grid(v1 = vars, v2 = vars, stringsAsFactors = FALSE))
  cl$r <- as.vector(cc[cbind(match(cl$v1, vars), match(cl$v2, vars))])
  cl$txt <- ifelse(abs(cl$r) > 0.7, "white", "#0b0b0b")
  g5 <- ggplot(cl, aes(v2, v1, fill = r)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", r), colour = txt), size = 3.4) +
    scale_colour_identity() +
    scale_fill_gradient2(low = amarillo, mid = "#f1f0ee", high = morado, midpoint = 0,
                         limits = c(-1, 1), name = "r") +
    scale_x_discrete(position = "top") + coord_equal() +
    labs(title = "Correlación entre los resultados del GLM de cada variable", x = NULL, y = NULL) +
    tema + theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = "right")
  ggsave(file.path(carpeta_figuras, "correlacion_resultados_glm.png"), g5, width = tam(length(vars))[1], height = tam(length(vars))[2], dpi = 150, bg = "white")
  
  #  SNPs significativos (Bonferroni) compartidos entre variables 
  hits <- lapply(res, function(d) d$SNP[d$hit])
  tr <- ejes(expand.grid(v1 = vars, v2 = vars, stringsAsFactors = FALSE))
  tr$n <- mapply(function(a, b) length(intersect(hits[[a]], hits[[b]])), as.character(tr$v1), as.character(tr$v2))
  tr$frac <- mapply(function(a, b) length(intersect(hits[[a]], hits[[b]])) / max(1, length(hits[[a]])), as.character(tr$v1), as.character(tr$v2))
  tr$txt <- ifelse(tr$frac > 0.55, "white", "#0b0b0b")
  g6 <- ggplot(tr, aes(v2, v1, fill = frac)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = n, colour = txt), size = 3.2) +
    scale_colour_identity() +
    scale_fill_gradient(low = "#eef7f1", high = "#1f7a52", limits = c(0, 1), name = "fracción\nde la fila") +
    scale_x_discrete(position = "top") + coord_equal() +
    labs(title = "SNPs significativos compartidos entre variables", x = NULL, y = NULL) +
    tema + theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = "right")
  ggsave(file.path(carpeta_figuras, "snps_compartidos.png"), g6, width = tam(length(vars))[1], height = tam(length(vars))[2], dpi = 150, bg = "white")
  
  # tabla: traslape observado vs el que se esperaría por puro azar
  # (si los hits de dos variables fueran independientes: n_A * n_B / N)
  N <- length(snps_comunes)
  hits_c <- lapply(hits, function(h) intersect(h, snps_comunes))
  pares <- expand.grid(var_A = vars, var_B = vars, stringsAsFactors = FALSE)
  pares <- pares[match(pares$var_A, vars) < match(pares$var_B, vars), ]
  pares$hits_A <- sapply(pares$var_A, function(v) length(hits_c[[v]]))
  pares$hits_B <- sapply(pares$var_B, function(v) length(hits_c[[v]]))
  pares$compartidos <- mapply(function(a, b) length(intersect(hits_c[[a]], hits_c[[b]])), pares$var_A, pares$var_B)
  pares$esperado_por_azar <- round(pares$hits_A * pares$hits_B / N, 1)
  pares$razon_obs_esp <- round(pares$compartidos / pmax(pares$esperado_por_azar, 0.1), 2)
  pares$r_resultados <- round(cc[cbind(match(pares$var_A, vars), match(pares$var_B, vars))], 2)
  write.csv(pares, file.path(carpeta_figuras, "traslape_observado_vs_azar.csv"), row.names = FALSE)
} else {
  cat("Solo hay 1 variable en 'archivos_glm': se omiten las figuras entre variables.\n")
}
cat("Figuras guardadas en:", carpeta_figuras, "\n")
