# ==============================================================================
# comparar_qq_binomial_quasi.r
# 
# compara el QQ-plot de la quasibinomial con el de un binomial puro (dispersión = 1)
# sin volver a correr ningún GLM: el binomial sale del mismo ajuste.
#
# quasibinomial y binomial dan el MISMO coeficiente; solo cambia el error estándar:
#   SE_quasi = SE_binomial * sqrt(dispersion)       (dispersion = chi2 de Pearson / grados de libertad residuales)
# así que  SE_binomial = error_estandar / sqrt(dispersion)  y  z = coef_bioclim / SE_binomial,  p = 2 * pnorm(-|z|)
# (el quasibinomial usa la distribución t con n - n_parametros grados de libertad; el binomial usa la normal)
# ==============================================================================

library(ggplot2)


# 1. rutas y parámetros

variable <- "bio1"
# salida de glm_frecuencia_bioclim.r CORRIDO CON usar_dispersion <- TRUE (quasibinomial)
ruta_glm <- paste0("/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/glm_10pc_", variable, ".csv")
carpeta_figuras <- "/home/sam/Documents/sur_ecoevo_lab/exp/nov_2026/03_nov/qq_binomial_vs_quasi"

n_parametros <- 12            # intercepto + bioclim + 10 PCs (para los grados de libertad del t: n_pop_usadas - n_parametros)
alfa <- 0.05
excluir_degenerados <- TRUE   # mismos SNPs excluidos que en visualizar_glm.r
ejes_iguales <- FALSE         # FALSE = cada panel con sus propios ejes (se lee mejor); TRUE = mismos ejes en los dos
mostrar_figuras <- TRUE       # TRUE = además de guardar, imprimir en RStudio
dir.create(carpeta_figuras, recursive = TRUE, showWarnings = FALSE)

verde <- "#2f9e6e"; morado <- "#6250d6"; col_texto <- "#52514e"
tema <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = col_texto),
        strip.text = element_text(face = "bold"))


# 2. leer y revisar

d <- read.csv(ruta_glm, stringsAsFactors = FALSE)
d$convergio <- toupper(as.character(d$convergio)) %in% c("TRUE", "T")
d$degenerado <- !d$convergio | d$error_estandar > 1 | d$dispersion > 50
m_pruebas <- sum(is.finite(d$p_valor))                      # número de pruebas hechas (para Bonferroni)
if (excluir_degenerados) d <- d[!d$degenerado, ]
d <- d[is.finite(d$p_valor) & is.finite(d$dispersion) & d$dispersion > 0, ]

if (all(abs(d$dispersion - 1) < 1e-9)) stop("La dispersión es 1 en todos los SNPs: este CSV parece de un binomial puro, no de un quasibinomial.")
cat(sprintf("SNPs usados: %d | dispersión: mediana %.2f, rango intercuartil %.2f a %.2f\n",
            nrow(d), median(d$dispersion), quantile(d$dispersion, .25), quantile(d$dispersion, .75)))

# comprobación: el p del quasibinomial recalculado debe coincidir con p_valor
t_val <- d$coef_bioclim / d$error_estandar
p_quasi_rec <- 2 * pt(-abs(t_val), df = d$n_pop_usadas - n_parametros)
desv <- max(abs(log10(p_quasi_rec) - log10(d$p_valor)), na.rm = TRUE)
if (desv > 1e-6) {
  cat(sprintf("el p recalculado no coincide con p_valor (diferencia máxima %.2g en log10). Revisar n_parametros (%d) y que el CSV sea del quasibinomial.\n", desv, n_parametros))
} else {
  cat("el p recalculado del quasibinomial coincide con p_valor; el binomial calculado con la misma fórmula es confiable.\n")
}


# 3. p-valor binomial (dispersión = 1) y lambdas
# se trabaja en -log10(p) calculado en escala logarítmica para no perder los p-valores diminutos

z_bin <- d$coef_bioclim / (d$error_estandar / sqrt(d$dispersion))
nlp_bin   <- -(log(2) + pnorm(-abs(z_bin), log.p = TRUE)) / log(10)
nlp_quasi <- -(log(2) + pt(-abs(t_val), df = d$n_pop_usadas - n_parametros, log.p = TRUE)) / log(10)

lambda_bin   <- median(z_bin^2) / qchisq(0.5, 1)          # misma definición: mediana de chi2 / mediana esperada (Devlin & Roeder 1999)
lambda_quasi <- {p <- d$p_valor[d$p_valor > 0]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)}
umbral_y <- -log10(alfa / m_pruebas)                      # umbral de Bonferroni en -log10(p)

resumen <- data.frame(
  modelo = c("quasibinomial", "binomial"),
  lambda = c(lambda_quasi, lambda_bin),
  prop_p05 = c(mean(nlp_quasi > -log10(0.05)), mean(nlp_bin > -log10(0.05))),
  hits_bonferroni = c(sum(nlp_quasi > umbral_y), sum(nlp_bin > umbral_y)))
cat(sprintf("\n(Bonferroni sobre %d pruebas: -log10(p) > %.2f)\n", m_pruebas, umbral_y))
print(resumen, row.names = FALSE, digits = 3)
write.csv(resumen, file.path(carpeta_figuras, paste0(variable, "_resumen_binomial_vs_quasi.csv")), row.names = FALSE)


# 4. QQ-plot lado a lado

cuantiles <- function(nlp, nombre, lam) {
  m <- length(nlp)
  o <- sort(nlp)                                             # menor a mayor; los esperados también van de menor a mayor
  data.frame(esperado = sort(-log10((seq_len(m) - 0.5) / m)), observado = o,
             hit = o > umbral_y,
             modelo = sprintf("%s (lambda = %.2f)", nombre, lam))
}
dq <- rbind(cuantiles(nlp_quasi, "quasibinomial", lambda_quasi),
            cuantiles(nlp_bin, "binomial (dispersión = 1)", lambda_bin))
dq$modelo <- factor(dq$modelo, levels = unique(dq$modelo))

g <- ggplot(dq, aes(esperado, observado, colour = hit)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed") +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_wrap(~ modelo, scales = if (ejes_iguales) "fixed" else "free") +
  scale_colour_manual(values = c(`FALSE` = verde, `TRUE` = morado), guide = "none") +
  labs(title = paste0(variable, ": QQ-plot de p-valores, quasibinomial vs binomial"),
       x = expression("-log"[10]*"(p) esperado por azar"), y = expression("-log"[10]*"(p) observado")) + tema
ggsave(file.path(carpeta_figuras, paste0(variable, "_qq_binomial_vs_quasi.png")), g, width = 11, height = 5.5, dpi = 150, bg = "white")
if (mostrar_figuras) print(g)
cat("\nFigura en:", carpeta_figuras, "\n")
