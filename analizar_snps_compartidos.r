# ==============================================================================
# analizar_snps_compartidos.r
# con las variables ya escogidas quiero ver
#   (1) cuántos hits tiene cada una
#   (2) qué SNPs son hit en más de una, y cuánto se traslapan frente a lo esperado por azar
#   (3) figuras 
# se leen de variables_escogidas.csv (salida de decidir_variables.r).
# 
# ==============================================================================

library(ggplot2)


# 1. rutas y parámetros

carpeta_glm <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/"
prefijo_glm <- "glm_10pc_"          # los archivos se llaman glm_10pc_bio1.csv, glm_10pc_bio2.csv, ...
todas_las_variables <- paste0("bio", 1:19)
archivos_glm <- setNames(file.path(carpeta_glm, paste0(prefijo_glm, todas_las_variables, ".csv")), todas_las_variables)

# lista de variables escogidas (salida de decidir_variables.r) y correlación entre variables
ruta_variables   <- "/home/sam/Documents/sur_ecoevo_lab/exp/nov_2026/03_nov/seleccion_variables_4/variables_escogidas.csv"
ruta_correlacion <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/bioclim/correlacion_bioclim_poblaciones.csv"

carpeta_salida  <- "/home/sam/Documents/sur_ecoevo_lab/exp/nov_2026/03_nov/snps_compartidos"
carpeta_figuras <- file.path(carpeta_salida, "figuras")

alfa <- 0.05                  # nivel de significancia antes de corregir por pruebas múltiples
excluir_degenerados <- TRUE   # quita SNPs con ajuste inestable (igual que en visualizar_glm.r)
k_resaltar <- 4               # los SNPs que son hit en al menos k variables se listan y se rotulan en la figura 4
mostrar_figuras <- TRUE       # TRUE = además de guardarlas, imprimirlas en RStudio
dir.create(carpeta_figuras, recursive = TRUE, showWarnings = FALSE)

# colores: paleta viridis (Smith & van der Walt; Garnier et al.), diseñada para ser legible con daltonismo
# y en escala de grises; scales ya viene instalado con ggplot2
v5 <- scales::viridis_pal(option = "D")(5)            # morado oscuro, azul, verde azulado, verde, amarillo
col_obs <- v5[2]; col_esp <- v5[4]; col_gris <- "#b8b8b4"; col_texto <- "#52514e"
tema <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = col_texto),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

mostrar_y_guardar <- function(g, nombre, ancho, alto) {
  ggsave(file.path(carpeta_figuras, nombre), g, width = ancho, height = alto, dpi = 150, bg = "white")
  if (mostrar_figuras) print(g)
}


# 2. leer archivos

lambda_gc <- function(p) {      # factor de inflación genómica (Devlin & Roeder 1999)
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

preparar <- function(ruta) {
  d <- read.csv(ruta, stringsAsFactors = FALSE)
  m_pruebas <- sum(is.finite(d$p_valor))          # número de pruebas hechas (para Bonferroni)
  d$convergio <- toupper(as.character(d$convergio)) %in% c("TRUE", "T")
  d$degenerado <- !d$convergio | d$error_estandar > 1 | d$dispersion > 50
  if (excluir_degenerados) d <- d[!d$degenerado, ]
  d <- d[is.finite(d$p_valor), ]
  d$p_bonf <- pmin(1, d$p_valor * m_pruebas)      # Bonferroni: p ajustado = min(1, p * m)
  d$hit <- d$p_bonf < alfa
  d[, c("SNP", "coef_bioclim", "error_estandar", "p_valor", "p_bonf", "hit")]
}

escogidas <- read.csv(ruta_variables, stringsAsFactors = FALSE)
variables_conservadas <- escogidas$variable
stopifnot(length(variables_conservadas) >= 2, all(variables_conservadas %in% todas_las_variables))
cat("Variables escogidas (de decidir_variables.r):", paste(variables_conservadas, collapse = ", "), "\n")

res <- lapply(archivos_glm[variables_conservadas], preparar)   # solo se leen las escogidas
r_clima <- as.matrix(read.csv(ruta_correlacion, row.names = 1, check.names = FALSE))


# 3. hits de cada variable escogida

tabla_var <- data.frame(
  variable = variables_conservadas,
  grupo = escogidas$grupo,
  hits_bonferroni = sapply(res, function(d) sum(d$hit)),
  prop_p05 = sapply(res, function(d) mean(d$p_valor < 0.05)),
  lambda = sapply(res, function(d) lambda_gc(d$p_valor)),
  row.names = NULL)
write.csv(tabla_var, file.path(carpeta_salida, "hits_variables_escogidas.csv"), row.names = FALSE)
cat("\n=== Variables escogidas: hits, proporción con p < 0.05 y lambda ===\n")
cat(sprintf("  %-6s (grupo %s) hits: %4d | p<0.05: %4.1f%% | lambda: %.2f\n",
            tabla_var$variable, tabla_var$grupo, tabla_var$hits_bonferroni, 100 * tabla_var$prop_p05, tabla_var$lambda), sep = "")
cat("(una lambda mucho mayor a 1 indica p-valores inflados; los hits de variables con lambda alta hay que leerlos con cautela)\n")


# 4. traslape entre las variables escogidas
# solo se comparan los SNPs que se pudieron analizar en TODAS las variables escogidas

cons <- res
snps_comunes <- Reduce(intersect, lapply(cons, function(d) d$SNP))
N <- length(snps_comunes)
hits <- lapply(cons, function(d) intersect(d$SNP[d$hit], snps_comunes))
cat(sprintf("\n=== Hits sobre %d SNPs analizados en todas las variables escogidas ===\n", N))
for (v in variables_conservadas) cat(sprintf("  %-6s %4d hits\n", v, length(hits[[v]])))

pares <- t(combn(variables_conservadas, 2))
pares <- data.frame(var_A = pares[, 1], var_B = pares[, 2], stringsAsFactors = FALSE)
pares$r_clima <- round(mapply(function(a, b) r_clima[a, b], pares$var_A, pares$var_B), 2)
pares$hits_A <- sapply(pares$var_A, function(v) length(hits[[v]]))
pares$hits_B <- sapply(pares$var_B, function(v) length(hits[[v]]))
pares$compartidos <- mapply(function(a, b) length(intersect(hits[[a]], hits[[b]])), pares$var_A, pares$var_B)
pares$esperado_por_azar <- round(pares$hits_A * pares$hits_B / N, 1)     # si los hits de A y B fueran independientes
pares$razon_obs_esp <- round(pares$compartidos / pmax(pares$esperado_por_azar, 0.1), 2)
write.csv(pares, file.path(carpeta_salida, "traslape_variables_escogidas.csv"), row.names = FALSE)


# 5. SNPs que son hit en más de una variable escogida

hit_matriz <- sapply(hits, function(h) snps_comunes %in% h)           # SNPs x variables (TRUE/FALSE)
n_variables <- rowSums(hit_matriz)
p_bonf_matriz <- sapply(cons, function(d) d$p_bonf[match(snps_comunes, d$SNP)])

# cuántos SNPs son hit en exactamente k variables, contra lo que se esperaría si las variables fueran independientes
# (distribución de una suma de Bernoullis con probabilidad hits_v / N, cada una distinta)
p_v <- sapply(hits, length) / N
dist_azar <- 1
for (p in p_v) dist_azar <- c(dist_azar * (1 - p), 0) + c(0, dist_azar * p)
distribucion <- data.frame(n_variables_hit = 0:length(p_v),
                           snps_observados = as.integer(table(factor(n_variables, levels = 0:length(p_v)))),
                           snps_esperados_por_azar = round(N * dist_azar, 1))
distribucion$razon_obs_esp <- round(distribucion$snps_observados / pmax(distribucion$snps_esperados_por_azar, 0.1), 2)
write.csv(distribucion, file.path(carpeta_salida, "distribucion_snps_por_numero_de_variables.csv"), row.names = FALSE)
cat("\n=== SNPs según en cuántas variables escogidas son hit ===\n")
print(distribucion, row.names = FALSE)

m <- regmatches(snps_comunes, regexec("^S([0-9]+)_([0-9]+)", snps_comunes))
tabla_snps <- data.frame(
  SNP = snps_comunes,
  cromosoma = as.integer(vapply(m, function(x) if (length(x)) x[2] else NA_character_, "")),
  posicion = as.numeric(vapply(m, function(x) if (length(x)) x[3] else NA_character_, "")),
  n_variables_hit = n_variables,
  variables_hit = apply(hit_matriz, 1, function(f) paste(variables_conservadas[f], collapse = ", ")),
  mejor_p_bonf = apply(p_bonf_matriz, 1, min, na.rm = TRUE),
  stringsAsFactors = FALSE)
tabla_snps <- tabla_snps[tabla_snps$n_variables_hit >= 1, ]
tabla_snps <- tabla_snps[order(-tabla_snps$n_variables_hit, tabla_snps$mejor_p_bonf), ]
write.csv(tabla_snps, file.path(carpeta_salida, "snps_hit_variables_escogidas.csv"), row.names = FALSE)

cat(sprintf("\nSNPs que son hit en al menos 1 variable escogida: %d (este sería el conjunto candidato más grande)\n", nrow(tabla_snps)))
cat(sprintf("SNPs que son hit en 2 o más: %d\n", sum(tabla_snps$n_variables_hit >= 2)))
muchos <- tabla_snps[tabla_snps$n_variables_hit >= k_resaltar, c("SNP", "cromosoma", "posicion", "n_variables_hit", "variables_hit", "mejor_p_bonf")]
muchos <- muchos[order(muchos$cromosoma, muchos$posicion), ]
cat(sprintf("\nSNPs que son hit en al menos %d variables (ordenados por cromosoma y posición):\n", k_resaltar))
print(muchos, row.names = FALSE)


# 6. figuras

# figura 1: hits por variable escogida
d1 <- tabla_var
d1$etiqueta <- sprintf("%d  (lambda %.1f)", d1$hits_bonferroni, d1$lambda)
d1$variable <- factor(d1$variable, levels = rev(d1$variable))
g1 <- ggplot(d1, aes(hits_bonferroni, variable)) +
  geom_col(width = 0.7, fill = col_obs) +
  geom_text(aes(label = etiqueta), hjust = -0.08, size = 3.2, colour = col_texto) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Hits por variable escogida", x = "hits (Bonferroni)", y = NULL) + tema
mostrar_y_guardar(g1, "fig1_hits_por_variable.png", 8, 1.5 + 0.45 * nrow(d1))

# figura 2: SNPs compartidos entre las variables escogidas (cuántos y cuántas veces lo esperado por azar)
d2 <- expand.grid(v1 = variables_conservadas, v2 = variables_conservadas, stringsAsFactors = FALSE)
d2$n <- NA_real_; d2$razon <- NA_real_
for (i in seq_len(nrow(pares))) {
  for (par in list(c(pares$var_A[i], pares$var_B[i]), c(pares$var_B[i], pares$var_A[i]))) {
    k <- d2$v1 == par[1] & d2$v2 == par[2]
    d2$n[k] <- pares$compartidos[i]; d2$razon[k] <- pares$razon_obs_esp[i]
  }
}
diag_k <- d2$v1 == d2$v2
d2$n[diag_k] <- sapply(d2$v1[diag_k], function(v) length(hits[[v]]))
d2$texto <- ifelse(diag_k, sprintf("%d\nhits", d2$n), sprintf("%d\n×%.1f", d2$n, d2$razon))
d2$v1 <- factor(d2$v1, levels = rev(variables_conservadas)); d2$v2 <- factor(d2$v2, levels = variables_conservadas)
d2$txt_col <- ifelse(!is.na(d2$razon) & d2$razon > 0.6 * max(d2$razon, na.rm = TRUE), "white", "#0b0b0b")
g2 <- ggplot(d2, aes(v2, v1, fill = razon)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = texto, colour = txt_col), size = 3, lineheight = 0.9) +
  scale_colour_identity() +
  scale_fill_viridis_c(direction = -1, na.value = "#f1f0ee", name = "veces lo esperado por azar") +
  scale_x_discrete(position = "top") + coord_equal() +
  labs(title = "SNPs compartidos entre variables escogidas", x = NULL, y = NULL) + tema +
  theme(panel.grid = element_blank(), legend.key.width = unit(1.4, "cm"))
mostrar_y_guardar(g2, "fig2_snps_compartidos_pares.png", 8, 7.5)

# figura 3: cuántos SNPs son hit en k variables (observado) contra lo esperado si las variables fueran independientes
# en escala log no se puede dibujar el 0: los puntos con 0 SNPs observados, o con menos de 0.1 esperados, se omiten
d3 <- distribucion[distribucion$n_variables_hit >= 1, ]
d3_obs <- d3[d3$snps_observados > 0, ]
d3_esp <- d3[d3$snps_esperados_por_azar >= 0.1, ]
d3_seg <- d3[d3$snps_observados > 0 & d3$snps_esperados_por_azar >= 0.1, ]
g3 <- ggplot() +
  geom_segment(data = d3_seg, aes(x = n_variables_hit, xend = n_variables_hit,
                                  y = snps_esperados_por_azar, yend = snps_observados), colour = "#c4cacf", linewidth = 0.8) +
  geom_point(data = d3_esp, aes(n_variables_hit, snps_esperados_por_azar, colour = "esperado por azar"), size = 3.5, shape = 18) +
  geom_point(data = d3_obs, aes(n_variables_hit, snps_observados, colour = "observado"), size = 3.2) +
  geom_text(data = d3_obs, aes(n_variables_hit, snps_observados, label = snps_observados), vjust = -1, size = 3, colour = col_texto) +
  scale_colour_manual(values = c(observado = col_obs, `esperado por azar` = col_esp), name = NULL) +
  scale_x_continuous(breaks = d3$n_variables_hit) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000), labels = scales::label_comma(accuracy = 0.1, drop0trailing = TRUE),
                expand = expansion(mult = c(0.1, 0.25))) +
  labs(title = "SNPs según el número de variables en que son hit",
       subtitle = "no se dibuja: 0 SNPs observados, ni esperados por azar < 0.1",
       x = "número de variables escogidas en que el SNP es hit", y = "número de SNPs (escala log)") + tema
mostrar_y_guardar(g3, "fig3_snps_por_numero_de_variables.png", 8, 5)

# figura 4: dónde están en el genoma los SNPs hit en 2 o más variables
d4 <- tabla_snps[tabla_snps$n_variables_hit >= 2 & !is.na(tabla_snps$cromosoma), ]
d4$posicion_mb <- d4$posicion / 1e6
d4$cromosoma_txt <- factor(paste("chr", d4$cromosoma), levels = paste("chr", sort(unique(d4$cromosoma))))
d4_resaltar <- d4[d4$n_variables_hit >= k_resaltar, ]
# el rótulo se ancla hacia adentro cuando el SNP está cerca de un borde del cromosoma, para que no se corte
max_cromosoma <- tapply(d4$posicion_mb, d4$cromosoma_txt, max)
frac <- d4_resaltar$posicion_mb / max_cromosoma[as.character(d4_resaltar$cromosoma_txt)]
d4_resaltar$hjust <- ifelse(frac > 0.8, 1, ifelse(frac < 0.2, 0, 0.5))
g4 <- ggplot(d4, aes(posicion_mb, n_variables_hit)) +
  geom_point(aes(colour = n_variables_hit >= k_resaltar), size = 1.6, alpha = 0.8) +
  geom_text(data = d4_resaltar, aes(label = SNP, hjust = hjust), size = 2.2, vjust = -0.8, colour = "#0b0b0b", check_overlap = TRUE) +
  facet_wrap(~ cromosoma_txt, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = c(`FALSE` = col_gris, `TRUE` = col_obs), guide = "none") +
  scale_y_continuous(breaks = 2:length(variables_conservadas), expand = expansion(mult = c(0.05, 0.2))) +
  labs(title = "Posición de los SNPs que son hit en 2 o más variables",
       x = "posición en el cromosoma (Mb)", y = "número de variables") + tema
mostrar_y_guardar(g4, "fig4_posicion_snps_compartidos.png", 13, 6)

cat("\nArchivos en:", carpeta_salida, "\nFiguras en:", carpeta_figuras, "\n")