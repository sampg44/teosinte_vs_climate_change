# ==============================================================================
# decidir_variables.r
# decidir qué variable bioclim conservar de cada grupo de variables correlacionadas.
# 
#   1) por grupo calculo centralidad = |r| promedio de la variable con las demás de su grupo
#      (Dormann et al. 2013 recomiendan quedarse con la variable central de cada grupo)
#   2) si hay empate (grupos de 2, o diferencias menores a tolerancia_empate), tendrÍa que buscar info. biolÓgica
#   3) con la lista variables_escogidas revisar que quede una por grupo, que ningún par tenga |r| >= umbral_r
#      y que el VIF del conjunto sea menor a umbral_vif
#
# correrlo una vez con variables_escogidas <- NULL para ver las candidatas; luego llenar variables_escogidas y volver a correrlo.
# figuras: 1) correlación por grupo, 2) VIF, 3) dendrograma (cómo se formaron los grupos y cuáles se escogieron)
# Guarda variables_escogidas.csv, que lee el sig script que se llama: analizar_snps_compartidos.r
#
# ==============================================================================

library(ggplot2)


# 1. rutas y parámetros

ruta_bioclim <- "/home/sam/Documents/sur_ecoevo_lab/data/worldclim/bioclim_por_pop.csv"   # el mismo de correlacion_bioclim.r
metodo       <- "pearson"                                                                  # el mismo de correlacion_bioclim.r
ruta_grupos  <- "/home/sam/Documents/sur_ecoevo_lab/exp/oct_2026/19_oct/bioclim/grupos_bioclim_correlacionadas.csv"

carpeta_salida  <- "/home/sam/Documents/sur_ecoevo_lab/exp/nov_2026/03_nov/seleccion_variables_4"
carpeta_figuras <- file.path(carpeta_salida, "figuras")

umbral_r   <- 0.7    # el mismo de correlacion_bioclim.r
umbral_vif <- 10     # convención (Dormann et al. 2013): VIF > 10 = colinealidad fuerte; es una regla práctica, no un teorema

# dos variables de un grupo cuentan como empatadas si sus centralidades difieren menos que esto
# (0.02 es un valor que puse yo para tus datos, donde todas valen ~0.9: cambiarlo si hace falta)
tolerancia_empate <- 0.02

# lista final (una variable por grupo). Criterio de Sam: variables que se repiten en otros estudios en teocintle.
# Dejarla en NULL para ver solo las candidatas por centralidad
variables_escogidas <- c("bio2", "bio3", "bio5", "bio12", "bio15", "bio17", "bio18", "bio19")

mostrar_figuras <- TRUE       # TRUE = además de guardarlas, imprimirlas en RStudio

etiquetas <- c(bio1  = "Temperatura media anual",
               bio2  = "Rango diurno medio",
               bio3  = "Isotermalidad (bio2 / bio7)",
               bio4  = "Estacionalidad de la temperatura",
               bio5  = "Temperatura máxima del mes más cálido",
               bio6  = "Temperatura mínima del mes más frío",
               bio7  = "Rango anual de temperatura",
               bio8  = "Temperatura media del trimestre más húmedo",
               bio9  = "Temperatura media del trimestre más seco",
               bio10 = "Temperatura media del trimestre más cálido",
               bio11 = "Temperatura media del trimestre más frío",
               bio12 = "Precipitación anual",
               bio13 = "Precipitación del mes más húmedo",
               bio14 = "Precipitación del mes más seco",
               bio15 = "Estacionalidad de la precipitación",
               bio16 = "Precipitación del trimestre más húmedo",
               bio17 = "Precipitación del trimestre más seco",
               bio18 = "Precipitación del trimestre más cálido",
               bio19 = "Precipitación del trimestre más frío")

dir.create(carpeta_figuras, recursive = TRUE, showWarnings = FALSE)

# colores: paleta viridis (Smith & van der Walt; Garnier et al.), diseñada para ser legible con daltonismo
# y en escala de grises; scales ya viene instalado con ggplot2
v5 <- scales::viridis_pal(option = "D")(5)            # morado oscuro, azul, verde azulado, verde, amarillo
col_dentro <- v5[2]                                    # ramas de uniones dentro de un grupo (dendrograma)
col_bien <- v5[3]; col_mal <- v5[5]                    # VIF aceptable / VIF > 10
col_texto <- "#52514e"

tema <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = col_texto),
        legend.position = "top")

mostrar_y_guardar <- function(g, nombre, ancho, alto) {
  ggsave(file.path(carpeta_figuras, nombre), g, width = ancho, height = alto, dpi = 150, bg = "white")
  if (mostrar_figuras) print(g)
}


# 2. leer archivos y ordenar las variables por grupo

bio <- read.csv(ruta_bioclim, stringsAsFactors = FALSE)
cols_bio <- grep("^bio[0-9]+$", colnames(bio), value = TRUE)
cols_bio <- cols_bio[order(as.integer(sub("bio", "", cols_bio)))]
X <- bio[complete.cases(bio[, cols_bio]), cols_bio]       # solo poblaciones con las 19 variables (así la matriz es válida para el VIF)
cat(sprintf("Poblaciones usadas: %d (descartadas por NA: %d)\n", nrow(X), nrow(bio) - nrow(X)))
R <- cor(X, method = metodo)
if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-10) cat("[AVISO] la matriz de correlación es singular: hay variables perfectamente correlacionadas\n")

grupos <- read.csv(ruta_grupos, stringsAsFactors = FALSE)
stopifnot(setequal(grupos$variable, colnames(R)))
grupos <- grupos[order(grupos$grupo, as.integer(sub("bio", "", grupos$variable))), ]
R <- R[grupos$variable, grupos$variable]


# 3. por variable: centralidad y VIF
# centralidad_j = promedio de |r_jk| sobre las otras variables k de su grupo (NA si el grupo tiene 1 sola)
# VIF_j = 1 / (1 - R2_j) = diagonal de la inversa de la matriz de correlación
#   R2_j = qué parte de la variable j se puede predecir con las demás; VIF = 1 si no se parece a ninguna otra

tabla <- data.frame(variable = grupos$variable, grupo = grupos$grupo,
                    etiqueta = unname(etiquetas[grupos$variable]), stringsAsFactors = FALSE)
tabla$centralidad <- sapply(tabla$variable, function(v) {
  otras <- setdiff(tabla$variable[tabla$grupo == tabla$grupo[tabla$variable == v]], v)
  if (length(otras) == 0) NA_real_ else mean(abs(R[v, otras]))
})
tabla$vif_19_variables <- diag(solve(R))[tabla$variable]
if (any(tabla$vif_19_variables < 1)) cat("[AVISO] hay VIF menores a 1: la matriz está mal condicionada, no confiar en el VIF de las 19 variables\n")

# candidata = variable que gana por centralidad (o empata dentro de la tolerancia)
tabla$candidata <- FALSE
for (g in unique(tabla$grupo)) {
  i <- tabla$grupo == g
  tabla$candidata[i] <- if (sum(i) == 1) TRUE else tabla$centralidad[i] >= max(tabla$centralidad[i]) - tolerancia_empate
}
write.csv(tabla, file.path(carpeta_salida, "tabla_variables.csv"), row.names = FALSE)

cat("\n=== Variables por grupo (* = gana o empata por centralidad) ===\n")
for (g in unique(tabla$grupo)) {
  t <- tabla[tabla$grupo == g, ]
  cat(sprintf("\ngrupo %d\n", g))
  cat(sprintf("  %s %-6s centralidad: %s | VIF (19 variables): %8.1f | %s\n",
              ifelse(t$candidata, "*", " "), t$variable,
              ifelse(is.na(t$centralidad), "  -- ", sprintf("%.3f", t$centralidad)),
              t$vif_19_variables, t$etiqueta), sep = "")
}

n_candidatas <- tapply(tabla$candidata, tabla$grupo, sum)
resueltos <- names(n_candidatas)[n_candidatas == 1]
empatados <- names(n_candidatas)[n_candidatas > 1]
cat("\n=== Resumen ===\n")
cat("Grupos con una sola ganadora (se pueden dejar como están):\n")
for (g in resueltos) cat(sprintf("  grupo %s: %s\n", g, tabla$variable[tabla$grupo == as.integer(g) & tabla$candidata]))
cat("Grupos EMPATADOS (tienes que decidir con criterio biológico):\n")
for (g in empatados) cat(sprintf("  grupo %s: %s\n", g, paste(tabla$variable[tabla$grupo == as.integer(g) & tabla$candidata], collapse = " o ")))


# 4. revisar la lista que tú escogiste

esc <- NULL
if (is.null(variables_escogidas)) {
  cat("\nvariables_escogidas está en NULL: llénala con una variable por grupo y vuelve a correr.\n")
} else {
  stopifnot(all(variables_escogidas %in% tabla$variable), !anyDuplicated(variables_escogidas))
  esc <- tabla[match(variables_escogidas, tabla$variable), ]
  
  cat("\n=== Revisión de variables_escogidas ===\n")
  conteo <- table(factor(esc$grupo, levels = sort(unique(tabla$grupo))))
  if (any(conteo == 0)) cat("[AVISO] grupos sin variable escogida:", paste(names(conteo)[conteo == 0], collapse = ", "), "\n")
  if (any(conteo > 1)) cat("[AVISO] grupos con más de una variable escogida:", paste(names(conteo)[conteo > 1], collapse = ", "), "\n")
  if (all(conteo == 1)) cat("OK: una variable por grupo.\n")
  
  Rsub <- R[variables_escogidas, variables_escogidas]
  altas <- which(abs(Rsub) >= umbral_r & upper.tri(Rsub), arr.ind = TRUE)
  if (nrow(altas) == 0) {
    cat(sprintf("OK: ningún par de variables escogidas tiene |r| >= %.2f (máximo = %.2f).\n", umbral_r, max(abs(Rsub[upper.tri(Rsub)]))))
  } else {
    cat(sprintf("[AVISO] pares con |r| >= %.2f:\n", umbral_r))
    for (k in seq_len(nrow(altas))) cat(sprintf("  %s - %s: r = %.2f\n", rownames(Rsub)[altas[k, 1]], colnames(Rsub)[altas[k, 2]], Rsub[altas[k, 1], altas[k, 2]]))
  }
  
  vif_conjunto <- tryCatch(diag(solve(Rsub))[variables_escogidas], error = function(e) rep(NA_real_, length(variables_escogidas)))
  esc$vif_en_conjunto <- unname(vif_conjunto)
  if (anyNA(vif_conjunto)) {
    cat("[AVISO] no se pudo calcular el VIF (matriz singular: hay variables perfectamente correlacionadas)\n")
  } else if (max(vif_conjunto) < umbral_vif) {
    cat(sprintf("OK: VIF del conjunto escogido entre %.2f y %.2f (umbral %d).\n", min(vif_conjunto), max(vif_conjunto), umbral_vif))
  } else {
    cat(sprintf("[AVISO] VIF por encima de %d en: %s\n", umbral_vif, paste(variables_escogidas[vif_conjunto >= umbral_vif], collapse = ", ")))
  }
  print(esc[, c("variable", "grupo", "centralidad", "vif_en_conjunto", "etiqueta")], row.names = FALSE)
  
  write.csv(esc[, c("variable", "grupo", "etiqueta", "centralidad", "vif_en_conjunto")],
            file.path(carpeta_salida, "variables_escogidas.csv"), row.names = FALSE)
  cat("\nGuardado:", file.path(carpeta_salida, "variables_escogidas.csv"), "\n")
}


# 5. figuras
# (si las ves en el panel de RStudio se ven más apretadas que el PNG guardado: abre el PNG o amplía el panel)

resaltar <- if (is.null(variables_escogidas)) tabla$variable[tabla$candidata] else variables_escogidas
nombre_resaltar <- if (is.null(variables_escogidas)) "candidata" else "escogida"
ord <- tabla$variable

# figura 1: correlacion entre variables, ordenadas por grupo (las lineas negras separan los grupos)
n_ord <- length(ord)
d1 <- expand.grid(v1 = ord, v2 = ord, stringsAsFactors = FALSE)
d1$r <- R[cbind(d1$v1, d1$v2)]
d1$v1 <- factor(d1$v1, levels = rev(ord)); d1$v2 <- factor(d1$v2, levels = ord)
d1$txt <- ifelse(d1$r < 0.2, "white", "#0b0b0b")   # texto blanco sobre los tonos oscuros de viridis (r bajo)
fronteras <- head(cumsum(table(tabla$grupo)), -1) + 0.5
# los nombres de las variables se dibujan como texto (para poder poner en negritas las escogidas)
etq <- data.frame(variable = ord, pos = seq_len(n_ord), resaltada = ord %in% resaltar, stringsAsFactors = FALSE)
etq$cara <- ifelse(etq$resaltada, "bold", "plain"); etq$color <- ifelse(etq$resaltada, "#0b0b0b", col_texto)
g1 <- ggplot(d1, aes(v2, v1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r), colour = txt), size = 2.6) +
  geom_vline(xintercept = fronteras, linewidth = 0.9, colour = "#0b0b0b") +
  geom_hline(yintercept = n_ord - fronteras + 1, linewidth = 0.9, colour = "#0b0b0b") +
  geom_text(data = etq, aes(x = 0.4, y = n_ord + 1 - pos, label = variable, fontface = cara, colour = color), hjust = 1, size = 3.6, inherit.aes = FALSE) +
  geom_text(data = etq, aes(x = pos, y = 0.4, label = variable, fontface = cara, colour = color), angle = 90, hjust = 1, size = 3.6, inherit.aes = FALSE) +
  scale_colour_identity() +
  scale_fill_viridis_c(option = "D", limits = c(-1, 1), name = "r (correlación de Pearson)",
                       guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  scale_x_discrete(expand = expansion(add = c(1.6, 0.5))) + scale_y_discrete(expand = expansion(add = c(1.6, 0.5))) +
  coord_equal() +
  labs(title = "Correlación entre variables bioclim, por grupo",
       subtitle = paste0("líneas negras = separan los grupos; negritas = variable ", nombre_resaltar), x = NULL, y = NULL) + tema +
  theme(panel.grid = element_blank(), axis.text = element_blank(), legend.key.width = unit(1.6, "cm"))
mostrar_y_guardar(g1, "fig1_correlacion_por_grupo.png", 10, 10)

# figura 2: VIF de las variables escogidas (solo tiene sentido con la lista final)
# VIF = 1 / (1 - R2), con R2 = parte de la variable que se puede predecir con las OTRAS variables escogidas
#   VIF 1 = nada (0%); VIF 2 = 50%; VIF 5 = 80%; VIF 10 = 90%
if (is.null(esc) || anyNA(esc$vif_en_conjunto)) {
  cat("\n(fig2 omitida: necesita variables_escogidas y un VIF calculable)\n")
} else {
  d2 <- data.frame(variable = esc$variable, vif = esc$vif_en_conjunto, stringsAsFactors = FALSE)
  d2 <- d2[order(d2$vif), ]
  d2$variable <- factor(d2$variable, levels = d2$variable)
  d2$pct <- round(100 * (1 - 1 / d2$vif))
  d2$nivel <- factor(ifelse(d2$vif > umbral_vif, "VIF > 10 (muy redundante)", "VIF <= 10 (aceptable)"),
                     levels = c("VIF <= 10 (aceptable)", "VIF > 10 (muy redundante)"))
  x_max <- max(umbral_vif * 1.15, max(d2$vif) * 1.35)
  marcas <- c(1, 2, 5, 10)
  g2 <- ggplot(d2, aes(vif, variable, fill = nivel)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = umbral_vif, linetype = "dashed", colour = "#0b0b0b") +
    geom_text(aes(label = sprintf("%.1f  (%d%%)", vif, pct)), hjust = -0.1, size = 3.4, colour = "#0b0b0b") +
    scale_fill_manual(values = setNames(c(col_bien, col_mal), levels(d2$nivel)), name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = marcas, name = "VIF",
                       sec.axis = sec_axis(~ ., breaks = marcas, labels = paste0(round(100 * (1 - 1 / marcas)), "%"),
                                           name = "% de la variable que las otras escogidas explican")) +
    coord_cartesian(xlim = c(1, x_max)) +
    labs(title = "Qué tan repetida está cada variable escogida",
         subtitle = paste(c("Cada barra es el VIF de una variable frente a las otras escogidas.",
                            "Junto a la barra: el VIF y, entre paréntesis, el % de esa variable",
                            "que se puede predecir con las demás. Línea punteada: VIF = 10."), collapse = "\n"),
         y = NULL) + tema +
    theme(axis.title.x.top = element_text(colour = col_texto, size = 10))
  mostrar_y_guardar(g2, "fig2_vif_escogidas.png", 8, 5)
}

# figura 3: dendrograma. Primero se ve cómo se formaron los grupos y despues cuáles se escogieron
# misma construccion que correlacion_bioclim.r: distancia 1 - |r|, enlace completo, corte en 1 - umbral_r
# (con enlace completo, la altura a la que se unen dos ramas = 1 - el |r| mas bajo entre las variables de esas ramas;
#  por eso, del corte hacia las hojas todos los pares de un grupo tienen |r| >= umbral_r)
hc <- hclust(as.dist(1 - abs(R)), method = "complete")
corte <- 1 - umbral_r
grupo_hc <- cutree(hc, h = corte)
mismos_grupos <- length(unique(grupo_hc)) == length(unique(tabla$grupo)) &&
  all(tapply(grupo_hc[tabla$variable], tabla$grupo, function(g) length(unique(g))) == 1)
if (!mismos_grupos) cat("[AVISO] el corte del dendrograma NO reproduce los grupos de grupos_bioclim_correlacionadas.csv (revisa umbral_r y metodo)\n")

# el eje usa la raiz cuadrada de la distancia (x = 1 - sqrt(distancia)) para separar las uniones cercanas a |r| = 1;
# los numeros del eje son |r| verdaderos
x_de <- function(h) 1 - sqrt(h)
n_v <- length(hc$labels)
y_hoja <- n_v - match(hc$labels, hc$labels[hc$order]) + 1       # posicion vertical de cada variable (la primera arriba)
grupo_hoja <- tabla$grupo[match(hc$labels, tabla$variable)]
y_nodo <- numeric(n_v - 1); g_nodo <- rep(NA_integer_, n_v - 1)
y_de <- function(k) if (k < 0) y_hoja[-k] else y_nodo[k]
h_de <- function(k) if (k < 0) 0 else hc$height[k]
g_de <- function(k) if (k < 0) grupo_hoja[-k] else g_nodo[k]
ramas <- vector("list", 3 * (n_v - 1))
for (i in seq_len(n_v - 1)) {
  a <- hc$merge[i, 1]; b <- hc$merge[i, 2]; h <- hc$height[i]
  ya <- y_de(a); yb <- y_de(b)
  y_nodo[i] <- (ya + yb) / 2
  ga <- g_de(a); gb <- g_de(b)
  g_nodo[i] <- if (!is.na(ga) && !is.na(gb) && ga == gb) ga else NA_integer_
  dentro <- h <= corte + 1e-9 && !is.na(g_nodo[i])
  col_rama <- if (dentro) col_dentro else "#bdbdb8"      # color = union dentro de un grupo; gris = union entre grupos
  grosor <- if (dentro) 1 else 0.6
  ramas[[3 * i - 2]] <- data.frame(x = x_de(h_de(a)), xend = x_de(h), y = ya, yend = ya, col = col_rama, lw = grosor)
  ramas[[3 * i - 1]] <- data.frame(x = x_de(h_de(b)), xend = x_de(h), y = yb, yend = yb, col = col_rama, lw = grosor)
  ramas[[3 * i]]     <- data.frame(x = x_de(h), xend = x_de(h), y = ya, yend = yb, col = col_rama, lw = grosor)
}
ramas <- do.call(rbind, ramas)

hojas <- data.frame(variable = hc$labels, y = y_hoja, grupo = grupo_hoja, stringsAsFactors = FALSE)
hojas$escogida <- hojas$variable %in% resaltar
hojas$cara <- ifelse(hojas$escogida, "bold", "plain")
# un corchete por grupo, a la derecha de los nombres, con el numero del grupo escrito una sola vez
corchetes <- do.call(rbind, lapply(sort(unique(hojas$grupo)), function(g) {
  yy <- hojas$y[hojas$grupo == g]; data.frame(grupo = g, y1 = min(yy) - 0.3, y2 = max(yy) + 0.3, ym = mean(yy))
}))

x_raiz <- x_de(max(hc$height)); x_corte <- x_de(corte)
r_marcas <- c(1, 0.9, 0.8, 0.7, 0.5, 0.3, 0)
g3 <- ggplot() +
  geom_segment(data = ramas, aes(x = x, xend = xend, y = y, yend = yend, colour = col, linewidth = lw), lineend = "round") +
  geom_segment(aes(x = x_corte, xend = x_corte, y = 0.4, yend = n_v + 0.6), linetype = "dashed", colour = "#0b0b0b") +
  annotate("text", x = x_corte, y = n_v + 1.1, label = sprintf("corte: |r| = %.1f", umbral_r), size = 3.4, fontface = "bold") +
  geom_point(data = hojas[!hojas$escogida, ], aes(x = 1, y = y), shape = 21, fill = "white", colour = col_texto, size = 2.2) +
  geom_point(data = hojas[hojas$escogida, ], aes(x = 1, y = y), shape = 21, fill = "#0b0b0b", colour = "#0b0b0b", size = 4) +
  geom_text(data = hojas, aes(x = 1.04, y = y, label = variable, fontface = cara), hjust = 0, size = 3.8, colour = "#0b0b0b") +
  geom_segment(data = corchetes, aes(x = 1.2, xend = 1.2, y = y1, yend = y2), colour = col_texto, linewidth = 0.6) +
  geom_text(data = corchetes, aes(x = 1.23, y = ym, label = paste("grupo", grupo)), hjust = 0, size = 3.4, colour = col_texto) +
  scale_colour_identity() + scale_linewidth_identity() +
  scale_x_continuous(breaks = x_de(1 - r_marcas), labels = sprintf("%.1f", r_marcas), name = "|r| mínimo entre las variables unidas (1 = idénticas)") +
  coord_cartesian(xlim = c(x_raiz - 0.02, 1.4), ylim = c(0.3, n_v + 1.5), clip = "off") +
  labs(title = "Cómo se formaron los grupos de variables",
       subtitle = paste0("Las variables unidas con línea de color (a la derecha del corte) forman un grupo.\n",
                         "Círculo negro y negritas: variable ", nombre_resaltar, " de cada grupo."),
       caption = "Eje no lineal: hay más espacio cerca de 1.", y = NULL) + tema +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(),
        plot.caption = element_text(hjust = 0, colour = col_texto))
mostrar_y_guardar(g3, "fig3_dendrograma_grupos_y_escogidas.png", 9, 6.5)

cat("\nFiguras en:", carpeta_figuras, "\n")