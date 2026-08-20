# ============================================================
# diagnóstico: ¿hay individuos con genotipo pero sin metadata,
# o con metadata pero sin genotipo? (full join, a diferencia del
# left join que usa convertir_a_bayescenv.R)
# ============================================================

Data <- read.csv2("data_teosinte.csv", header = TRUE, sep = ",", as.is = FALSE)
fam  <- read.table("T3604_33929_all.fam", header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

# all = TRUE es lo que hace que sea FULL join: conserva TODO de ambos lados,
# aunque no encuentren pareja (a diferencia de all.x = TRUE, que solo
# garantiza conservar el lado izquierdo).
full_check <- merge(fam, Data, by.x = "IID", by.y = "Sample_name", all = TRUE)

sin_metadata <- full_check[is.na(full_check$Accession) & !is.na(full_check$FID), "IID"]
sin_genotipo <- full_check[is.na(full_check$FID), "IID"]

cat("Individuos CON genotipo pero SIN metadata:", length(sin_metadata), "\n")
print(sin_metadata)

cat("\nFilas de metadata SIN genotipo correspondiente:", length(sin_genotipo), "\n")
print(sin_genotipo)

# Si alguna lista no está vacía, busca coincidencias aproximadas (typos,
# espacios de más/menos) contra el otro lado -- agrep() es una función base
# de R para "aproximate grep" (no necesita paquetes extra).
if (length(sin_metadata) > 0) {
  cat("\nPosibles coincidencias aproximadas (revisa a mano, no es automático):\n")
  for (id in sin_metadata) {
    candidatos <- agrep(id, Data$Sample_name, max.distance = 0.1, value = TRUE)
    if (length(candidatos) > 0) cat(id, "-> posible:", paste(candidatos, collapse = ", "), "\n")
  }
}

