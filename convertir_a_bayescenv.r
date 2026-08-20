# ============================================================
# Convertir genotipos PLINK (.bed/.bim/.fam) al formato de
# entrada de BayeScEnv
#
# formato de salida esperado:
#   [loci]=100
#   [populations]=16
#   [pop]=1
#   1 40 2 29 11
#   2 40 2 4 36
#   ...
# columnas por línea: locus  cromosomas  n_alelos  conteo_alelo1  conteo_alelo2
#
# Usa Accession (276 niveles) como población, NO POB_CODE (74),
# consistente con la corrección hecha en basic_analysis_teosinte.pdf
#
#
# ============================================================

library(SNPRelate)
library(dplyr)

# --- 1. Rutas ---------
ruta_metadata <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/data_teosinte.csv"
ruta_fam      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.fam"
ruta_bim      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.bim"
ruta_bed      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/T3604_33929_all.bed"
ruta_gds      <- "/home/sam/Documents/sur_ecoevo_lab/data/teosinte/archivos/teosinte.gds"

ruta_salida_geno  <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/276_pob_bayesformato.txt"  # input para bayescenv
ruta_salida_orden <- "/home/sam/Documents/sur_ecoevo_lab/data/bayescenv/orden_pop.csv"             # mapping pop_index -> Accession

filtrar_monomorficos <- TRUE  # quita loci sin variación en TODO el dataset (recomendado, no obligatorio)

# Para un ejemplo juguete: pon aquí un vector de códigos Accession (p. ej. los
# que salgan de seleccionar_poblaciones_prueba.R) y solo esas poblaciones se
# escribirán. Deja NULL para correr con las 276 poblaciones completas.
poblaciones_incluir <- NULL  # ej: c("CIM10003", "JLHNM-661", ...)

# --- 2. Metadata y asignación de población (Accession, no POB_CODE) ------
Data <- read.csv2(ruta_metadata, header = TRUE, sep = ",", as.is = FALSE)

fam <- read.table(ruta_fam, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

meta_fam <- merge(fam, Data, by.x = "IID", by.y = "Sample_name", all.x = TRUE)
stopifnot(sum(is.na(meta_fam$Accession)) == 0)  # todos los individuos deben tener población asignada

# --- 3. Abrir (o crear) el GDS --------------------------------------------
if (!file.exists(ruta_gds)) {
  snpgdsBED2GDS(bed.fn = ruta_bed, bim.fn = ruta_bim, fam.fn = ruta_fam,
                out.gdsfn = ruta_gds)
}
genofile <- snpgdsOpen(ruta_gds)

samp   <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

poblacion <- meta_fam$Accession[match(samp, meta_fam$IID)]
stopifnot(sum(is.na(poblacion)) == 0)

# Orden fijo de poblaciones -- MUY IMPORTANTE: este mismo orden es el que
# vas a tener que usar después para el archivo de variables ambientales de
# BayeScEnv (una línea por población, en el mismo orden que aquí).
poblaciones <- sort(unique(poblacion))

if (!is.null(poblaciones_incluir)) {
  faltantes <- setdiff(poblaciones_incluir, poblaciones)
  if (length(faltantes) > 0) warning("Estos Accession no existen en los datos: ",
                                     paste(faltantes, collapse = ", "))
  poblaciones <- poblaciones[poblaciones %in% poblaciones_incluir]
}
n_pops <- length(poblaciones)
if (is.null(poblaciones_incluir)) stopifnot(n_pops == 276)

write.csv(data.frame(pop_index = seq_along(poblaciones), Accession = poblaciones),
          ruta_salida_orden, row.names = FALSE)
message("Orden de poblaciones guardado en: ", ruta_salida_orden,
        "  usar mismo orden para el archivo ambiental")

# --- 4. Extraer matriz de genotipos (individuos x SNPs) -------------------
# Valores: 0/1/2 = copias del alelo fijado por SNPRelate para ese locus
# NA = dato faltante.
# Memoria aproximada: n_individuos x n_snps x 8 bytes (~1 GB para 3604 x 33929).
# Si tu clúster tiene poca RAM disponible en el nodo, pide más memoria en el
# job (sbatch --mem=...) en vez de intentar reducir esta matriz.
geno <- snpgdsGetGeno(genofile, snp.id = snp_id, sample.id = samp,
                      with.id = FALSE, verbose = TRUE)
# dim(geno) debe ser: length(samp) x length(snp_id)

# --- 5. (Opcional) filtrar loci monomórficos en TODO el dataset -----------
if (filtrar_monomorficos) {
  af_total <- snpgdsSNPRateFreq(genofile)
  loci_polimorficos <- af_total$MinorFreq > 0
  message(sum(!loci_polimorficos), " loci monomórficos removidos de ",
          length(loci_polimorficos), " totales.")
  geno   <- geno[, loci_polimorficos, drop = FALSE]
  snp_id <- snp_id[loci_polimorficos]
}
n_loci <- length(snp_id)

# --- 5b. Diccionario de identidad real de los alelos -----------------------
# "conteo_alelo1"/"conteo_alelo2" en el archivo de salida NO son A1/A2 del
# .bim directamente -- hay ambigüedad documentada en cómo SNPRelate asigna
# "A"/"B" al convertir de PLINK a GDS (ver github.com/zhengxwen/SNPRelate,
# issues #4 y #84). Para no depender de adivinar esa conversión, leemos el
# campo snp.allele DEL MISMO GDS que usamos para contar -- así es
# garantizadamente consistente con los conteos 0/1/2, sin importar qué pasó
# internamente en la conversión.
bim <- read.table(ruta_bim, header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim")

alelos_gds <- read.gdsn(index.gdsn(genofile, "snp.allele"))  # formato "A/B", p.ej. "C/T"
alelos_split <- strsplit(alelos_gds, "/")

diccionario_alelos <- data.frame(
  locus_index_bayescenv = seq_len(n_loci),        # el índice que aparece en el archivo de BayeScEnv
  snp_id_gds            = snp_id,
  CHR                    = bim$CHR[match(snp_id, seq_len(nrow(bim)))],
  SNP                    = bim$SNP[match(snp_id, seq_len(nrow(bim)))],
  A1_bim                 = bim$A1_bim[match(snp_id, seq_len(nrow(bim)))],
  A2_bim                 = bim$A2_bim[match(snp_id, seq_len(nrow(bim)))],
  alelo1_conteo          = vapply(alelos_split, `[`, character(1), 1),  # el que se cuenta en "conteo_alelo1"
  alelo2_conteo          = vapply(alelos_split, `[`, character(1), 2)   # el que se cuenta en "conteo_alelo2"
)
write.csv(diccionario_alelos, "diccionario_alelos.csv", row.names = FALSE)
message("Diccionario de alelos (qué nucleótido es alelo1/alelo2 por locus) guardado en: diccionario_alelos.csv")
# Chequeo rápido de cordura: alelo1_conteo/alelo2_conteo deben coincidir con
# A1_bim/A2_bim en la mayoría de los loci -- si NO coinciden en muchos casos,
# eso confirma que SNPRelate sí hizo un swap respecto al bim, y hay que usar
# SIEMPRE diccionario_alelos.csv (nunca el bim original) para interpretar
# después cuáles nucleótidos están detrás de cada locus outlier.
mismatch <- mean(diccionario_alelos$alelo1_conteo != diccionario_alelos$A1_bim)
message(round(100 * mismatch, 1), "% de loci tienen alelo1_conteo distinto de A1_bim -- ",
        "usa siempre diccionario_alelos.csv, no el bim, para interpretar resultados.")

# --- 6. Escribir el archivo en formato BayeScan/BayeScEnv -----------------
con <- file(ruta_salida_geno, "w")
writeLines(paste0("[loci]=", n_loci), con)
writeLines(paste0("[populations]=", n_pops), con)

for (p_index in seq_along(poblaciones)) {
  p   <- poblaciones[p_index]
  idx <- which(poblacion == p)
  sub_geno <- geno[idx, , drop = FALSE]
  
  n_ind_por_locus <- colSums(!is.na(sub_geno))       # individuos genotipados en esta pob, por locus
  sample_size     <- 2 * n_ind_por_locus              # copias de gen = tamaño de muestra que pide el formato
  conteo_alelo1   <- colSums(sub_geno, na.rm = TRUE)  # suma de dosis = copias del "alelo 1"
  conteo_alelo2   <- sample_size - conteo_alelo1
  
  writeLines(paste0("[pop]=", p_index), con)
  bloque <- paste(seq_len(n_loci), sample_size, 2L, conteo_alelo1, conteo_alelo2)
  writeLines(bloque, con)
}
close(con)
snpgdsClose(genofile)

message("Archivo de genotipos para BayeScEnv escrito en: ", ruta_salida_geno)

# --- 7. Chequeo de cordura (sanity check) ---------------------------------
# Revisa a mano un locus/población al azar y confirma que los conteos suman
# el tamaño de muestra esperado, y que no quedaron poblaciones con sample
# size 0 en TODOS sus loci (indicaría que esa Accession no tiene genotipos).
chequeo <- read.table(ruta_salida_geno, skip = 2, fill = TRUE,
                      col.names = c("marca", paste0("V", 1:5)))
# no olvidar checarlo manualmente con una revisión rápida


































