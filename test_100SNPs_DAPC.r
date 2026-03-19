################################################################################
# test rápido: 100 SNPs DAPC 
# Sam Pacheco
# 18 marzo 2026 
# (Basado en DAPC_juguete_50PCs_v3.r )
################################################################################

library(dplyr)
library(terra)
library(geodata)
library(adegenet)
library(SNPRelate)
library(ggplot2)

# --- 1. CONFIGURACIÓN DE RUTAS (CORREGIDAS) ---
path_teosinte = "/mnt/data/sur/users/spacheco/data/teosinte/"
# CAMBIO AQUÍ: de 'resultados' a 'results'
path_resultados = "/mnt/data/sur/users/spacheco/results/" 
path_clima = "/mnt/data/sur/users/spacheco/data/worldclim/mesoamerica_30s.tif"


cat("Rutas configuradas y verificadas\n")

# --- 2. CARGA DE DATA (Igual que el original) ---
Data = read.csv2(paste0(path_teosinte, "data_teosinte.csv"), header = TRUE, sep=",", as.is = FALSE)
Data$Latitude  = as.numeric(as.character(Data$Latitude))
Data$Longitude = as.numeric(as.character(Data$Longitude))

# Cargamos el archivo único (automáticamente detecta las 19 capas)

clima_raw = rast(path_clima)

nombres_chidos = c("Temp_Media_Anual", "Rango_Diurno", "Isotermalidad", 
                   "Estacionalidad_Temp", "Max_Temp_Mes_Calido", "Min_Temp_Mes_Frio", 
                   "Rango_Anual_Temp", "Temp_Cuartil_Humedo", "Temp_Cuartil_Seco", 
                   "Temp_Cuartil_Calido", "Temp_Cuartil_Frio", "Prec_Anual", 
                   "Prec_Mes_Humedo", "Prec_Mes_Seco", "Estacionalidad_Prec", 
                   "Prec_Cuartil_Humedo", "Prec_Cuartil_Seco", "Prec_Cuartil_Calido", 
                   "Prec_Cuartil_Frio")
names(clima_raw) = nombres_chidos

cat("Data de clima mesoamerica cargada\n")

# --- 3. PROCESAMIENTO GEOGRÁFICO ---
# Primero calculamos los promedios por población

geo_pop = Data %>%
  group_by(POB_CODE) %>%
  summarise(lat = mean(Latitude, na.rm=TRUE),
            lon = mean(Longitude, na.rm=TRUE),
            .groups="drop") %>%
  filter(!is.na(lat) & !is.na(lon))

# Ahora extraemos los valores (esto será instantáneo ahora)
valores_climaticos = extract(clima_raw, geo_pop[, c("lon", "lat")])

# Unimos todo
geo_pop_clim = cbind(geo_pop, valores_climaticos)
Data_Final = left_join(Data, geo_pop_clim %>% select(-lat, -lon, -ID), by = "POB_CODE")

cat("procesamiento geográfico terminado\n")


# --- 4. LECTURA GENÉTICA Y SUBSET DE 100 SNPs ---
archivo_gds = paste0(path_teosinte, "teosinte.gds")
genofile = snpgdsOpen(archivo_gds)
# Obtenemos todos los IDs
all_snp_ids = read.gdsn(index.gdsn(genofile, "snp.id"))

# SELECCIÓN ALEATORIA DE 100 SNPs
set.seed(123)
subset_snps = sample(all_snp_ids, 100)

# Leemos solo esos 100 SNPs
geno_mat = snpgdsGetGeno(genofile, snp.id = subset_snps)
samp_ids = read.gdsn(index.gdsn(genofile, "sample.id"))
snpgdsClose(genofile)

# Creamos el genlight pequeño
gl_test = new("genlight", geno_mat)
indNames(gl_test) = samp_ids
locNames(gl_test) = subset_snps

cat("Submuestreo de 100 SNPs completado\n")

# --- 5. DAPC EXPRÉS ---
set.seed(147)
grupos_test = find.clusters(gl_test, n.pca = 50, n.clust = 5, choose.n.clust = FALSE)
dapc_res = dapc(gl_test, pop = grupos_test$grp, n.pca = 50, n.da = 4)

# --- 6. GUARDAR (LA PRUEBA REAL) ---
pdf_name = paste0(path_resultados, "test_100SNPs_DAPC.pdf")
pdf(pdf_name, width = 10, height = 8)
# Tu plot de ggplot y el scatter...
scatter(dapc_res, main="TEST 100 SNPs - Si ves esto, el path funciona")
dev.off()

cat("¡Prueba terminada! El archivo debería estar en:", pdf_name, "\n")


