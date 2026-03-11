################################################################################
# Samantha M. Pacheco Gómez
# 11 marzo 2026
# DAPC_juguete_50PCs_v2.r
# formato script para división por clusters k= 5
# data teosinte
# versión cluster PRUEBA con 50 PCs
# cambio: usar data del clima desde mi carpeta local
################################################################################

# cargar librerías

library(dplyr)
library(terra)
library(geodata)
library(adegenet)
library(SNPRelate)
library(ggplot2)

# --- 1. CONFIGURACIÓN DE RUTAS ---
# Definimos las rutas como variables para que sea fácil cambiarlas después
path_teosinte = "/mnt/data/sur/users/spacheco/data/teosinte/"
path_clima    = "/mnt/data/sur/users/spacheco/data/worldclim/climdata_30s/"
path_resultados = "/mnt/data/sur/users/spacheco/resultados/"

# --- 2. CARGA DE METADATA Y CLIMA ---
Data = read.csv2(paste0(path_teosinte, "data_teosinte.csv"), 
                 header = TRUE, sep=",", as.is = FALSE)

# Convertimos coordenadas a numérico para que terra no de problemas
Data$Latitude  = as.numeric(as.character(Data$Latitude))
Data$Longitude = as.numeric(as.character(Data$Longitude))

# --- Lectura directa de archivos locales ---
# archivos .tif en mi carpeta climdata_30s
listado_tifs = list.files(path_clima, pattern = ".tif$", full.names = TRUE)

# Es vital ordenarlos (1, 2, 3...) para que coincidan con tus 'nombres_chidos'
# Esta línea extrae el número del nombre del archivo y los ordena numéricamente
listado_tifs = listado_tifs[order(as.numeric(gsub("\\D", "", basename(listado_tifs))))]

# Cargamos todos los TIFs como un solo objeto 'SpatRaster'
clima_raw = rast(listado_tifs)

nombres_chidos = c("Temp_Media_Anual", "Rango_Diurno", "Isotermalidad", 
                   "Estacionalidad_Temp", "Max_Temp_Mes_Calido", "Min_Temp_Mes_Frio", 
                   "Rango_Anual_Temp", "Temp_Cuartil_Humedo", "Temp_Cuartil_Seco", 
                   "Temp_Cuartil_Calido", "Temp_Cuartil_Frio", "Prec_Anual", 
                   "Prec_Mes_Humedo", "Prec_Mes_Seco", "Estacionalidad_Prec", 
                   "Prec_Cuartil_Humedo", "Prec_Cuartil_Seco", "Prec_Cuartil_Calido", 
                   "Prec_Cuartil_Frio")
names(clima_raw) = nombres_chidos

# Extraemos el clima para nuestras poblaciones
geo_pop = Data %>%
  group_by(POB_CODE) %>%
  summarise(lat = mean(Latitude, na.rm=TRUE),
            lon = mean(Longitude, na.rm=TRUE),
            .groups="drop") %>%
  filter(!is.na(lat) & !is.na(lon))

valores_climaticos = extract(clima_raw, geo_pop[, c("lon", "lat")])
geo_pop_clim = cbind(geo_pop, valores_climaticos)
Data_Final = left_join(Data, geo_pop_clim %>% select(-lat, -lon, -ID), by = "POB_CODE")

# --- 3. CONVERSIÓN Y LECTURA GENÉTICA ---
# Solo crea el GDS si no existe todavía en la carpeta teosinte
archivo_gds = paste0(path_teosinte, "teosinte.gds")

if (!file.exists(archivo_gds)) {
  snpgdsBED2GDS(bed.fn = paste0(path_teosinte, "T3604_33929_all.bed"), 
                fam.fn = paste0(path_teosinte, "T3604_33929_all.fam"), 
                bim.fn = paste0(path_teosinte, "T3604_33929_all.bim"), 
                out.gdsfn = archivo_gds)
}

genofile = snpgdsOpen(archivo_gds)
geno_mat = snpgdsGetGeno(genofile)
samp_ids = read.gdsn(index.gdsn(genofile, "sample.id"))
snp_ids  = read.gdsn(index.gdsn(genofile, "snp.id"))
snpgdsClose(genofile)

gl_teosinte = new("genlight", geno_mat)
indNames(gl_teosinte) = samp_ids
locNames(gl_teosinte) = snp_ids

# --- 4. DAPC DIAGNÓSTICO (PRUEBA DE 50 PCs) ---
set.seed(147) # Tu semilla favorita

# Buscamos 5 clusters usando solo 50 PCs para que sea rápido
grupos_test = find.clusters(gl_teosinte, n.pca = 50, n.clust = 5, choose.n.clust = FALSE)

# Corremos el DAPC también con 50 PCs
dapc_test = dapc(gl_teosinte, pop = grupos_test$grp, n.pca = 50, n.da = 4)

# --- 5. GUARDAR RESULTADOS ---
# pdf con nombre significativo segun yo jiji
pdf(paste0(path_resultados, "DAPC_juguete_50PCs.pdf"), width = 10, height = 8)

# Gráfica del nicho térmico por taxón
print(ggplot(Data_Final, aes(x=Taxon, y=Temp_Media_Anual, fill=Taxon)) +
        geom_boxplot() + theme_minimal() +
        labs(title="Diferencias de Temperatura Media (30s res)", y="Temp (°C)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Gráfica del DAPC (La constelación de puntos)
scatter(dapc_test, col = spectral(5), scree.da = FALSE, bg = "white", pch = 20, 
        main = "DAPC Diagnóstico (50 PCs - Semilla 147)")

dev.off()

# Guardamos la tabla final que une Genética + Clima de 30s
pertenencia = as.data.frame(dapc_test$posterior)
pertenencia$Sample_name = rownames(pertenencia)
Data_Master = left_join(Data_Final, pertenencia, by = "Sample_name")

# Guardamos el .RData también en la carpeta de resultados
save(Data_Master, dapc_test, file = paste0(path_resultados, "DAPC_juguete_50PCs.RData"))


