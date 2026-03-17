################################################################################
# Samantha M. Pacheco Gómez
# 17 marzo 2026
# LD_juguete_10SNPs.r
# este script se supone que toma 10 SNPs aleatorios y busca su LD entre ellos
# para ver que funciona y cómo lo hace
################################################################################

# cargar librerías

#solo necesito esta 
library(SNPRelate)

#pero para no jugarle al vivo y que algo falle, no vaya a ser
library(dplyr)
library(terra)
library(geodata)
library(adegenet)
library(ggplot2)

# 1. Abrir el archivO gds
genofile = snpgdsOpen("/mnt/data/sur/users/spacheco/data/teosinte/teosinte.gds")

# 2. Elegimos 10 SNPs al azar
set.seed(123)
snps_disponibles = read.gdsn(index.gdsn(genofile, "snp.id"))
diez_snps = sample(snps_disponibles, 10)

# 3. Calculamos la matriz de LD (usando r como medida de correlación)
matriz_ld = snpgdsLDMat(genofile, snp.id = diez_snps, method = "r")

# 4. Elevamos al cuadrado para obtener r^2 (0 = nada de LD, 1 = clones)
r_cuadrado = matriz_ld$LD^2

# Le ponemos nombres para entender qué vemos
colnames(r_cuadrado) = diez_snps
rownames(r_cuadrado) = diez_snps

print(round(r_cuadrado, 3))

snpgdsClose(genofile)


#para guardarlo en imagen

#abrimos un "túnel" hacia un archivo PDF
  pdf("/mnt/data/sur/users/spacheco/results/Mapa_LD_juguete_10SNPs.pdf")

# corremos la función de la imagen (la misma que ya usaste)
image(t(r_cuadrado[nrow(r_cuadrado):1,]), 
      main="Mapa de Calor de LD (10 SNPs)", 
      col=rev(heat.colors(20)))

# CERRAMOS el túnel (¡SUPER IMPORTANTE!)
dev.off()
