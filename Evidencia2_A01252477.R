#Jordan Barba
#A01252477
#Evidencia 2
#Analisis de Biologia Computacional
#Profesor Heriberto Garcia Coronado

#Librerias
library(ggplot2)
library(seqinr)
library(gridExtra)
library(ape)
library(Biostrings)
library(viridis)
library(DECIPHER)
library(phangorn)
library(phytools)
library(geiger)
library(RSQLite)
library(ade4)

#Directorio en el que estaremos trabajando
setwd("C:/Users/jorda/OneDrive/Documents/Rstudio codes")

#Leer los archivos en .fasta
MT324062 <- read.fasta("MT324062.fasta")[[1]] #South Africa
MN908947 <- read.fasta("MN908947.fasta")[[1]] #China
MT020781 <- read.fasta("MT020781.fasta")[[1]] #Finland
MT039888 <- read.fasta("MT039888.fasta")[[1]] #USA
OM737994 <- read.fasta("OM737994.fasta")[[1]] #Palau
ON148315 <- read.fasta("ON148315.fasta")[[1]] #Brazil
ON248600 <- read.fasta("ON248600.fasta")[[1]] #Canis lupus familiaris/India
MZ914594 <- read.fasta("MZ914594.fasta")[[1]] #Canis lupus familiaris/Spain
OL790397 <- read.fasta("OL790397.fasta")[[1]] #Odocoileus virginianus/USA
OK555092 <- read.fasta("OK555092.fasta")[[1]] #Felis catus/Thailand

#Lista de secuencias
coronavirus_Sequences <- c(MT324062, MN908947, MT020781, MT039888,
                          OM737994, ON148315, ON248600, MZ914594,
                          OL790397, OK555092)

#Funcion que genera una lista con las bases de nucleotidos, la longitud, el contenido GC, y parte de la secuencia complementaria
Descripcion <- function(Secuencia){
  HEAD <- head(comp(Secuencia))
  HEAD1 <- paste(HEAD, collapse = "")
  #lista con los datos de la secuencia, 1. count, 2. length, 3. GC, 4. head
  lista_componentes <- list(count(Secuencia, wordsize=1), length(Secuencia), GC(Secuencia), HEAD1)
  return(lista_componentes)
}

#lista con los nombres de las variantes (Para facilitar el uso de dataframes)
NAMES <- c("MT324062", "MN908947", "MT020781", "MT039888",
           "OM737994", "ON148315", "ON248600", "MZ914594",
           "OL790397", "OK555092")

#Se aplica la funcion 'Descrpicion' y se guarda una variable para cada variante
MT324062_list <- Descripcion(MT324062)
MN908947_list <- Descripcion(MN908947)
MT020781_list <- Descripcion(MT020781)
MT039888_list <- Descripcion(MT039888)
OM737994_list <- Descripcion(OM737994)
ON148315_list <- Descripcion(ON148315)
ON248600_list <- Descripcion(ON248600)
MZ914594_list <- Descripcion(MZ914594)
OL790397_list <- Descripcion(OL790397)
OK555092_list <- Descripcion(OK555092)

#vector con los contenidos de Adenina de cada variante
contenidos_A <- c(MT324062_list[[1]][[1]], MN908947_list[[1]][[1]], MT020781_list[[1]][[1]], MT039888_list[[1]][[1]],
                  OM737994_list[[1]][[1]], ON148315_list[[1]][[1]], ON248600_list[[1]][[1]], MZ914594_list[[1]][[1]],
                  OL790397_list[[1]][[1]], OK555092_list[[1]][[1]])

#vector con los contenidos de Timina de cada variante
contenidos_T <- c(MT324062_list[[1]][[2]], MN908947_list[[1]][[2]], MT020781_list[[1]][[2]], MT039888_list[[1]][[2]],
                  OM737994_list[[1]][[2]], ON148315_list[[1]][[2]], ON248600_list[[1]][[2]], MZ914594_list[[1]][[2]],
                  OL790397_list[[1]][[2]], OK555092_list[[1]][[2]])

#vector con los contenidos de Guanina de cada variante
contenidos_G <- c(MT324062_list[[1]][[3]], MN908947_list[[1]][[3]], MT020781_list[[1]][[3]], MT039888_list[[1]][[3]],
                  OM737994_list[[1]][[3]], ON148315_list[[1]][[3]], ON248600_list[[1]][[3]], MZ914594_list[[1]][[3]],
                  OL790397_list[[1]][[3]], OK555092_list[[1]][[3]])

#vector con los contenidos de Citosina de cada variante
contenidos_C <- c(MT324062_list[[1]][[4]], MN908947_list[[1]][[4]], MT020781_list[[1]][[4]], MT039888_list[[1]][[4]],
                  OM737994_list[[1]][[4]], ON148315_list[[1]][[4]], ON248600_list[[1]][[4]], MZ914594_list[[1]][[4]],
                  OL790397_list[[1]][[4]], OK555092_list[[1]][[4]])

#vector con la longitud de cada variante
LENGTHS <- c(MT324062_list[[2]], MN908947_list[[2]], MT020781_list[[2]], MT039888_list[[2]],
             OM737994_list[[2]], ON148315_list[[2]], ON248600_list[[2]], MZ914594_list[[2]],
             OL790397_list[[2]], OK555092_list[[2]])

##vector con los contenidos de GC de cada variante
CONTENIDO_GC <- c(MT324062_list[[3]], MN908947_list[[3]], MT020781_list[[3]], MT039888_list[[3]],
                  OM737994_list[[3]], ON148315_list[[3]], ON248600_list[[3]], MZ914594_list[[3]],
                  OL790397_list[[3]], OK555092_list[[3]])

#vector con los primeros nucleotidos de la secuencia complementaria de cada variante
COMPLEMENTARIA <- c(MT324062_list[[4]], MN908947_list[[4]], MT020781_list[[4]], MT039888_list[[4]],
                    OM737994_list[[4]], ON148315_list[[4]], ON248600_list[[4]], MZ914594_list[[4]],
                    OL790397_list[[4]], OK555092_list[[4]])

#Tabla comparativa de las bases de cada variante
tabla_bases <- data.frame(Variante = NAMES, Adenina = contenidos_A, Timina = contenidos_T,
                          Guanina = contenidos_G, Citosina = contenidos_C)
tabla_bases

#Tabla comparativa con la longitud, el contenido GC, y los primeros nucleotidos del complementario de cada variante
tabla1_ev01 <- data.frame(Virus = NAMES, length = LENGTHS, Contenido_GC = CONTENIDO_GC, Complementaria = COMPLEMENTARIA)
tabla1_ev01

#Grafica de barras del contenido de GC y longitud de cada variante
plot_tabla_ev01_GC <- ggplot(tabla1_ev01, aes(x=Virus, y=Contenido_GC, fill = Virus)) + coord_cartesian(ylim = c(0.375, 0.385)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_bar(stat="identity") 
plot_tabla_ev01_length <- ggplot(tabla1_ev01, aes(x=Virus, y=length, fill = Virus)) + coord_cartesian(ylim = c(29700, 29950)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_bar(stat="identity")
plot_tabla_ev01_GC
plot_tabla_ev01_length

#Contenido de nucleotidos transformado a dataframe para facilitar su graficacion
MT324062_count <- as.data.frame(MT324062_list[[1]])
MN908947_count <- as.data.frame(MN908947_list[[1]])
MT020781_count <- as.data.frame(MT020781_list[[1]])
MT039888_count <- as.data.frame(MT039888_list[[1]])
OM737994_count <- as.data.frame(OM737994_list[[1]])
ON148315_count <- as.data.frame(ON148315_list[[1]])
ON248600_count <- as.data.frame(ON248600_list[[1]])
MZ914594_count <- as.data.frame(MZ914594_list[[1]])
OL790397_count <- as.data.frame(OL790397_list[[1]])
OK555092_count <- as.data.frame(OK555092_list[[1]]) 

#Definicion de las graficas de bases de cada variante
plot_MT324062_count <- ggplot(MT324062_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("MT324062")
plot_MN908947_count <- ggplot(MN908947_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("MN908947")
plot_MT020781_count <- ggplot(MT020781_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("MT020781")
plot_MT039888_count <- ggplot(MT039888_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("MT039888")
plot_OM737994_count <- ggplot(OM737994_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("OM737994")
plot_ON148315_count <- ggplot(ON148315_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("ON148315")
plot_ON248600_count <- ggplot(ON248600_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("ON248600")
plot_MZ914594_count <- ggplot(MZ914594_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("MZ914594")
plot_OL790397_count <- ggplot(OL790397_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("OL790397")
plot_OK555092_count <- ggplot(OK555092_count, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = 'identity') + xlab("OK555092")

#Muestreo de las graficas de bases de cada variante
grid.arrange(plot_MT324062_count, plot_MN908947_count, plot_MT020781_count, plot_MT039888_count, ncol = 2)
grid.arrange(plot_OM737994_count, plot_ON148315_count, plot_ON248600_count, plot_MZ914594_count, ncol = 2)
grid.arrange(plot_OL790397_count, plot_OK555092_count)

coronavirus_secuencias <- c("MT324062", "MN908947", "MT020781", "MT039888",
                "OM737994", "ON148315", "ON248600", "MZ914594",
                "OL790397", "OK555092")

#lectura de las secuencias utilizando read.GenBank
virus_sequences <- read.GenBank(coronavirus_secuencias)
length(virus_sequences)
class(virus_sequences)
typeof(virus_sequences)
str(virus_sequences)

#Creacion de documento en .fasta utilizando la secuencia leida en GenBank
write.dna(virus_sequences, file = 'coronavirus_secuencias.fasta', format = "fasta")

#Lectura de secuencias
virus_seq_not_align <- readDNAStringSet("coronavirus_secuencias.fasta", format = "fasta")

#Orientacion de nucleotidos de las secuencias
virus_seq_not_align <- OrientNucleotides(virus_seq_not_align)

#Alineacion de las secuencias
virus_seq_align <- AlignSeqs(virus_seq_not_align)

#Funcion para observar la secuencia formateado por colores, para compararlas
BrowseSeqs(virus_seq_align)

#Genera archivo de salida en nuestro directorio
writeXStringSet(virus_seq_align, file = "coronavirus_align.fasta")

#Lectura del archivo de salida generado anteriormente
virus_aligned <- read.alignment("coronavirus_align.fasta", format = "fasta")

#Generacion de matriz de distancia entre las 10 secuencias
matriz_distancia <- dist.alignment(virus_aligned, matrix = "similarity")
matriz_distancia

#Transformacion de la matriz de distancia a dataframe
temp <- as.data.frame(as.matrix(matriz_distancia))
class(temp)
typeof(temp)
dim(temp)
str(temp)

#Ilustracion de la matriz de distancia en forma de tabla
table.paint(temp, cleg=0, clabel.row=.5, clabel.col = .5) + scale_color_viridis()

#Creacion y ploteo de arbol filogenetico
virus_tree <- nj(matriz_distancia)
class(virus_tree)
virus_tree <- ladderize(virus_tree)
plot(virus_tree)
