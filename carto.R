#### carto ####

#fichier "Occurences_Alps.txt"

#A vérifier avec Julien mais la taxonomie c'est celle de Flora Europeae + Flora Mediterranea, P. cottia doit être sous le nom P. villosa et P. appenina sous le nom P. pedemontana

#récupération de données dans data
#install.packages("rgdal")
library(sp);library(raster)
library(rgdal)

WD = "~/Bureau/BEE/Stage/Pedemontana/data_carto"

plot.raster = function(rast,ext) {
  plot(crop(raster(rast), ext)) ;  lines(read.table("data_carto/WORLD_lowres.dat"))
}
plot.obj = function(rast) {
  plot(rast) ;  lines(read.table("data_carto/WORLD_lowres.dat"))
}


# création de cartes mensuelle de TMEAN et PRECMEAN (moy climatique 1970-2000) ####
LFprec = list.files(WD,pattern ="PRECm",full = T)
LFtm = list.files(WD,pattern = "TMEANm",full = T)
PRECm.01 = raster(LFprec[1]) #res is O.1666° = 10' ; coordinate in dec degres ; no projection
#graphics.off(); x11(15,10) #permet de fermer la fenetre graphique et afficher sur fenetre annexe
plot(raster(LFprec[2]), axes=T)
lines(read.table("data_carto/WORLD_lowres.dat")) #contour des continents

GRE = cbind(5.720,45.180) #coordonnées de grnoble
GRE
extract(PRECm.01,GRE) #extrait les prec moyenne de grenoble en janvier

PRECm.GRE = rep(NA,12)
for (i in 1:12) PRECm.GRE[i] = extract(raster(LFprec[i]),GRE)
barplot(PRECm.GRE, col = "darkblue", ylab = "Precipitation (mm)", xlab = "Mois",main = paste("Precipitation Grenoble (annuel :",sum(PRECm.GRE),"mm)",sep = " "),
        ylim = c(0,120), names = month.abb)

PREC = stack(LFprec) #empile les couches
extract(PREC,GRE) # meme chose qu'avant, mais sans boucle for

RF = raster("data_carto/RAINFOREST_10min.tif") #tif des rainforest
plot.obj(RF)
#mais rainforest en europe, donc peut etre besoin d'un seuil

RF80 = RF # même chose que précedement, mais en mode binaire
RF80[RF>=80] = 1
RF80[RF<80] = 0 #on peut aussi faire ça par raster::reclassify
plot.obj(RF80)

#quel est le profil thermique et hydrique des RF (a 80%) ####
FT80 = xyFromCell(RF,ID.RF80)
PREC = stack(LFprec) #empile les couches
TEMP = stack(LFtm)
FTprec = extract(PREC,FT80) # meme chose qu'avant, mais sans boucle for
FTtemp = extract(TEMP,FT80)

Pmean = apply(FTprec,2,mean, na.rm =T)
Tmean = apply(FTtemp,2,mean, na.rm =T)/10

par(mar=c(5,4,4,5))
barplot(Pmean, col = "cadetblue3", ylab = "Precipitation (mm)", xlab = "Mois",
        ylim = c(0,250), names = month.abb, las = 2)
par(new=T)
plot(Tmean, col = "darkred", xlab=NA, ylab=NA,main = paste("Diagramme Tropical (P.annuel :",sum(Pmean),"mm, T.mean :",mean(Tmean),"°)",sep = " "),
     ylim = c(0,125), type = "l", axes = F, lw = 3)
axis(side = 4) ; mtext(side = 4, line = 3, "Temperature (°)")

Alpine.ext = extent (-10,40,35,55)
plot.raster("data_carto/GRASSLANDS_10min.tif",Alpine.ext)

#plot(raster("data_carto/PIA03395.tif"),extent (0,21600,0,9049))
#plot.raster("data_carto/PIA03395.tif",extent (10800,12000,7350,8000))

