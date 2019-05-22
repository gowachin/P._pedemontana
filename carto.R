#### carto ####
# carto raster WC ####
library(readr)
library(raster)
Location <- read_delim("Location.csv", "&", escape_double = FALSE, trim_ws = TRUE)
loc = cbind(Location$Longitude,Location$Latitude)
loc

WD = "~/Bureau/BEE/Stage/Pedemontana/wc2.0_30s_bio"
wc_bio = list.files(WD,pattern ="wc2.0_bio_30s_",full = T)
wc_bio

var = c("AMT","MDR","IsoT","TSeas","MxTWm","MnTCom","TAR","MTWeQ","MTDrQ","MTWaQ","MTCoQ","APre","PWem","PDrm","PSeas","PWeQ","PDrQ","PWaQ","PCoQ")
title = c("Annual Mean Temperature",
          "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
          "Isothermality (BIO2/BIO7) (* 100)",
          "Temperature Seasonality (standard deviation *100)",
          "Max Temperature of Warmest Month",
          "Min Temperature of Coldest Month",
          "Temperature Annual Range (BIO5-BIO6)",
          "Mean Temperature of Wettest Quarter",
          "Mean Temperature of Driest Quarter",
          "Mean Temperature of Warmest Quarter",
          "Mean Temperature of Coldest Quarter",
          "Annual Precipitation",
          "Precipitation of Wettest Month",
          "Precipitation of Driest Month",
          "Precipitation Seasonality (Coefficient of Variation)",
          "Precipitation of Wettest Quarter",
          "Precipitation of Driest Quarter",
          "Precipitation of Warmest Quarter",
          "Precipitation of Coldest Quarter"
)

wc_bio.1 = raster(wc_bio[1])
plot(wc_bio.1)
WC_bio = stack(wc_bio)

Alpine.ext = extent (-10,40,35,55)
Europ = crop(WC_bio,Alpine.ext)
for(i in 1:dim(Europ)[3]) {
  writeRaster(Europ[[i]], filename=paste("Europ_clim/wc2.0_bio_30s_Europ_",i,".tif",sep=""), format="GTiff", overwrite=TRUE)
}
loc_wc =extract(Europ,loc)
colnames(loc_wc)=var ; row.names(loc_wc) = Location$Code
loc_wc
library(ape)

draw_ACP=function(ACP,Xaxe,Yaxe,ind=F){
  library(ade4)
  par(mfrow = c(1, 1),mar = c(5, 4, 4, 2)) ; s.arrow(4*ACP$co,xax=Xaxe,yax=Yaxe,sub=paste("axe",as.character(Xaxe),"/axe",as.character(Yaxe),sep =" "))
  if (ind == T) s.label(acpT$li,clab = 0.7,xax=Xaxe,yax=Yaxe, add.plot = T) else print("no individual plot")
  par(new=TRUE)
  # Set plot layout
  layout(mat = matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2),
         heights = c(2, 1),    # Heights of the two rows
         widths = c(1, 2))     # Widths of the two columns
  par(mar = c(2, 2, 2, 2))
  color=vector("character",length(ACP$eig)) ; color[c(Xaxe,Yaxe)]="black" ; color[-c(Xaxe,Yaxe)]="darkgrey"
  barplot(ACP$eig,col=color,main ="Eigen Values")
}
draw_Eigen=function(ACP){
  par(mfrow = c(1, 1),mar = c(5, 4, 4, 2))
  label=round((ACP$eig/sum(ACP$eig))*100,2) ; label=paste(label,c(rep ("%",length(ACP$eig))), sep = "") #ajout de légendes stylées ma gueuule
  bb=barplot(ACP$eig,col="Darkgrey",main ="Eigen Values",ylab="Inertia")
  end_point = 0.5 + length(ACP$eig) + length(ACP$eig)-1
  text(seq(1,end_point,1.2), par("usr")[3]-0.1,srt = 60, adj= 1, xpd = TRUE,labels = label, cex=0.8)
  #lines(x = bb, y = ACP$eig)
}
draw_Eigen(acpT)

table = loc_wc[,c(1,2,4,7,9,14)]
acpT=dudi.pca(table, center=T, scale=T, scannf=F, nf=5)
inertia.dudi(acpT, col=T, row=T)$tot.inertia #garder 3 axes 83,6%
#contributions des variables aux axes (quelle variable et en quelles proportions)
inertia.dudi(acpT, col=T, row=T)$col.abs #axe 1= morphoG et IMC   axe2 = tailleG   axe3 = physio et age
#interpretation des graphiques (tracer graphes)
draw_Eigen(acpT)
draw_ACP(acpT,1,2,T)
draw_ACP(acpT,1,3,T)
draw_ACP(acpT,3,2,T)

Location
plot(Europ[[1]])

plot(r, breaks=cuts, col = pal(7))
plot.raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif",ext = extent (5.6,16,44.2,47.3) ,line = F)

donensis = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

villosa = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

hirsuta = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

cottia = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

apennina = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

pedemontana = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

ecrins = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
polygon(x=f$x,y = f$y,col="red",density = 50,border ="red", lwd = 0.5)

text(x = Location$Longitude, y = Location$Latitude, Location$Code, cex = 0.7)


# TD choler ####

#fichier "Occurences_Alps.txt"

#A vérifier avec Julien mais la taxonomie c'est celle de Flora Europeae + Flora Mediterranea, P. cottia doit être sous le nom P. villosa et P. appenina sous le nom P. pedemontana

#récupération de données dans data
#install.packages("rgdal")
library(sp);library(raster)
library(rgdal)

WD = "~/Bureau/BEE/Stage/Pedemontana/data_carto"

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
plot.raster(RF, raw=F,line =T)
#mais rainforest en europe, donc peut etre besoin d'un seuil

RF80 = RF # même chose que précedement, mais en mode binaire
RF80[RF>=80] = 1
RF80[RF<80] = 0 #on peut aussi faire ça par raster::reclassify
plot.raster(RF80, raw = F, line = T)

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
plot.raster("data_carto/GRASSLANDS_10min.tif",ext = Alpine.ext,line = T)

#plot(raster("data_carto/PIA03395.tif"),extent (0,21600,0,9049))
plot.raster("data_carto/STRM.tif",ext =extent (10800,12000,7350,8000))
