# plot de l'évolution du pourcentage de données en fonction des seuils ####
par(mfrow = c(2,2))
#seuil de lignes
perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov10,n.r=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,2], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par lignes")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,2], res[,3], pos = 3)

perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov20,n.r=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,2], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par lignes")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,2], res[,3], pos = 3)

#seuil d'ind
perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov10,n.c=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,4], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par ind")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,4], res[,5], pos = 3)

perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov20,n.c=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,4], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par ind")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,4], res[,5], pos = 3)
par(mfrow = c(1,1))


# plot carto ####
plot.raster = function(rast,ext=extent (0,0,0,0),line = F,raw =T) {
  if(raw == T) {rast=raster(rast)}
  if (ext != extent (0,0,0,0)) {plot(crop(rast, ext), col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
  } else {plot(rast, col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)) }
  if(line==T){lines(read.table("data_carto/WORLD_lowres.dat"))}
}

graphics.off()
Location <- read_delim("Location.csv", "&", escape_double = FALSE, trim_ws = TRUE)
colors = c("darkolivegreen3","lightsalmon3","goldenrod2","steelblue3","slateblue2","darkorchid3","violetred2")
par(mar = c(4.85, 2.825, 4, 4.56))
barplot(100,width=1,col= "lightblue", axes = F)
par(new=TRUE, mar = c(5,4,4,2))
layout(mat = matrix(c(1), nrow = 1, ncol = 1),
       heights = c(1),
       widths = c(1))
plot(crop(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), extent (3.95,16,43,47.5), scalebar=F),
     col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), legend=F,
     xlab = "Longitude", ylab = "Latitude")
#Occurences_Alps <- read_csv("data_carto/Occurences_Alps.csv")[,-7]
#Repart = list(Repart.pedemontana = Occurences_Alps[Occurences_Alps$Taxon == "Primula pedemontana",],
#              Repart.villosa = Occurences_Alps[Occurences_Alps$Taxon == "Primula villosa",],
#              Repart.hirsuta = Occurences_Alps[Occurences_Alps$Taxon == "Primula hirsuta",],
#              Repart.daonensis = Occurences_Alps[Occurences_Alps$Taxon == "Primula daonensis",])
#tail = c(0.1,0.1,0.1,0.1) ; color = c("red","orange","green","purple")
#for (i in 1:length(Repart)) {points(x = Repart[[i]]$lon , y = Repart[[i]]$lat, cex = tail[i], col=color[i]) }

#hirsuta = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=30, id=FALSE, xy=T, type="n")
hirsuta = data.frame(x = c(9.629167, 9.070833, 8.395833, 8.179167, 7.762500, 7.562500, 7.287500, 7.070833, 6.687500, 6.412500, 6.270833, 6.037500, 5.837500,5.820833, 6.079167, 6.437500, 6.612500,
                           6.712500, 7.029167, 7.045833, 6.970833, 7.104167, 7.462500, 7.837500, 8.120833, 8.379167, 8.512500, 8.870833, 9.345833,10.020833,10.29583, 10.51250, 10.51250, 10.00417),
                     y= c(47.08750,47.08750,46.99583,46.97917,46.77083,46.32083,46.23750,46.21250,45.84583,45.65417,45.51250,45.37083,45.16250,44.91250,44.56250,44.66250,44.95417,45.15417,45.37083,
                          45.55417,45.69583,45.83750,45.86250,45.94583,46.12917,46.15417,46.05417,46.01250,45.98750,46.11250,46.25417, 46.58750, 46.95417, 47.07917))
polygon(x=hirsuta$x,y = hirsuta$y,col=colors[1],density = 75,border =colors[1], lwd = 0.5)
#hirsutap = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=15, id=FALSE, xy=T, type="n")
hirsutap = data.frame(x = c(0.60416667, -0.12916667, -0.80416667, -0.52916667, -0.10416667,  0.98750000,  1.47916667,  1.94583333,  2.04583333,  1.35416667,  1.03750000,0.77083333),
                      y= c(42.89583, 42.98750, 42.87917, 42.68750, 42.60417, 42.44583, 42.41250, 42.40417, 42.58750, 42.83750, 42.86250, 42.87917))
polygon(x=hirsutap$x,y = hirsutap$y,col=colors[1],density = 75,border =colors[1], lwd = 0.5)
#daonensis = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=12, id=FALSE, xy=T, type="n")
daonensis = data.frame(x = c(10.58750, 10.43750, 10.33750, 10.36250, 10.39583, 10.36250, 10.26250, 10.39583, 10.67917, 10.87083, 10.97917, 10.81250),
                       y= c(46.61250, 46.61250, 46.53750, 46.41250, 46.29583, 46.07083, 45.81250, 45.87083, 46.07083, 46.25417, 46.42917, 46.58750))
polygon(x=daonensis$x,y = daonensis$y,col=colors[2] ,density = 75,border =colors[2], lwd = 0.5)
#villosa = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=12, id=FALSE, xy=T, type="n")
villosa = data.frame(x = c(14.42083, 14.13750, 13.89583, 13.57083, 13.55417, 13.87917, 14.13750, 14.37083, 14.62917, 15.02917, 15.20417, 15.07083),
                     y= c(47.07917, 47.07917, 47.16250, 47.00417, 46.82917, 46.67917, 46.81250, 46.98750, 46.87917, 46.89583, 47.07917, 47.22917))
polygon(x=villosa$x,y = villosa$y,col=colors[3],density = 75,border =colors[3], lwd = 0.5)
#cottia = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
cottia = data.frame(x = c(6.979167, 6.845833, 6.720833, 6.720833, 6.870833, 7.087500, 7.095833, 7.079167),
                    y= c(45.05417, 44.94583, 44.86250, 44.71250, 44.57917, 44.56250, 44.73750, 44.88750))
polygon(x=cottia$x,y = cottia$y,col=colors[4],density = 75,border =colors[4], lwd = 0.5)
#apennina = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=10, id=FALSE, xy=T, type="n")
apennina = data.frame(x = c(9.820833,  9.629167,  9.470833,  9.537500,  9.787500, 10.070833, 10.362500, 10.587500, 10.404167, 10.070833),
                      y= c(44.58750, 44.60417, 44.57083, 44.42917, 44.32083, 44.26250, 44.14583, 44.17083, 44.30417, 44.47083))
polygon(x=apennina$x,y = apennina$y,col=colors[5],density = 75,border =colors[5], lwd = 0.5)
#pedemontana = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=10, id=FALSE, xy=T, type="n")
pedemontana = data.frame(x = c(7.470833, 7.054167, 6.704167, 6.737500, 6.870833, 7.054167, 7.279167, 7.495833, 7.629167, 7.629167),
                         y= c(45.67917, 45.67917, 45.49583, 45.20417, 45.11250, 45.22083, 45.37083, 45.45417, 45.56250, 45.65417))
polygon(x=pedemontana$x,y = pedemontana$y,col=colors[6],density = 75,border =colors[6], lwd = 0.5)
#ecrins = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
ecrins = data.frame(x = c(6.437500, 6.312500, 6.129167, 6.062500, 6.104167, 6.237500, 6.404167, 6.470833),
                    y= c(44.87917, 44.94583, 44.97083, 44.89583, 44.79583, 44.72917, 44.74583, 44.80417))
polygon(x=ecrins$x,y = ecrins$y,col=colors[7],density = 75,border =colors[7], lwd = 0.5)

text(x = Location$Longitude, y = Location$Latitude, Location$Code, cex = 0.65)

points(x = c(6.1667,9.1859243,5.724524), y = c(46.2,45.4654219,45.188529), pch = 4, lwd= 2)
text(x = c(6.1667,9.1859243,5.724524), y = c(46.2,45.4654219,45.188529), c("Geneva","Milano","Grenoble"), font=c(2), pos = c(3,3,2), cex =0.8)

text(x = c(12.5), y = c(44), c("N"), font=c(2), pos = 3, cex =1.5)
arrows(12.5, 43.1, 12.5, y1 = 44, length = 0.2, angle = 30, code = 2, lwd = 5)

legend("bottomright", c("P. hirsuta (6)","P. daonensis (2)","P. villosa (4)","P. cottia (4)","P. apennina (2)","P. pedemontana (2)","P. pedemontana Ecrins (2)"),
       fill = colors, text.font=c(3), cex = 0.7)

legend("bottomleft", c("                                  "," "," "," "," "), cex = 0.86)
par(new=TRUE)
layout(mat = matrix(c( 0, 0,0,0,1,0, 0, 0,0), nrow = 3, ncol = 3),
       heights = c(5, 2,1.7),
       widths = c(0.75,1.85, 5))
barplot(100,width=1.2,col= "lightblue", axes = F)
par(new=TRUE)
layout(mat = matrix(c( 0, 0,0,0,1,0, 0, 0,0), nrow = 3, ncol = 3),
       heights = c(5, 2,1.7),
       widths = c(0.8,1.7, 5))
plot(crop(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), extent (-2.1,4,41.5,44),
          scalebar=F), col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), box=F, legend=F)
#hirsutap = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=15, id=FALSE, xy=T, type="n")
#hirsutap = data.frame(x = c(0.60416667, -0.12916667, -0.80416667, -0.52916667, -0.10416667,  0.98750000,  1.47916667,  1.94583333,  2.04583333,  1.35416667,  1.03750000,0.77083333),
#                      y= c(42.89583, 42.98750, 42.87917, 42.68750, 42.60417, 42.44583, 42.41250, 42.40417, 42.58750, 42.83750, 42.86250, 42.87917))
#polygon(x=hirsutap$x,y = hirsutap$y,col="red",density = 50,border ="red", lwd = 0.5)
#polygon(x=hirsutap$x,y = hirsutap$y,col="red",density = 50,border ="red", lwd = 0.5)
par(new=TRUE)
layout(mat = matrix(c(1), nrow = 1, ncol = 1),
       heights = c(1),
       widths = c(1))
par(mar = c(4.85, 2.825, 4, 4.56))
plot(100,100,type="n", axes = F, ylab="", xlab="")

text(71,80,"Pyrénées", font= 2, cex= 0.8)

hirsutap = data.frame(x = c(71.72795, 69.90683, 68.34588, 67.82556, 68.60604, 70.42715, 73.15883, 73.93931, 73.02875, 71.85803, 70.94747, 68.99628),
                      y= c(74.09680, 75.14700, 75.14700, 73.30915, 71.99640, 70.94620, 69.10835, 71.99640, 73.83425, 73.83425, 74.09680, 75.67210))
polygon(x=hirsutap$x,y = hirsutap$y,col=colors[1],density = 75,border =colors[1], lwd = 0.5)
text(70,75,"HP1", cex= 0.65)


# phrapl tree ####

library(phrapl)
library(rgl)

data(TestData)
migrationArray[[4]]

#rajout des growth map car le jeu de base ne les contient pas!!!
for (i in 1:6) {
  migrationArray[[i]]$growthMap = matrix(c(0,0,0,0,NA,0),ncol =2)
}

PlotModel2D( migrationArray[[2]])

ab_plot = migrationArray[[2]]

ab_plot$collapseMatrix = matrix(c(1,1,0,0,
                                  1,NA,1,0,
                                  1,NA,NA,1), ncol = 3)
ab_plot$n0multiplierMap = matrix(c(1,1,1,1,
                                   1,NA,1,1,
                                   1,NA,NA,1), ncol = 3)
ab_plot$growthMap = matrix(c(0,0,0,0,
                             0,NA,0,0,
                             0,NA,NA,0), ncol = 3)
d_plot = ba_plot = ab_plot
d_plot$migrationArray = array(c(NA,0,0,0,
                                0,NA,0,0,
                                0,0,NA,0,
                                0,0,0,NA,
                                NA,NA,0,0,
                                NA,NA,NA,NA,
                                0,NA,NA,0,
                                0,NA,0,NA,
                                NA,NA,NA,0,
                                NA,NA,NA,NA,
                                NA,NA,NA,NA,
                                0,NA,NA,NA),c(4,4,3))
ab_plot$migrationArray = array(c(NA,0,0,0,
                                 0,NA,1,0,
                                 0,1,NA,0,
                                 0,0,0,NA,
                                 NA,NA,0,0,
                                 NA,NA,NA,NA,
                                 0,NA,NA,0,
                                 0,NA,0,NA,
                                 NA,NA,NA,0,
                                 NA,NA,NA,NA,
                                 NA,NA,NA,NA,
                                 0,NA,NA,NA),c(4,4,3))

ba_plot$migrationArray = array(c(NA,0,1,0,
                                 0,NA,0,0,
                                 1,0,NA,0,
                                 0,0,0,NA,
                                 NA,NA,0,0,
                                 NA,NA,NA,NA,
                                 0,NA,NA,0,
                                 0,NA,0,NA,
                                 NA,NA,NA,0,
                                 NA,NA,NA,NA,
                                 NA,NA,NA,NA,
                                 0,NA,NA,NA),c(4,4,3))

par(mfrow=c(1,3))
PlotModel2D(d_plot, taxonNames = c("A","B","B","A"), tree.col = c("grey"),tree.lwd = 10)
PlotModel2D(d_plot, taxonNames = c("B","A","B","A"), tree.col = c("grey"),tree.lwd = 10)
PlotModel2D(ab_plot, taxonNames = c("A","B","B","A"), tree.col = c("grey"),tree.lwd = 10)
PlotModel2D(ba_plot, taxonNames = c("B","A","B","A"), tree.col = c("grey"),tree.lwd = 10)



