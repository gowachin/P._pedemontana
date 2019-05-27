manual()

library(pbapply) # uberclass but maybe not usefull

# phrapl ####

# installer gmp https://gmplib.org/ avec dll, tar puis in repertory
#	./configure
#make
#make check		<= VERY IMPORTANT!!
#	make install

#sudo apt-get install libfreetype6-dev
devtools::install_github('bomeara/phrapl')


if(!require("ape")){install.packages("ape",repos="http://cran.rstudio.com")}
if(!require("partitions")){install.packages("partitions",repos="http://cran.rstudio.com")}
if(!require("lattice")){install.packages("lattice",repos="http://cran.rstudio.com")}
if(!require("polynom")){install.packages("polynom",repos="http://cran.rstudio.com")}
if(!require("gmp")){install.packages("gmp",repos="http://cran.rstudio.com")}
if(!require("rgenoud")){install.packages("rgenoud",repos="http://cran.rstudio.com")}
if(!require("parallel")){install.packages("parallel",repos="http://cran.rstudio.com")}
if(!require("optimx")){install.packages("optimx",repos="http://cran.rstudio.com")}
if(!require("igraph")){install.packages("igraph",repos="http://cran.rstudio.com")}
if(!require("numDeriv")){install.packages("numDeriv",repos="http://cran.rstudio.com")}
if(!require("nloptr")){install.packages("nloptr",repos="http://cran.rstudio.com")}
if(!require("Matrix")){install.packages("Matrix",repos="http://cran.rstudio.com")}
if(!require("rgl")){install.packages("rgl",repos="http://cran.rstudio.com")}
if(!require("RColorBrewer")){install.packages("RColorBrewer",repos="http://cran.rstudio.com")}
if(!require("igraph")){install.packages("igraph",repos="http://cran.rstudio.com")}
if(!require("diagram")){install.packages("diagram",repos="http://cran.rstudio.com")}
if(!require("binom")){install.packages("binom",repos="http://cran.rstudio.com")}
if(!require("P2C2M")){install.packages("P2C2M",repos="http://cran.rstudio.com")}

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

