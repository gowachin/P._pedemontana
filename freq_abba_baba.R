table = data.frame(AML = c("A","A","B","A"),
                   PT1 = c("A","B","A","B"),
                   VL2 = c("A","B","B","A"),
                   HS3 = c("B","A","A","B"))
table

Data=data.frame()
P1=c("apennina","cottia","pedemontana")
P2=c("apennina","cottia","ecrins")
P3=c("hirsuta")
Root=c("daonensis")
i=1
pb <- txtProgressBar(min = 1, max = length(P1)*length(P2)*length(P3)*length(Root), style = 3)
for (k in P1){
  for (j in P2) {
    for (l in P3) {
      for (m in Root) {
        #cat("\n",c(k,j,l,m),"\n")
        Data[i,1]=k
        Data[i,2]=j
        Data[i,3]=l
        Data[i,4]=m
        D = durand_freq_D(file,P1 =k,P2=j,P3 =l,Root =m,n=1e4)
        d = D[[1]]
        inf = D[[2]]
        Data[i,5]= mean(d)
        Data[i,6] = sum(d<0)/length(d)
        Data[i,7] = sum(d>0)/length(d)

        z <- abs(d[1]/sd(d[-1]))
        new.pval <- 2 * (1 - pnorm(z))
        Data[i,8] = new.pval
        Data[i,9] = z
        Data[i,10] = inf

        setTxtProgressBar(pb, i)
        i=i+1
      }
    }
  }
}
close(pb)

colnames(Data) = c("P1","P2","P3","Root",
                   "d","p<0","p>0","pval","z","site inf")
View(Data)
beep(3)

Data = rbind(Data[7,-c(6,7,9)],Data[8,-c(6,7,9)],Data[9,-c(6,7,9)],Data[3,-c(6,7,9)],Data[6,-c(6,7,9)])

write.table(Data,"Rendu/fig/ABBA.csv",sep=" & ")

length(summary(file$.pop))
P1 = "pedemontana" ; P2 = "hirsuta" ; P3 = "daonensis" ; Root = "lutea"

x =durand_freq_D(file,P1,P2,P3,Root,n=1e4)
hist(x)
z <- abs(x[1]/sd(x[-1]))
2 * (1 - pnorm(z))
sum(x<0)/length(x)
sum(x>0)/length(x)



z <- abs(d[1]/sd(d[-1]))
new.pval <- 2 * (1 - pnorm(z))
new.pval


durand_freq_D = function(file,P1,P2,P3,Root,n) {

  table = suppressMessages(substr(as.matrix(readr::read_delim(file$.csv,"\t", escape_double = FALSE, trim_ws = TRUE)[,-c(1:9)]),1,3))
  table = rbind(as.character(file$.pop),table)
  #verif
  if (P1 %in% file$.pop == F ) {message("P1 is not in file$.pop")}
  if (P2 %in% file$.pop == F ) {message("P2 is not in file$.pop")}
  if (P3 %in% file$.pop == F ) {message("P3 is not in file$.pop")}
  if (Root %in% file$.pop == F ) {message("Root is not in file$.pop")}

  #if ( length(unique.default(c(P1,P2,P3,Root)))<4 ) {
  #  forward = menu(c("yes","no"), title = "2 populations are the same, wish to continue?")
  #  if (forward == 1) {} else {stop("Assign other populations")}
  #  }

  echanti = c(P1,P2,P3,Root)
  manus = matrix(ncol = 4, nrow = dim(table)[1])
  for(i in 1:4) {
    matrix = table[,which(table[1,]==echanti[i])]
    H. = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) != ".")*2 )
    H1 = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) == "1")+sum(substr(as.character(matrix),3,3) == "1") )/H.
    manus[,i] = H1
  }

  x = apply(manus,1,sum) == 0
  manus = manus[which(x == F),]

  abba = (1-manus[,1])*manus[,2]*manus[,3]*(1-manus[,4])
  baba = manus[,1]*(1-manus[,2])*manus[,3]*(1-manus[,4])

  a_b = abba != 0 ; b_a = baba != 0
  inf = sum(a_b == T & b_a == T)
  #print(c("sites inform" = inf))
  D= c(sum(abba-baba)/sum(abba+baba))

  #pb <- txtProgressBar(min = 1, max = n, style = 3)
  for (i in 1: n+1){
    tmp.manus = manus[sample(c(1:dim(manus)[1]),dim(manus)[1], replace = T ),] ; tmp.manus

    abba = (1-tmp.manus[,1])*tmp.manus[,2]*tmp.manus[,3]*(1-tmp.manus[,4])
    baba = tmp.manus[,1]*(1-tmp.manus[,2])*tmp.manus[,3]*(1-tmp.manus[,4])

    D[i] = sum(abba-baba)/sum(abba+baba)

   # setTxtProgressBar(pb, i)
  }
  #close(pb)

  return(list(D,inf))
}




sample(c(1:10),10, replace = T)

