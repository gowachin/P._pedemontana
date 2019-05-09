table = data.frame(AML = c("A","A","B","A"),
                   PT1 = c("A","B","A","B"),
                   VL2 = c("A","B","B","A"),
                   HS3 = c("B","A","A","B"))
table

ab_ba = function(table) {
  r = "no.inf"
  if ( table[1] == table[4] & table[2] == table[3] & table[1] != table[2] ) {r = "abba"
  } else if ( table[1] == table[3] & table[2] == table[4] & table[1] != table[2] ) {r = "baba"
  } #else {r = NA}
  return(r)
}

x = as.factor(unlist(apply(table, 1, function(table) ab_ba(table))))
abba = sum(x == "abba") ; abba
baba = sum(x == "baba") ; baba
x
d <- (abba - baba) / (abba + baba) ; d


file$.pop
P1 = "pedemontana" ; P2 = "valgau" ; P3 = "hirsuta" ; Root = "daonensis"

x =durand_D(file,P1,P2,P3,Root,1e3)
x[[2]]
hist(x[[1]],nclass = 24)
summary(x[[1]])

durand_D = function(file,P1,P2,P3,Root,n) {

  table = substr(as.matrix(readr::read_delim(file$.csv,"\t", escape_double = FALSE, trim_ws = TRUE)[,-c(1:9)]),1,3)
  table = rbind(as.character(file$.pop),table)
#verif
if (P1 %in% file$.pop == F ) {message("P1 is not in file$.pop")}
if (P2 %in% file$.pop == F ) {message("P2 is not in file$.pop")}
if (P3 %in% file$.pop == F ) {message("P3 is not in file$.pop")}
if (Root %in% file$.pop == F ) {message("Root is not in file$.pop")}

if ( length(unique.default(c(P1,P2,P3,Root)))<4 ) {
  forward = menu(c("yes","no"), title = "2 populations are the same, wish to continue?")
  if (forward == 1) {} else {stop("Assign other populations")}
  }

D =c()
pb <- txtProgressBar(min = 1, max = n, style = 3)
for (i in 1:n) {
setTxtProgressBar(pb, i)

c1 = sample(colnames(table[,which(table[1,]==P1)]),1) ;c1
c2 = sample(colnames(table[,which(table[1,]==P2)]),1) ;c2
c3 = sample(colnames(table[,which(table[1,]==P3)]),1) ;c3
root = sample(colnames(table[,which(table[1,]==Root)]),1) ; root

manus = cbind(table[,which(colnames(table) == c1 )],
              table[,which(colnames(table) == c2 )],
              table[,which(colnames(table) == c3 )],
              table[,which(colnames(table) == root )] )
colnames(manus) = c(c1,c2,c3,root)

hap1 = substr(as.matrix(manus),1,1)
hap2 = substr(as.matrix(manus),3,3)
homo = hap1 == hap2
homo
manus = manus[]
manus[which(homo==F)] = "."
manus = as.data.frame(manus)
manus$inf = as.factor(unlist(apply(manus, 1, function(manus) ab_ba(manus))))
summary(manus$inf)
abba = sum(manus$inf == "abba") ; abba
baba = sum(manus$inf == "baba") ; baba
d <- (abba - baba) / (abba + baba) ; d
D[i] =d
}
close(pb)
return(list(D,summary(manus$inf)))
}




sample(c(1:10),10, replace = T)

