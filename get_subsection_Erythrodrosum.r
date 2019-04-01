library(ape)
mat=read.dna(file='data/freebayes_-F0.3-n10-m30_-q20_mincov20_90samples_SNPs_only_FASTA.fasta',format='fasta')
mat
row.names(mat)


mat2=mat[which(row.names(mat)%in%c('AP1','AMB','AML','AOL','CS1','CP1','CP4','DGB','DRL','DMB','HC1','HGL','HS2','HP1','HPB','PT1','PV1','GA2','GA4','VR3','VR1','VL2','VB1')),] # with one outgroup: AP1 (P. lutea)
mat2

# subset to keep SNPs only
as.character(mat2[,1])

N_bases=rep(NA,dim(mat2)[2])
for (i in 1:length(N_bases)){
  temp=c(unique(as.character(mat2[,i])))
  N_bases[i]=sum(temp!='?')
  print(i)
}

table(N_bases)
as.character(mat2[,which(N_bases==7)])
mat3=mat2[,which(N_bases>1)]
mat3
write.dna(mat3,file='data/freebayes_-F0.3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_only_FASTA.fasta',format='fasta')


# sous-groupe 'pedemontana s.l.'
pedemontana=c('AMB','AML','AOL','CS1','CP1','CP4','PT1','PV1','GA2','GA4')

