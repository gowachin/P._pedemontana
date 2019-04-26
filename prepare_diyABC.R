# preparation pour DIYABC ####

# Erythro ####
#.geno as a table
test.geno = read.geno("data_vcf/Eryth10_r5q20_t8i8_pos1e4.geno")
test.geno= rbind(rep("A",dim(test.geno)[2]),test.geno)
pop = c("POP","apennina", "apennina"#,"apennina"
        ,"pedemontana","pedemontana"
        ,"valgau","valgau"
        ,"cottia","cottia","cottia"
        ,"villosa","villosa","villosa","villosa"
        ,"hirsuta","hirsuta","hirsuta","hirsuta","hirsuta","hirsuta"
        ,"daonensis","daonensis"
)
ind=c("IND","AMB", "AOL",#"AML",
      "PT1", "PV1", "GA2", "GA4",
      "CS1", "CP1", "CP4",
      "VR3", "VR1", "VL2", "VB1",
      "DMB", "HC1", "HGL", "HS2", "HP1", "HPB",
      "DGB", "DRL")
test.geno = cbind(ind,c("SEX",rep("9",dim(test.geno)[1]-1)),pop,test.geno)
colnames(test.geno) = NULL
name = "tryhard_Eryth10" ; .snp = paste(name,".snp",sep="")
write.table(test.geno, .snp ,sep = "\t", quote = F, row.names=F, col.names = F)
# coller la première
system(paste("sed -i '1 i\ <NM=1NF>' ",.snp,sep=""))

