#!/usr/bin/python2
#usage: python2 convert_vcf2introgress_loci.py input.vcf > input.introgress

# creates loci input file for introgress (R package) from vcf file

from sys import argv

read_file=open(argv[1],"r")

print "locus"+"\t"+"type"+"\t"+"chromosome"+"\t"+"location"
for line in read_file:
  if line.startswith("#"):
    pass
  else:
    chrom=line.split("\t")[0]
    pos=line.split("\t")[1]
    print chrom+"_"+pos+"\t"+"C"+"\t"+chrom+"\t"+chrom+"."+pos
    
read_file.close()


      
      
