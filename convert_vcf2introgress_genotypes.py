#!/usr/bin/python2
#usage: python2 convert_vcf2introgress.py input.vcf > input.introgress

from sys import argv

read_file=open(argv[1],"r")

for line in read_file:
  if line.startswith("#"):
    pass
  else:
    inds="\t".join(i.split(":")[0] for i in line.split("\t")[9:])
    inds=inds.replace(".","NA")
    print inds
read_file.close()


      
      
