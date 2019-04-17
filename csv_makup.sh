#!/bin/bash
# tryhard :

name=$1
head=$2
csv=$3

echo $name $head $csv

grep  '##' $name > $head
sed '/^##/ d' $name > $csv
sed -i '1 s/^.\{1\}//' $csv


#sed -i '1s/.*/#&/' $name.csv # rajoute une # devant la premiere ligne
#cat vcf_head.txt $name.csv > tryhard.vcf # recolle l'entete supprimé avant





#mv tryhard.csv tryhard.vcf ; ls

