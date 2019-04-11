#!/bin/bash
# tryhard :

name=$1

sed -i '1 i\ ##fileformat=VCFv4.1'   $name
sed -i '2 i\ ##fileDate=20151126'   $name
sed -i '3 i\ ##source=freeBayes v0.9.21-26-gbfd9832'   $name
sed -i '4 i\ ##reference=./Pveris_genome/Pveris_pbjelly_updated_assembly_edited_all_ambig_turned_to_N.fasta'   $name
sed -i '5 i\ ##phasing=none'   $name
sed -i '6 i\ ##commandline="./freebayes/bin/freebayes --stdin -f ./Pveris_genome/Pveris_pbjelly_updated_assembly_edited_all_ambig_turned_to_N.fasta --min-coverage 20 -n 10 -m 30 -q 20 -F 0.3"'   $name
sed -i '7 i\ ##filter="TYPE = snp"'   $name
sed -i '8 i\ ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">'   $name
sed -i '9 i\ ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">'   $name
sed -i '10 i\ ##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">'   $name
sed -i '11 i\ ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">'   $name
sed -i '12 i\ ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'   $name
sed -i '13 i\ ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">'   $name
sed -i '14 i\ ##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count, with partial observations recorded fractionally">'   $name
sed -i '15 i\ ##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">'   $name
sed -i '16 i\ ##INFO=<ID=PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">'   $name
sed -i '17 i\ ##INFO=<ID=PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">'   $name
sed -i '18 i\ ##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">'   $name
sed -i '19 i\ ##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">'   $name
sed -i '20 i\ ##INFO=<ID=PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">'   $name
sed -i '21 i\ ##INFO=<ID=PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">'   $name
sed -i '22 i\ ##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">'   $name
sed -i '23 i\ ##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">'   $name
sed -i '24 i\ ##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">'   $name
sed -i '25 i\ ##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">'   $name
sed -i '26 i\ ##INFO=<ID=SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '27 i\ ##INFO=<ID=SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '28 i\ ##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">'   $name
sed -i '29 i\ ##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '30 i\ ##INFO=<ID=RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">'   $name
sed -i '31 i\ ##INFO=<ID=RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '32 i\ ##INFO=<ID=RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '33 i\ ##INFO=<ID=RPL,Number=A,Type=Float,Description="Reads Placed Left: number of reads supporting the alternate balanced to the left (5) of the alternate allele">'   $name
sed -i '34 i\ ##INFO=<ID=RPR,Number=A,Type=Float,Description="Reads Placed Right: number of reads supporting the alternate balanced to the right (3) of the alternate allele">'   $name
sed -i '35 i\ ##INFO=<ID=EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '36 i\ ##INFO=<ID=EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffdings inequality">'   $name
sed -i '37 i\ ##INFO=<ID=DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">'   $name
sed -i '38 i\ ##INFO=<ID=ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">'   $name
sed -i '39 i\ ##INFO=<ID=GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">'   $name
sed -i '40 i\ ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">'   $name
sed -i '41 i\ ##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that = is replaced by M to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">'   $name
sed -i '42 i\ ##INFO=<ID=NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">'   $name
sed -i '43 i\ ##INFO=<ID=MEANALT,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">'   $name
sed -i '44 i\ ##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">'   $name
sed -i '45 i\ ##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">'   $name
sed -i '46 i\ ##INFO=<ID=MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">'   $name
sed -i '47 i\ ##INFO=<ID=PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">'   $name
sed -i '48 i\ ##INFO=<ID=PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">'   $name
sed -i '49 i\ ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'   $name
sed -i '50 i\ ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">'   $name
sed -i '51 i\ ##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">'   $name
sed -i '52 i\ ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'   $name
sed -i '53 i\ ##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">'   $name
sed -i '54 i\ ##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">'   $name
sed -i '55 i\ ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">'   $name
sed -i '56 i\ ##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">'   $name
sed -i '57s/.*/#&/' $name

sed -i '1 s/^.\{1\}//'   $name
sed -i '2 s/^.\{1\}//'   $name
sed -i '3 s/^.\{1\}//'   $name
sed -i '4 s/^.\{1\}//'   $name
sed -i '5 s/^.\{1\}//'   $name
sed -i '6 s/^.\{1\}//'   $name
sed -i '7 s/^.\{1\}//'   $name
sed -i '8 s/^.\{1\}//'   $name
sed -i '9 s/^.\{1\}//'   $name
sed -i '10 s/^.\{1\}//'   $name
sed -i '11 s/^.\{1\}//'   $name
sed -i '12 s/^.\{1\}//'   $name
sed -i '13 s/^.\{1\}//'   $name
sed -i '14 s/^.\{1\}//'   $name
sed -i '15 s/^.\{1\}//'   $name
sed -i '16 s/^.\{1\}//'   $name
sed -i '17 s/^.\{1\}//'   $name
sed -i '18 s/^.\{1\}//'   $name
sed -i '19 s/^.\{1\}//'   $name
sed -i '20 s/^.\{1\}//'   $name
sed -i '21 s/^.\{1\}//'   $name
sed -i '22 s/^.\{1\}//'   $name
sed -i '23 s/^.\{1\}//'   $name
sed -i '24 s/^.\{1\}//'   $name
sed -i '25 s/^.\{1\}//'   $name
sed -i '26 s/^.\{1\}//'   $name
sed -i '27 s/^.\{1\}//'   $name
sed -i '28 s/^.\{1\}//'   $name
sed -i '29 s/^.\{1\}//'   $name
sed -i '30 s/^.\{1\}//'   $name
sed -i '31 s/^.\{1\}//'   $name
sed -i '32 s/^.\{1\}//'   $name
sed -i '33 s/^.\{1\}//'   $name
sed -i '34 s/^.\{1\}//'   $name
sed -i '35 s/^.\{1\}//'   $name
sed -i '36 s/^.\{1\}//'   $name
sed -i '37 s/^.\{1\}//'   $name
sed -i '38 s/^.\{1\}//'   $name
sed -i '39 s/^.\{1\}//'   $name
sed -i '40 s/^.\{1\}//'   $name
sed -i '41 s/^.\{1\}//'   $name
sed -i '42 s/^.\{1\}//'   $name
sed -i '43 s/^.\{1\}//'   $name
sed -i '44 s/^.\{1\}//'   $name
sed -i '45 s/^.\{1\}//'   $name
sed -i '46 s/^.\{1\}//'   $name
sed -i '47 s/^.\{1\}//'   $name
sed -i '48 s/^.\{1\}//'   $name
sed -i '49 s/^.\{1\}//'   $name
sed -i '50 s/^.\{1\}//'   $name
sed -i '51 s/^.\{1\}//'   $name
sed -i '52 s/^.\{1\}//'   $name
sed -i '53 s/^.\{1\}//'   $name
sed -i '54 s/^.\{1\}//'   $name
sed -i '55 s/^.\{1\}//'   $name
sed -i '56 s/^.\{1\}//'   $name

#mv tryhard.csv tryhard.vcf ; ls

