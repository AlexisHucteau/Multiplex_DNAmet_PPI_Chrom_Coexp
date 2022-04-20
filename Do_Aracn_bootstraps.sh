#! /bin/bash

read -p 'Phenotype to analyse: ' Pheno
read -p 'pvalue threshold (ex 1E-8): ' pvalue

java -Xmx5G -jar ~/ARACNe-AP/dist/aracne.jar -e DATA/RNAseq4Aracn_$Pheno.tsv \
-o ./output_$Pheno/ \
--tfs DATA/TF.txt \
--pvalue $pvalue --seed 1 --calculateThreshold

for i in {1..100}
do
  date
  echo $i
  java -Xmx5G -jar ~/ARACNe-AP/dist/aracne.jar \
  -e ./DATA/RNAseq4Aracn_$Pheno.tsv \
  -o ./output_$Pheno \
  --tfs DATA/TF.txt \
  --pvalue $pvalue --seed $i
done

java -Xmx5G -jar ~/ARACNe-AP/dist/aracne.jar -o output_$Pheno/ --consolidate
