# reQTL-mapping input files organization

Four types of input files are required for reQTL-mapping

- `File_A`
 The file includes top selected cis-eQTL and its relevant regulated gene. The format is as below:

 | gene | Chrom | Position | SNP |
 | :---: | :---: |:---: |:---: |
 |Zm00001d036337|Chr6| 84748285 | S6_84748285 |
 |Zm00001d036337|Chr6| 84748320 | S6_84748320 |

- `File_B`
 Genotype file includes variants in the studied panel. The format is as below:

 | id | GenoA | GenoB | GenoC | GenoC |
 | :---: | :---: |:---: |:---: |:--: |
 |S1_43699| 0 | 2 | 0 | 0 |
 |S1_43964| 0 | 2 | 2 | 2 |

- `File_C`
 Two covariate files derived from two conditions were employed for controlling confounding variables (File_C1 and File_C2). The format is as below:

 | id | GenoA | GenoB | GenoC | GenoC |
 | :---: | :---: |:---: |:---: |:--: |
 |PC1| -123 | -551 | 157 | 102 |
 |PC2| 466 | -44 | -111 | -48 |

- `File_D`
 Two normalized expression files of genes employed in the study (File_D1 and File_D2). The format is as below:

 | id | GenoA | GenoB | GenoC | GenoC |
 | :---: | :---: |:---: |:---: |:--: |
 |ENSRNA049458448| -0.398 | 0.133 | 0.477 | 2.068 |
 |Zm00001d000098| 2.150 | -0.834 | 1.747 | -1.51 |
 
# reQTL mapping
`Rscript cis-eQTL-detection.R File_A File_B File_C1 File_C2 File_D1 File_D2 > reQTLmapping-result.txt`
