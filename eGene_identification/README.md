# reQTL-mapping in crop species

## install required R packages for running the model
`install.packages(c("dplyr","tidyr","tibble","lmerTest"))`

## file organization
Four types of input files are required for reQTL-mapping

- `File_A`
 The file includes top selected cis-eQTL and its relevant regulated gene. The format is as below:

 | gene | Chrom | Position | SNP |
 | :---: | :---: |:---: |:---: |
 |Zm00001d036337|Chr6| 84748285 | S6_84748285 |
 |Zm00001d036337|Chr6| 84748320 | S6_84748320 |

- `File_B`
 Genotype file includes variants in the studied panel. The format is as below:

 | id | GenoA | GenoB | GenoC | GenoD |
 | :---: | :---: |:---: |:---: |:--: |
 |S1_43699| 0 | 2 | 0 | 0 |
 |S1_43964| 0 | 2 | 2 | 2 |

- `File_C`
 Two covariate files derived from two conditions were employed for controlling confounding variables (File_C1 and File_C2). Principal components were estimated using R function [prcomp](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp) and PEER covariates were estimated using [PEER](https://github.com/PMBio/peer). The format is as below:

 | id | GenoA | GenoB | GenoC | GenoD |
 | :---: | :---: |:---: |:---: |:--: |
 |PC1| -123 | -551 | 157 | 102 |
 |PC2| 466 | -44 | -111 | -48 |

- `File_D`
 Two normalized expression files of genes employed in the study (File_D1 and File_D2). Expression data was normalized using R function [bestNormalize](https://github.com/petersonR/bestNormalize) The format is as below:

 | id | GenoA | GenoB | GenoC | GenoD |
 | :---: | :---: |:---: |:---: |:--: |
 |ENSRNA049458448| -0.398 | 0.133 | 0.477 | 2.068 |
 |Zm00001d000098| 2.150 | -0.834 | 1.747 | -1.51 |

## running the command for reQTL mapping
Using the command below for reQTL mapping

`Rscript cis-eQTL-detection.R File_A File_B File_C1 File_C2 File_D1 File_D2 > reQTLmapping-result.txt`
