library(dplyr)
library(tidyr)
library(tibble)
library(lmerTest)

args <- commandArgs(trailingOnly = TRUE)
# args[1]: gene-SNP pair
# args[2]: genotype file
# args[3]: 5PCs + 25 PEER covariate file in control
# args[4]: 5PCs + 25 PEER covariate file in treatment
# args[5]: normalized expression file in control
# args[6]: normalized expression file in treatment
# args[7]: output file
pair = read.table(args[1],sep='\t',head=T)
head(pair)
pgeno = read.table(args[2],sep='\t',head=T, check.names = F)
rownames(pgeno) = pgeno$id
pgeno = pgeno[,-1]
geno = data.frame(t(pgeno))

pc_cov = read.table(args[3],head=T,sep='\t',check.names = F)
ph_cov = read.table(args[4],head=T,sep='\t',check.names = F)
rownames(pc_cov) = pc_cov$id
rownames(ph_cov) = ph_cov$id
pc_cov = pc_cov[,-1]
ph_cov = ph_cov[,-1]
ccov = data.frame(t(pc_cov)) %>% rownames_to_column(var="id")
hcov = data.frame(t(ph_cov)) %>% rownames_to_column(var="id")

cpheno = read.table(args[5],sep='\t',head=T,check.names = F)
colnames(cpheno) = c("id",colnames(pc_cov))
rownames(cpheno) = cpheno$id
cpheno = cpheno[,-1]
hpheno = read.table(args[6],sep='\t',head=T, check.names = F)
colnames(hpheno) = c("id",colnames(ph_cov))
rownames(hpheno) = hpheno$id
hpheno = hpheno[,-1]
cp = data.frame(t(cpheno))
hp = data.frame(t(hpheno))

results = data.frame()
total = nrow(pair)
#start = args[1]
#stop = args[2]
start = 1
stop = total
for (i in c(start:stop)){
    data = pair[i,]
    snp = c(as.character(data[,4]))
    gene = c(as.character(data[,1]))
    name = paste(gene, snp, sep=":")
    print (paste(gene, snp, sep=';'))
    if (gene %in% colnames(cp) & gene %in% colnames(hp)){
        csubpheno = cp %>% rownames_to_column(var="id") %>% select(id, gene)
        hsubpheno = hp %>% rownames_to_column(var="id") %>% select(id, gene)
        cdata = geno %>% select(snp) %>% rownames_to_column(var="id") %>%
            left_join(csubpheno, by="id") %>% left_join(ccov, by="id") %>% mutate(Condition="control")
        hdata = geno %>% select(snp) %>% rownames_to_column(var="id") %>%
            left_join(hsubpheno, by="id") %>% left_join(hcov, by="id") %>% mutate(Condition="heat")
        data1 = rbind(cdata, hdata)
        variables = c(snp, "Condition", "PC1", "PC2", "PC3", "PC4", "PC5", "V1", "V2", "V3",
                      "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15",
                      "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25",
                      "PC1:Condition","PC2:Condition","PC3:Condition","PC4:Condition","PC5:Condition","V1:Condition", "V2:Condition", "V3:Condition",
                      "V4:Condition", "V5:Condition", "V6:Condition", "V7:Condition", "V8:Condition",
                      "V9:Condition", "V10:Condition", "V11:Condition", "V12:Condition", "V13:Condition",
                      "V14:Condition", "V15:Condition", "V16:Condition", "V17:Condition", "V18:Condition",
                      "V19:Condition", "V20:Condition", "V21:Condition", "V22:Condition", "V23:Condition",
                      "V24:Condition", "V25:Condition", paste(snp, "Condition", sep=":"), "(1|id)")
        fm = as.formula(paste(gene, paste(variables, collapse = "+"),sep="~"))
        data2 = data1
        data2[,c(4:8)] = scale(data1[,c(4:8)])
        OLSexamp = lmerTest::lmer(fm, data = data2)
        result = data.frame(anova(OLSexamp)[63,])
        result[,7] = name
        colnames(result)[7] = "name"
        results = rbind(results, result)
    }
}
write.table(results, file=args[7], sep='\t')
