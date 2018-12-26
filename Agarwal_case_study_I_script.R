#Loading the required libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)

#Reading the input files and storing it in data table
file_path = "~/Desktop/TUM Course/Data Analysis and Visualization in R/Case Study/Data_sets"
gene <- fread(file.path(file_path,"gene.txt"))
genotype <- fread(file.path(file_path,"genotype.txt"))
growth <- fread(file.path(file_path,"growth.txt"))
marker <- fread(file.path(file_path,"marker.txt"))
expression <- fread(file.path(file_path,"expression.txt"))

glimpse(gene)
dim(genotype)
glimpse(genotype[,1:5])
glimpse(growth)
glimpse(marker)
dim(expression)
glimpse(expression[,1:5])

#Data Analysis I
growth_gather <- gather(growth,key="env",value="growth_rate",-1)

ggplot(growth_gather,aes(x=env,y=growth_rate)) + geom_boxplot() + 
  geom_jitter(width = 0.4,alpha=0.2)+labs(x="Environment",y="Growth Rate")

ggplot(growth_gather,aes(growth_rate)) + geom_density(aes(fill=env,alpha=0.25))+ labs(x="Growth Rate",y="Density")

#Data Analysis - II
genotype_mrk <- gather(genotype,"marker","value",-1)
genotype_mrk <- as.data.table(genotype_mrk)
genotype_strain_count <- genotype_mrk[,.(count=.N),by=.(value,strain)]
genotype_strain_count <- genotype_strain_count[,c(2,3,1)]
genotype_strain_count <- spread(genotype_strain_count,key=value,value=count)

#Differentiate a segment based on dominant strain
genotype_strain_count$genotype <- ifelse(genotype_strain_count$`Lab strain`>genotype_strain_count$`Wild isolate`,"Lab strain","Wild isolate")
growth_merged <- merge(growth_gather,genotype_strain_count)

ggplot(growth_merged,aes(x=env,y=growth_rate,color=genotype)) +labs(x="Environment",y="Growth Rate")+ geom_boxplot()

#Data Analysis - III
marker <- read.delim(file.path(file_path,"marker.txt"))
growth <- read.delim(file.path(file_path,"growth.txt"))
genotype <- read.delim(file.path(file_path,"genotype.txt"))
mygeno <- genotype[, which(marker$chrom=="chr07" & marker$start== 1069229)]
names(mygeno) <- genotype$strain
plot(YPMalt ~ mygeno[strain], data=growth)

#Finding the index of the desired marker
which(marker$chrom=="chr07" & marker$start== 1069229)
#Find the marker id for the obtained index
marker[which(marker$chrom=="chr07" & marker$start== 1069229),]
#Column name of the marker id used
colnames(genotype)[which(marker$chrom=="chr07" & marker$start== 1069229)]
#Desired column name of the marker id
colnames(genotype)[which(marker$chrom=="chr07" & marker$start== 1069229)+1]

#Desired Boxplot graph
mygeno <- genotype[, which(marker$chrom=="chr07" & marker$start== 1069229)+1]
names(mygeno) <- genotype$strain
plot(YPMalt ~ mygeno[strain], data=growth)