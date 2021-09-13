library(tidyverse)
library(vcfR)
library(ape)
library(beeswarm)
library(vroom)
library(RColorBrewer)

# to do
	# visualize mapping rate data
	# visualize PCA
	# visualize admixture/cross validation

	# read in populations.sumstats.tsv.gz, look at HWE, etc

	# dig into VCF
		# depth
		# heterozygosity
		# hdplot

##############################
# visualize mapping rate data
##############################

# these steps aggregate
	# 1) metadata
	# 2) mapping data
	# 3) barcode/pool data
	
# read in metadata table
meta <- read.table("../meta/Urban_ddRAD_FishIDs_Bioinformatics_2018.tsv",header=TRUE,sep="\t",quote="",comment.char="")

# read in sam stats table
sn <- read.table("../results/align_stats/SN.txt",sep="\t") %>% t()
	sn <- cbind(Sample=rownames(sn),data.frame(sn))
	colnames(sn) <- colnames(sn) %>% str_replace(regex("\\.*$"),"")

	sn <- sn[,c(1,2,8,9,14,20,24,39)]

# add mapping data to metadata table
meta <- left_join(x=meta,y=sn,by=c("bioinformatics.id"="Sample")) 

# edit column names to replace periods with underscores
colnames(meta) <-  gsub("\\.", "_",colnames(meta))

# read in barcode and pool information
	# loop over files, bind them into one table
# get file list
f <- list.files("../meta",pattern="barcode*",full.names=TRUE)

# initialize dataframe
bc <- data.frame()
# for loop
for(i in f){
	seqpool <- str_extract(i,regex("[0-9][0-9][0-9]"))
	print(seqpool)
	bc <- rbind(bc,cbind(pool=seqpool,read.table(i)))
	}

	colnames(bc) <- c("pool","barcode","Sample")

# merge barcode and pool info with metadata and mapping data
meta <- left_join(x=meta,y=bc,by=c("bioinformatics_id"="Sample")) 
	rownames(meta) <- meta$bioinformatics_id

# have a look
	# a bunch of individuals in pools 1 and 2 have major problems. these will get dropped when filtering the VCF
plot(meta$reads_mapped/meta$raw_total_sequences,col=factor(meta$pool))


##############################
# plot PCA
##############################

# get the percent explained by each component
	# first get eigenvalues, scale them by sum to get pve
eval <- read.table("../results/plink_pca/fb.eigenval")
pve <- eval/sum(eval)*100

# get the PCA eignvectors
evec <- read.table("../results/plink_pca/fb.eigenvec")

# add population and "region" to the table
evec <- cbind(evec, population=meta[evec[,2],"popmap_initial"],region=meta[evec[,2],"region"])

# plot PC1 and PC2
ggplot(evec,aes(x=V3,y=V4,color=population,shape=region)) +
	geom_point() + 
	xlab(paste0("PC1 (", signif(pve[1,1], 3), "%)")) + 
	ylab(paste0("PC2 (", signif(pve[2,1], 3), "%)"))


#######################
# visualize admixture
#######################
 # adapted from https://github.com/mishaploid/Bo-demography/blob/master/src/plot_admixture_results.R

# get cross validation errors
logs <- list()
for(i in 1:10){
	f <- paste0("../results/admixture/log",i,".out")
	logs[[i]] <- scan(f,what="character",sep="\n")
	logs[[i]] <- logs[[i]][grep("CV error",logs[[i]])] %>% str_extract_all(regex("[0-9.]+")) %>% unlist() %>% as.numeric()
}

cv <- do.call(rbind,logs)
	colnames(cv) <- c("K","cv_error")

# k with min CV error is a good choice to plot, we'll plot them all though
plot(cv)

# list admixture output files 
# assumes files are located in current working directory
afiles <- list.files(path = "../results/admixture", pattern = "*.Q$", full.names = TRUE) 
# exclude k=1,k=10
afiles <- afiles[-(1:2)]

# population ids for each sample 
# text file with two columns (sample and population)
groups <- evec[,c(2,13)] %>% rename(sample=V2)

# read in admixture results to a single df 
q_df <- map_df(afiles, ~vroom(.x, col_names = FALSE, id = "k", delim = " ")) %>%   
  mutate(k = word(k,4,sep="[.]")) %>% # may need to adjust this - wanted to extract k value from file name
  cbind(groups, .) %>%  # combine results with population ids 
  pivot_longer(cols=starts_with("X"), names_to="cluster", values_to="proportion") %>% # convert to long format 
  drop_na()

# specify color palette (replace number in parentheses with max k value)
cols <- colorRampPalette(brewer.pal(8, "Dark2"))(9)

# find a factor order for populations
popmeanprop <- q_df %>% 
  group_by(population,k,cluster) %>% 
  summarize(proportion=mean(proportion)) %>% 
  pivot_wider(id_cols=population,names_from=c(k,cluster),values_from=proportion) %>%
  data.frame()
rownames(popmeanprop) <- popmeanprop[,1]

popord <- (dist(as.matrix(popmeanprop[,-1])) %>% hclust())$order
popord <- popmeanprop[popord,1]

# find a factor order for individuals based on a single k

kvalue <- 6
kk <- q_df %>%
  filter(k==kvalue) %>% 
  pivot_wider(id_cols=sample,names_from=cluster,values_from=proportion) %>% 
  mutate_at(vars(starts_with("X")),function(x){round(x,digits=2)}) %>%
  arrange(across(starts_with("X")))

sampleord <- kk[[1]]

# sort samples by proportion ----------------------------------------------
# https://stackoverflow.com/questions/41679888/how-to-group-order-data-in-r-for-a-barplot

q_ordered <- q_df %>%
  group_by(sample) %>%
  mutate(likely_assignment = cluster[which.max(proportion)], assignment_prob = max(proportion)) %>% # determine most likely assignment 
  arrange(likely_assignment, desc(assignment_prob),.by_group=TRUE) %>% # sort samples by assignment 
  ungroup() %>% # make sure samples are ordered
  mutate(sample = factor(sample,levels=sampleord), population = factor(population,levels=popord)) 

# plot the results
q_ordered %>% 
  ggplot(., aes(sample, proportion, fill = cluster)) +
  geom_bar(stat = "identity",width=1) + # stacked barplot
  scale_fill_manual(values = cols) + # use custom colors
  scale_x_discrete(guide = guide_axis(angle = 90)) + # rotate x-axis labels 
  facet_grid(k ~ population, scales = "free_x", space = "free") + # facet by k value and population id
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90, hjust = 0), axis.text.x = element_text(size = 4), panel.spacing.x = unit(0.001, "lines"), legend.position = "none")   # rotate strip text for easier reading 


######################################
# read in populations.sumstats.tsv.gz
######################################

# gzipped in advance, not by stacks
sumstats <- read.table("../results/stacks/refmap/populations.sumstats.tsv.gz",skip=28,header=TRUE,comment.char="",sep="\t")
	colnames(sumstats) <- colnames(sumstats) %>%
		str_replace(.,"^X..","") %>% 
		str_replace_all(.,"\\.","_")

# plot observed and expected heterozygosity, color points by violations of HWE

filter(sumstats, Pop_ID=="Kup6" & N > 15) %>% 
	ggplot(.,aes(x=Exp_Het,y=Obs_Het,col=-log(HWE_P_value,10) > 3)) + 
	geom_point()

# count up how many populations have HWE deviations for each site

ooHWE <- group_by(sumstats,Locus_ID,Chr,BP) %>% 
	summarize(pops_oo_hwe=sum(HWE_P_value < 0.005))

# how many sites have N deviations from HWE at p < 0.005
table(ooHWE[[4]])


######################################
# dig into VCF, run HDplot
######################################

# source McKinney et al 2017 hdplot function
	# https://github.com/gjmckinney/HDplot
source("hdplot.R")

vcf<-read.vcfR("../results/filtered_vcfs/refmap_final.vcf.gz")

HDres <- HDplot(vcf)

# plot of heterozygosity against read ratio deviation
HDres %>% ggplot() + geom_point(aes(x=H,y=D),size=0.5)

# again but zoom in on Y
	HDres %>% ggplot() + geom_point(aes(x=H,y=D),size=0.5) + ylim(-10,10)


# plot of heterozygosity against read ratio
HDres %>% ggplot()+geom_point(aes(x=H,y=ratio),size=0.5)

# extract the depth of coverage for each individual genotype
dp<-extract.gt(vcf, element = "DP", mask = FALSE, as.numeric = TRUE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)

# plot depth of each snp along the genome, colored by chromosome
plot(rowSums(dp),pch=20,cex=.2,col=factor(vcf@fix[,1]))
	# again but zoom in on the Y
	plot(rowSums(dp),pch=20,cex=.2,col=factor(vcf@fix[,1]),ylim=c(0,100000))

# plot a histogram of depth of coverage for each genotype
hist(rowSums(dp),breaks=1000,xlim=c(0,50000))

hist((dp),breaks=1000,xlim=c(0,300))


# plot total site depth as a function of heterozygosity
plot(HDres$H,rowSums(dp),ylim=c(0,100000),pch=20,cex=.4)


# look more closely at individuals---------------------------

# extract genotypes, recode as 0,1,2
genos<-extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
	genos[is.na(genos)]<-"./."
	genos <- gsub("\\|","/",genos)
	
	genos<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=2,'0/1'=1,'1/0'=1,.default=NA_real_))

# create a matrix of distances between genotypes
gdist <- genos %>% t() %>% dist()

# plot an ugly neighbor-joining tree
	# it's not too hard to add color to the labels
as.matrix(gdist) %>% nj() %>% plot(.,"unrooted",show.tip.label=FALSE)

# make an mds plot (similar to a PCA)
mds <- cmdscale(gdist) 
	colnames(mds) <- c("mds1","mds2")
	mds <- cbind(meta[rownames(mds),c("bioinformatics_id","popmap_initial","region")],mds)
	
ggplot(mds,aes(x=mds1,y=mds2,color=popmap_initial,shape=region)) +
    geom_point()


# again, but zoom in on region "Kup"
# 
gdist <- genos[,meta[colnames(genos),"region"]=="Kup"] %>% t() %>% dist()

nj(gdist) %>% plot(.,"unrooted")

# make an mds plot (similar to a PCA)
mds <- cmdscale(gdist) 
	colnames(mds) <- c("mds1","mds2")
	mds <- cbind(meta[rownames(mds),c("bioinformatics_id","popmap_initial","region")],mds)
	
ggplot(mds,aes(x=mds1,y=mds2,color=popmap_initial,shape=region)) +
    geom_point()

ihz <- colMeans(genos==1,na.rm=TRUE)

beeswarm(ihz ~ meta[colnames(genos),"region"],pch=20,ylab="individual heterozygosity",xlab="region")

#########################################
# how to create a site list for filtering
#########################################


# for use in vcftools: tab separated, two columns, SEQUENCE POSITION

	ex <- rowSums(dp) > 40000 

	vcf@fix[ex,1:2]

	# then use write.table() to write the file

# stacks wants a list of its locus IDs to use as a blacklist

	# this is pretty easy if you're using the sumstats.tsv file
		# grab a set of markers with, say > 5 failed HWE tests and write out the locus IDs

	ooHWE <- group_by(sumstats,Locus_ID,Chr,BP) %>% 
		summarize(pops_oo_hwe=sum(HWE_P_value < 0.005))

	# if you want to filter on information you extracted from the VCF, then you'll need to do a "join" operation


