#The script is to analysis QA_nrt1.1a-data sequenced by SongLY in 2022/11.
# Date: 2024/04/14
# By: hupeiyao.


#### Part I. ASV construction ##########
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)

path<- "~/Documents/Research/LabMember/YR/PLCP_16S/"
setwd(path)
rawseq.dir<-paste0(path,"/trim")
fnFs <- sort(list.files(rawseq.dir, pattern="-R1.fq.gz.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(rawseq.dir, pattern="-R2.fq.gz.trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)
QC.F<-plotQualityProfile(fnFs[1:6])
QC.R<-plotQualityProfile(fnRs[1:6])
ggsave("2_Process/Figure/QC_F.pdf",QC.F,width=20,height=20,units = "cm")
ggsave("2_Process/Figure/QC_R.pdf",QC.R,width=20,height=20,units = "cm")

filtFs <- file.path(path, "2_Process/0_dada2/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "2_Process/0_dada2/filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
write.csv(out,"2_Process/0_dada2/filt_statout.csv")


passid<- !sample.names %in% c("YR_2_14","YR_2_15","YR_2_16")

errF <- learnErrors(filtFs[passid], multithread=TRUE)
errR <- learnErrors(filtRs[passid], multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs[passid], err=errF, multithread=TRUE)
dadaRs <- dada(filtRs[passid], err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs[passid], dadaRs, filtRs[passid], verbose=TRUE)
head(mergers[[1]])
saveRDS(mergers,"2_Process/0_dada2/mergers.RData")
saveRDS(dadaFs,"2_Process/0_dada2/dadaFs.RData")
saveRDS(dadaRs,"2_Process/0_dada2/dadaRs.RData")
saveRDS(errF,"2_Process/0_dada2/errF.RData")
saveRDS(errR,"2_Process/0_dada2/errR.RData")
saveRDS(filtFs,"2_Process/0_dada2/filtFs.RData")
saveRDS(filtRs,"2_Process/0_dada2/filtRs.RData")

seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)#12 6835

sum(seqtab.nochim)/sum(seqtab) #0.7949234

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, "2_Process/0_dada2/dada2_track.csv")

head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/Research/github_16S/silva138/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa2 <- addSpecies(taxa,"~/Documents/Research/github_16S/silva138/silva_species_assignment_v138.1.fa" )
saveRDS(taxa,"2_Process/0_dada2/taxa.RData")

#rm(errF,errR,mergers,seqtab,seqtab.nochim,taxa,dadaFs,dadaRs,filtFs,filtRs,fnFs,fnRs)
save.image(file = "2_Process/0_dada2/dada2.RData")
#### Part II. Phyloseq-object construct ############
library(Biostrings)
library(dplyr)
library(phyloseq)
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)

path<- "~/Documents/Research/LabMember/YR/PLCP_16S/"
setwd(path)
load("2_Process/0_dada2/dada2.RData")
asv.seq<-colnames(seqtab.nochim)
asv.seq<-DNAStringSet(asv.seq)
names(asv.seq)<-paste0("ASV_",1:length(asv.seq))
dir.create("2_Process/0_dada2/Table")
writeXStringSet(asv.seq,"2_Process/0_dada2/Table/Table0.asvseq.fasta")

seqtab.nochim.bk<-seqtab.nochim
colnames(seqtab.nochim.bk)<- paste0("ASV_",1:length(asv.seq))
asv.tab<-t(seqtab.nochim.bk)
rm(seqtab.nochim.bk)
write.csv(asv.tab,"2_Process/0_dada2/Table/Table0.asvtab.csv")

asv.tax<-readRDS("2_Process/0_dada2/taxa.RData")
rownames(asv.tax)<-paste0("ASV_",1:nrow(asv.tax))
write.csv(asv.tax,"2_Process/0_dada2/Table/Table0.asvtax.csv")

asv.tax<- read.csv("2_Process/0_dada2/Table/Table0.asvtax.csv",header = T,row.names = 1)
asv.tab<- read.csv("2_Process/0_dada2/Table/Table0.asvtab.csv",header = T,row.names = 1)
# metadata<- read.table("metadata.txt",sep="\t",header = T)
metadata<- read.table("metadata_v2.txt", sep="\t", header = T)
rownames(metadata)<- metadata$Seq_ID
table(colnames(asv.tab) %in% rownames(metadata))
metadata<- metadata[rownames(metadata) %in% colnames(asv.tab),] # 3 samples were removed due to low depth in version1
asv.tab<- asv.tab[,rownames(metadata)]


raw.ps<-phyloseq(otu_table(as.matrix(asv.tab),taxa_are_rows = T),
                 tax_table(as.matrix(asv.tax)),
                 sample_data(metadata)) 
# version1: 27858 taxa and 46 samples 
# version2: 27858 taxa and 31 samples
sample_sums(raw.ps)
summary(sample_sums(raw.ps))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 53156   88483  117878  125639  162645  239731  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 53156   84153  118707  126894  164660  239731 

source("~/Documents/Research/mGWAS/mGWAS2024/999_R_function/ps_filter.R")
ps_filter<- ps_filter(raw.ps,prev_n = 3)

ps.filter<- ps_filter$ps.filter 
#3735 taxa and 46 samples 
#3119 taxa and 31 samples
# write.csv(otu_table(ps.filter), "Table/Table0.asvtab_filt3735.csv")
# write.csv(tax_table(ps.filter), "Table/Table0.asvtax_filt3735.csv")
write.csv(otu_table(ps.filter), "Table/Table0.asvtab_filt3119.csv")
write.csv(tax_table(ps.filter), "Table/Table0.asvtax_filt3119.csv")

library(tibble)
ps_filter$discard_log %>%  
  rownames_to_column("Seq_ID") %>%
  left_join(metadata) 

source("~/Documents/Research/mGWAS/mGWAS2024/999_R_function/evalution_retained_otu.R")
retain_prevalence(raw.ps,
                  filtN = 0, RA = 0.0001, prevN = 0.8, var_name = "raw.ps")
#### for LEFse ####
rarefy.ps.filter <- rarefy_even_depth(ps.filter, sample.size = min(sample_sums(ps.filter))) #3119 taxa and 31 samples; depth:50632
rarefy.ps.filter %>% otu_table() %>% data.frame() %>% write.csv("Table/Table0.asvtab_filt3119_rarefydepth50632.csv")
rarefy.ps.filter %>% tax_table() %>% data.frame() %>% write.csv("Table/Table0.asvtax_filt3119_rarefydepth50632.csv")
#### END ####

rarefy.root.ps<- prune_samples(data.frame(sample_data(ps.filter))$compartment == "Root",ps.filter) 
rarefy.root.ps<-   rarefy_even_depth(rarefy.root.ps, sample.size = min(sample_sums(rarefy.root.ps)))
rarefy.ps<- rarefy_even_depth(ps.filter, sample.size = min(sample_sums(ps.filter)))
sample_sums(rarefy.root.ps) 
#50585 in version1
#50632 in version2
########## CPCoA #######
library(amplicon)
beta_cpcoa(data.frame(otu_table(rarefy.root.ps)),
           data.frame(sample_data(rarefy.root.ps)),
           dis="bray",
           groupID= "Line")+
  scale_fill_manual(values = color_subspecies2)

ggsave("Figure/Figure2.CPCoA.pdf", width=5, height = 4)


##### alpha #######
method = c("Observed","Shannon","Simpson","Chao1")
# group_level<- c("BS","WT","Line5-7",  "Line5-2",  "Line22-1", "Line22-3", "Line49-3", "Line9-1" )
group_level<- c("BS","WT","5-7",  "22-1", "9-1" )

# set  TukeyHSD_Df
TukeyHSD_df<- data.frame()
#color_subspecies2<- c("#8A4B08","#FFCF48","#A3C1AD","#1E88E5","#F1A9A0")

 sub_rps<- rarefy.ps
  sub_alpha<- estimate_richness(sub_rps,measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
  sub_design<- sub_rps%>%
    sample_data() %>%
    data.frame() %>% 
    cbind(sub_alpha[rownames(.),])
  
  color_subspecies2<- c("#8A4B08","#A3C1AD","#756bb1","#bcbddc","#3182bd","#9ecae1","#F1A9A0","#f7fcb9")
  # if(i!=1){color_subspecies2<- color_subspecies2[c(1,5:7)]}
  library(agricolae)
  alpha_plot_list<- list()
  i=1
  for(mi in 1:4){
    #m = "shannon"
    m=method[mi]
    model = aov(sub_design[[m]] ~ Line, data=sub_design)
    Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
    Tukey_HSD_table = as.data.frame(Tukey_HSD$Line) 
    Tukey_HSD_table$group<- m
    
    TukeyHSD_df<- rbind(TukeyHSD_df, Tukey_HSD_table)
    suppressWarnings(write.table(TukeyHSD_df, file=paste("./Table/Table2.alpha_TukeyHSD",m,".txt",sep=""), 
                                 append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", 
                                 row.names = T, col.names = T))
    
    
    # LSD test for stat label
    out = LSD.test(model,"Line", p.adj="fdr") # alternative fdr
    stat = out$groups
    sub_design$stat=stat[as.character(sub_design$Line),]$groups
    max=max(sub_design[,c(m)])
    min=min(sub_design[,c(m)])
    x = sub_design[,c("Line",m)]
    y = x %>% group_by(Line) %>% summarise_(Max=paste('max(',m,')',sep=""))
    y=as.data.frame(y)
    rownames(y)=y$Line
    sub_design$y=y[as.character(sub_design$Line),]$Max + (max-min)*0.1
    
    # set factor levels:
    level1<-group_level 
    sub_design$Group=factor(sub_design$Line,levels = level1)
    
    # to plot
    if(m == "Shannon"){
      p<- ggplot(sub_design,aes(x=Group, y=.data[[m]]))+theme_minimal()+
        geom_boxplot(aes(fill=Group),alpha=.75,lwd=.2,width=.5,color="black",outlier.shape = NA)+
        # geom_boxplot(aes(fill=Group),alpha=.75,outlier.size = .5,linewidth=.3)+
        geom_jitter(aes(color=Group),position=position_jitter(0.17), size=.5, alpha=0.7)+scale_y_continuous(limits = c(0,1.2*max))+
        labs(title="",x="Compartment",y=m)+geom_text(data=sub_design, aes(x=Group, y=y, color=Group, label= stat))+
        theme(
          legend.position="none",
          axis.title.y=element_text(hjust=.5),
          plot.title=element_text(size=11,hjust=.5),
        )+ylim(c(0,8))+
        scale_fill_manual(values = color_subspecies2)+
        scale_color_manual(values = color_subspecies2)+
        theme(axis.ticks = element_line(linewidth = .5),
              axis.text.x = element_text(size=12,angle = 45,hjust = 1,colour = "black"),
              axis.text.y = element_text(size=12,colour = "black"),
              axis.title = element_text(size = 15,colour = "black")
        )+theme(legend.position = "none")+
        ggtitle(m)+labs(x="")
    }else{
      p<- ggplot(sub_design,aes(x=Group, y=.data[[m]]))+theme_minimal()+
        geom_boxplot(aes(fill=Group),alpha=.75,lwd=.2,width=.5,color="black",outlier.shape = NA)+
        # geom_boxplot(aes(fill=Group),alpha=.75,outlier.size = .5,linewidth=.3)+
        geom_jitter(aes(color=Group),position=position_jitter(0.17), size=.5, alpha=0.7)+scale_y_continuous(limits = c(0,1.2*max))+
        labs(title="",x="Compartment",y=m)+geom_text(data=sub_design, aes(x=Group, y=y, color=Group, label= stat))+
        theme(
          legend.position="none",
          axis.title.y=element_text(hjust=.5),
          plot.title=element_text(size=11,hjust=.5),
        )+
        scale_fill_manual(values = color_subspecies2)+
        scale_color_manual(values = color_subspecies2)+
        theme(axis.ticks = element_line(linewidth = .5),
              axis.text.x = element_text(size=12,angle = 45,hjust = 1,colour = "black"),
              axis.text.y = element_text(size=12,colour = "black"),
              axis.title = element_text(size = 15,colour = "black")
        )+theme(legend.position = "none")+
        ggtitle(m)+labs(x="")
    }
    
    ggsave(paste("./Figure/Figure2A-",m, ".pdf", sep=""), p, width = 8, height = 9,units = "cm")
  }


##### 4) Compositional ####
comp_list<- list()
stack_comp_list<- list()
# for(i in 1:3){
  # source data (otutab and taxonomy files) preparation
  sub_rps<- rarefy.ps
  rarefytab<- sub_rps %>% otu_table() %>% data.frame()
  taxa_rarefy<- sub_rps %>% tax_table() %>% data.frame()
  sub_design<- sub_rps%>% sample_data() %>% data.frame()
  colnames(taxa_rarefy) = c("Kingdom","Phylum","Class","Order","Family","Genus")
  
  # select Proteobacteria line
  idx = taxa_rarefy$Phylum == "Proteobacteria"
  idx[is.na(idx)]<- FALSE
  taxa_rarefy$full=as.character(taxa_rarefy$Phylum) 
  taxa_rarefy[idx,]$full=as.character(taxa_rarefy[idx,]$Class)
  tax_count = merge(taxa_rarefy, rarefytab, by="row.names")
  
  tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[8], FUN=sum) # mean
  rownames(tax_count_sum) = tax_count_sum$full
  tax_count_sum = tax_count_sum[,-1]
  per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100
  dim(per)
  
  mean_sort = per[(order(-rowSums(per))), ] # decrease sort
  colSums(mean_sort)
  write.csv(per,paste("./Table/Table3.taxonomy_sum_Phyla_PLCP",".csv",sep=""))
  
  phyla_display<- c("Gammaproteobacteria" ,"Alphaproteobacteria", "Acidobacteriota"  ,   "Actinobacteriota" ,
                    "Bacteroidota"    ,    "Chloroflexi"    ,   "Firmicutes",
                    "Myxococcota"  ,        "Spirochaetota" )
  mean_sort=as.data.frame(mean_sort)
  other = colSums(mean_sort[!rownames(mean_sort) %in% phyla_display,])
  mean_sort = mean_sort[phyla_display, ]
  mean_sort = rbind(mean_sort,other)
  rownames(mean_sort)[length(phyla_display)+1] = c("Others")
  #mean_sort<-mean_sort[c(sort(rownames(mean_sort)[1:num_top]),rownames(mean_sort)[num_top+1]),]
  write.table(mean_sort, paste("./Table/Table3.2.Top10phylum_ProClass_","PLCP.txt",sep=""), append = F, sep="\t", quote=F, row.names=T, col.names=T)
  
  
  topN=rownames(mean_sort)
  topN<- c(topN[grep("proteo",topN)],  topN[-grep("proteo",topN)])
  #select group
  sub_design 
  sub_design$Group  = factor(sub_design$Line, levels=group_level)
  
  
  mean_sort = mean_sort[, rownames(sub_design)] # reorder according to design
  mean_sort$phylumpro = rownames(mean_sort)
  mean_sort$phylumpro =factor(mean_sort$phylumpro,levels = rev(topN))
  colotax<-c("#20B2AA","#9AD9CA","#E9AFA3","#FA8072","#F1E684",
             "#466983","#D8BFD8","#5DB1DD", "#F59ABE","#EADEBD")
  colotax2<- rev(colotax)
  
  # if(i!=2){
  sample_fact=t(mean_sort) %>% data.frame() %>%
    rownames_to_column("SampleID") %>%
    arrange(Gammaproteobacteria) %>% pull(SampleID)
  sample_fact= sample_fact[-length(sample_fact)]
  # }else{
  #   sample_fact=t(mean_sort) %>% data.frame() %>%
  #     rownames_to_column("SampleID") %>% 
  #     arrange(Firmicutes) %>% pull(SampleID)
  #   sample_fact= sample_fact[-length(sample_fact)]
  #   mean_sort$phylumpro =factor(mean_sort$phylumpro,levels = c("Others", "Spirochaetota","Myxococcota" , "Chloroflexi","Bacteroidota" ,"Actinobacteriota",    "Acidobacteriota", "Gammaproteobacteria" ,"Alphaproteobacteria","Firmicutes" ))
  #   colotax2<- c("#EADEBD", "#F59ABE", "#5DB1DD", "#466983", "#F1E684","#FA8072" , "#E9AFA3", "#9AD9CA", "#20B2AA","#D8BFD8")
  # }
  
  
  data_all = as.data.frame(melt(mean_sort, id.vars=c("phylumpro")))
  data_all = merge(data_all, sub_design[c("Group")], by.x="variable", by.y = "row.names")
  data_all$variable<- factor(data_all$variable, levels = sample_fact)
  write.csv(data_all,paste("./Table/Table3.allPro_","PLCP",".csv",sep=""))
  
a<- data_all %>%
  mutate(variable=as.character(data_all$variable)) %>%
  arrange(Group, variable) 
  
data_all<- data_all %>%
  mutate(variable = factor(data_all$variable, levels= unique(a$variable)))
  #"#F39C12","#802268"
  #if(i==1){
  #  colotax2<- colotax[-8]
  #}else{colotax2<- colotax}
  p = ggplot(data_all, aes(x=variable, y = value, fill = phylumpro )) + 
    geom_bar(stat = "identity",position="fill", width=1)+ 
    scale_y_continuous(labels = scales::percent) + 
    # facet_grid( ~ Group, scales = "free_x", switch = "x")  +theme_classic()+
    theme(axis.ticks.x = element_blank(), 
          # legend.position="top", 
          axis.text.x = element_text(angle=90) ,
          strip.background = element_blank())+
    xlab("Groups")+ylab("Percentage (%)")+
    scale_fill_manual(values = colotax2)+
    ggtitle("PLCP")
  p
  ggsave(paste("./Figure/Figure3_AllTaxonomy_PLCP","v3.pdf", sep=""), p, width = 10, height = 6)
  
  
  mat = mean_sort[,1:(dim(mean_sort)[2]-1)]
  mat_t = t(mat)
  mat_t2 = merge(sub_design[c("Line")], mat_t, by="row.names")
  mat_t2 = mat_t2[,-1]
  mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
  mat_mean_final = do.call(rbind, mat_mean)[-1,]
  geno = mat_mean$Line
  colnames(mat_mean_final) = geno
  mat_mean_final = as.data.frame(mat_mean_final)
  mat_mean_final$phylumpro = rownames(mat_mean_final)
  data_all = as.data.frame(melt(mat_mean_final, id.vars=c("phylumpro")))
  data_all$variable<-factor(data_all$variable,levels=group_level)
  phyla_order<- c(unique(data_all$phylumpro)[grep("proteo",unique(data_all$phylumpro))],
                  unique(data_all$phylumpro)[-grep("proteo",unique(data_all$phylumpro))])
  data_all$phylumpro<-factor(data_all$phylumpro,levels = rev(phyla_order))
  data_all$value<- as.numeric(data_all$value)
  p=ggplot(data_all, aes(x=variable, y = value, fill = phylumpro )) + 
    geom_bar(stat = "identity",position="fill", width=0.7)+ 
    scale_y_continuous(labels = scales::percent) + 
    xlab("Compartment")+ylab("Relative Abundance (%)")+theme_classic()+
    labs(fill="Taxonomy")+
    theme(axis.text.x = element_text(size=12,angle=45,hjust = 1,vjust = 1),
          axis.line.x=element_line(size=.5),
          axis.line.y=element_line(size=.5),
          axis.text=element_text(size=12, color="black"),
          axis.title = element_text(size=12, color="black"),
          legend.text= element_text(size=6),
          legend.title = element_text(size=7),
          legend.key.size = unit(.3, 'cm'),
          legend.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))+
    scale_fill_manual(values = colotax2)+
    ggtitle("PLCP")
  p
  ggsave(paste("Figure/Figure3B_","PLCP_all",".pdf",sep=""),width = 12,height = 10,units = "cm")
  write.csv(data_all,paste("Table/Figure3B.taxonomyplot_all","PLCP",".csv",sep=""))

