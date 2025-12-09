setwd("~/Desktop/EXPORTS_2018_Updated")
#### Load Libraries ####
library("devtools")
library("dada2")
library("ggplot2")
library("phyloseq")
library("vegan")
library("DESeq2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("ape")
library("ggrepel")
library("tibble")
library("ggsci")
library("plyr")
library("dplyr")
library("treemap")
library(wesanderson)
library(patchwork)
library("ggbreak")
library("ggpubr")
library("gridExtra")
library(microViz)

#### Load Data and Make Physeq Object ####
rev18_count_tab <- read.table("revelle_pr2_bradley_ASV_counts.tsv", header=T, row.names=1,check.names=F, sep="\t")
rev18_tax_tab <- as.matrix(read.table("revelle_pr2_bradley_taxonomy.tsv",header=T, row.names=1, check.names=F, sep="\t"))
rev18_sample_info_tab <- read.table("exports2018_revelle_18S_samples.tsv", header=T, row.names=1,check.names=F, sep="\t")

rev18_sample_info_tab[is.na(rev18_sample_info_tab)] <- 0
rev18_sample_info_tab$position <- cut(rev18_sample_info_tab$light,
                                      breaks=c(0.0, 0.10, 0.65),
                                      labels=c('1-10%','20-65%'))
rev18_sample_info_tab$light_level[rev18_sample_info_tab$light_level == "40%" |rev18_sample_info_tab$light_level == "65%" ] <- "40-65%"

count_tab <- otu_table(rev18_count_tab, taxa_are_rows=T)
tax_tab <- tax_table(rev18_tax_tab)
sample_info <- sample_data(rev18_sample_info_tab)

rev18_sample_info_tab <- rev18_sample_info_tab %>% mutate(depth_position = case_when(
  depth < MLD ~ "above",
  depth > MLD ~ "below"
))

ASV_physeq_npac <- phyloseq(count_tab, tax_tab, sample_info)

classification_counts <- read.table("Classification_Mode_Counts.tsv", header=T, check.names=F, sep="\t")

#### Basic Pruning of Unnecessary Samples and Zero sum rows ####
ASV_exports <- subset_samples(ASV_physeq_npac, !type=="Fe_addition" & !type=="deep_grazing" & !type=="nutr_addition")
# 54 Samples, 2813 ASVs 
ASV_exports_pruned1 <- prune_taxa(taxa_sums(ASV_exports) > 0, ASV_exports)
# 54 samples, 2401 ASVs
# 731 syndiniales
ASV_exports_mz <- subset_taxa(ASV_exports_pruned1, Division == "Dinoflagellata" | Division == "Ciliophora")
# 54 Samples, 1351 ASVs
# 731 Syndiniales ASVs ; 585 ASVs dino & ciliate remain
# 409 Dinoflagellate ASVs
# 176 Ciliate ASVs
ASV_exports_pruned2 <- prune_taxa(taxa_sums(ASV_exports_mz) > 50, ASV_exports_mz)
# 54 Samples, 712 ASVs
# Means that roughly half of the ID'd ASVs had < 50 total read across all 54 samples (super low read count)
# 50 estimates ~ 0.3% (maximum) of a sample
ASV_exports_phyto_sub <- subset_taxa(ASV_exports_pruned1, Division == "Chlorophyta" | Division=="Cryptophyta"| Division=="Haptophyta"| Division=="Ochrophyta")
ASV_exports_phyto_pruned <- prune_taxa(taxa_sums(ASV_exports_phyto_sub) > 50, ASV_exports_phyto_sub)

#### FIGURE 1 ASV Totals by Trophic Mode  ####
# Of total euk dataset by order or family
# Show # of ASVs
counts_for_supp <- as.data.frame(otu_table(ASV_exports_pruned2))
write.table(counts_for_supp, file="counts_for_supp.tsv", sep="\t", quote=F)

taxa <- as.data.frame(tax_table(ASV_exports_pruned2))
write.table(taxa, file="modes_added_tax_table.tsv", sep="\t", quote=F)
asv_totals_by_level <- taxa %>% group_by(Division, Class, Order, Family, Genus) %>% count

#not shown -- load classification_counts dataframe 
classification_counts$Level <- factor(classification_counts$Level, levels=c("Phylum","Class","Order","Family","Genus"))
classification_counts$Division <- as.factor(classification_counts$Division)
classification_counts$mode <- as.factor(classification_counts$mode)

ggplot(classification_counts, aes(x=mode, y=n, fill=Level)) +
  geom_bar(position="stack",width=1, stat="identity", colour="black") +
  facet_wrap(~Division) +
  scale_y_continuous(limits=c(0,375),expand = expansion(mult = c(0, 0.05))) +
  scale_y_break(breaks=c(110, 325), ticklabels=c(325,350,375))+
  scale_fill_manual(values = c('#117733','#999933','#DDCC77', '#CC6677','#882255'))+
  #c('#0077BB', '#88CCEE', '#228833', '#CCBB44','#882255'))+
  theme_bw()+
  theme(axis.text.y.right=element_blank(), axis.ticks.y.right=element_blank(),axis.title.y.left = element_text(size=12), axis.text.y.left = element_text(size=12, color="black"), 
        legend.key.size = unit(1,"line"), text=element_text(size=12, color = "black"),legend.title=element_text())+
  labs(x= " ", y="Total ASVs")

#### FIGURE 2 ASV TOTAL BY ORDER ####
asv_totals_byorder <- taxa %>% group_by(Division, Order) %>% count
test1 <- asv_totals_byorder %>%  mutate(Order=replace(Order, Division == "Dinoflagellata" & is.na(Order), "D_Unknown"))
test2 <- test1 %>%  mutate(Order=replace(Order, Division == "Ciliophora" & is.na(Order), "C_Unknown"))

# add trophic mode for orders 
#Heterotrophic (red)
#Mixotrophic (green)
#Multiple (Yellow)
#Parasitic (blue)
#Unknown (grey)
modes <- c("Heterotrophic","Heterotrophic","Heterotrophic","Heterotrophic","Heterotrophic",
           "Heterotrophic","Heterotrophic","Heterotrophic","Mixotrophic","Mixotrophic","Mixotrophic",
           "Mixotrophic","Heterotrophic","Heterotrophic","Unclassified","Parasitic","Parasitic",
           "Parasitic","Parasitic","Parasitic","Unclassified","Heterotrophic","Mixotrophic","Multiple","Multiple",
           "Mixotrophic","Mixotrophic","Mixotrophic","Unclassified")
orders_mode <-cbind(test2,modes)
orders_mode <- rename(orders_mode, mode = ...4)
orders_mode$Order <- as.factor(orders_mode$Order)
orders_mode$Division <- as.factor(orders_mode$Division)
orders_mode$mode <- as.factor(orders_mode$mode)
orders_mode <- orders_mode[order(orders_mode$n),]

dinos <- subset(orders_mode, Division == "Dinoflagellata")
cil <- subset(orders_mode, Division == "Ciliophora")

# All Orders combined ciliate / dinoflagellates
figure3_asvs_byorder <- ggplot(orders_mode, aes(x=reorder(Order,-n), y=n, fill=mode)) +
  geom_bar(width=1, stat="identity", color="black") +
  scale_y_continuous(limits=c(0,285),expand = expansion(mult = c(0, 0.05))) +
  scale_y_break(breaks=c(75, 275), ticklabels=c(275,280,285))+
  scale_fill_manual(values=c("#800000FF","#8A9045FF","#FFA319FF","#155F83FF", "#767676FF"))+
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45,hjust=1, color="black"),axis.text.y.right=element_blank(), axis.ticks.y.right=element_blank(),axis.title.y.left = element_text(size=12), axis.text.y.left = element_text(size=12, color="black"), 
        legend.key.size = unit(1,"line"), text=element_text(size=12, family="Times", color = "black"),legend.title=element_text())+
  labs(x= "Order", y="Total ASVs", fill = "Feeding Modes") 
ggsave(plot=figure3_asvs_byorder,"figure3_asvs_byorder.pdf", device="pdf", width=9.8, height=7.5, units = c("in"))

#### Presence/Absence Table Write Out ####
Species_PA <- ASV_exports_pruned2 %>% tax_transform("binary")
#write.table(otu_table(Species_PA), file="Species_PA.tsv", sep="\t", quote=F, col.names=NA)
#write.table(sample_data(ASV_exports_pruned2), file="ASV_exports_pruned2_samples.tsv", sep="\t", quote=F, col.names=NA)

#### Data CLR Transformation ####
exports_clr <- microbiome::transform(ASV_exports_pruned2, 'clr')

#### Trophic Richness Data ####
#ASV_exports_pruned2
#Heterotrophic (red) #800000FF
#Mixotrophic (blue) #155F83FF
#Multiple (Yellow) #FFA319FF
#Parasitic (green) #8A9045FF
#Unknown (grey) #767676FF

richness <- read.table("MS1_Richness_ByMode.tsv", header=T, check.names=F, sep="\t")
richness_long <- gather(richness, key="mode",value="asvs",Heterotroph_ASV,Mixotroph_ASV, Parasite_ASV, Unclassified)

#### FIGURE 3 Richness over time & statistical testing ####
# TUKEY-KRAMER TEST 
# days with n < 3:  229, 232, 237, 240 (n=2), 242, 248
# days with n >= 3: 228, 230, 236, 241, 244, 245, 246, 249, 250
n_subset_high <- subset(richness, day =="228" | day == "230" | day == "236" | 
                          day == "241" | day == "244" | day == "245" | day== "246" | day == "249" | day=="250")

par(mfrow = c(2, 2))
het_lm <- lm(Heterotroph_ASV~as.factor(day), data=n_subset_high)
het_anova <- anova(het_lm)
# Analysis of Variance Table
# Response: Heterotroph_ASV
#                 Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(day)  8  245.19  30.649  0.9647 0.4775
# Residuals      38 1207.23  31.769
plot(TukeyHSD(aov(het_lm)))
het_lm <- lm(Heterotroph_ASV~as.factor(day), data=richness)
het_anova <- anova(het_lm)
# Analysis of Variance Table
# Response: Heterotroph_ASV
#                 Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(day) 14  407.6  29.114   0.879 0.5859
# Residuals      39 1291.7  33.121 
plot(TukeyHSD(aov(het_lm)))
hist(richness$Heterotroph_ASV)
hist(n_subset_high$Heterotroph_ASV)
# No difference in mean heterotroph ASVs, all CIs contain zero including those with days n<3
het_kruskal <- kruskal.test(Heterotroph_ASV~as.factor(day), data=n_subset_high)
# Kruskal-Wallis rank sum test
# A Kruskal-Wallis test is a non-parametric analog of a one-way ANOVA. 
#It does not assume that the variable has a normal distribution. 
#(Instead, it tests whether the variable has the same distribution with the same mean in each group.)
# data:  Heterotroph_ASV by as.factor(day)
# Kruskal-Wallis chi-squared = 7.164, df = 8, p-value = 0.519
het_kruskal <- kruskal.test(Heterotroph_ASV~as.factor(day), data=richness)
# Kruskal-Wallis rank sum test
# data:  Heterotroph_ASV by as.factor(day)
# Kruskal-Wallis chi-squared = 13.545, df = 14, p-value = 0.4842

mix_lm <- lm(Mixotroph_ASV~as.factor(day), data=n_subset_high)
mix_anova <- anova(mix_lm)
# Analysis of Variance Table
# Response: Mixotroph_ASV
#                Df Sum Sq Mean Sq F value   Pr(>F)   
# as.factor(day)  8 1611.1 201.391    3.73 0.002629 **
# Residuals      38 2051.7  53.992   
plot(TukeyHSD(aov(mix_lm)))
# YD 244-236, 245-236, 250-236 SIG DIFF p < 0.05
mix_lm <- lm(Mixotroph_ASV~as.factor(day), data=richness)
mix_anova <- anova(mix_lm)
# Analysis of Variance Table
# Response: Mixotroph_ASV
#                Df Sum Sq Mean Sq F value   Pr(>F)   
# as.factor(day) 14 2209.1  157.79  2.8627 0.004817 **
# Residuals      39 2149.7   55.12                   
plot(TukeyHSD(aov(mix_lm)))
# YD 244-236, 245-236, 250-236 SIG DIFF p < 0.05
hist(richness$Mixotroph_ASV)
hist(n_subset_high$Mixotroph_ASV)
# Significant difference 
mix_kruskal <- kruskal.test(Mixotroph_ASV~as.factor(day), data=n_subset_high)
# Kruskal-Wallis rank sum test
# data:  Mixotroph_ASV by as.factor(day)
# Kruskal-Wallis chi-squared = 17.424, df = 8, p-value = 0.02598
mix_kruskal <- kruskal.test(Mixotroph_ASV~as.factor(day), data=richness)
# Kruskal-Wallis rank sum test
# data:  Mixotroph_ASV by as.factor(day)
# Kruskal-Wallis chi-squared = 23.939, df = 14, p-value = 0.0466

par_lm <- lm(Parasite_ASV~as.factor(day), data=n_subset_high)
par_anova <- anova(par_lm)
# Analysis of Variance Table
# Response: Parasite_ASV
#               Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(day)  8   4228  528.50   1.451 0.2077
# Residuals      38  13840  364.22  
plot(TukeyHSD(aov(par_lm)))
# No difference in mean parasite ASVs, all CIs contain zero including those with days n<3
par_lm <- lm(Parasite_ASV~as.factor(day), data=richness)
par_anova <- anova(par_lm)
# Analysis of Variance Table
# Response: Parasite_ASV
#                 Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(day) 14  6231.7  445.12  1.1652 0.3383
# Residuals      39 14898.4  382.01
hist(richness$Parasite_ASV)
par_kruskal <- kruskal.test(Parasite_ASV~as.factor(day), data=n_subset_high)
# Kruskal-Wallis rank sum test
# data:  Parasite_ASV by as.factor(day)
# Kruskal-Wallis chi-squared = 11.535, df = 8, p-value = 0.1732
par_kruskal <- kruskal.test(Parasite_ASV~as.factor(day), data=richness)
# Kruskal-Wallis rank sum test
# data:  Parasite_ASV by as.factor(day)
# Kruskal-Wallis chi-squared = 18.837, df = 14, p-value = 0.1713

unc_lm <- lm(Unclassified~as.factor(day), data=n_subset_high) #only days with n >= 3
unc_anova <- anova(unc_lm)
# Analysis of Variance Table
# Response: Unclassified
# Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(day)  8 325.65  40.707  1.5561 0.1709
# Residuals      38 994.05  26.159               
plot(TukeyHSD(aov(unc_lm)))
# No difference in mean unclassified ASVs, all CIs contain zero including those with days n<3
unc_lm <- lm(Unclassified~as.factor(day), data=richness) 
unc_anova <- anova(unc_lm)
# Analysis of Variance Table
# Response: Unclassified
#                 Df Sum Sq Mean Sq F value Pr(>F)
# as.factor(day) 14 577.88  41.277  1.6162 0.1179
# Residuals      39 996.05  25.540     
hist(richness$Unclassified)
hist(n_subset_high$Unclassified)
uncl_kruskal <- kruskal.test(Unclassified~as.factor(day), data=n_subset_high)
# Kruskal-Wallis rank sum test
# data:  Unclassified by as.factor(day)
# Kruskal-Wallis chi-squared = 10.845, df = 8, p-value = 0.2106
uncl_kruskal <- kruskal.test(Unclassified~as.factor(day), data=richness)
# Kruskal-Wallis rank sum test
# data:  Unclassified by as.factor(day)
# Kruskal-Wallis chi-squared = 19.787, df = 14, p-value = 0.137

###
heterotroph_time <- ggplot(richness, aes(x=as.factor(day), y=Heterotroph_ASV))+
  geom_boxplot(fill="#800000FF", color='black',outlier.shape = 21)+
  theme_bw()+
  #ggtitle("Heterotroph")+
  theme(plot.title = element_text(hjust=1),axis.text.x=element_blank(),text=element_text(size=14,  color = "black")) +
  ylim(20,50)+
  ylab("")+
  xlab("")
mixotroph_time <- ggplot(richness, aes(x=as.factor(day), y=Mixotroph_ASV))+
  geom_boxplot(fill="#155F83FF",color='black',outlier.shape = 21)+
  theme_bw()+
  #ggtitle("Mixotroph")+
  theme(plot.title = element_text(hjust=1),axis.text.x=element_blank(),text=element_text(size=14, color = "black")) +
  ylim(45,90)+
  ylab("")+
  xlab("")
parasite_time <- ggplot(richness, aes(x=as.factor(day), y=Parasite_ASV))+
  geom_boxplot(fill="#8A9045FF",color='black', outlier.shape = 21)+
  theme_bw()+
  #ggtitle("Parasite")+
  theme(plot.title = element_text(hjust=1),axis.text.x=element_blank(),text=element_text(size=14, color = "black")) +
  ylim(60,185)+
  ylab("Observed ASV Richness")+
  xlab("")
unclassified_time <- ggplot(richness, aes(x=as.factor(day), y=Unclassified))+
  geom_boxplot(fill="#767676FF",color='black',outlier.shape = 21)+
  theme_bw()+
  #ggtitle("Unclassified")+
  theme(plot.title = element_text(hjust=1),text=element_text(size=14, color = "black")) +
  ylim(15,50)+
  ylab("")+
  xlab("Year-Day")
figure3_asvs_bymode_time <- plot_spacer() / heterotroph_time / plot_spacer() / mixotroph_time / plot_spacer() / parasite_time / plot_spacer() / unclassified_time / plot_spacer() + plot_layout(heights = c(2,5, -1.1 ,5,-1.1,5,-1.1, 5,2))
ggsave("figure3_asvs_bymode_time.png")
ggsave(plot=figure3_asvs_bymode_time,"figure3_asvs_bymode_time_v2.pdf", device="pdf")

#### SUPP FIGURE S2 - SAMPLING DEPTHS OVER TIME ####
sampling_time_z <- ggplot(richness, aes(x=as.factor(day)))+
  geom_boxplot(aes(y=depth),fill="lightgrey", color='black',outlier.shape = 21)+
  scale_y_reverse()+
  theme_bw()+
  theme(plot.title = element_text(hjust=1),text=element_text(size=14,  color = "black")) +
  ylab("Depth (m)")+
  xlab("")

sampling_time_par <- ggplot(richness, aes(x=as.factor(day),y=light))+
  geom_boxplot(fill="lightgrey", color='black',outlier.shape = 21)+
  theme_bw()+
  theme(plot.title = element_text(hjust=1),axis.text.x=element_blank(),text=element_text(size=14,  color = "black")) +
  ylab("")+
  xlab("")
sampling_time_MLD <- ggplot(richness, aes(x=as.factor(day), y=MLD))+
  geom_point(color='black')+
  scale_y_reverse(limits=c(40,20))+
  geom_hline(yintercept=29.23)+
  theme_bw()+
  theme(plot.title = element_text(hjust=1),text=element_text(size=14,  color = "black")) +
  ylab("")+
  xlab("")
sampling_time_par /sampling_time_MLD
sampling_variation <- sampling_time_par / plot_spacer() /sampling_time_MLD + plot_layout(heights = c(4, -1.1 ,4))
ggsave(plot=sampling_time_MLD,"sampling_time_MLD.pdf", device="pdf", width=8, height=4, units = c("in")) 

#### FIGURE 6: ASV Richness ENV Correlates  ####
library("ggpubr")
zxn <- ggscatter(richness, x = "depth", y = "NO3", color= "black",
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Depth (m)", ylab = "Nitrate")
zxp <- ggscatter(richness, x = "depth", y = "PO4", color= "black",
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Depth (m)", ylab = "Phosphate")
zxt <- ggscatter(richness, x = "depth", y = "Temp", color= "black",
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Depth (m)", ylab = "Temperature")
zxl <- ggscatter(richness, x = "depth", y = "light", color= "black",
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab = "Depth (m)", ylab = "Light (% PAR)")
Depth_correlates <- zxn / zxp | zxt /zxl
ggsave(plot=Depth_correlates,"Depth_correlates.pdf", device="pdf", width=7, height=7, units = c("in")) 


Het_z <- ggscatter(richness, x = "depth", y = "Heterotroph_ASV", color= "#800000FF",
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = " ", ylab = "Heterotroph Richness")
Het_par <- ggscatter(richness, x = "light", y = "Heterotroph_ASV", color= "#800000FF",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = " ", ylab = " ")
Het_n <- ggscatter(richness, x = "NO3", y = "Heterotroph_ASV", color= "#800000FF",
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = " ", ylab = " ")
Het_p <- ggscatter(richness, x = "PO4", y = "Heterotroph_ASV", color= "#800000FF",
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = " ", ylab = " ")
Het_temp <- ggscatter(richness, x = "Temp", y = "Heterotroph_ASV", color= "#800000FF",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = " ", ylab = " ")
Het_chl <- ggscatter(richness_v2, x = "chl_a_T0_gff", y = "Heterotroph_ASV", color= "#800000FF",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = " ", ylab = " ")
Het_grazing <- ggscatter(richness_v2, x = "g_gff", y = "Heterotroph_ASV", color= "#800000FF",
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = " ", ylab = " ")
Het_correlates <- Het_z | Het_chl | Het_grazing
ggsave(plot=Het_correlates,"Het_correlates.pdf", device="pdf", width=7, height=6, units = c("in")) 

Mix_z <- ggscatter(richness, x = "depth", y = "Mixotroph_ASV", color="#8A9045FF",
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = " ", ylab = "Mixotroph Richness")
Mix_par <- ggscatter(richness, x = "light", y = "Mixotroph_ASV", color="#8A9045FF",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = " ", ylab = " ")
Mix_n <- ggscatter(richness, x = "NO3", y = "Mixotroph_ASV", color="#8A9045FF",
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = " ", ylab = " ")
Mix_p <- ggscatter(richness, x = "PO4", y = "Mixotroph_ASV", color="#8A9045FF",
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = " ", ylab = " ")
Mix_temp <- ggscatter(richness, x = "Temp", y = "Mixotroph_ASV", color="#8A9045FF",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = " ", ylab = " ")
Mix_chl <- ggscatter(richness_v2, x = "chl_a_T0_gff", y = "Mixotroph_ASV", color="#8A9045FF",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = " ", ylab = " ")
Mix_grazing <- ggscatter(richness_v2, x = "g_gff", y = "Mixotroph_ASV", color="#8A9045FF",
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = " ", ylab = " ")

Mix_correlates <- Mix_z | Mix_chl | Mix_grazing  
ggsave(plot=Mix_correlates,"Mix_correlates.pdf", device="pdf", width=7, height=6, units = c("in")) 

Para_z <- ggscatter(richness, x = "depth", y = "Parasite_ASV", color="#155F83FF",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = " ", ylab = "Parasite Richness")
Para_par <- ggscatter(richness, x = "light", y = "Parasite_ASV", color="#155F83FF",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = " ", ylab = " ")
Para_n <- ggscatter(richness, x = "NO3", y = "Parasite_ASV", color="#155F83FF",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = " ", ylab = " ")
Para_p <- ggscatter(richness, x = "PO4", y = "Parasite_ASV", color="#155F83FF",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = " ", ylab = " ")
Para_temp <- ggscatter(richness, x = "Temp", y = "Parasite_ASV", color="#155F83FF",
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = " ", ylab = " ")
Para_chl <- ggscatter(richness_v2, x = "chl_a_T0_gff", y = "Parasite_ASV", color="#155F83FF",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = " ", ylab = " ")
Para_grazing <- ggscatter(richness_v2, x = "g_gff", y = "Parasite_ASV", color="#155F83FF",
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = " ", ylab = " ")
Para_correlates <- Para_z |Para_chl | Para_grazing
ggsave(plot=Para_correlates,"Para_correlates.pdf", device="pdf", width=7, height=6, units = c("in")) 

Uncl_z <- ggscatter(richness, x = "depth", y = "Unclassified", color="#767676FF",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = "Depth (m)", ylab = "Unclassified Richness")
Uncl_par <- ggscatter(richness, x = "light", y = "Unclassified", color="#767676FF",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Light Level (PAR)", ylab = " ")
Uncl_n <- ggscatter(richness, x = "NO3", y = "Unclassified", color="#767676FF",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = "[NO3]", ylab = " ")
Uncl_p <- ggscatter(richness, x = "PO4", y = "Unclassified", color="#767676FF",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = "[PO4]", ylab = " ")
Uncl_temp <- ggscatter(richness, x = "Temp", y = "Unclassified", color="#767676FF",
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = "Temperature", ylab = " ")
Uncl_correlates <- Uncl_z | Uncl_par |Uncl_n | Uncl_p | Uncl_temp 
ggsave(plot=Uncl_correlates,"Uncl_correlates.pdf", device="pdf", width=7, height=6, units = c("in")) 

correlates <- Het_correlates / Mix_correlates / Para_correlates / Uncl_correlates
correlates2 <- Het_correlates / Mix_correlates / Para_correlates
ggsave(plot=correlates,"all_correlates_basics.pdf", device="pdf", width=16, height=12, units = c("in")) 
ggsave(plot=correlates2,"fig6_correlates.pdf", device="pdf", width=9, height=9, units = c("in")) 

sample_bind_info <- as.data.frame(sample_data(ASV_exports_pruned2))
richness_v2 <- cbind(richness,sample_bind_info$chl_a_T0_gff,sample_bind_info$g_gff)
colnames(richness_v2)[colnames(richness_v2) == 'sample_bind_info$chl_a_T0_gff'] <- 'chl_a_T0_gff'
colnames(richness_v2)[colnames(richness_v2) == 'sample_bind_info$g_gff'] <- 'g_gff'
#### Testing statistical difference in ASV richness across modes and light levels ####
#richness_long
#Two Sample t Test to compare mean richness in each light region 
richness_long$position <- cut(richness_long$light,
                              breaks=c(0.0, 0.10, 0.65),
                              labels=c('1-10%','20-65%'))
summary(richness_long %>% filter(mode == "Heterotroph_ASV", position=="1-10%") %>% .$asvs)
summary(richness_long %>% filter(mode == "Heterotroph_ASV", position=="20-65%") %>% .$asvs)
summary(richness_long %>% filter(mode == "Mixotroph_ASV", position=="1-10%") %>% .$asvs)
summary(richness_long %>% filter(mode == "Mixotroph_ASV", position=="20-65%") %>% .$asvs)
summary(richness_long %>% filter(mode == "Parasite_ASV", position=="1-10%") %>% .$asvs)
summary(richness_long %>% filter(mode == "Parasite_ASV", position=="20-65%") %>% .$asvs)
richness_long %>% filter(mode == "Heterotroph_ASV") %>% ggplot(aes(position, asvs)) +
  geom_boxplot()
richness_long %>% filter(mode == "Mixotroph_ASV") %>% ggplot(aes(position, asvs)) +
  geom_boxplot()
richness_long %>% filter(mode == "Parasite_ASV") %>% ggplot(aes(position, asvs)) +
  geom_boxplot()

# het asv distribution within light regionn
library(gridExtra) 
h1 <- richness_long %>% filter(mode == "Heterotroph_ASV") %>% ggplot(aes(asvs)) +
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ position)

h2 <- richness_long %>% filter(mode == "Heterotroph_ASV") %>%  ggplot(aes(asvs)) +
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ position) +
  scale_x_log10()
grid.arrange(h1, h2, nrow = 2)

# mixo asv distribution within light regionn
library(gridExtra) 
m1 <- richness_long %>% filter(mode == "Mixotroph_ASV") %>% ggplot(aes(asvs)) +
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ position)

m2 <- richness_long %>% filter(mode == "Mixotroph_ASV") %>%  ggplot(aes(asvs)) +
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ position) +
  scale_x_log10()
grid.arrange(m1, m2, nrow = 2)

# para asv distribution within light regionn
library(gridExtra) 
p1 <- richness_long %>% filter(mode == "Parasite_ASV") %>% ggplot(aes(asvs)) +
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ position)

p2 <- richness_long %>% filter(mode == "Parasite_ASV") %>%  ggplot(aes(asvs)) +
  geom_histogram(fill = "white", color = "grey30") +
  facet_wrap(~ position) +
  scale_x_log10()
grid.arrange(p1, p2, nrow = 2)

het_df <- richness_long %>% filter(mode == "Heterotroph_ASV")
t.test(asvs ~ position, data = het_df)

t.test(log(asvs) ~ position, data = het_df)

wilcox.test(asvs ~ position, data = het_df)

mix_df <- richness_long %>% filter(mode == "Mixotroph_ASV")
t.test(asvs ~ position, data = mix_df)

t.test(log(asvs) ~ position, data = mix_df)

wilcox.test(asvs ~ position, data = mix_df)

par_df <- richness_long %>% filter(mode == "Parasite_ASV")
t.test(asvs ~ position, data = par_df)

t.test(log(asvs) ~ position, data = par_df)

wilcox.test(asvs ~ position, data = par_df)

aov.df.low <- richness_long %>% filter(position == "1-10%")
low.aov <- aov(asvs ~ mode, data = aov.df.low)
summary(low.aov)

aov.df.high <- richness_long %>% filter(position == "20-65%")
high.aov <- aov(asvs ~ mode, data = aov.df.high)
summary(high.aov)

total.aov <- aov(asvs ~ position*mode, data = richness_long)
summary(total.aov)

sub_rich <- subset(richness_long, !mode=="Unclassified")
total.aov <- aov(asvs ~ mode, data = sub_rich)
summary(total.aov)

tukey.one.way<-TukeyHSD(total.aov)
tukey.one.way
tukey.one.way$mode[,"p adj"]
plot(tukey.one.way)

#### Ordination File Prep ####

updated_tax_tab <- as.matrix(read.table("ASV_exports_pruned2.tsv",header=T, row.names=1, check.names=F, sep="\t"))
modes_tax_tab <- tax_table(updated_tax_tab)
updated_count_tab <- read.table("ASV_exports_pruned2_otus.tsv", header=T, row.names=1,check.names=F, sep="\t")
modes_count_tab <- otu_table(updated_count_tab, taxa_are_rows=T)

ASV_physeq_modes <- phyloseq(modes_count_tab, modes_tax_tab, sample_info)
sample_data(ASV_physeq_modes)$light<-factor(sample_data(ASV_physeq_modes)$light, levels=c("1%","5%","10%","20%","40-65%"))
levels(sample_data(ASV_physeq_modes)$light)

parasites_physeq <- subset_taxa(ASV_physeq_modes, Mode == "Parasite")
sample_data(parasites_physeq)$light_level<-factor(sample_data(parasites_physeq)$light_level, levels=c("1%","5%","10%","20%","40-65%"))
levels(sample_data(parasites_physeq)$light_level)

heterotrophs_physeq <- subset_taxa(ASV_physeq_modes, Mode == "Heterotroph")
sample_data(heterotrophs_physeq)$light_level<-factor(sample_data(heterotrophs_physeq)$light_level, levels=c("1%","5%","10%","20%","40-65%"))
levels(sample_data(heterotrophs_physeq)$light_level)

mixotrophs_physeq <- subset_taxa(ASV_physeq_modes, Mode == "Mixotroph")
sample_data(mixotrophs_physeq)$light_level<-factor(sample_data(mixotrophs_physeq)$light_level, levels=c("1%","5%","10%","20%","40-65%"))
levels(sample_data(mixotrophs_physeq)$light_level)

unclassified_mode_physeq <- subset_taxa(ASV_physeq_modes, Mode == "Unclassified")
sample_data(unclassified_mode_physeq)$light_level<-factor(sample_data(unclassified_mode_physeq)$light_level, levels=c("1%","5%","10%","20%","40-65%"))
levels(sample_data(unclassified_mode_physeq)$light_level)

# mixo_uncl_physeq <- subset_taxa(ASV_physeq_modes, Mode == "Unclassified" |Mode == "Mixotroph")
# meyer_phyto_physeq <- merge_phyloseq(mixo_uncl_physeq, ASV_exports_phyto_pruned)
# write.table(tax_table(meyer_phyto_physeq), file="meyer_phyto_npac_tax_table.tsv", sep="\t", quote=F, col.names=NA)
# write.table(otu_table(meyer_phyto_physeq), file="meyer_phyto_npac_counts_table.tsv", sep="\t", quote=F, col.names=NA)

#### FIGURE 7a Heterotrophic RDA ####
heterotrophs_physeq_binary <- heterotrophs_physeq %>% tax_transform("binary", rank = NA)
het_counts_tab <- otu_table(heterotrophs_physeq_binary)
het_tax_tab <- tax_table(heterotrophs_physeq_binary)
het_species_tmp_vec <- het_tax_tab[,8]
het_asv_tmp_vec <- row.names(het_tax_tab)
het_tax_tab_prep1 <- data.frame("Species"=het_species_tmp_vec, row.names = het_asv_tmp_vec)
het_tax_tab_prep1 <- het_tax_tab_prep1 %>% rownames_to_column("ASV")
het_tax_tab_prep1$Species_ASV <- paste(het_tax_tab_prep1$Species, het_tax_tab_prep1$ASV, sep="_")
het_tax_tab_prep1 %>% column_to_rownames("Species_ASV")
rownames(het_counts_tab) <- as.vector(het_tax_tab_prep1$Species_ASV)
head(het_counts_tab)
het_df_counts_tax <- as.data.frame(het_counts_tab)
het_df2 <- het_df_counts_tax %>% rownames_to_column("Species_ASV")
het_df2.g  <- gather(het_df2, Sample, Presence, -Species_ASV)
sample_info_for_merge <- data.frame("Sample"=row.names(rev18_sample_info_tab), "light"=rev18_sample_info_tab$light,
                                    "cruiseID"=rev18_sample_info_tab$cruiseID,"date"=rev18_sample_info_tab$date, 
                                    "MLD"=rev18_sample_info_tab$MLD, "julian"=rev18_sample_info_tab$julian, 
                                    "depth"=rev18_sample_info_tab$depth, "position"=rev18_sample_info_tab$position,
                                    "light_level"=rev18_sample_info_tab$light_level, "epoch"=rev18_sample_info_tab$epoch, 
                                    "PO4"=rev18_sample_info_tab$PO4,"SE_po4"=rev18_sample_info_tab$SE_po4,
                                    "NO3"=rev18_sample_info_tab$NO3,"SE_no3"=rev18_sample_info_tab$SE_no3,
                                    "chl_big"=rev18_sample_info_tab$chl_a_T0_5um,"chl_tot"=rev18_sample_info_tab$chl_a_T0_gff,
                                    "chl_big_sd"=rev18_sample_info_tab$chl_a_T0_sd_5um,"chl_tot_sd"=rev18_sample_info_tab$chl_a_T0_sd_gff,
                                    "g_big"=rev18_sample_info_tab$g_5um,"g_big_se"=rev18_sample_info_tab$g_se_5um,
                                    "g_tot"=rev18_sample_info_tab$g_gff,"g_tot_se"=rev18_sample_info_tab$g_se_gff,
                                    "temp"=rev18_sample_info_tab$temp, stringsAsFactors=F)
het_npac_final_df <- merge(het_df2.g, sample_info_for_merge)
het_combined_wide <- het_npac_final_df %>%
  select(Sample, Species_ASV, Presence, light, cruiseID,date,MLD,julian,depth,position, light_level, epoch,PO4,NO3,chl_big,chl_tot,g_big,g_tot,temp) %>%
  pivot_wider(names_from = Species_ASV,
              values_from = Presence)
het_dbRDA_df <- column_to_rownames(het_combined_wide, var = "Sample")

#This is my dataframe with env AND species presence info so I am subsetting it to get just the species columns
het_species_only <- het_dbRDA_df[,17:105] 
het_species_only <- as.data.frame(het_species_only)

#transform using hellinger method
het_species_transformed <- decostand(het_species_only, method = "hellinger") 

#subsetting just envi data
het_environment_only <- het_dbRDA_df[,1:16] 

#colinearity
het_environment_numeric_only <- het_environment_only[,-c(2:4,7:9)]

# We can visually look for correlations between variables:
heatmap(abs(cor(het_environment_numeric_only)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(10)), 
        Colv = NA, Rowv = NA)
legend("bottom",title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(10)))

#function for normalizing env data
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#Normalize env data
env_norm <- as.data.frame(lapply(het_environment_only[10:16], min_max_norm))

#Add date and other categorical columns back in 
env_norm$light <- het_environment_only$light
env_norm$cruiseID <- het_environment_only$cruiseID
env_norm$date <- het_environment_only$date
env_norm$MLD <- het_environment_only$MLD
env_norm$julian <- het_environment_only$julian
env_norm$depth <- het_environment_only$depth
env_norm$position <- het_environment_only$position
env_norm$light_level <- het_environment_only$light_level
env_norm$epoch <- het_environment_only$epoch
row.names(env_norm) <- row.names(het_environment_only)

#make distance matrix using jaccard 
het_species_dist_jaccard <- vegdist(het_species_only, method = "jaccard")
#make it a matrix
het_species_dist_jaccard.mat <- as.matrix(het_species_dist_jaccard)

#run the RDA -- envi parameters are listed and then the env dataframe is listed 
het_dbRDA <- rda(het_species_transformed ~ julian +depth + PO4 + NO3 + chl_tot + g_tot + temp, env_norm, dist = "jaccard")
#het_dbRDA <- rda(het_species_transformed ~ julian +depth, env_norm, dist = "jaccard")

summary(het_dbRDA)
# Constrained Proportion: variance of Y explained by X (35%)
# Unconstrained Proportion: unexplained variance in Y (65%)
# The included environmental variables explain 35% of the variation in community composition across sites

# Accumulated constrained eigenvalues
# Importance of components:
#                       RDA1    RDA2
# Eigenvalue            0.08079 0.02217
# Proportion Explained  0.78469 0.21531
# Cumulative Proportion 0.78469 1.00000

# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(het_dbRDA)
# $r.squared
# [1] 0.2811991
# 
# $adj.r.squared
# [1] 0.2530108
screeplot(het_dbRDA)
#Get R2 values
het_dbRDA_R2 <- RsquareAdj(het_dbRDA)$r.squared
# 0.2811991

### Monte Carlo Permutation test (999) to test tb-RDA significance 
#R2.obs <- RsquareAdj (het_dbRDA)$r.squared
## Calculate Variance explained bt randomized env variables
env.rand <- env_norm[sample (1:54),]  # the function "sample" will reshuffle the rows with environmental variables
tbRDA.rand <- rda (het_species_transformed ~ julian +  depth, data=env.rand, dist = "jaccard")
RsquareAdj(tbRDA.rand)$r.squared 
# This value represents the variance explained by two random explanatory variables. 
# (Note that we randomize the rows of the environmental matrix, and not each variable independently; 
# in this way, if environmental variables are correlated with each other, this correlation will be preserved.) 

n.perm <- 999  # set the number of permutations
R2.rand <- replicate (n = n.perm, expr = {
  env.rand <- env_norm[sample (1:54),]
  tbRDA.rand <- rda (het_species_transformed ~ julian +  depth, data=env.rand, dist = "jaccard")
  RsquareAdj (tbRDA.rand)$r.squared
})

# The vector R2.rand contains 999 values of variance explained by random variables. 
# In the next step, we will merge them with the observed R2 (R2.obs), since this is considered to be also 
# the part of null distribution (and this is also the reason why the numbers of permutation are usually 
# ending with 9):
R2 <- c (R2.rand, het_dbRDA_R2)
# Draw the histogram of values and highlight the observed R2 there:
hist (R2, nclass = 100)  # ; argument "nclass" separates the bars into 100 categories (compare with hist without this argument)
abline (v = het_dbRDA_R2, col = 'red')  # red line to indicate where in the histogram is the observed value

hist(R2)
abline (v = het_dbRDA_R2, col = 'red')
# Already now we can see that our observed value (0.35, red vertical line) is far higher than all values 
# generated by randomized env. variables (range (R2.rand) shows that these values are between 0.05 and 0.20, 
# but again, this range will change in another run of the analysis, and will likely increase 
# if the number of permutations increases, since there will be a higher chance to get values more different 
# from the mean).

# To calculate the significance, we need to know what is the number of cases in which the variance explained 
# by random explanatory variables was higher or equal to the observed one (explained by real variables). 
# Remember that R.obs is the part of our null distribution (we added it there), so we always obtain at least one 
# value in such comparison (the R.obs itself).

P <- sum (R2 >= het_dbRDA_R2)/(n.perm + 1)  # 0.001
# You can see that our observed R2 (R2.obs) is far higher than any R2 generated by a random variable, 
# so the calculation of P-value contains only one value equal to or higher than observed R2 itself. 
# The resulting P-value is 0.001, which is the lowest P-value we can get with the Monte Carlo permutation 
# test based on 999 permutations:
1/(999 + 1)  # 0.001

#basic plot
plot(het_dbRDA)

anova(het_dbRDA) #Overall test of the significance of the analysis
anova(het_dbRDA, by = "axis", perm.max = 5000) #test axes for sig
anova(het_dbRDA, by = "terms", permu = 5000) #test for significance of env. variables
# light, julian, temp, and depth sig for heterotrophs

#Getting the scores which will be helpful for ordination plot
scores(het_dbRDA) 
het_scores_dbRDA <- scores(het_dbRDA)
#site scores (like sample scores)
het_site_scores <- data.frame(het_scores_dbRDA$sites)
het_site_scores <- rownames_to_column(het_site_scores, var = "Sample")
#species scores
het_species_scores <- het_scores_dbRDA$species

#This is joining scores DF with env data so that I can get env information for each row
env1 <- rownames_to_column(env_norm, var = "Sample")
env1_updated <- env1 %>% 
  mutate(depth_position = case_when(
    depth < MLD ~ "above",
    depth > MLD ~ "below"
  ))
het_site_scores <- het_site_scores %>%
  left_join(env1_updated, by ="Sample")

library(factoextra)
# Compute K-means clustering
fviz_nbclust(x=het_species_dist_jaccard.mat, kmeans, method = "wss", k.max=5, nboot=999,print.summary = TRUE)
# 2 appears to be the optimal number of clusters
km.res <- kmeans(het_species_dist_jaccard.mat, 2, nstart = 25)
print(km.res)
km.res$centers
cluster_wide <- data.frame("cluster"=km.res$cluster)
cluster <- cluster_wide %>% rownames_to_column("Sample")
het_site_scores_cluster <- het_site_scores %>%
  left_join(cluster, by ="Sample")

het_site_scores <- column_to_rownames(het_site_scores, var = "Sample")

het_species_scores_forname <- as.data.frame(het_scores_dbRDA$species)
#subsetting for only a given number of species since many were clustered around the origin
library(psychTools)
het_species_scores_reordered <- dfOrder(het_species_scores_forname, c(-1,-2),absolute=TRUE,ascending=FALSE)
het_species_scores_select <- het_species_scores_reordered[1:30,]
het_species_scores_select <- het_species_scores_select %>% rownames_to_column("Species")
# Model: rda(formula = het_species_transformed ~ julian + depth + PO4 + NO3 + chl_tot + g_tot + temp, data = env_norm, dist = "jaccard")
# Df Variance       F Pr(>F)    
# julian    1 0.022448  5.0585  0.001 ***
#   depth     1 0.080504 18.1408  0.001 ***
#   PO4       1 0.022527  5.0761  0.001 ***
#   NO3       1 0.007756  1.7477  0.052 .  
# chl_tot   1 0.004729  1.0657  0.352    
# g_tot     1 0.005043  1.1364  0.288    
# temp      1 0.018976  4.2761  0.002 ** 
#   Residual 46 0.204136                   
#write.csv(het_species_scores_reordered, file="heterotroph_asv_rda_scores.csv")

#making another dataframe for env data used for envi vectors of RDA
het_env_forfit <- env_norm %>%
  select(julian,depth)

#the actual env fit for RDA
het_envfit <- envfit(het_dbRDA, het_env_forfit, permutations = 999)
#getting env scores
het_env_scores <- data.frame((het_envfit$vectors)$arrows, (het_envfit$vectors)$r, (het_envfit$vectors)$pvals)
het_env_scores <- rownames_to_column(het_env_scores, var = "Variable")

#This is basic ordination with arrows for species, envi variables and envi labels. It also has a box to highlight the species arrows
het.rda.plot_segment <- ggplot()+
  geom_segment(data = het_species_scores_select, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_segment(data = het_env_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = 0.01, 1.0, label = "Year-Day", fontface =2)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_rect(aes(xmin = -0.25, xmax = 0.3, ymin = -0.25, ymax = 0.25),
            fill = "transparent", color = "grey", linewidth = 1.0) +
  annotate("Text", x = -0.13, y = 0.32, label = "Species", color="black", family="Times", fontface = 2, size=4)
het.rda.plot_segment

#this is the complete RDA with points colored by light_level 
complete_rda_hets <- het.rda.plot_segment +
  geom_point(het_site_scores, mapping = aes(x=RDA1, y=RDA2, fill = depth_position), pch =21, color = "black", size = 3)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")+
  labs(x = "tb-RDA1 [78%]", y = "tb-RDA2 [22%]")+
  scale_fill_viridis(option = "B", discrete=TRUE, begin=0.25, end=0.85, direction=-1)+
  stat_ellipse(data=het_site_scores_cluster, aes(x =RDA1, y = RDA2, group=cluster),
               geom = "path", alpha=0.2, level=0.95)+
  scale_color_manual(values = c("gray27", "gray80"))+
  theme_bw()+
  theme(legend.position = "right",axis.title = element_text(size = 18, face = "bold"), axis.text = element_text(size = 16, face = "bold"))
complete_rda_hets
ggsave(plot=complete_rda_hets,"HETS_dbRDA_plot_complete_no_species.pdf", device="pdf",width=10, height=10, units = c("in")) 

#This is just the species box with manually added species labels
species_rda <- ggplot()+
  geom_segment(data = het_species_scores_select, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text_repel(data = het_species_scores_select, aes(x=RDA1,y=RDA2,label=Species), max.overlaps = Inf,
                  colour="red") +
  labs(x = "tb-RDA1 [78%]", y = "tb-RDA2 [22%]")+
  geom_vline(xintercept=c(-0,0), linetype="dotted")+
  geom_hline(yintercept=c(-0,0), linetype="dotted")+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1), axis.text = element_text(size = 16, face = "bold"),axis.title = element_text(size = 18, face = "bold"))+
  theme_bw()
species_rda
ggsave(plot=species_rda,"HETS_dbRDA_plot_species_April2024.pdf", device="pdf",width=10, height=10, units = c("in")) 

#### FIGURE 7c Parasite RDA ####
# Parasite Community
parasites_physeq_binary <- parasites_physeq %>% tax_transform("binary", rank = NA)
par_counts_tab <- otu_table(parasites_physeq_binary)
par_tax_tab <- tax_table(parasites_physeq_binary)
par_species_tmp_vec <- par_tax_tab[,8]
par_asv_tmp_vec <- row.names(par_tax_tab)
par_tax_tab_prep1 <- data.frame("Species"=par_species_tmp_vec, row.names = par_asv_tmp_vec)
par_tax_tab_prep1 <- par_tax_tab_prep1 %>% rownames_to_column("ASV")
par_tax_tab_prep1$Species_ASV <- paste(par_tax_tab_prep1$Species, par_tax_tab_prep1$ASV, sep="_")
par_tax_tab_prep1 %>% column_to_rownames("Species_ASV")
rownames(par_counts_tab) <- as.vector(par_tax_tab_prep1$Species_ASV)
head(par_counts_tab)
par_df_counts_tax <- as.data.frame(par_counts_tab)
par_df2 <- par_df_counts_tax %>% rownames_to_column("Species_ASV")
par_df2.g  <- gather(par_df2, Sample, Presence, -Species_ASV)
sample_info_for_merge <- data.frame("Sample"=row.names(rev18_sample_info_tab), "light"=rev18_sample_info_tab$light,
                                    "cruiseID"=rev18_sample_info_tab$cruiseID,"date"=rev18_sample_info_tab$date, 
                                    "MLD"=rev18_sample_info_tab$MLD, "julian"=rev18_sample_info_tab$julian, 
                                    "depth"=rev18_sample_info_tab$depth, "position"=rev18_sample_info_tab$position,
                                    "light_level"=rev18_sample_info_tab$light_level, "epoch"=rev18_sample_info_tab$epoch, 
                                    "PO4"=rev18_sample_info_tab$PO4,"SE_po4"=rev18_sample_info_tab$SE_po4,
                                    "NO3"=rev18_sample_info_tab$NO3,"SE_no3"=rev18_sample_info_tab$SE_no3,
                                    "chl_big"=rev18_sample_info_tab$chl_a_T0_5um,"chl_tot"=rev18_sample_info_tab$chl_a_T0_gff,
                                    "chl_big_sd"=rev18_sample_info_tab$chl_a_T0_sd_5um,"chl_tot_sd"=rev18_sample_info_tab$chl_a_T0_sd_gff,
                                    "g_big"=rev18_sample_info_tab$g_5um,"g_big_se"=rev18_sample_info_tab$g_se_5um,
                                    "g_tot"=rev18_sample_info_tab$g_gff,"g_tot_se"=rev18_sample_info_tab$g_se_gff,
                                    "temp"=rev18_sample_info_tab$temp, stringsAsFactors=F)
par_npac_final_df <- merge(par_df2.g, sample_info_for_merge)
par_combined_wide <- par_npac_final_df %>% 
  select(Sample, Species_ASV, Presence, light, cruiseID,date,MLD,julian,depth,position, light_level, epoch,PO4,NO3,chl_big,chl_tot,g_big,g_tot,temp) %>%
  pivot_wider(names_from = Species_ASV,
              values_from = Presence)
par_dbRDA_df <- column_to_rownames(par_combined_wide, var = "Sample")
ncol(par_dbRDA_df) #change the column range to the max here
#This is my dataframe with env AND species presence info so I am subsetting it to get just the species columns
par_species_only <- par_dbRDA_df[,17:386] 
par_species_only <- as.data.frame(par_species_only)

#transform using hellinger method
par_species_transformed <- decostand(par_species_only, method = "hellinger") 

#subsetting just envi data
#par_environment_only <- par_dbRDA_df[,1:9] #change this
# use env_norm from heterotroph section, all the same env data

#make distance matrix using jaccard 
par_species_dist_jaccard <- vegdist(par_species_only, method = "jaccard")
#make it a matrix
par_species_dist_jaccard.mat <- as.matrix(par_species_dist_jaccard)

#run the RDA -- envi parameters are listed and then the env dataframe is listed 
#change this
#par_dbRDA <- rda(par_species_transformed ~ light  + depth + julian + chl_big + chl_tot + g_big +g_tot +temp+ NO3+ PO4, env_norm, dist = "jaccard")
# test diff rda models
#par_dbRDA2 <- rda(par_species_transformed ~ light  + julian + temp+ NO3, env_norm, dist = "jaccard")
# all terms significant
par_dbRDA3 <- rda(par_species_transformed ~ julian + depth +PO4 +NO3 +chl_tot +g_tot +temp, env_norm, dist = "jaccard")
#par_dbRDA4 <- rda(par_species_transformed ~ light  +  julian  +chl_big +chl_tot +g_big +g_tot+temp+ NO3, env_norm, dist = "jaccard")

RsquareAdj(par_dbRDA3)
# $r.squared
# [1] 0.2904978
# 
# $adj.r.squared
# [1] 0.2626742
screeplot(par_dbRDA3)

### Monte Carlo Permutation test (999) to test tb-RDA significance 
R2.obs <- RsquareAdj (par_dbRDA3)$r.squared
## Calculate Variance explained bt randomized env variables
env.rand <- env_norm[sample (1:54),]  # the function "sample" will reshuffle the rows with environmental variables
tbRDA.rand <- rda (par_species_transformed ~ julian + depth, env.rand, dist = "jaccard")
RsquareAdj (tbRDA.rand)$r.squared #0.08699748
n.perm <- 999  # set the number of permutations
R2.rand <- replicate (n = n.perm, expr = {
  env.rand <- env_norm[sample (1:54),]
  tbRDA.rand <- rda (par_species_transformed ~ julian + depth, env.rand, dist = "jaccard")
  RsquareAdj (tbRDA.rand)$r.squared
})
R2 <- c (R2.rand, R2.obs)
hist (R2, nclass = 100)  
abline (v = R2.obs, col = 'red')  # red line to indicate where in the histogram is the observed value
hist(R2)
abline (v = R2.obs, col = 'red')
P <- sum (R2 >= R2.obs)/(n.perm + 1)  # 0.001


#basic plot
plot(par_dbRDA3)
anova(par_dbRDA3) #Overall test of the significance of the analysis
anova(par_dbRDA3, by = "axis", perm.max = 5000) #test axes for sig
anova(par_dbRDA3, by = "terms", permu = 5000) #test for significance of env. variables
# Model: rda(formula = par_species_transformed ~ julian + depth + PO4 + NO3 + chl_tot + g_tot + temp, data = env_norm, dist = "jaccard")
# Df Variance       F Pr(>F)    
# julian    1 0.030260  6.0836  0.001 ***
#   depth     1 0.093926 18.8834  0.001 ***
#   PO4       1 0.030960  6.2243  0.001 ***
#   NO3       1 0.009053  1.8201  0.069 .  
# chl_tot   1 0.006278  1.2623  0.180    
# g_tot     1 0.004967  0.9986  0.373    
# temp      1 0.023245  4.6733  0.001 ***
#   Residual 46 0.228805  
summary(par_dbRDA3)
# Accumulated constrained eigenvalues
# Importance of components:
#                        RDA1    RDA2
# Eigenvalue            0.0960 0.02818
# Proportion Explained  0.7731 0.22693
# Cumulative Proportion 0.7731 1.00000

#Getting the scores which will be helpful for ordination plot
scores(par_dbRDA3) 
par_scores_dbRDA <- scores(par_dbRDA3)
#site scores (like sample scores)
par_site_scores <- data.frame(par_scores_dbRDA$sites)
par_site_scores <- rownames_to_column(par_site_scores, var = "Sample")
#species scores
par_species_scores <- par_scores_dbRDA$species

#This is joining scores DF with env data so that I can get env information for each row
par_env1 <- rownames_to_column(env_norm_updated, var = "Sample")
par_site_scores <- par_site_scores %>%
  left_join(par_env1, by ="Sample")

# Compute K-means clustering
parasite_clust <- fviz_nbclust(x=par_species_dist_jaccard.mat, kmeans, method = "wss", k.max=5, nboot=999,print.summary = TRUE)
# 2 appears to be the optimal number of clusters
km.p <- kmeans(par_species_dist_jaccard.mat, 2, nstart = 25)
print(km.p)
km.p$centers
pcluster_wide <- data.frame("cluster"=km.p$cluster)
pcluster <- pcluster_wide %>% rownames_to_column("Sample")
par_site_scores_cluster <- par_site_scores %>%
  left_join(pcluster, by ="Sample")

par_site_scores <- column_to_rownames(par_site_scores, var = "Sample")

par_species_scores_forname <- as.data.frame(par_scores_dbRDA$species)
#subsetting for only a given number of species since many were clustered around the origin
par_species_scores_reordered <- dfOrder(par_species_scores_forname, c(-1,-2),absolute=TRUE,ascending=FALSE)
par_species_scores_select <- par_species_scores_reordered[1:30,]
#par_species_scores_select <- par_species_scores_forname[1:20,]
par_species_scores_select <- rownames_to_column(par_species_scores_select, var = "Species")

#write.csv(par_species_scores_reordered, file="parasite_asv_rda_scores.csv")

#making another dataframe for env data used for envi vectors of RDA
par_env_forfit <- env_norm %>%
  select(julian, depth)

#the actual env fit for RDA
par_envfit <- envfit(par_dbRDA3, par_env_forfit, permutations = 999)
#getting env scores
par_env_scores <- data.frame((par_envfit$vectors)$arrows, (par_envfit$vectors)$r, (par_envfit$vectors)$pvals)
par_env_scores <- rownames_to_column(par_env_scores, var = "Variable")

#This is basic ordination with arrows for species, envi variables and envi labels. It also has a box to highlight the species arrows
par.rda.plot_segment <- ggplot()+
  geom_segment(data = par_species_scores_select, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_segment(data = par_env_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = 0.25, 1.0, label = "Year-Day", fontface =2)+
  annotate(geom = "text", x = -0.85,y = 0.15, label = "Depth (m)", fontface =2)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_rect(aes(xmin = -0.25, xmax = 0.3, ymin = -0.25, ymax = 0.25),
            fill = "transparent", color = "grey", size = 1.0) +
  annotate("Text", x = -0.13, y = 0.32, label = "Species", color="black", family="Times", fontface = 2, size=4)
par.rda.plot_segment

#this is the complete RDA with points colored by light_level 
par_complete_rda <- par.rda.plot_segment +
  geom_point(par_site_scores, mapping = aes(x=RDA1, y=RDA2, fill = depth_position), pch =21, color = "black", size = 3)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")+
  labs(x = "tb-RDA1 [77%]", y = "tb-RDA2 [23%]")+
  scale_fill_viridis(option = "B", discrete=TRUE, begin=0.25, end=0.85, direction=-1)+
  stat_ellipse(data=par_site_scores_cluster, aes(x =RDA1, y = RDA2, group=cluster),
               geom = "polygon", alpha=0.2, level=0.95)+
  scale_color_manual(values = c("gray27", "gray80"))+
  theme(axis.text = element_text(size = 16, face = "bold"), axis.title = element_text(size = 18, face = "bold"),legend.position = "right")+
  theme_bw()
par_complete_rda
ggsave(plot=par_complete_rda,"PARA_dbRDA_plot_complete_no_species.pdf", device="pdf",width=10, height=10, units = c("in")) 

#This is just the species box with manually added species labels
para_species_rda <- ggplot()+
  geom_segment(data = par_species_scores_select, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text_repel(data = par_species_scores_select, aes(x=RDA1,y=RDA2,label=Species), max.overlaps = Inf,colour="red") +
  labs(x = "tb-RDA1 [77%]", y = "tb-RDA2 [23%]")+
  geom_vline(xintercept=c(-0,0), linetype="dotted")+
  geom_hline(yintercept=c(-0,0), linetype="dotted")+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),axis.text = element_text(size = 16, face = "bold"), axis.title = element_text(size = 18, face = "bold"))+
  theme_bw()
para_species_rda
ggsave(plot=para_species_rda,"PARA_dbRDA_plot_species.pdf", device="pdf",width=10, height=10, units = c("in")) 

#### FIGURE 7b Mixotroph RDA ####
# Mixotroph Community
mixotrophs_physeq_binary <- mixotrophs_physeq %>% tax_transform("binary", rank = NA)
mix_counts_tab <- otu_table(mixotrophs_physeq_binary)
mix_tax_tab <- tax_table(mixotrophs_physeq_binary)
mix_species_tmp_vec <- mix_tax_tab[,8]
mix_asv_tmp_vec <- row.names(mix_tax_tab)
mix_tax_tab_prep1 <- data.frame("Species"=mix_species_tmp_vec, row.names = mix_asv_tmp_vec)
mix_tax_tab_prep1 <- mix_tax_tab_prep1 %>% rownames_to_column("ASV")
mix_tax_tab_prep1$Species_ASV <- paste(mix_tax_tab_prep1$Species, mix_tax_tab_prep1$ASV, sep="_")
mix_tax_tab_prep1 %>% column_to_rownames("Species_ASV")
rownames(mix_counts_tab) <- as.vector(mix_tax_tab_prep1$Species_ASV)
head(mix_counts_tab)
mix_df_counts_tax <- as.data.frame(mix_counts_tab)
mix_df2 <- mix_df_counts_tax %>% rownames_to_column("Species_ASV")
mix_df2.g  <- gather(mix_df2, Sample, Presence, -Species_ASV)
sample_info_for_merge <- data.frame("Sample"=row.names(rev18_sample_info_tab), "light"=rev18_sample_info_tab$light,
                                    "cruiseID"=rev18_sample_info_tab$cruiseID,"date"=rev18_sample_info_tab$date, 
                                    "MLD"=rev18_sample_info_tab$MLD, "julian"=rev18_sample_info_tab$julian, 
                                    "depth"=rev18_sample_info_tab$depth, "position"=rev18_sample_info_tab$position,
                                    "light_level"=rev18_sample_info_tab$light_level, "epoch"=rev18_sample_info_tab$epoch, 
                                    "PO4"=rev18_sample_info_tab$PO4,"SE_po4"=rev18_sample_info_tab$SE_po4,
                                    "NO3"=rev18_sample_info_tab$NO3,"SE_no3"=rev18_sample_info_tab$SE_no3,
                                    "chl_big"=rev18_sample_info_tab$chl_a_T0_5um,"chl_tot"=rev18_sample_info_tab$chl_a_T0_gff,
                                    "chl_big_sd"=rev18_sample_info_tab$chl_a_T0_sd_5um,"chl_tot_sd"=rev18_sample_info_tab$chl_a_T0_sd_gff,
                                    "g_big"=rev18_sample_info_tab$g_5um,"g_big_se"=rev18_sample_info_tab$g_se_5um,
                                    "g_tot"=rev18_sample_info_tab$g_gff,"g_tot_se"=rev18_sample_info_tab$g_se_gff,
                                    "temp"=rev18_sample_info_tab$temp, stringsAsFactors=F)
mix_npac_final_df <- merge(mix_df2.g, sample_info_for_merge)
mix_combined_wide <- mix_npac_final_df %>%
  select(Sample, Species_ASV, Presence, light, cruiseID,date,MLD,julian,depth,position, light_level, epoch,PO4,NO3,chl_big,chl_tot,g_big,g_tot,temp) %>%
  pivot_wider(names_from = Species_ASV,
              values_from = Presence)
mix_dbRDA_df <- column_to_rownames(mix_combined_wide, var = "Sample")

#This is my dataframe with env AND species presence info so I am subsetting it to get just the species columns
mix_species_only <- mix_dbRDA_df[,17:167] #change this
mix_species_only <- as.data.frame(mix_species_only)

#transform using hellinger method
mix_species_transformed <- decostand(mix_species_only, method = "hellinger") 

#subsetting just envi data
#mix_environment_only <- mix_dbRDA_df[,1:9] 
# use env_norm for transformed env data 

#make distance matrix using jaccard 
mix_species_dist_jaccard <- vegdist(mix_species_only, method = "jaccard")
#make it a matrix
mix_species_dist_jaccard.mat <- as.matrix(mix_species_dist_jaccard)

#run the RDA -- envi parameters are listed and then the env dataframe is listed 
#change this too
mix_dbRDA <- rda(mix_species_transformed ~ julian + depth +NO3 +PO4 +temp +chl_tot +g_tot, env1_updated, dist = "jaccard")
summary(mix_dbRDA)
# Accumulated constrained eigenvalues
# Importance of components:
#                       RDA1    RDA2
# Eigenvalue            0.07006 0.02078
# Proportion Explained  0.77123 0.22877
# Cumulative Proportion 0.77123 1.00000

# Partitioning of variance:
#                Inertia Proportion
# Total          0.3240     1.0000
# Constrained    0.1335     0.4121
# Unconstrained  0.1905     0.5879

screeplot(mix_dbRDA)
#Get R2 values
mix_dbRDA_R2 <- RsquareAdj(mix_dbRDA)$r.squared
# 0.2803976

### Monte Carlo Permutation test (999) to test tb-RDA significance 
R2.obs <- RsquareAdj (mix_dbRDA)$r.squared
## Calculate Variance explained bt randomized env variables
env.rand <- env_norm[sample (1:54),]  # the function "sample" will reshuffle the rows with environmental variables
tbRDA.rand <- rda (mix_species_transformed ~ julian + depth, env.rand, dist = "jaccard")
RsquareAdj (tbRDA.rand)$r.squared #0.09776623
n.perm <- 999  # set the number of permutations
R2.rand <- replicate (n = n.perm, expr = {
  env.rand <- env_norm[sample (1:54),]
  tbRDA.rand <- rda (mix_species_transformed ~ julian + depth, env.rand, dist = "jaccard")
  RsquareAdj (tbRDA.rand)$r.squared
})
R2 <- c (R2.rand, R2.obs)
hist (R2, nclass = 100)  
abline (v = R2.obs, col = 'red')  # red line to indicate where in the histogram is the observed value
hist(R2)
abline (v = R2.obs, col = 'red')
P <- sum (R2 >= R2.obs)/(n.perm + 1)  # 0.001

#basic plot
plot(mix_dbRDA)
anova(mix_dbRDA) #Overall test of the significance of the analysis
anova(mix_dbRDA, by = "axis", perm.max = 5000) #test axes for sig
anova(mix_dbRDA, by = "terms", permu = 5000) #test for significance of env. variables
# Model: rda(formula = mix_species_transformed ~ julian + depth + NO3 + PO4 + temp + chl_tot + g_tot, data = env1_updated, dist = "jaccard")
#             Df Variance       F Pr(>F)    
#   julian    1 0.023018  6.0325  0.001 ***
#   depth     1 0.067823 17.7747  0.001 ***
#   NO3       1 0.020985  5.4995  0.001 ***
#   PO4       1 0.007106  1.8622  0.034 *  
#   temp      1 0.019224  5.0382  0.001 ***
#   chl_tot   1 0.005741  1.5045  0.101    
#   g_tot     1 0.004555  1.1936  0.222    
# Residual 46 0.175522 

#Getting the scores which will be helpful for ordination plot
scores(mix_dbRDA) 
mix_scores_dbRDA <- scores(mix_dbRDA)
#site scores (like sample scores)
mix_site_scores <- data.frame(mix_scores_dbRDA$sites)
mix_site_scores <- rownames_to_column(mix_site_scores, var = "Sample")
#species scores
mix_species_scores <- mix_scores_dbRDA$species

#This is joining scores DF with env data so that I can get env information for each row
mix_env1 <- rownames_to_column(env_norm_updated, var = "Sample")
mix_site_scores <- mix_site_scores %>%
  left_join(mix_env1, by ="Sample")

# Compute K-means clustering
mixotroph_clust <- fviz_nbclust(x=mix_species_dist_jaccard.mat, kmeans, method = "wss", k.max=5, nboot=999,print.summary = TRUE)
# 2 appears to be the optimal number of clusters
km.m <- kmeans(mix_species_dist_jaccard.mat, 2, nstart = 25)
print(km.m)
km.m$centers
mcluster_wide <- data.frame("cluster"=km.m$cluster)
mcluster <- mcluster_wide %>% rownames_to_column("Sample")
mix_site_scores_cluster <- mix_site_scores %>%
  left_join(mcluster, by ="Sample")

mix_site_scores <- column_to_rownames(mix_site_scores, var = "Sample")

mix_species_scores_forname <- as.data.frame(mix_scores_dbRDA$species)
#subsetting for only a given number of species since many were clustered around the origin
library(psychTools)
mix_species_scores_reordered <- dfOrder(mix_species_scores_forname, c(-1,-2),absolute=TRUE,ascending=FALSE)
mix_species_scores_select <- mix_species_scores_reordered[1:30,]
mix_species_scores_select <- rownames_to_column(mix_species_scores_select, var = "Species")

write.csv(mix_species_scores_reordered, file="mixotroph_asv_rda_scores.csv")

#making another dataframe for env data used for envi vectors of RDA
mix_env_forfit <- env_norm %>%
  select(julian, depth)

#the actual env fit for RDA
mix_envfit <- envfit(mix_dbRDA, mix_env_forfit, permutations = 999)
#getting env scores
mix_env_scores <- data.frame((mix_envfit$vectors)$arrows, (mix_envfit$vectors)$r, (mix_envfit$vectors)$pvals)
mix_env_scores <- rownames_to_column(mix_env_scores, var = "Variable")


#This is basic ordination with arrows for species, envi variables and envi labels. It also has a box to highlight the species arrows
mix.rda.plot_segment <- ggplot()+
  geom_segment(data = mix_species_scores_select, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_segment(data = mix_env_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = 0.01, 1.0, label = "Year-Day", fontface =2)+
  annotate(geom = "text", x = -0.75,y = 0.25, label = "Depth (m)", fontface =2)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  geom_rect(aes(xmin = -0.25, xmax = 0.3, ymin = -0.25, ymax = 0.25),
            fill = "transparent", color = "grey", size = 1.0) +
  annotate("Text", x = -0.13, y = 0.32, label = "Species", color="black", family="Times", fontface = 2, size=4)
mix.rda.plot_segment

#this is the complete RDA with points colored by light_level 
mix_complete_rda <- mix.rda.plot_segment +
  geom_point(mix_site_scores, mapping = aes(x=RDA1, y=RDA2, fill = depth_position), pch =21, color = "black", size = 4)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")+
  labs(x = "tb-RDA1 [77%]", y = "tb-RDA2 [23%]")+
  scale_fill_viridis(option = "B", discrete=TRUE, begin=0.25, end=0.85, direction=-1)+
  stat_ellipse(data=mix_site_scores_cluster, aes(x =RDA1, y = RDA2, group=cluster),
               geom = "path", alpha=0.2, level=0.95)+
  scale_color_manual(values = c("gray27", "gray80"))+
  theme(legend.position = "right",axis.title = element_text(size = 18, face = "bold"),axis.text = element_text(size = 16, face = "bold"))+
  theme_bw()
mix_complete_rda
ggsave(plot=mix_complete_rda,"MIX_dbRDA_plot_complete_no_species.pdf", device="pdf",width=10, height=10, units = c("in")) 

#This is just the species box with manually added species labels
mix_species_rda <- ggplot()+
  geom_segment(data = mix_species_scores_select, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), color = "grey", arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text_repel(data = mix_species_scores_select, aes(x=RDA1,y=RDA2,label=Species), max.overlaps = Inf,colour="red") +
  labs(x = "tb-RDA1 [77%]", y = "tb-RDA2 [23%]")+
  geom_vline(xintercept=c(-0,0), linetype="dotted")+
  geom_hline(yintercept=c(-0,0), linetype="dotted")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1),axis.title = element_text(size = 18, face = "bold"), axis.text = element_text(size = 16, face = "bold"))
mix_species_rda
ggsave(plot=mix_species_rda,"MIX_dbRDA_plot_species.pdf", device="pdf",width=10, height=10, units = c("in")) 



parasite_fixed <- tax_fix(parasites_physeq)
p3 <- parasite_fixed %>% tax_transform("clr", rank = NA) %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color="light_level", size=4) +
  scale_color_manual(values = wes_palette("Cavalcanti1"))+
  stat_ellipse(aes(lty=position))+
  ggtitle("Parasites")+
  labs(color = "Light Level")+
  theme(legend.position = "none",plot.title = element_text(hjust=1),text=element_text(size=14, color = "black"))

# PerMANOVA - partitioning the euclidean distance matrix by species
ps_parasite_bray <- phyloseq::distance(parasite_fixed, method = "bray")
sampledf <- data.frame(sample_data(parasite_fixed))
adonis2(ps_parasite_bray ~ position, data = sampledf)

heterotroph_fixed <- tax_fix(heterotrophs_physeq)
p1 <- heterotroph_fixed %>% tax_transform("clr", rank = NA) %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color="light_level", size=4) +
  scale_color_manual(values = wes_palette("Cavalcanti1"))+
  stat_ellipse(aes(lty=position))+
  ggtitle("Heterotrophs")+
  labs(color = "Light Level")+
  theme(legend.position = "none",plot.title = element_text(hjust=1),text=element_text(size=14,  color = "black"))

# PerMANOVA - partitioning the euclidean distance matrix by species
ps_heterotroph_bray <- phyloseq::distance(heterotroph_fixed, method = "bray")
sampledf <- data.frame(sample_data(heterotroph_fixed))
adonis2(ps_heterotroph_bray ~ position, data = sampledf)

mixotroph_fixed <- tax_fix(mixotrophs_physeq)
p2 <- mixotroph_fixed %>%  tax_transform("clr", rank = NA) %>%
  #dist_calc(dist="aitchison") %>%
  ord_calc() %>% 
  ord_plot(color="light_level", size=4) +
  scale_color_manual(values = wes_palette("Cavalcanti1"))+
  stat_ellipse(aes(lty=position))+
  ggtitle("Mixotrophs")+
  labs(color = "Light Level")+
  theme(plot.title = element_text(hjust=1),text=element_text(size=14, color = "black"))

# PerMANOVA - partitioning the euclidean distance matrix by species
ps_mixotroph_bray <- phyloseq::distance(mixotroph_fixed, method = "bray")
sampledf <- data.frame(sample_data(mixotroph_fixed))
adonis2(ps_mixotroph_bray ~ position, data = sampledf)

unclass_fixed <- tax_fix(unclassified_mode_physeq)
p4 <- unclass_fixed %>% tax_transform("clr", rank = NA) %>% 
  ord_calc(method = "PCA") %>% 
  ord_plot(color="light_level", size=4) +
  scale_color_manual(values = wes_palette("Cavalcanti1"))+
  stat_ellipse(aes(lty=position))+
  ggtitle("Unclassified Mode")+
  labs(color = "Light Level")+
  theme(legend.position = "none",plot.title = element_text(hjust=1),text=element_text(size=14, color = "black"))
pca_plot <- wrap_plots(p1, p2, p3, p4)
ggsave(plot=pca_plot,"Figure5_PCA_CLR.pdf", device="pdf",width=10, height=10, units = c("in")) 

# PerMANOVA - partitioning the euclidean distance matrix by species
ps_unclass_bray <- phyloseq::distance(unclass_fixed, method = "bray")
sampledf <- data.frame(sample_data(unclass_fixed))
adonis2(ps_unclass_bray ~ position, data = sampledf)

#### Presence/Absence Data Processing #####
# whole_comm_fixed -- phyloseq object with taxa, modes, sample info, and asv counts
whole_comm_fixed_binary <- whole_comm_fixed %>% tax_transform("binary", rank  = NA)
presence_df <- as.data.frame(otu_table(whole_comm_fixed_binary))
presence_df2 <- presence_df %>% rownames_to_column("ASV")
presence_df_long <- gather(presence_df2, key="sample",value="presence",-ASV)
taxa_mode <- as.data.frame(tax_table(whole_comm_fixed_binary))
genus_species_mode <- taxa_mode %>% dplyr::select(Genus, Species, Mode, unique)
genus_species_mode$ASV <- genus_species_mode$unique
test_merge <- merge(presence_df_long,genus_species_mode, by="ASV")
sample_for_merging <- rev18_sample_info_tab %>% rownames_to_column("sample")

presence_full_long <- merge(test_merge, sample_for_merging, by="sample")
#write.table(presence_full_long, file="presence_full_long.tsv", sep="\t", quote=F, col.names=NA)
presence_full_long <- presence_full_long %>% mutate(depth_position = case_when(
  depth < MLD ~ "above",
  depth > MLD ~ "below"
))
presence_full_long <- read.table("presence_full_long.tsv", header=T, row.names=1,check.names=F, sep="\t")

#### FIGURE 4 VENN DIAGRAMS #####
library(ggVennDiagram)
#combine rows with same value and aggregate remaining columns
test_wider <- presence_full_long %>%
  group_by(Mode, sample, depth_position) %>%
  summarise(presence=sum(presence))
sub_wider <- subset(test_wider,!Mode=="Unclassified")

anova_depth <- lm(presence~Mode*depth_position, data=sub_wider)
summary(anova_depth)
anova(anova_depth)

wide_light <- test_wider %>% spread(light_level, presence)
wide_light$`1%`[wide_light$`1%` > 0] <- 1
wide_light$`5%`[wide_light$`5%` > 0] <- 1
wide_light$`10%`[wide_light$`10%` > 0] <- 1
wide_light$`20%`[wide_light$`20%` > 0] <- 1
wide_light$`40-65%`[wide_light$`40-65%` > 0] <- 1
wider <- wide_light %>% column_to_rownames("ASV")
#write.table(wide_light, file="wide_light.tsv", sep="\t", quote=F, col.names=NA)

Het_venn <- wider %>% subset(Mode=="Heterotroph")
Mix_venn <- wider %>% subset(Mode=="Mixotroph")
Para_venn <- wider %>% subset(Mode=="Parasite")
Uncl_venn <- wider %>% subset(Mode=="Unclassified")

het_light_list2 <- list('high' = c(rownames(Het_venn)[Het_venn$`20%` == 1 | Het_venn$`40-65%` == 1]),
                        'low' = c(rownames(Het_venn)[Het_venn$`1%` == 1 |Het_venn$`5%` == 1  | Het_venn$`10%` == 1 ]))

hets2 <- ggVennDiagram(het_light_list2, category.names = c("65-20% PAR", "10-1% PAR"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_distiller(palette = "Reds", direction = 1)+
  scale_color_grey()+
  theme(legend.position = "none")

hets2 <- ggvenn(het_light_list2, c("high","low"), auto_scale = TRUE, fill_color = c("pink","darkred"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  theme(legend.position = "none")

mix_light_list2 <- list('low' = c(rownames(Mix_venn)[Mix_venn$`1%` == 1 |Mix_venn$`5%` == 1  | Mix_venn$`10%` == 1 ]),
                        'high' = c(rownames(Mix_venn)[Mix_venn$`20%` == 1 | Mix_venn$`40-65%` == 1]))

mixo2 <- ggvenn(mix_light_list2, c("high","low"), auto_scale = TRUE, fill_color = c("yellowgreen","darkgreen"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  theme(legend.position = "none")

mixo2 <- ggVennDiagram(mix_light_list2, category.names = c("1-10% PAR","20-65% PAR"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_distiller(palette = "Blues", direction = 1)+
  scale_color_grey()+
  theme(legend.position = "none")

para_light_list2 <- list('low' = c(rownames(Para_venn)[Para_venn$`1%` == 1 |Para_venn$`5%` == 1  | Para_venn$`10%` == 1 ]),
                         'high' = c(rownames(Para_venn)[Para_venn$`20%` == 1 | Para_venn$`40-65%` == 1]))

para2 <- ggvenn(para_light_list2, c("high","low"), auto_scale = TRUE, fill_color = c("lightblue","darkblue"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  theme(legend.position = "none")

parasite2 <- ggVennDiagram(para_light_list2, category.names = c("1-10% PAR","20-65% PAR"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_distiller(palette = "Greens", direction = 1)+
  scale_color_grey()+
  theme(legend.position = "none")

uncl_light_list2 <- list('low' = c(rownames(Uncl_venn)[Uncl_venn$`1%` == 1 |Uncl_venn$`5%` == 1  | Uncl_venn$`10%` == 1 ]),
                         'high' = c(rownames(Uncl_venn)[Uncl_venn$`20%` == 1 | Uncl_venn$`40-65%` == 1]))

uncl2 <- ggvenn(uncl_light_list2, c("high","low"), auto_scale = TRUE, fill_color = c("grey","black"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  theme(legend.position = "none")

unclassified2 <- ggVennDiagram(uncl_light_list2, category.names = c("1-10% PAR","20-65% PAR"))+
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_distiller(palette = "Greys", direction = 1)+
  scale_color_grey()+
  theme(legend.position = "none")

venns_two <- grid.arrange(hets2, mixo2, para2, uncl2, nrow=4, ncol=1)
ggsave(plot=venns_two,"venns_two_levels.pdf", device="pdf",width=8, height=11, units = c("in")) 

#### Frequency of Occurrence ####
Microzoop_frequency_position <- presence_full_long %>% group_by(Genus, Species,ASV, depth_position, Mode) %>% 
  summarise(absolute = sum(presence > 0),Freq = absolute / length(ASV),odds = Freq / (1-Freq))

# write.table(Microzoop_frequency_epochs, file="Microzoop_frequency_epochs.tsv", sep="\t", quote=F, col.names=NA)
# write.table(Microzoop_total_frequency_cruise, file="Microzoop_total_frequency_cruise.tsv", sep="\t", quote=F, col.names=NA)
# write.table(Microzoop_frequency_light_position, file="Microzoop_frequency_light_position.tsv", sep="\t", quote=F, col.names=NA)

# Load Files
odds_test <- read.table("log_odds.tsv", header=T,check.names=F, sep="\t")
odds_test

#### Distribution of log odds ratios of grazer occurrence ####
# Distribution of log odds ratio for heterotrophs and mixotrophs -- negative is low light (10-1%), positive high light (65-20%), zero is all light levels -- but close to zero is also weak preference
odds <- read.table("mean_log_odds_light.tsv", header=T)
odds$mode <- as.factor(odds$mode)

mu <- ddply(odds_test,"Mode", summarise, grp.mean=mean(LOG))
head(mu)
#         Mode   grp.mean
# 1  Heterotroph -0.7349886
# 2    Mixotroph  0.0936381
# 3     Parasite -0.5593703
# 4 Unclassified -0.2899389
# 
median_light <- ddply(odds_test,"Mode", summarise, grp.median=median(LOG))
head(median_light)
#           Mode  grp.median
# 1  Heterotroph -0.60205999
# 2    Mixotroph  0.07463362
# 3     Parasite -0.58733673
# 4 Unclassified -0.57526857

ggplot(odds_test, aes(x=LOG, color=Mode, fill=Mode))+
  geom_histogram(aes(y=..density..),position="dodge", binwidth=1, alpha=0.4)+
  geom_density(alpha=0.6)+
  facet_wrap(~Mode)+
  geom_vline(data=median_light, aes(xintercept=grp.median, color=Mode),linetype="dashed")+
  scale_color_manual(values=c("#800000FF","#155F83FF","#8A9045FF","#767676FF"))+
  scale_fill_manual(values=c("#800000FF", "#155F83FF","#8A9045FF","#767676FF"))+
  labs(x="Log Ratio", y = "Density")+
  theme_classic()+
  theme(text=element_text(size=14, family="Times New Roman", color = "black"),axis.text = element_text(size = 12), legend.text = element_text(size=12),
        axis.title = element_text(size = 12), legend.title = element_text(size=12))

#### Chi Squared Test ####
#ctable2 <- table(odds_test$Mode, odds_test$Light)
chi2 <- chisq.test(ctable2)
chi2

chisq.test(ctable2)$expected
# All      High       Low
# Heterotroph   2.750000  36.62500  49.62500
# Mixotroph     4.665730  62.13904  84.19522
# Parasite     11.432584 152.26124 206.30618
# Unclassified  3.151685  41.97472  56.87360

#### Fisher's Exact Test ####
#ctable2 <- table(odds_test$Mode, odds_test$Light)
fish <- fisher.test(ctable2, hybrid=TRUE,simulate.p.value = TRUE)

#### FIGURE 5 - Proportion of Grazer Occurrence per Zone ####
#of proportion of grazers more likely to occur in High, Low or All light levels
# Define the patterns for each level of "Light"

# Proportion of grazers per zone 
# Light:                  65-20%  10-1%  All
# Heterotrophs (n=89 ASV)   29      56    4
# Mixotrophs (n=151)        79      61    11
# Parasites (n=370)         141     222   7
# Unclassified (n=102)      44      58    0

Light <- c('65-20%','All','10-1%','65-20%','All','10-1%','65-20%','All','10-1%','65-20%','All','10-1%')
Proportion_light <-c(29/89, 4/89, 56/89, 79/151, 11/151,61/151, 141/370, 7/370, 222/370,44/102,0/102,58/102)
Mode <- c('Heterotroph','Heterotroph','Heterotroph','Mixotroph','Mixotroph','Mixotroph','Parasite','Parasite','Parasite','Unclassified','Unclassified','Unclassified')
prop_taxa_par <- data.frame(Light,Proportion_light,Mode)
prop_taxa_par$Light <- factor(prop_taxa_par$Light, levels = c('65-20%','10-1%','All'))

fig5 <- ggplot(prop_taxa_par, aes(x=Mode, y=Proportion_light, fill=Mode, alpha=Light)) +
  geom_bar(position=position_dodge(),stat="identity",) +
  scale_fill_manual(" ", values = c( "Heterotroph"="#800000FF","Mixotroph"="#8A9045FF","Parasite"="#155F83FF", "Unclassified"="#767676FF"))+
  scale_y_continuous(limits = c(0,0.8), expand = c(0, 0)) +
  labs(x="Feeding Mode",y="Proportion of ASVs",fill = "Relative Position")+
  theme_classic()+
  theme(text=element_text(size=14, color = "black"),legend.position = c(1, 1),legend.justification = c("right", "top"), 
        axis.text = element_text(size = 12), legend.text = element_text(size=12),
        axis.title = element_text(size = 12), legend.title = element_text(size=12))
ggsave(plot=fig5,"figure5_propasvs_bymode_light_occur.pdf", device="pdf", width=7, height=8, units = c("in")) 

fig5_alt <- ggplot(prop_taxa_par, aes(x=Light, y=Proportion_light, fill=Mode)) +
  geom_bar(position=position_dodge(),stat="identity",) +
  scale_fill_manual(" ", values = c( "Heterotroph"="#800000FF","Mixotroph"="#8A9045FF","Parasite"="#155F83FF", "Unclassified"="#767676FF"))+
  scale_y_continuous(limits = c(0,0.8), expand = c(0, 0)) +
  labs(x="Depth Region",y="Proportion of ASVs",fill = "Trophic Mode")+
  theme_classic()+
  theme(text=element_text(size=14, color = "black"),legend.position = c(1, 1),legend.justification = c("right", "top"), 
        axis.text = element_text(size = 12), legend.text = element_text(size=12),
        axis.title = element_text(size = 12), legend.title = element_text(size=12))
ggsave(plot=fig5_alt,"figure5_propasvs_bylight_occur_alt.pdf", device="pdf", width=7, height=8, units = c("in")) 

#### SUPP FIGURES S3-5 Heat Maps ####
hets <- subset(Microzoop_frequency_position,Mode=="Heterotroph", select=c(Species,ASV,depth_position,Freq))
mixo <- subset(Microzoop_frequency_position,Mode=="Mixotroph" & !ASV=="ASV_22", select=c(Species,ASV,depth_position,Freq))
para <- subset(Microzoop_frequency_position,Mode=="Parasite", select=c(Species,ASV,depth_position,Freq))
# "#800000FF","#155F83FF","#8A9045FF",
library(reshape)

## HETEROTROPH
# convert light level to factor and order of levels
levels(hets$light_level) <- c("1%", "10%","20%","5%","40-65%")
hets$light_level <- factor(hets$light_level, levels = c("1%", "5%", "10%","20%","40-65%"))
#hets <- hets[order(hets$Freq),]

heat_het <- ggplot(hets, aes(reorder(Species, -Freq), depth_position, fill=Freq)) +  
  geom_tile()+
  scale_y_discrete(expand=c(0, 0))+
  labs(x=" ", y=" ")+
  coord_fixed(ratio = 2)+
  coord_flip()+
  scale_x_discrete(expand=c(0, 0))+
  theme_grey(base_size=12)+
  scale_fill_gradient(high = "#800000FF",low="lightcoral")+
  guides(fill=guide_legend(title="Frequency of\nOccurrence",reverse=T))+
  theme(text=element_text(size=14, color = "black"),
        legend.position="right", legend.direction="vertical",
        legend.title=element_text(colour='black'),
        legend.margin=margin(grid::unit(0, "cm")),
        legend.text=element_text(colour='black', size=7, face="bold"),
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.6, "cm"),
        axis.text.x=element_text(size=10, colour='black',angle = 45, hjust=1),
        axis.text.y=element_text(vjust=0.2, colour='black'),
        axis.title = element_text(colour='black'),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm"),
        plot.title=element_text(colour='black', hjust=0, size=14, face="bold"))
ggsave("figure7a_het_heat.png")
ggsave("figure7a_het_heat.pdf") 

hets_wide <- pivot_wider(hets, names_from = depth_position, values_from = Freq)
hets_wider <- unite(hets_wide, ID, c(Species, ASV))
hets_m <- column_to_rownames(hets_wider, var="ID")
hets_mat <- as.matrix(hets_m)
palette <- colorRampPalette(c("white","#800000FF"))(50)

# 50% filter 38 ASVs
library(heatmaply)
hets_short <- hets_m %>% 
  #mutate(fileSum = select(., 'above':'below') %>% rowSums()) %>%
  #filter(fileSum >= 0.50) %>%
  filter(above >= 0.50 | below >=0.50)
#select(-c(fileSum))
hets_mat <-  as.matrix(hets_short)
het_heat_new <- ggheatmap(hets_mat,dendrogram=c("row"), k_row = 4, colors=palette, xlab="Depth Position", ylab="Species",
                          key.title="Frequency of\nOccurrence",margin = c(100, 90))
ggsave(plot=het_heat_new,"figure_het_heat_updated_APRIL.pdf",device="pdf", width=8.5, height=11, units = c("in")) 

## MIXOTROPH
# convert light level to factor and order of levels
levels(mixo$light_level) <- c("1%", "10%","20%","5%","40-65%")
mixo$light_level <- factor(mixo$light_level, levels = c("1%", "10%","20%","5%","40-65%"))

heat_mix <- ggplot(mixo, aes(reorder(Species, -Freq), depth_position, fill=Freq)) +  
  geom_tile()+
  scale_y_discrete(expand=c(0, 0))+
  labs(x="", y="Light Level (% PAR)")+
  coord_flip()+
  #coord_fixed(ratio = 2)+
  scale_x_discrete(expand=c(0, 0))+
  theme_grey(base_size=12)+
  scale_fill_gradient(high = "#155F83FF",low="lightblue")+
  guides(fill=guide_legend(title="Frequency of\nOccurrence",reverse=T))+
  theme(text=element_text(size=14,  color = "black"),
        legend.position="right", legend.direction="vertical",
        legend.title=element_text(colour='black'),
        legend.margin=margin(grid::unit(0, "cm")),
        legend.text=element_text(colour='black', size=7, face="bold"),
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.6, "cm"),
        axis.text.x=element_text(size=10, colour='black',angle = 45, hjust=1),
        axis.text.y=element_text(vjust=0.2, colour='black'),
        axis.title = element_text(colour='black'),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm"),
        plot.title=element_text(colour='black', hjust=0, size=14, face="bold"))
ggsave("figure7b_mixo_heat.png")
ggsave("figure7b_mixo_heat.pdf") 


mixo_wide <- pivot_wider(mixo, names_from = depth_position, values_from = Freq)
mixo_wider <- unite(mixo_wide, ID, c(Species, ASV))
mixo_m <- column_to_rownames(mixo_wider, var="ID")
mixo_mat <- as.matrix(mixo_m)
palette_m <- colorRampPalette(c("white","darkolivegreen"))(50)
mixo_heat_new_long <- ggheatmap(mixo_mat,dendrogram=c("row"), colors=palette_m, xlab="Light Level (% PAR)", ylab="Species",
                                key.title="Frequency of\nOccurrence",margin = c(100, 90))
# with 50% filter there are 85 ASVs
# with 75% filter there are 60 ASVs
mixo_short <- mixo_m %>% 
  #mutate(fileSum = select(., 'above':'below') %>% rowSums()) %>%
  filter(above >= 0.75 | below >=0.75)
#select(-c(fileSum))
mixo_mat <-  as.matrix(mixo_short)
mixo_heat_new <- ggheatmap(mixo_mat, dendrogram=c("row"), k_row = 3,colors=palette_m, xlab="Depth Position", ylab="Species",
                           key.title="Frequency of\nOccurrence", margins=c(50,50))
ggsave(plot=mixo_heat_new,"figure_mixo_heat_APRIL.pdf",device="pdf", width=8.5, height=11, units = c("in")) 

# convert light level to factor and order of levels
levels(para$light_level) <- c("1%", "10%","20%","5%","40-65%")
para$light_level <- factor(para$light_level, levels = c("1%", "5%", "10%","20%","40-65%"))

# assign text colour
heat_par <- ggplot(para, aes(Species, light_level, fill=Freq)) +  
  geom_tile()+
  scale_y_discrete(expand=c(0, 0))+
  labs(x="", y="Light Level (% PAR)")+
  coord_flip()+
  scale_x_discrete(expand=c(0, 0))+
  theme_grey(base_size=12)+
  scale_fill_gradient(high = "#8A9045FF",low="white")+
  guides(fill=guide_legend(title="Frequency of\nOccurrence",reverse=T))+
  theme(text=element_text(size=14, color = "black"),
        legend.position="right", legend.direction="vertical",
        legend.title=element_text(colour='black'),
        legend.margin=margin(grid::unit(0, "cm")),
        legend.text=element_text(colour='black', size=7, face="bold"),
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.6, "cm"),
        axis.text.x=element_text(size=10, colour='black',angle = 45, hjust=1),
        axis.text.y=element_text(vjust=0.2, colour='black'),
        axis.title = element_text(colour='black'),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm"),
        plot.title=element_text(colour='black', hjust=0, size=14, face="bold"))
ggsave("figure7c_para_heat.png")
ggsave("figure7c_para_heat.pdf") 

para_wide <- pivot_wider(para, names_from = depth_position, values_from = Freq)
para_wider <- unite(para_wide, ID, c(Species, ASV))
para_m <- column_to_rownames(para_wider, var="ID")

para_mat <- as.matrix(para_m)
palette_p <- colorRampPalette(c("white","#155F83FF"))(100)
para_heat_new <- ggheatmap(para_mat,dendrogram=c("row"), colors=palette_p, xlab="Light Level (% PAR)", ylab="Species",
                           key.title="Frequency of\nOccurrence",margin = c(100, 90))

avg_fo_para <- para %>% dplyr::group_by(light_level,Species) %>% summarize(mean(Freq))
para_wide2 <- pivot_wider(avg_fo_para, names_from = light_level, values_from = "mean(Freq)")
para_m2 <- column_to_rownames(para_wide2, var="Species")
para_m2_reordered <- para_m2 %>%  relocate("40-65%",.before="1%") %>% 
  relocate("10%",.before="5%") %>%
  relocate("20%",.before="1%") %>%
  relocate("1%",.after="5%")
para_mat2 <- as.matrix(para_m2_reordered)
palette_p <- colorRampPalette(c("white","#155F83FF"))(100)
para_heat_new2 <- ggheatmap(para_mat2, dendrogram=c("row"),k_row=3, colors=palette_p, xlab="Light Level (% PAR)", ylab="Species",
                            key.title="Frequency of\nOccurrence",margin = c(100, 90))
# 95% filter 39 ASVs
# 90% filter 55 ASVs
# 85% filter 63 ASVs
# 80% filter 82 ASVs
# 75% filter 93 ASVs
# 50% filter 148 ASVs
para_short <- para_m %>% 
  #mutate(fileSum = select(., 'above':'below') %>% rowSums()) %>%
  #filter(fileSum >= 0.95) %>%
  filter(above >= 0.90 | below >=0.90)
#select(-c(fileSum))
para_mat <-  as.matrix(para_short)
para_heat_new <- ggheatmap(para_mat, dendrogram=c("row"), colors=palette_p, xlab="Depth Position", ylab="Species",
                           key.title="Frequency of\nOccurrence",margin = c(100, 90))
ggsave(plot=para_heat_new,"figure_para_heat_APRIL.pdf",device="pdf", width=8.5, height=11, units = c("in")) 



heat_het / heat_mix / heat_par

para_wide <- pivot_wider(para, names_from= c(Species,light_level), values_from = Freq)
