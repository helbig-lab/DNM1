#Load packages ----
library(tidyverse)
library(viridis)
library(Hmisc)

##Transcript Analysis ----
#Read transcript quantities and sample info
#Transcript quantification per sample from bulk RNAseq using Kallisto, doi:10.1038/nbt.3519
#Subsetted for DNM1 Ensembl transcripts and merged across samples
raw_quant <- read_csv("DNM1.kallisto.csv")
samples <- readxl::read_excel("./Documents/DNM1_supplemental.xlsx", skip=1) %>% 
  select(ID, AGE=AGE_AT_COLLECTION) 

#Assign full-length coding transcripts to corresponding exon 10
dnm1b <- c("ENST00000628346.2","ENST00000341179.11","ENST00000372923.7")
dnm1a <- c("ENST00000393594.7","ENST00000486160.3")

#Build comparison tables: exon with age, age-collapsed, exon-collapsed
exon_compare <- raw_quant %>% 
  arrange(pect_est %>% desc) %>% 
  mutate(exon=case_when(target_id %in% dnm1b ~ "DNM1b", 
                        target_id %in% dnm1a ~ "DNM1a", 
                        T ~ "NA")) %>% 
  filter(length>3000) %>% 
  group_by(ID) %>% 
  mutate(total=sum(pect_est)) %>% 
  mutate(pect_est_fig = pect_est/total*100) %>% 
  mutate(pect_est_fig = case_when(exon == "DNM1b" ~ -pect_est_fig, 
                              exon == "DNM1a" ~ pect_est_fig, 
                              T ~ -1000)) %>% 
  filter(pect_est_fig != -1000) %>% 
  ungroup() %>% 
  group_by(ID, exon) %>% 
  dplyr::summarize(pect_est_fig=sum(pect_est_fig)) %>% 
  mutate(pect_est = abs(pect_est_fig)) %>% 
  arrange(-pect_est_fig) %>% 
  left_join(samples) %>% 
  mutate(ID = factor(ID, levels=unique(exon_compare$ID))) 

age_compare <- exon_compare %>% 
  pivot_wider(id_cols=c(ID, AGE), names_from=exon, values_from=pect_est) 

exon_expand <- exon_compare %>% pivot_wider(id_cols=ID, names_from=exon, values_from=pect_est)

#Statistics: Correlation test of transcript expression with age; t-test of transcript ratio
cor.test(age_compare$AGE, age_compare$DNM1a)
cor.test(age_compare$AGE, age_compare$DNM1b)
t.test(exon_expand$DNM1a/exon_expand$DNM1b, mu=1)

#Figure 2A
ggplot(data=exon_compare, aes(y=ID, x = pect_est_fig, fill = exon)) +
  geom_col() +
  theme_classic() + 
  scale_fill_viridis(option="mako", discrete=T, begin=.25) + 
  geom_vline(xintercept=0, size=.1) + 
  theme(axis.text.y = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none") + 
  scale_x_continuous(expand=c(0,0), breaks=c(-40,-20,0,20,40,60), limits=c(-50,70))

#Figure 2B
ggplot(data=age_compare, aes(x=AGE, y=DNM1b)) + 
  geom_smooth(size=1, method="lm", se=F, linetype="dashed", color="gray80") + 
  geom_point(size=3, color="#d1ffb8") + 
  geom_smooth(data=age_compare, method="lm", aes(y=DNM1a), size=1, se=F, linetype="dashed", color="gray80") + 
  geom_point(data=age_compare, aes(y=DNM1a), size=3, color="#0e148d") + 
  theme_classic() + 
  scale_y_continuous(limits=c(0,80), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,25), expand=c(0,0)) + 
  labs(x="Age", y="Exon Expression (%)") + 
  theme(legend.position="none", 
        axis.line=element_blank())

#Figure 2C
ggplot(data=exon_compare, aes(x=exon, y=pect_est)) + 
  geom_violin(aes(fill=exon)) + 
  geom_boxplot(width=.1, aes(color=exon)) + 
  scale_fill_viridis(option="mako", begin=.25, end=1, discrete=T) + 
  scale_color_viridis(option="mako",begin=1,end=.25,discrete=T) + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,80)) + 
  theme(legend.position = "none", 
        text=element_blank(), 
        axis.line=element_blank(), 
        axis.ticks=element_blank())

##Splicing Variant Analysis ----
#Build table with SpliceAI outputs for gnomAD variants
#https://spliceailookup.broadinstitute.org/
splicevars <- tibble(
  variant = c(rep("NM_004408.3:c.589+1G>C",4), 
              rep("NM_004408.3:c.1197-7C>T",4), 
              rep("NM_004408.3:c.1197-7C>A",4), 
              rep("NM_004408.3:c.1197-7C>G",4), 
              rep("NM_001288739.1:c.1335+1G>A",4), 
              rep("NM_004408.3:c.1557+1G>T",4), 
              rep("NM_001288739.1:c.1197-8G>A",4)), 
  type = rep(c("Acceptor Loss", "Donor Loss", 
               "Acceptor Gain", "Donor Gain"), 7),
  score = c(0, .99, 0, .75, 
            .01, 0, 0, 0, 
            .03, 0, .22, 0, 
            .01, 0, .01, 0, 
            0, .88, 0, 0, 
            .83, .93, 0, 0, 
            .41, 0, .93, 0)
  ) %>% 
  mutate(variant=factor(variant, levels=unique(variant)), 
         score = case_when(grepl("Loss", type) ~ -score, T ~ score), 
         type = case_when(grepl("Acceptor", type) ~ "Acceptor", T ~ "Donor"))

#Figure 3A
ggplot(splicevars, aes(x=variant, y=score, fill=type)) + 
  geom_col(position="dodge", width=.5) + 
  theme_classic() + 
  geom_rect(aes(xmin=6.5, xmax=7.5, ymin=-1, ymax=1), fill="white", color="black", linetype="dashed") + 
  scale_fill_manual(values=c("Donor" = "steelblue", "Acceptor" = "red")) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = .2, linetype = "dashed", color = "darkgreen") + 
  geom_hline(yintercept = .5, linetype = "dashed", color = "goldenrod") + 
  geom_hline(yintercept = .8, linetype = "dashed", color = "darkred") + 
  geom_hline(yintercept = -.2, linetype = "dashed", color = "darkgreen") + 
  geom_hline(yintercept = -.5, linetype = "dashed", color = "goldenrod") + 
  geom_hline(yintercept = -.8, linetype = "dashed", color = "darkred") + 
  geom_col(position="dodge", width=.5) + 
  scale_y_continuous(limits=c(-1,1), expand=c(0,0)) + 
  theme(legend.position = "none", text = element_blank(), 
        axis.line.x=element_blank(), axis.ticks=element_blank())

##Population Missense Variant Analysis ----
#Read in gnomAD transcript-specific constraint data and sort transcripts
#Constraint data from gnomAD: doi:10.1038/s41586-020-2308-7, Supplementary Dataset 11 
gnomad <- read_tsv("gnomad.v2.1.1.lof_metrics.by_transcript.txt") %>% 
  filter(gene=="DNM1") %>% 
  mutate(exon = case_when(
    transcript %in% c("ENST00000341179","ENST00000372923") ~ "10b",
    transcript %in% c("ENST00000393594","ENST00000486160","ENST00000475805") ~ "10a"), 
    .before = canonical)

#Figure S2A
ggplot(gnomad, aes(x=transcript, y=oe_mis, fill=exon)) + 
  geom_col() + 
  geom_errorbar(aes(ymin=oe_mis_lower, ymax=oe_mis_upper), width=.2) + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c("10b"="#d1ffb8", "10a"="#0e148d")) + 
  theme(text=element_blank(), axis.ticks=element_blank(), legend.position="none")

#Figure S2B
ggplot(gnomad, aes(x=transcript, y=mis_z, fill=exon)) + 
  geom_col() + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c("10b"="#d1ffb8", "10a"="#0e148d")) + 
  theme(text=element_blank(), axis.ticks=element_blank(), legend.position="none")

#Read in gnomAD DNM1 missense variants
#Exported from gnomAD browser: https://gnomad.broadinstitute.org/gene/ENSG00000106976
misvar <- read_csv("dnm1_gnomad.csv") %>% 
  select(transcript=Transcript, cvar=`Transcript Consequence`, pvar=`Protein Consequence`, 
         maf=`Allele Frequency`) %>% 
  mutate(pos = as.numeric(str_extract(cvar, "\\d+"))) %>% 
  mutate(isoform = case_when(
    (pos > 1336 | pos < 1196) ~ "10a_10b", 
    transcript == "ENST00000372923.3" ~ "10b", 
    transcript == "ENST00000393594.3" ~ "10a"),
    maf = log10(maf)) %>% 
  select(-transcript) %>% separate_rows(isoform) 

#Figure S2C
ggplot(misvar %>% filter(isoform=="10a"), aes(x=pos, y=maf)) + 
  geom_vline(xintercept=1196, linetype="dashed") + 
  geom_vline(xintercept=1335, linetype="dashed") + 
  geom_step(color="#0e148d") + 
  theme_classic() + 
  scale_y_continuous(limits=c(-6,-1),expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,2592), expand=c(0,0)) + 
  theme(text=element_blank(), axis.ticks=element_blank()) 
ggplot(misvar %>% filter(isoform=="10b"), aes(x=pos, y=maf)) + 
  geom_vline(xintercept=1196, linetype="dashed") + 
  geom_vline(xintercept=1335, linetype="dashed") + 
  geom_step(color="#d1ffb8") + 
  theme_classic() + 
  scale_y_continuous(limits=c(-6,-1),expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,2592), expand=c(0,0)) + 
  theme(text=element_blank(), axis.ticks=element_blank()) 

#Subset gnomAD missense variants to exon 10
misvar_10 <- misvar %>% 
  bind_rows(tibble(cvar=c("fakelb","fakeub"), pvar=c("fakelb","fakeub"),
                   maf=-5.400531, pos=c(1197,1334), isoform="10a")) %>% 
  bind_rows(tibble(cvar=c("fakelb","fakeub"), pvar=c("fakelb","fakeub"),
                   maf=-5.400531, pos=c(1197,1334), isoform="10b"))

#Figure S2D
ggplot(misvar_10 %>% filter(pos %in% (1197:1334), isoform=="10a"), aes(x=pos, y=maf)) + 
  geom_step(color="#0e148d", size=2) + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0), limits=c(-5.5,-4)) + 
  scale_x_continuous(limits=c(1197,1334),expand=c(0,0)) + 
  theme(text=element_blank(), axis.ticks=element_blank(), 
        axis.line = element_line(size=2))  
ggplot(misvar_10 %>% filter(pos %in% (1197:1334), isoform=="10b"), aes(x=pos, y=maf)) + 
  geom_step(color="#d1ffb8", size=2) + 
  theme_classic() + 
  scale_y_continuous(expand=c(0,0), limits=c(-5.5,-4)) + 
  scale_x_continuous(limits=c(1197,1334),expand=c(0,0)) + 
  theme(text=element_blank(), axis.ticks=element_blank(), 
        axis.line = element_line(size=2)) 

##Vesicle Diameter Plot ----
#Build summary table for vesicle measurements
vesicle <- tibble(source = c(rep("Th",2),rep("SpC",2),rep("GP",2),rep("Cb",2),rep("Skin",2)),
                  sample = rep(c("Control","Affected"),5), 
                  n = c(34,258,43,53,24,69,94,102,373,272), 
                  mean = c(28.35,44.78,33.34,41.07,27.42,43.88,31.31,42.38,28.43,51.17), 
                  sem = c(.76,.59,.87,1.06,.93,1.01,.69,.86,.29,.45)) %>% 
  unite("unid", c(source, sample), remove=F) %>% 
  mutate(unid=factor(unid,levels=unid), 
         sample=factor(sample, levels=c("Control","Affected"))) %>% 
  mutate(sourcelab=case_when(sample=="Affected" ~ "", T ~ paste("     ",source)))

#Figure 6N
ggplot(vesicle, aes(x=unid, y=mean, color=sample)) + 
  geom_rect(xmin=2.5, xmax=4.5, ymin=26, ymax=52, color=NA, fill="gray80") + 
  geom_rect(xmin=6.5, xmax=8.5, ymin=26, ymax=52, color=NA, fill="gray80") + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem), width=.5) + 
  theme_classic() + 
  scale_x_discrete(labels=vesicle$sourcelab) + 
  scale_y_continuous(limits=c(26,53), expand=c(0,0)) + 
  labs(y = "Vesicle Diameter (mm)", color = "Sample") + 
  scale_color_manual(values=c("Control"="steelblue","Affected"="red")) + 
  guides(color = guide_legend(title.hjust=0.5)) + 
  theme(axis.title.x = element_blank(), 
        axis.title.y=element_text(size=10), 
        legend.position="none")

