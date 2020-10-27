# tonsa multi-generational study


### map to reference genome

```bash

my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=/data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fna.gz
my_samblstr=~/bin/samblaster/samblaster

cd ~/tonsa_genomics/data/trimmed/

for sample in `ls ~/tonsa_genomics/data/trimmed | grep '.fq.gz' | cut -f 1 -d "."| uniq`

do

    echo "starting sample ${sample}"
    #starting with only name of rep. need to pull out files

    rep_1=$(ls ~/tonsa_genomics/data/trimmed | grep ${sample} | grep 'R1')
    rep_2=$(ls ~/tonsa_genomics/data/trimmed | grep ${sample} | grep 'R2')

    echo $rep_1
    echo $rep_2

    rg=$(echo \@RG\\tID:$sample\\tPL:Illumina\\tPU:x\\tLB:$sample\\tSM:$sample)

    $my_bwa mem -t 4 -R $rg $bwagenind $rep_1 $rep_2 | \
    $my_samblstr |\
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o ~/tonsa_genomics/data/aligned/${sample}.bam

done

```


### run popoolation for genetic diversity

```bash

cd ~/tonsa_multigen/analysis/popoolation

for i in $(ls ~/tonsa_genomics/data/aligned/merged | cut -f 1 -d '.'); do

    echo "Status: starting $i"

    samtools mpileup -Q 20 -B --max-depth 2000 --skip-indels \
    -f /data/copepods/tonsa_genome/GCA_900241095.1_Aton1.0_genomic.fa \
    ~/tonsa_genomics/data/aligned/merged/${i}.bam \
    > ~/tonsa_multigen/analysis/popoolation/${i}.mpileup

    echo "Status: $i pileup done; starting popoolation"

    perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --pool-size 50 --min-qual 20 \
    --min-coverage 30 --min-count 2 --max-coverage 1000 --min-covered-fraction 0.5 \
    --input ~/tonsa_multigen/analysis/popoolation/${i}.mpileup \
    --window-size 100 --step-size 100 \
    --output ~/tonsa_multigen/analysis/popoolation/${i}.l2.pi --measure pi

    rm ~/tonsa_multigen/analysis/popoolation/${i}.mpileup

    echo "Status: $i popoolation done, pileup removed"

done

```


### analyze pi

```r

# pi and selection.
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)

#locate the directory containing the files.
dir <- "~/tonsa_multigen/analysis/popoolation"
files <- file.path(dir, list.files(dir))
files <- files[grep("params", files, invert=TRUE)]
files <- files[grep("l2.pi", files, invert=FALSE)]

# read in files
d <- lapply(files, read.table)

# rename without path, etc.
names(d) <- (str_split_fixed(files, "/", 5)[,5] %>% str_split_fixed( "[.]", 3))[,1]

# assign column names
d <- lapply(d, function(x) {
  colnames(x) <- c("gene", "position", "n_variants", "prop_covered", "pi")
  x
})

d <- lapply(d, function(x) {
  x$pi <- as.numeric(as.character(x$pi))
  x
})

d <- lapply(d, function(x) {
  x$snp <- paste(x$gene, x$position, sep="_")
  x
})

for(i in 1:length(d)){
 d[[i]]$gp <- names(d)[i]

}

pops <- names(d)

for(i in 1:length(d)){
 print(sum(!is.na(d[[i]]$pi)))

}

d2 <- vector(mode = "list", length = 16)

for(i in 1:length(d)){
 d2[[i]] <- d[[i]][!is.na(d[[i]]$pi),]
}


###############################################
########
######## reformat for stats
########
################################################


v1 <- vector(mode = "list", length = 16)

for(i in 1:length(d)){
 v1[[i]] <- d2[[i]]$snp
}

overlap <- Reduce(intersect, v1)

# then parse down to only the intersect
d3 <- vector(mode = "list", length = 16)
for(i in 1:length(d2)){
 d3[[i]] <- d2[[i]][d2[[i]]$snp %in% overlap,]
}


#merge entire dist.
mdat <- bind_rows(d3)
nrow(mdat)

mdat$treatment <- substr(mdat$gp,1,2)

forstats <- cbind(mdat$gp, mdat$treatment, mdat$pi)
colnames(forstats) <- c("sample", "treatment", "pi")

write.table(forstats,"~/tonsa_multigen/forstats.txt", quote=FALSE, sep="\t", row.names=F)


###############################################
########
######## run stats
########
################################################


dat <- read.table("~/tonsa_multigen/forstats.txt", header=T)

dat$rep <- substr(dat$sample, 8,11)


out.p <- pairwise.wilcox.test(dat$pi,
                          dat$sample,
                          p.adjust.method="fdr")$p.value

write.table(out.p,"~/tonsa_multigen/pi_stats.txt", quote=FALSE, sep="\t", row.names=T)

dat %>% group_by(treatment) %>%
    summarise(
        n = n(),
        mean = round(mean(pi, na.rm=T),4),
        sd = sd(pi, na.rm=T))

dat %>% group_by(treatment, rep) %>%
    summarise(
        n = n(),
        mean = round(mean(pi, na.rm=T),4),
        sd = sd(pi, na.rm=T))



################
#
# figure for hans
#
################

ylim1 = boxplot.stats(mdat$pi)$stats[c(1, 5)]

pc <- ggplot(data=mdat, aes(x=gp, y=pi, fill=treatment)) +
  geom_boxplot(outlier.shape=NA, show.legend = FALSE) +
 # scale_fill_manual(values=c("gray96", "gray57"))+
  theme_bw() +
  coord_cartesian(ylim = ylim1*1.05)+
  scale_fill_manual(values=c("royalblue3", "forestgreen", "orange2", "red3"))+
  scale_x_discrete(breaks=unique(mdat$gp),
        labels=c("", "", "AA", "",
                 "", "", "AH", "",
                 "", "", "HA", "",
                 "", "", "HH", ""))


ggsave("~/tonsa_multigen/figures/pi_fig.pdf", pc, width = 6, height = 4)


```
