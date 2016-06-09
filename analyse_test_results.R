## Analyse test results
where <- "~/relatives/v13.1/v13.1_phibd_"

data <- read.table(paste0(where, "results.txt"), as.is=TRUE, header=TRUE)
data$mean <- ifelse(data$count==0, 0, data$mean)

## Remove SGDP
sgdp <- grepl("^[ABST]R?_", data$ID1)|grepl("^[ABST]R?_", data$ID2)
data <- data[!sgdp,]
## Remove failed
data <- data[!is.na(data$mean),]

dd1 <- ifelse(data$ID1<data$ID2, data$ID1, data$ID2)
dd2 <-  ifelse(data$ID1<data$ID2, data$ID2, data$ID1)
data$ID1 <- dd1
data$ID2 <- dd2

N_cutoff<-1e5

rels <- FALSE
## rel <- read.table("~/phibd/Mathieson_etal_rel", as.is=TRUE)
## rels <- paste0(data$ID1, data$ID2) %in% paste0(rel[,1], rel[,2])

pdf("~/phibd/p_vs_10mblength.pdf")
plot(data$p, data$IBD_filtered, pch=ifelse(data$N_SNP<N_cutoff, 1, 16), col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Pseudo-haploid heterozygosity", ylab="Proportion of genome in IBD chunks >10mb")
legend("topright", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
abline(h=0.75, col="grey", lty=2)
abline(h=0.5, col="grey", lty=2)
abline(h=0.25, col="grey", lty=2)
abline(h=0.125, col="grey", lty=2)
dev.off()

pdf(paste0(where, "SNPS_v_length.pdf"), width=6, height=6)
plot(data$N_SNP, data$IBD, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Number of SNPs", ylab="Proportion of genome in IBD chunks", pch=ifelse(data$N_SNP<N_cutoff, 1, 16))
points(data$N_SNP[rels], data$IBD[rels], col="#E41A1C80")
legend("topleft", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()
pdf(paste0(where, "SNPS_v_10mblength.pdf"), width=6, height=6)
plot(data$N_SNP, data$IBD_filtered, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Number of SNPs", ylab="Proportion of genome in IBD chunks > 10mb", pch=ifelse(data$N_SNP<N_cutoff, 1, 16))
points(data$N_SNP[rels], data$IBD_filtered[rels], col="#E41A1C80")

legend("topright", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()

pdf(paste0(where, "ids.pdf"), width=6, height=6)
plot(data$p, data$IBD, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Sharing probability", ylab="Proportion of genome in IBD chunks")
abline(v=0.7, lty=2, col="grey")
abline(v=0.87, lty=2, col="grey")
abline(h=0.75, lty=2, col="grey")
abline(v=0.8, lty=2, col="grey")
abline(v=0.775, col="grey")
abline(v=0.825, col="grey")
abline(h=0.5, lty=2, col="grey")
abline(h=0.25, lty=2, col="grey")
abline(h=0.675, col="grey")
abline(h=0.375, col="grey")
dev.off()

## Make list
data$Status <- rep("Unknown", NROW(data))
duplicate <- data$p>0.85 & data$IBD<0.25
parent.child <- data$p>0.80&data$p<0.825&data$IBD<0.25
unrelated <- data$p>0.7 & data$p<=0.775 & data$IBD<0.25
deg3 <-  data$p>0.7 & data$p<=0.775 & data$IBD>0.25 &data$IBD<0.375
deg2 <-  data$p>0.7 & data$p<=0.775 & data$IBD<0.675 &data$IBD>0.375
sib <-  data$p>0.7 & data$p<=0.775 & data$IBD>0.675
data$Status[duplicate] <- "Duplicate"
data$Status[parent.child] <- "Parent_child"
data$Status[unrelated] <- "Unrelated"
data$Status[deg3] <- "Second_degree"
data$Status[deg2] <- "Third_degree"
data$Status[sib] <- "Siblings"

write.table(data, paste0(where, "results_annotated.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
