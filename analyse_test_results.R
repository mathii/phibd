## Analyse test results
data <- read.table("~/phibd/test_results", as.is=TRUE, header=TRUE)
data$mean <- ifelse(data$count==0, 0, data$mean)

dd1 <- ifelse(data$ID1<data$ID2, data$ID1, data$ID2)
dd2 <-  ifelse(data$ID1<data$ID2, data$ID2, data$ID1)
data$ID1 <- dd1
data$ID2 <- dd2

rel <- read.table("~/phibd/Mathieson_etal_rel", as.is=TRUE)

rels <- paste0(data$ID1, data$ID2) %in% paste0(rel[,1], rel[,2])

pdf("~/phibd/p_vs_10mblength.pdf")
plot(data$p, data$IBD_filtered, pch=16, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Pseudo-haploid heterozygosity", ylab="Proportion of genome in IBD chunks >10mb")
legend("topright", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
abline(h=0.75, col="grey", lty=2)
abline(h=0.5, col="grey", lty=2)
abline(h=0.25, col="grey", lty=2)
abline(h=0.125, col="grey", lty=2)
dev.off()

pdf("~/phibd/SNPS_v_length.pdf", width=6, height=6)
plot(data$N_SNP, data$IBD, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Number of SNPs", ylab="Proportion of genome in IBD chunks", pch=16)
points(data$N_SNP[rels], data$IBD[rels], col="#E41A1C80")
legend("topleft", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()
pdf("~/phibd/SNPS_v_10mblength.pdf", width=6, height=6)
plot(data$N_SNP, data$IBD_filtered, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Number of SNPs", ylab="Proportion of genome in IBD chunks > 10mb", pch=16)
points(data$N_SNP[rels], data$IBD_filtered[rels], col="#E41A1C80")

legend("topright", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()
