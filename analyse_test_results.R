## Analyse test results
data <- read.table("~/phibd/test_results", as.is=TRUE, header=FALSE)
names(data) <- c("ID1", "ID2", "SNPS", "p", "mean", "chunks", "length", "length.10mb")
data$mean <- ifelse(chunks==0, 0, data$mean)

dd1 <- ifelse(data$ID1<data$ID2, data$ID1, data$ID2)
dd2 <-  ifelse(data$ID1<data$ID2, data$ID2, data$ID1)
data$ID1 <- dd1
data$ID2 <- dd2

rel <- read.table("~/phibd/Mathieson_etal_rel", as.is=TRUE)
rels <- paste0(data$ID1, data$ID2) %in% paste0(rel[,1], rel[,2])

pdf("~/phibd/p_vs_10mblength.pdf")
plot(data$p, data$length.10mb, pch=16, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Pseudo-haploid heterozygosity", ylab="Proportion of genome in IBD chunks >10mb")
legend("topleft", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()

pdf("~/phibd/SNPS_v_length.pdf", width=6, height=6)
plot(data$SNPS, data$length, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Number of SNPs", ylab="Proportion of genome in IBD chunks", pch=16)
points(data$SNPS[rels], data$length[rels], col="#E41A1C80")
legend("topleft", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()
pdf("~/phibd/SNPS_v_10mblength.pdf", width=6, height=6)
plot(data$SNPS, data$length.10mb, col=ifelse(rels, "#E41A1C80" , "#377EB840"), xlab="Number of SNPs", ylab="Proportion of genome in IBD chunks > 10mb", pch=16)
points(data$SNPS[rels], data$length.10mb[rels], col="#E41A1C80")

legend("topleft", c("Related", "Unrelated"), col=c("#E41A1C80" , "#377EB840"), pch=16, bty="n")
dev.off()
