library(VariantAnnotation)
library(GenomicFeatures)
library(doParallel)
library(BiocParallel)
library(data.table)

change.nm <- function(temp.mat,te.mat,id.nm,temp.nm,te.nm,te1="[.]"){
    rownames(temp.mat) <- temp.mat[,temp.nm]
    colnames(te.mat) <- gsub(te1,"-",colnames(te.mat))
    over.nm <- intersect(colnames(te.mat)[te.nm],temp.mat[,temp.nm])
    cn <- c(id.nm,over.nm)
    if (class(te.mat)[1] == "data.table"){
        f.mat <- te.mat[,..cn]
    }
    else {
        f.mat <- te.mat[,cn]
    }
    colnames(f.mat)[colnames(f.mat) != id.nm] <- temp.mat[over.nm,"wgs_id"]
    return (f.mat)
}


vcf.names <- c("NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_",".recalibrated_variants.Broad_Rush.vcf.gz")

exp.data <- gsub(" ","",as.matrix(read.table("./5.rsem_result/ROSMAP_exp",sep='\t',header=T)))
medi.exp <- apply(exp.data[,-1],1,function(x)    median(as.double(x)))
exp.data <- exp.data[which(medi.exp > 0.1),]

methyl.exp <- fread("./ROSMAP_arrayMethylation_imputed.tsv")
methyl.info <- gsub(" ","",as.matrix(read.csv("./ROSMAP_arrayMethylation_metaData.tsv",sep='\t',header=T,comment.char="")))
methyl.info <- methyl.info[methyl.info[,"CHR"] != "" & methyl.info[,"MAPINFO"] != "",]
methyl.ran <- GRanges(Rle(methyl.info[,"CHR"]),IRanges(as.integer(methyl.info[,"MAPINFO"]),as.integer(methyl.info[,"MAPINFO"])))

colnames(exp.data) <- gsub("ROSMAP_","",colnames(exp.data))
ROSMAP.cli <- gsub(" ","",as.matrix(read.table("./ROSMAP_clinical.csv",sep=',',header=T)))
ROSMAP.ID <- gsub(" ","",as.matrix(read.table("./ROSMAP_IDkey.csv",sep=',',header=T)))
ROSMAP.mer.cli <- as.matrix(merge(ROSMAP.cli,ROSMAP.ID,by.x="projid",by.y="projid"))
ROSMAP.mer.cli <- ROSMAP.mer.cli[ROSMAP.mer.cli[,"rnaseq_id"] != "" & ROSMAP.mer.cli[,"wgs_id"] != "" & ROSMAP.mer.cli[,"mwas_id"] != "",]
f.cli.mat <- ROSMAP.mer.cli

AD.sam <- f.cli.mat[as.integer(f.cli.mat[,"ceradsc"]) < 4,"wgs_id"]
CN.sam <- f.cli.mat[as.integer(f.cli.mat[,"ceradsc"]) == 4,"wgs_id"]

# change sample ids for expression and methylation into WGS ids.
te.exp.data <- change.nm(f.cli.mat,exp.data,"TXID","rnaseq_id",c(2:637))
te.methyl.exp <- change.nm(f.cli.mat,methyl.exp,"TargetID","mwas_id",c(2:741))

GTFdb <- loadDb("~/txDB_GTF75")
promo.re <- promoters(GTFdb,2000,200)
promo.re <- promo.re[exp.data[,"TXID"],]
promo.mat <- cbind(as.matrix(elementMetadata(promo.re))[,2],as.character(seqnames(promo.re)),start(promo.re),end(promo.re))
over.re <- as.matrix(findOverlaps(promo.re,methyl.ran))
u.over.re <- unique(over.re[,1])

registerDoParallel(cores=12)
paral.re <- foreach(i=1:length(u.over.re),.combine=rbind ,.errorhandling = "pass") %dopar% {
    ea.re <- rbind(over.re[over.re[,1] == u.over.re[i],])
    tx.id <- promo.mat[as.integer(u.over.re[i]),1]
    methyl.id <- methyl.info[as.integer(ea.re[,2]),"TargetID"]

    tx.exp <- rbind(te.exp.data[te.exp.data[,1] == tx.id,])
    ea.pro.re <- promo.mat[promo.mat[,1] == tx.id,]

    params <- ScanVcfParam(which=GRanges(Rle(ea.pro.re[2]),IRanges(as.integer(ea.pro.re[3]),as.integer(ea.pro.re[4]))))
    vcf.file.path <- paste("./6.WGS/",vcf.names[1],ea.pro.re[2],vcf.names[2],sep="")
    geno.data <- rbind(geno(readVcf(TabixFile(vcf.file.path), "hg19", params))$GT)
    te.geno.data <- rbind(geno.data[,is.element(colnames(geno.data),intersect(colnames(te.methyl.exp),colnames(tx.exp)))])
    AD.geno.dt <- rbind(te.geno.data[,is.element(colnames(te.geno.data),AD.sam)])
    maf.f <- apply(AD.geno.dt,1,function(x){
        t.x <- table(unlist(strsplit(x,"/")))
        sum(na.omit(t.x[c("1","2")]))/(length(x[x != "./."]) * 2)
    })
    com.mu <- names(maf.f[which(maf.f > 0.05 & maf.f < 0.95)])
    te.snp.nm <- com.mu[grep("_[A,C,G,T]/[A,C,G,T]$",com.mu)]
    te.geno <- rbind(geno.data[te.snp.nm,])
    rownames(te.geno) <- te.snp.nm
    if (length(te.geno)){
        pre.re <- NULL
        for (j in 1:length(methyl.id)){
            tx.met <- as.matrix(te.methyl.exp[which(te.methyl.exp[,1] == methyl.id[j]),])
            colnames(tx.met) <- c("id","Met")
            for (z in 1:length(te.geno[,1])){
                ea.te.geno <- te.geno[z,]
                ea.te.geno.num <- do.call(rbind,strsplit(ea.te.geno,"/"))
                ea.te.geno.num <- as.integer(ea.te.geno.num[,1]) + as.integer(ea.te.geno.num[,2])
                ea.te.geno.mat <- cbind(ids=names(ea.te.geno),geno=ea.te.geno.num)
                tx.exp.mat <- t(rbind(ids=colnames(tx.exp),geno=tx.exp))
                colnames(tx.exp.mat) <- c("ids","exp")

                tx.exp.geno.cli <- as.matrix(merge(f.cli.mat,ea.te.geno.mat,by.x="wgs_id",by.x="ids"))
                tx.exp.geno.cli <- as.matrix(merge(tx.exp.geno.cli,tx.exp.mat,by.x="wgs_id",by.y="ids"))
                tx.exp.geno.cli <- unique(merge(tx.exp.geno.cli,tx.met,by.x="wgs_id",by.y="id"))
                AD.data <- tx.exp.geno.cli[as.integer(tx.exp.geno.cli[,"ceradsc"]) < 4,]
                AD.data[,"msex"] <- gsub("0","female",gsub("1","male",AD.data[,"msex"]))
                
                merged.data[,"exp"] <- as.double(merged.data[,"exp"])
                merged.data[,"geno"] <- as.integer(merged.data[,"geno"])
                merged.data[,"Met"] <- as.double(merged.data[,"Met"])
                merged.data[,"age_death"] <- as.integer(merged.data[,"age_death"])
                f.lm <- lm(exp ~ geno + Met + geno*Met + msex + age_death, data=merged.data)
                r.lm <- lm(exp ~ geno + Met + msex + age_death, data=merged.data)

                snp.p.v <- anova(snp.reduced.lm,snp.full.lm,test="LRT")$"Pr(>Chi)"[2]
                snp.OR <- summary(r.lm)$coefficient["geno:Met","Estimate"]
                snp.lm <- summary(r.lm)$coefficient["geno","Pr(>|t|)"]
                met.lm <- summary(r.lm)$coefficient["Met","Pr(>|t|)"]

                # Revision codes
                #snp.f.lm <- lm(exp ~ geno + Met + geno*Met + msex + age_death + braaksc + ceradsc + apoe_genotype, data=merged.data)
                #snp.r.lm <- lm(exp ~ geno + Met + msex + age_death + braaksc + ceradsc + apoe_genotype, data=merged.data)
                #snp.p.v <- anova(snp.f.lm,snp.r.lm,test="LRT")$"Pr(>Chi)"[2]

                #snp.f.lm.2 <- lm(exp ~ geno + Met + geno*Met, data=merged.data)
                #snp.r.lm.2 <- lm(exp ~ geno + Met, data=merged.data)
                #snp.p.v.2 <- anova(snp.f.lm.2,snp.r.lm.2,test="LRT")$"Pr(>Chi)"[2]
                #pre.re <- rbind(pre.re,c(tx.id,methyl.id[j],rownames(te.geno)[z],snp.p.v,snp.p.v.2))
                
                pre.re <- rbind(pre.re,c(tx.id,methyl.id[j],rownames(te.geno)[z],snp.p.v,snp.OR,snp.lm,met.lm))
                }
            }
        pre.re
        }
    }
                  
save(paral.re,file="./results/interaction.re")
