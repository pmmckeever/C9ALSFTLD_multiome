library(Seurat)
set.seed(1234)

##merger
seurC1_RNA<-readRDS(file="objects/seurC1_RNA.RDS")
seurC2_RNA<-readRDS(file="objects/seurC2_RNA.RDS")
seurC3_RNA<-readRDS(file="objects/seurC3_RNA.RDS")
seurC4_RNA<-readRDS(file="objects/seurC4_RNA.RDS")
seurC5_RNA<-readRDS(file="objects/seurC5_RNA.RDS")
seurC6_RNA<-readRDS(file="objects/seurC6_RNA.RDS")
seur1_RNA<-readRDS(file="objects/seur1_RNA.RDS")
seur2_RNA<-readRDS(file="objects/seur2_RNA.RDS")
seur3_RNA<-readRDS(file="objects/seur3_RNA.RDS")
seur4_RNA<-readRDS(file="objects/seur4_RNA.RDS")
seur5_RNA<-readRDS(file="objects/seur5_RNA.RDS")
seur6_RNA<-readRDS(file="objects/seur6_RNA.RDS")
seur7_RNA<-readRDS(file="objects/seur7_RNA.RDS")
seur8_RNA<-readRDS(file="objects/seur8_RNA.RDS")
seur11_RNA<-readRDS(file="objects/seur11_RNA.RDS")
seur12_RNA<-readRDS(file="objects/seur12_RNA.RDS")
seur13_RNA<-readRDS(file="objects/seur13_RNA.RDS")
seur14_RNA<-readRDS(file="objects/seur14_RNA.RDS")
seur15_RNA<-readRDS(file="objects/seur15_RNA.RDS")
seur16_RNA<-readRDS(file="objects/seur16_RNA.RDS")
seur21_RNA<-readRDS(file="objects/seur21_RNA.RDS")
seur22_RNA<-readRDS(file="objects/seur22_RNA.RDS")
seur23_RNA<-readRDS(file="objects/seur23_RNA.RDS")
seur70_RNA<-readRDS(file="objects/seur70_RNA.RDS")

##annotate metadata
seurC1_RNA$sample <- "CTRL1"
seurC2_RNA$sample <- "CTRL2"
seurC3_RNA$sample <- "CTRL3"
seurC4_RNA$sample <- "CTRL4"
seurC5_RNA$sample <- "CTRL5"
seurC6_RNA$sample <- "CTRL6"
seur1_RNA$sample <- "sALSnoFTLD1"
seur2_RNA$sample <- "sALSnoFTLD2"
seur3_RNA$sample <- "sALSnoFTLD3"
seur4_RNA$sample <- "sALSnoFTLD4"
seur5_RNA$sample <- "sALSnoFTLD5"
seur6_RNA$sample <- "sALSnoFTLD6"
seur7_RNA$sample <- "sALSnoFTLD7"
seur8_RNA$sample <- "sALSnoFTLD8"
seur11_RNA$sample <- "C9ALSFTLD1"
seur12_RNA$sample <- "C9ALSFTLD2"
seur13_RNA$sample <- "C9ALSFTLD3"
seur14_RNA$sample <- "C9ALSFTLD4"
seur15_RNA$sample <- "C9ALSFTLD5"
seur16_RNA$sample <- "C9ALSFTLD6"
seur21_RNA$sample <- "C9ALSnoFTLD1"
seur22_RNA$sample <- "C9ALSnoFTLD2"
seur23_RNA$sample <- "C9ALSnoFTLD3"
seur70_RNA$sample <- "C9noALSnoFTLD"

##annotate metadata
seurC1_RNA$diagnosis <- "control"
seurC2_RNA$diagnosis <- "control"
seurC3_RNA$diagnosis <- "control"
seurC4_RNA$diagnosis <- "control"
seurC5_RNA$diagnosis <- "control"
seurC6_RNA$diagnosis <- "control"
seur1_RNA$diagnosis <- "sALSnoFTLD"
seur2_RNA$diagnosis <- "sALSnoFTLD"
seur3_RNA$diagnosis <- "sALSnoFTLD"
seur4_RNA$diagnosis <- "sALSnoFTLD"
seur5_RNA$diagnosis <- "sALSnoFTLD"
seur6_RNA$diagnosis <- "sALSnoFTLD"
seur7_RNA$diagnosis <- "sALSnoFTLD"
seur8_RNA$diagnosis <- "sALSnoFTLD"
seur11_RNA$diagnosis <- "C9ALSFTLD"
seur12_RNA$diagnosis <- "C9ALSFTLD"
seur13_RNA$diagnosis <- "C9ALSFTLD"
seur14_RNA$diagnosis <- "C9ALSFTLD"
seur15_RNA$diagnosis <- "C9ALSFTLD"
seur16_RNA$diagnosis <- "C9ALSFTLD"
#seur17_RNA$diagnosis <- "C9ALSFTLD"
seur21_RNA$diagnosis <- "C9ALSnoFTLD"
seur22_RNA$diagnosis <- "C9ALSnoFTLD"
seur23_RNA$diagnosis <- "C9ALSnoFTLD"
#seur24_RNA$diagnosis <- "C9ALSnoFTLD"
seur70_RNA$diagnosis <- "C9noALSnoFTLD"

seurC1_RNA$sex <- "female"
seurC2_RNA$sex <- "male"
seurC3_RNA$sex <- "female"
seurC4_RNA$sex <- "male"
seurC5_RNA$sex <- "female"
seurC6_RNA$sex <- "female"
seur1_RNA$sex <- "male"
seur2_RNA$sex <- "female"
seur3_RNA$sex <- "male"
seur4_RNA$sex <- "female"
seur5_RNA$sex <- "male"
seur6_RNA$sex <- "female"
seur7_RNA$sex <- "male"
seur8_RNA$sex <- "male"
seur11_RNA$sex <- "male"
seur12_RNA$sex <- "male"
seur13_RNA$sex <- "male"
seur14_RNA$sex <- "female"
seur15_RNA$sex <- "female"
seur16_RNA$sex <- "male"
seur21_RNA$sex <- "female"
seur22_RNA$sex <- "female"
seur23_RNA$sex <- "female"
seur70_RNA$sex <- "male"

seurC1_RNA$chemistry <- "V3"
seurC2_RNA$chemistry <- "V3"
seurC3_RNA$chemistry <- "V3"
seurC4_RNA$chemistry <- "V3"
seurC5_RNA$chemistry <- "V3"
seurC6_RNA$chemistry <- "V3"
seur1_RNA$chemistry <- "V2"
seur2_RNA$chemistry <- "V2"
seur3_RNA$chemistry <- "V3"
seur4_RNA$chemistry <- "V3"
seur5_RNA$chemistry <- "V3"
seur6_RNA$chemistry <- "V2"
seur7_RNA$chemistry <- "V2"
seur8_RNA$chemistry <- "V3"
seur11_RNA$chemistry <- "V3"
seur12_RNA$chemistry <- "V3"
seur13_RNA$chemistry <- "V3"
seur14_RNA$chemistry <- "V3"
seur15_RNA$chemistry <- "V2"
seur16_RNA$chemistry <- "V2"
seur21_RNA$chemistry <- "V3"
seur22_RNA$chemistry <- "V2"
seur23_RNA$chemistry <- "V3"
seur70_RNA$chemistry <- "V3"

seurC1_RNA[["percent.rps"]] <- PercentageFeatureSet(seurC1_RNA, pattern = "^RPS")
seurC2_RNA[["percent.rps"]] <- PercentageFeatureSet(seurC2_RNA, pattern = "^RPS")
seurC3_RNA[["percent.rps"]] <- PercentageFeatureSet(seurC3_RNA, pattern = "^RPS")
seurC4_RNA[["percent.rps"]] <- PercentageFeatureSet(seurC4_RNA, pattern = "^RPS")
seurC5_RNA[["percent.rps"]] <- PercentageFeatureSet(seurC5_RNA, pattern = "^RPS")
seurC6_RNA[["percent.rps"]] <- PercentageFeatureSet(seurC6_RNA, pattern = "^RPS")

seurC1_RNA[["percent.rpl"]] <- PercentageFeatureSet(seurC1_RNA, pattern = "^RPL")
seurC2_RNA[["percent.rpl"]] <- PercentageFeatureSet(seurC2_RNA, pattern = "^RPL")
seurC3_RNA[["percent.rpl"]] <- PercentageFeatureSet(seurC3_RNA, pattern = "^RPL")
seurC4_RNA[["percent.rpl"]] <- PercentageFeatureSet(seurC4_RNA, pattern = "^RPL")
seurC5_RNA[["percent.rpl"]] <- PercentageFeatureSet(seurC5_RNA, pattern = "^RPL")
seurC6_RNA[["percent.rpl"]] <- PercentageFeatureSet(seurC6_RNA, pattern = "^RPL")

seur1_RNA[["percent.rps"]] <- PercentageFeatureSet(seur1_RNA, pattern = "^RPS")
seur2_RNA[["percent.rps"]] <- PercentageFeatureSet(seur2_RNA, pattern = "^RPS")
seur3_RNA[["percent.rps"]] <- PercentageFeatureSet(seur3_RNA, pattern = "^RPS")
seur4_RNA[["percent.rps"]] <- PercentageFeatureSet(seur4_RNA, pattern = "^RPS")
seur5_RNA[["percent.rps"]] <- PercentageFeatureSet(seur5_RNA, pattern = "^RPS")
seur6_RNA[["percent.rps"]] <- PercentageFeatureSet(seur6_RNA, pattern = "^RPS")
seur7_RNA[["percent.rps"]] <- PercentageFeatureSet(seur7_RNA, pattern = "^RPS")
seur8_RNA[["percent.rps"]] <- PercentageFeatureSet(seur8_RNA, pattern = "^RPS")

seur1_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur1_RNA, pattern = "^RPL")
seur2_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur2_RNA, pattern = "^RPL")
seur3_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur3_RNA, pattern = "^RPL")
seur4_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur4_RNA, pattern = "^RPL")
seur5_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur5_RNA, pattern = "^RPL")
seur6_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur6_RNA, pattern = "^RPL")
seur7_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur7_RNA, pattern = "^RPL")
seur8_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur8_RNA, pattern = "^RPL")

seur11_RNA[["percent.rps"]] <- PercentageFeatureSet(seur11_RNA, pattern = "^RPS")
seur12_RNA[["percent.rps"]] <- PercentageFeatureSet(seur12_RNA, pattern = "^RPS")
seur13_RNA[["percent.rps"]] <- PercentageFeatureSet(seur13_RNA, pattern = "^RPS")
seur14_RNA[["percent.rps"]] <- PercentageFeatureSet(seur14_RNA, pattern = "^RPS")
seur15_RNA[["percent.rps"]] <- PercentageFeatureSet(seur15_RNA, pattern = "^RPS")
seur16_RNA[["percent.rps"]] <- PercentageFeatureSet(seur16_RNA, pattern = "^RPS")

seur11_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur11_RNA, pattern = "^RPL")
seur12_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur12_RNA, pattern = "^RPL")
seur13_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur13_RNA, pattern = "^RPL")
seur14_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur14_RNA, pattern = "^RPL")
seur15_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur15_RNA, pattern = "^RPL")
seur16_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur16_RNA, pattern = "^RPL")

seur21_RNA[["percent.rps"]] <- PercentageFeatureSet(seur21_RNA, pattern = "^RPS")
seur22_RNA[["percent.rps"]] <- PercentageFeatureSet(seur22_RNA, pattern = "^RPS")
seur23_RNA[["percent.rps"]] <- PercentageFeatureSet(seur23_RNA, pattern = "^RPS")

seur21_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur21_RNA, pattern = "^RPL")
seur22_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur22_RNA, pattern = "^RPL")
seur23_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur23_RNA, pattern = "^RPL")

#merge prior to integration
merger<-merge(x=seurC1_RNA,
    y=c(seurC2_RNA, seurC3_RNA, seurC4_RNA, seurC5_RNA, seurC6_RNA,
        seur1_RNA, seur2_RNA, seur3_RNA, seur4_RNA, seur5_RNA, seur6_RNA, seur7_RNA, seur8_RNA,
        seur11_RNA, seur12_RNA, seur13_RNA, seur14_RNA, seur15_RNA, seur16_RNA,
        seur21_RNA, seur22_RNA, seur23_RNA, seur70_RNA),
    add.cell.ids=c("CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6",
                    "sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8",
                    "C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD4","C9ALSFTLD5","C9ALSFTLD6",
                    "C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3","C9noALSnoFTLD"), project="C9ALS")
DefaultAssay(merger)<-"RNA"

saveRDS(merger,file="objects/unint_RNA.RDS")