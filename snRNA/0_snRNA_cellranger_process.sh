#directory structure removed for simplicity. Add home directory on your end

cellranger count --id=CTRL1_snRNA \
                 --transcriptome=/ALS_/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL1_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL2_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL2_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL3_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL3_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL4_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL4_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL5_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL5_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL6_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL5_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=C9noALSnoFTLD_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=C9noALSnoFTLD_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL1_snRNA \
                 --transcriptome=/ALS_/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL1_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL2_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL2_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL3_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL3_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL4_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL4_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL5_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL5_snRNA \
                 --expect-cells=6000 &
wait;

cellranger count --id=CTRL6_snRNA \
                 --transcriptome=/ALS_snRNA/refdata-gex-GRCh38-2020-A \
                 --fastqs=/ALS_snRNA/fastq/ \
                 --include-introns \
                 --sample=CTRL5_snRNA \
                 --expect-cells=6000 &
wait;
