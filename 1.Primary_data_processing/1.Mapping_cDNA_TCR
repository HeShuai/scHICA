###-------------------------cDNA mapping---------------------------------------###
####-------human cDNA mapping--------###
#BSUB -n THREAD
#BSUB -o %J.our
#BSUB -e %J.err
#BSUB -R span[hosts=1]
##BSUB -q smp
#######BSUB -m s001

#--force-cells=20000 \

cellranger6 count \
--force-cells=20000 \
--id=SAMPLENAME  \
--localcores=THREAD \
--transcriptome=/data/home/heshuai/reference_data/refdata-cellranger6.0-GRCh38_18Y_HBV \
--fastqs=ABSPATH \
--sample=SAMPLENAME

##---submit to cluster
for i in $(ls /data4/heshuai/RAW_data/1-SingleCell/24.Others/1.Rawdata|grep fastq.gz|grep -E R2|grep -E "cDNA"|cut -d '_' -f1|sort -u|grep -E "F17-cDNA")
do
ABSPATH="/data4/heshuai/RAW_data/1-SingleCell/24.Others/1.Rawdata"
THREAD=16
SAMPLENAME=$(echo $i)
sed "s#ABSPATH#$ABSPATH#g;s#SAMPLENAME#$SAMPLENAME#g;s#THREAD#$THREAD#g" /data4/heshuai/RAW_data/1-SingleCell/3-HCA/HCA2.0/1.Adult_20210729/4.Script/1.cDNA_mapping.lsf|bsub
done


###-------------------------TCR mapping---------------------------------------###
####-------human VDJ mapping--------###
#BSUB -n THREAD
#BSUB -o %J.our
#BSUB -e %J.err
#BSUB -R span[hosts=1]
##BSUB -q smp
#######BSUB -m s001

cellranger6 vdj \
--id=SAMPLENAME \
--localcores=THREAD  \
--reference=/data/home/heshuai/reference_data/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--fastqs=ABSPATH \
--sample=SAMPLENAME

##---submit to cluster
for i in $(ls /data4/heshuai/RAW_data/1-SingleCell/3-HCA/HCA2.0/4.20211101_Fetal3/1.Rawdata|grep fastq.gz|grep -E R2|grep -E "TCR|BCR"|cut -d '_' -f1-2|sort -u)
do
ABSPATH="/data4/heshuai/RAW_data/1-SingleCell/3-HCA/HCA2.0/4.20211101_Fetal3/1.Rawdata"
THREAD=8
SAMPLENAME=$(echo $i)
sed "s#ABSPATH#$ABSPATH#g;s#SAMPLENAME#$SAMPLENAME#g;s#THREAD#$THREAD#g" /data4/heshuai/RAW_data/1-SingleCell/3-HCA/HCA2.0/1.Adult_20210729/4.Script/2.TCR_BCR_mapping.lsf|bsub
done
