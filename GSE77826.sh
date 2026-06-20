#!/usr/bin/bash


# File: RNA seq analysis
# downloading SRA file (GSE77826)


var=$(tail -n +2 SraRunTable.csv | cut -d ',' -f 1)

for i in ${var}
        do
        if [[ -f ${i}.fastq.gz ]]
        then
                echo "the ${i} is present in the directory"
        else
                echo "downloading sra entry: ${i}"
                # Downloading Sra entry from GEO dataset
                fastq-dump --gzip --defline-qual '+' ${i}
                echo "completed downloading  ${i}"
        fi

        done

# unzip reference file
gunzip Solanum_tuberosum.SolTub_3.0.cds.all.fa.gz
# indexing salmon file
salmon index -t Solanum_tuberosum.SolTub_3.0.cds.all.fa -p 12 -i solanum_index_1 --gencode

## Run Salmon
var=$(tail -n +2 SraRunTable.csv | cut -d ',' -f 1)
for i in ${var}
        do
                SRX=${i}.fastq.gz
                echo ${SRX}

salmon quant -i solanum_index_1 -l A -r ${SRX}  -o ${i}
        done
# pull transcript and gene id from the reference gene

grep -P -o "PGSC0003DMT\d{9}" Solanum_tuberosum.SolTub_3.0.cds.all.fa > transcript.txt
grep -P -o "PGSC0003DMG\d{9}" Solanum_tuberosum.SolTub_3.0.cds.all.fa > gene.txt
paste -d ',' transcript.txt gene.txt > gene_id.txt




















