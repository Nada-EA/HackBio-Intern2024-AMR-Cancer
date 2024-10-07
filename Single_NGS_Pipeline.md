**Step 1: Make directories (mkdir)**

mkdir raw\_data QC\_reports aligned bam\_files bcf\_vcf

   
**Step 2: Download samples and reference (wget)**

cd raw\_data

wget \-O ERR8774458\_1.fastq.gz [https://zenodo.org/records/10426436/files/ERR8774458\_1.fastq.gz?download=1](https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz?download=1)  
wget \-O ERR8774458\_2.fastq.gz [https://zenodo.org/records/10426436/files/ERR8774458\_2.fastq.gz?download=1](https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz?download=1)  
wget \-O ../reference/reference.fasta [https://zenodo.org/records/10886725/files/Reference.fasta?download=1](https://zenodo.org/records/10886725/files/Reference.fasta?download=1)

**Step 3: Quality Control (fastqc)**

fastqc ERR8774458\_1.fastq.gz ERR8774458\_2.fastq.gz \-o ../QC\_reports/

**Step 4: Trimming of low-quality bases and adapter sequences (fastp)**

fastp \-i ERR8774458\_1.fastq.gz \-I ERR8774458\_2.fastq.gz \\  
\-o ../trimmed/trimmed\_1.fastq.gz \-O ../trimmed/trimmed\_2.fastq.gz

**Step 5: Genome Mapping (bwa)**

\#\#Indexing reference genome  
bwa index ../reference/reference.fasta

\#\#Align sample reads to reference genome  
bwa mem ../reference/reference.fasta \\   
../trimmed/trimmed\_1.fastq.gz ../trimmed/trimmed\_2.fastq.gz \> ../aligned/maped\_reads.sam  
**Step 6: Convert SAM to BAM & sort the BAM file (samtools)**

\#\#convert SAM to BAM  
samtools view \-bS ../aligned/maped\_reads.sam \\   
\-o ../bam\_files/sample.bam

\#\#Check the BAM file  
samtools view ../bam\_files/sample.bam | head		\#Not sorted

\#\#sort BAM file  
samtools sort ../bam\_files/sample.bam \\   
\-o ../bam\_files/sorted\_sample.bam

\#\#View the sorted BAM file  
samtools tview ../bam\_files/sorted\_sample.bam \\   
../reference/reference.fasta	  \#\#example: variants at POS 38 & 133

**Step 7: Variant Calling (bcftools)**

\#\#Indexing sample sorted BAM file  
samtools index ../bam\_files/sorted\_sample.bam

\#\#Alignment pileup  
bcftools mpileup \-f ../reference/reference.fasta \\   
../bam\_files/sorted\_sample.bam \> ../bcf\_vcf/sample.bcf

\#\#Variant calling  
bcftools call \--ploidy 1 \-m \-v \-o ../bcf\_vcf/sample\_variants.vcf ../bcf\_vcf/sample.bcf

