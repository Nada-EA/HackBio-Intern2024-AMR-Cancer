<div align="center">

![image](https://github.com/user-attachments/assets/aa8c0e18-f0a4-4b38-b04f-6eaca1b9f004)


   # **NGS Analysis Pipeline** 
   ## **(Sample: ERR8774458)**
   
</div>



## **Step 1: Make directories (mkdir)**

<pre>
mkdir raw\_data QC\_reports aligned bam\_files bcf\_vcf  
</pre>

---

## **Step 2: Download samples and reference (wget)**
<pre>
cd raw\_data

wget \-O ERR8774458\_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458\_1.fastq.gz?download=1
wget \-O ERR8774458\_2.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458\_2.fastq.gz?download=1
wget \-O ../reference/reference.fasta https://zenodo.org/records/10886725/files/Reference.fasta?download=1
</pre>
---

## **Step 3: Quality Control (fastqc)**
<pre>
fastqc ERR8774458\_1.fastq.gz ERR8774458\_2.fastq.gz \-o ../QC\_reports/
</pre>
---

## **Step 4: Trimming of low-quality bases and adapter sequences (fastp)**
<pre>
fastp \-i ERR8774458\_1.fastq.gz \-I ERR8774458\_2.fastq.gz \-o ../trimmed/trimmed\_1.fastq.gz \-O ../trimmed/trimmed\_2.fastq.gz
</pre>
---

## **Step 5: Genome Mapping (bwa)**
<pre>
\#\#Indexing reference genome  
bwa index ../reference/reference.fasta

\#\#Align sample reads to reference genome  
bwa mem ../reference/reference.fasta ../trimmed/trimmed\_1.fastq.gz ../trimmed/trimmed\_2.fastq.gz \> ../aligned/maped\_reads.sam  
</pre>
---

## **Step 6: Convert SAM to BAM & sort the BAM file (samtools)**
<pre>
\#\#convert SAM to BAM  
samtools view \-bS ../aligned/maped\_reads.sam \-o ../bam\_files/sample.bam
</pre>
<pre>
\#\#Check the BAM file  
samtools view ../bam\_files/sample.bam | head		\#Not sorted

\#\#sort BAM file  
samtools sort ../bam\_files/sample.bam \-o ../bam\_files/sorted\_sample.bam

\#\#View the sorted BAM file  
samtools tview ../bam\_files/sorted\_sample.bam ../reference/reference.fasta
</pre>
---

## **Step 7: Variant Calling (bcftools)**
<pre>
\#\#Indexing sample sorted BAM file  
samtools index ../bam\_files/sorted\_sample.bam
</pre>
<pre>
\#\#Alignment pileup  
bcftools mpileup \-f ../reference/reference.fasta ../bam\_files/sorted\_sample.bam \> ../bcf\_vcf/sample.bcf
</pre>
<pre>
\#\#Variant calling  
bcftools call \--ploidy 1 \-m \-v \-o ../bcf\_vcf/sample\_variants.vcf ../bcf\_vcf/sample.bcf
</pre>
---

<div align="center">

  ### **END**
</div>
