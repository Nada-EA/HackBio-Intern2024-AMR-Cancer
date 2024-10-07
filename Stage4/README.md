<div align="center">
  
![image](https://github.com/user-attachments/assets/e55501e8-1448-47df-99a7-19cd510cbd01)
</div>

<div align="center"> 
  
# **Stage 4: NGS Analysis Pipeline**
</div>

![image](https://github.com/user-attachments/assets/ad348f80-7067-4d6e-bc2b-2e973cc4f692)

---

## This stage of the HackBio internship included two tasks. 

In **task one**, we were required to create a simple NGS analysis pipeline for a single sample.
This included downloading the datasets (forward and reverse sample reads and reference genome), applying quality control check, trimming adapter content as well as poor-quality bases, aligning sample reads to the reference genome, and variant calling. 
You can find the step by step process of this analysis in the "Single_NGS_Pipeline.md" file.

---

In **task two**, we were asked to expand the analysis we did to have an automated pipeline that can be applied to five samples.
We created three files:

## **"1. script.sh"**: This is the main pipeline script that automates the entire NGS analysis process.

It performs the following tasks:

1- Downloads samples read files and reference genome file.

2- Performs quality control on the raw data.

3- Trims the reads for adaptor content and poor-quality bases.

4- Maps the reads to the reference genome.

5- Converts the SAM files to BAM files, then sorts the produced BAM files.

6- Performs variant calling and generates VCF files.

## **"2. setup.sh"**: This script installs all the necessary tools required for the pipeline.
These include wget, fastqc, fastp, bwa, samtools, and bcftools. It ensures that your environment is ready for executing the script.sh.

## **"3. requirements.txt"**: 
This file lists all the software tools used in the project.

---

<div align="center">
  
## **End**
</div>
