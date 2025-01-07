# 评估饱和度的一些软件

## RSeQC
`junction_saturation.py -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam -r hg19.refseq.bed12 -o output`

[RSeQC](https://rseqc.sourceforge.net/#usage-information)  基于比对数据重新采样

![image](https://github.com/user-attachments/assets/6a5cb595-2bb7-4d32-9e24-53e345977199)

## Qualimap
`qualimap counts -d condition.txt -s mouse`

[Qualimap](http://qualimap.conesalab.org/doc_html/command_line.html#counts-qc)  基于表达矩阵下采样，可以用于UMI的数据

![image](https://github.com/user-attachments/assets/21ddd1c5-e566-4b12-ba28-dae2f9a6b5ca)


## bbmap
`bbcountunique.sh in=reads.fq out=histogram.txt`

[bbmap](https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/43538-how-to-plot-the-saturation-curves-to-assess-the-sequencing-depth)  基于kmer相似性

![image](https://github.com/user-attachments/assets/1c391598-d987-47f3-af9c-249477473992)
