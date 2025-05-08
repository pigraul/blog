- model: `onecompartment` or `twocompartment`
  
`onecompartment`是全细胞提取物的数据

`twocompartment`是分别提取细胞核(nucleus)和细胞质的数据(cytosol)，要一一对应（样本，时间点等），在samples_sheet.tsv文件中要先写nucleus部分的数据，再写cytosol的数据，
  参考[demo](https://github.com/IMSBCompBio/Halfpipe/blob/main/docs/mocksamplesheet_singleend_twocompartment.tsv)

</br>

- filterstrategy: `pseudosingleend`, `singleend`, `pairedend`

`pseudosingleend`在分析突变阶段，虽然是PE数据，但是只使用read1的数据，根据文章的建库类型（bulk链特异性文库RF），R1更靠近3UTR区

比对时仍然使用PE比对；其他两个类型字面意思
