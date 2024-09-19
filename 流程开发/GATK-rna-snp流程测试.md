# GATK使用相关问题

</br>

## Gatk VariantFiltration 字符报错的问题
`java.lang.NumberFormatException: For input string: "3.89"`

输入的过滤参数需要用浮点型而不是整型，例如

`--filter-name "QD" --filter "QD < 2.0"`

一定要写2.0，如果写2会报错

https://gatk.broadinstitute.org/hc/en-us/community/posts/360072157591-NumberFormatException-Error-in-VariantFiltration

</br>

## STAR two pass run可以分开跑
https://github.com/alexdobin/STAR/issues/733

https://groups.google.com/g/rna-star/c/rBQK-ujtSh8

这样可以接着单细胞的分析步骤跑，不用重新从fastq开始第一次比对

</br>
