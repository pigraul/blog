
- [切分bam](#切分bam)
- [bulk slamseq流程](#bulk-slamseq流程)
- [Starsolo Statistics Explanation](#starsolo-statistics-explanation)
- [使用命令行的输出作为一个变量在script里使用](#使用命令行的输出作为一个变量在script里使用)
- [删除临时文件](#删除临时文件)
- [不识别python文件](#不识别python文件)
- [空channel的传递](#空channel的传递)
- [一个process的多个输出自定义](#一个process的多个输出自定义)
- [Optional input](#optional-input)
- [nextflow不输出到最终路径](#nextflow不输出到最终路径)
- [scanpy docker权限问题](#scanpy-docker权限问题)

<br>

## 切分bam
切分bam，并行计算，减少时间

 https://nf-co.re/subworkflows/bam_split_by_region
 
切分等级考虑（目前选择4）：

	1. 按chr：单个文件大，运行时间过长
 
	2. 按gene：文件过多，io问题，后续合并文件数目限制问题
 
	3. 按block：将gene分成多个block切分；需要设置block原则，份数或者大小，不同基因组考虑较多；index时设置
 
	4. 按reads number且，目前选择1M reads一个文件（参数可调整，大文件调大，时间换空间），大多在300-500个文件


## Bulk slamseq流程
https://github.com/nf-core/slamseq


</br>

## Starsolo Statistics Explanation
https://github.com/alexdobin/STAR/issues/1887


</br>

## 使用命令行的输出作为一个变量在script里使用

https://stackoverflow.com/questions/66568781/how-to-call-a-variable-created-in-the-script-in-nextflow


</br>

## 删除临时文件
`process.scratch = true`

https://github.com/nextflow-io/nextflow/issues/165


</br>

## 不识别python文件
- 文件需要可执行
- 文件编码，换行文件需要LF格式
![7f2571f755deca7c60359adfa72ba478](https://github.com/user-attachments/assets/d5921500-a745-4dcf-94ea-b6019cbec6b4)

</br>

## 空channel的传递
https://nextflow-io.github.io/patterns/process-when-empty/

</br>

## 一个process的多个输出自定义
publishDir使用数组定义
```
publishDir = [
                    [
                        path: { "${params.outdir}/${params.trimmer}/fastqc" },
                        mode: params.publish_dir_mode,
                        pattern: "*.{html,zip}"
                    ],
                    [
                        path: { "${params.outdir}/${params.trimmer}" },
                        mode: params.publish_dir_mode,
                        pattern: "*.fq.gz",
                        enabled: params.save_trimmed
                    ],
                    [
                        path: { "${params.outdir}/${params.trimmer}" },
                        mode: params.publish_dir_mode,
                        pattern: "*.txt"
                    ]
]
```
https://github.com/nf-core/rnaseq/blob/3bec2331cac2b5ff88a1dc71a21fab6529b57a0f/conf/modules.config#L237-L254

</br>

## Optional input

https://nextflow-io.github.io/patterns/optional-input/

</br>

## nextflow不输出到最终路径
```
    withName: GATK4_BEDTOINTERVALLIST {
        publishDir  = [ enabled: false ]
    }
```

</br>

## scanpy docker权限问题

遇到numba等包的权限问题，可以设置临时路径

```
process LABEL {
    container "raulee/sgr-scanpy"
    containerOptions '--env HOME=/tmp'
}
```

