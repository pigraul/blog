
---

# 主要优化内容

- 增加断点执行

- 基因组相关配置写入配置文件

- 模块化

- 自动计算基因组大小


# 新增参数

- `--startPoint` 开始执行的步骤，默认0。步骤：

    0:create

    1:basic analysis

    2:Dimensionality Reduction

    3:scRNA-seq based annotation

    4:Identify peaks and plot marker peak

    5:Motif analysis

    6:Co-accessibility and Peak2GeneLinkage


- `archrdir` 当步骤不为0时，需提供ArchRProject所在路径，该路径下需包含`Save-ArchR-Project.rds`文件；
    如文件已被删除，可以去0.rds/路径下查找相应步骤的rds文件，拷贝至指定路径下并重命名。

- `genomeSize` 有效基因组大小，若不提供，流程根据自动genome_annotation的chromSizes自动计算。

- `extraLib` 额外的库文件路径，即BSgenome/TxDb/OrgDb的存放路径，建议新增物种放在统一路径下方便管理。

- `fdr` 如果从第5步开始执行，需提供步骤4中的FDR最大值，可以在.o文件中找到，如`[1] "FDR <= 0.65 & Log2FC >= 1"`


# 参考脚本

```
Rscript "./archr_singlesample_V1.R" \
    --fragment  ./fragments_corrected_dedup_count.tsv.gz \
    --spname GF_LWF_43W_1 \
    --threads 1 \
    --prefix GF_LWF_43W_1 \
    --species chicken \
    --outdir ./  \
    --resolution 0.8 \
    --dim 8 \
    --scRNA no \
    --testMethod wilcoxon \
    --peakcells 75 \
    --startPoint 4 --archrdir ./
```
