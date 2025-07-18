```warning
genome大于4G，谨慎！
http://github.com/Bioconductor/BSgenomeForge/issues/45
```

</br>

# BSgenome

- 查询是否已有BSgenome

  网页查询：https://bioconductor.org/packages/release/BiocViews.html#___BSgenome 或在R中代码查询

  ```
  BiocManager::available("BSgenome")
  ```

- 通过[BSgenomeForge](https://github.com/Bioconductor/BSgenomeForge)自建BSgenome
  
  - 对于NCBI或者UCSC来源的基因组数据，直接使用`forgeBSgenomeDataPkgFromNCBI()`或`forgeBSgenomeDataPkgFromUCSC()`即可，使用方法参考BSgenomeForge的[quick introduction](https://bioconductor.org/packages/release/bioc/vignettes/BSgenomeForge/inst/doc/QuickBSgenomeForge.html)

  - 对于Ensembl来源的基因组数据，通过[Ensembl release note](https://ftp.ensembl.org/pub/release-114/species_EnsemblVertebrates.txt)查找Ensembl与NCBI的对应关系，然后使用`forgeBSgenomeDataPkgFromNCBI()`

  - 对于其他来源的基因组数据，根据[Advanced BSgenomeForge usage](https://bioconductor.org/packages/release/bioc/vignettes/BSgenomeForge/inst/doc/AdvancedBSgenomeForge.pdf)来构建。如果序列较多，建议用2bit文件；如果基因组大于4G，建议用fasta来构建。

    ```
    # fasta
    mkdir seqs_srcdir
    faSplit byname Sta.fasta seqs_srcdir/
    cd seqs_srcdir/
    gzip *

    # 2bit
    library("BSgenomeForge")
    fastaTo2bit("./Sta.fasta", "./Sta.2bit", assembly_accession=NA)

    # 配置seed文件后，构建package
    forgeBSgenomeDataPkg("./BSgenome.Stamariscina.UN.stav1-seed")
    ```

    seed文件可以参考：https://github.com/Bioconductor/BSgenome/issues/62

  **我们的系统目前无法建BSgenome, https://github.com/Bioconductor/BSgenome/issues/54 ；通过在gitpot上测试也许是权限问题**

  在gitpod上拉取bioconductor_docker，启动R环境构建

  ```
  docker pull bioconductor/bioconductor_docker:3.21-R-4.5.1 
  docker images 
  # bioconductor/bioconductor_docker   3.21-R-4.5.1   07b54e5b8cf3   45 hours ago    4.65GB
  # 启动R环境，由于权限问题只能在根目录下操作，不能挂载根目录，所以不执行--rm，保留容器用于cp文件
  docker run -it --name my_bioconductor_app bioconductor/bioconductor_docker:3.21-R-4.5.1 R

  ## 在启动的R环境中执行BSgenomeForge的操作

  # 关闭容器，将结果cp出来，本例结果路径为/BSgenome.Stamariscina.NCBI.ASM302478v1
  docker cp bioconductor/bioconductor_docker:/BSgenome.Stamariscina.NCBI.ASM302478v1 /workspace/voting-app/test/
  ```

  最后用R命令构建包

  ```
  R CMD build ./BSgenome.Stamariscina.NCBI.ASM302478v1
  R CMD check BSgenome.Stamariscina.NCBI.ASM302478v1_1.0.0.tar.gz  # check
  R CMD INSTALL BSgenome.Stamariscina.NCBI.ASM302478v1_1.0.0.tar.gz   # install

  # 在R中调用
  > library(BSgenome.Stamariscina.NCBI.ASM302478v1)
  ```

</br>

# TxDb

- 查询是否已有TxDb

  网页查询：https://bioconductor.org/packages/devel/BiocViews.html#___TxDb 或在R中代码查询

  ```
  BiocManager::available("TxDb")
  ```

- 通过[txdbmaker](https://github.com/Bioconductor/txdbmaker)构建

  构建参考[Manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/txdbmaker/inst/doc/txdbmaker.html)，可以直接构建UCSC/Ensembl/Biomart来源的数据，非常规数据库来源的可以使用gtf或gff3文件来构建（要确保符合相应的文件格式要求）

  ```
  > library(txdbmaker)
  > txdb_gtf <- makeTxDbFromGFF(file = "Sta.gtf", format = "gtf", dataSource = "NCBI assembly ASM302478v1", organism = "Selaginella tamariscina")
  > saveDb(txdb_gtf, file="Sta_txdb.sqlite")
  TxDb object:
  # Db type: TxDb
  # Supporting package: GenomicFeatures
  # Data source: NCBI assembly ASM302478v1
  # Organism: Selaginella tamariscina
  # Taxonomy ID: 137178
  # miRBase build ID: NA
  # Genome: NA
  # Nb of transcripts: 27761
  # Db created by: txdbmaker package from Bioconductor
  # Creation time: 2025-07-17 10:13:19 +0800 (Thu, 17 Jul 2025)
  # txdbmaker version at creation time: 1.2.0
  # RSQLite version at creation time: 2.4.1
  # DBSCHEMAVERSION: 1.2
  
  > sta_txdb <- loadDb("Sta_txdb.sqlite")
  ```

</br>

# OrgDb

- 查询是否已有OrgDb

  网页查询：https://bioconductor.org/packages/devel/BiocViews.html#___OrgDb 或在R中代码查询

  ```
  BiocManager::available("OrgDb")
  ```

- 通过[OrganismDbi](https://github.com/Bioconductor/OrganismDbi)构建

  构建参考[Manual](https://bioconductor.org/packages/release/bioc/vignettes/OrganismDbi/inst/doc/OrganismDbi.html)。
  简单的根据TxDb构建使用[makeOrganismDbFromTxDb](https://rdrr.io/bioc/OrganismDbi/man/makeOrganismDbFromTxDb.html)（仅限Bioconductor中已有的OrgDb package）。

  ```
  library(OrganismDbi)
  library(AnnotationHub)
  makeOrganismDbFromTxDb(txdb)
  ```

</br>

# ensembldb and GenomicRanges

当无法制作TxDb或OrgDb时，可以提供三个`GRanges` object: TSS/exons/genes来构建`geneAnnotation`，参考ArchR[官方文档](https://www.archrproject.com/bookdown/getting-set-up.html#)和[issues#112](https://github.com/GreenleafLab/ArchR/issues/112#issuecomment-630380660)。

```
library(ensembldb)
library(GenomicRanges)
db <- ensDbFromGtf(gtf="gene.gtf",organism = "Gallus_gallus", genomeVersion = "GRCg6a", version = 111)
edb <- EnsDb(db)
rm(db)
edb = addFilter(edb, filter = GeneBiotypeFilter('protein_coding',"=="))
gene.ranges <- genes(edb)
tss.ranges <- resize(gene.ranges, 1, "start")
exon.ranges = exons(edb)

# 可以将构建好的object保存供后续使用
saveRDS(tss.ranges,file="tss.rds")
saveRDS(gene.ranges,file="gene.rds")
saveRDS(exon.ranges,file="exon.rds")
```

非Ensembl来源的gtf文件需要满足以下条件：
- transcript部分，需要有`transcript_id`和`gene_id`
- exon部分，需要有`transcript_id`和`exon_number`
- version = "NA"

最简的gtf格式，少于这些信息，会无法构建db
```
scaffold329     maker   transcript      12736   14066   2394    -       .       transcript_id "ST329G00.1"; gene_id "ST329G00";
scaffold329     maker   exon    12736   12996   .       -       .       transcript_id "ST329G00.1"; gene_id "ST329G00"; exon_number "3";
scaffold329     maker   exon    13622   13724   .       -       .       transcript_id "ST329G00.1"; gene_id "ST329G00"; exon_number "2";
scaffold329     maker   exon    13914   14066   .       -       .       transcript_id "ST329G00.1"; gene_id "ST329G00"; exon_number "1";
scaffold329     maker   transcript      15125   15584   266     +       .       transcript_id "ST329G02.1"; gene_id "ST329G02";
scaffold329     maker   exon    15125   15126   .       +       .       transcript_id "ST329G02.1"; gene_id "ST329G02"; exon_number "1";
scaffold329     maker   exon    15205   15306   .       +       .       transcript_id "ST329G02.1"; gene_id "ST329G02"; exon_number "2";
scaffold329     maker   exon    15423   15584   .       +       .       transcript_id "ST329G02.1"; gene_id "ST329G02"; exon_number "3";
```

</br>

# 有效基因组大小[faCount](https://github.com/ENCODE-DCC/kentUtils?tab=readme-ov-file)

通过`faCount genome.fa  -summary > gemone.stat`计算基因组各碱基含量后，用非N长度(len - N)作为有效基因组大小。

```
#seq     len           A            C            G            T            N          cpg
total    1053332251    303321718    221538399    221562717    303525499    3383918    12649701
prcnt    1.0           0.2880       0.2103       0.2103       0.2882       0.0032     0.0120
```



**以上都可以使用集群ucsc_kents环境**