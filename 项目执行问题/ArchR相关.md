# 软件流程使用参考

[官方教程](https://www.archrproject.com/bookdown/index.html)

[HuBMAP Consortium流程代码](https://github.com/hubmapconsortium/sc-atac-seq-pipeline/blob/pennycuda/archr-steps/bin/run_ArchR_analysis_pt2.R)

</br>

# 使用`loadArchRProject`进行读取再操作

情形：读取别人的ArchRObject(rds文件)，虽然用`readRDS`可以读取，也可以操作，但是涉及到写入的时候，会默认写入到原来（别人）的路径，涉及读写权限问题，即使设置`outputDirectory = "xx"`也不行。

更进一步，用读取后再保存到指定路径的方法，仍然没有改写输出路径，如下图所示，即使已经保存到当前目录，查看输出路径仍然是别人的路径
![image](https://github.com/user-attachments/assets/9a868ecb-6332-40cb-83f0-3517df3e5e3f)

通过`loadArchRProject`读取，则可以解决这个问题

![image](https://github.com/user-attachments/assets/e9ba1674-dd0d-4a53-bd46-73f42e1f4e39)


参考：https://www.archrproject.com/bookdown/saving-and-loading-an-archrproject.html

**流程优化建议：因为每一步执行完会保存一次数据，所以流程可以添加断点续跑参数。**

</br>

# `addGroupCoverages`报错问题解决

关于找不到基因组的报错

```
Error in .requirePackage(genome) :
  Required package : BSgenome.Ggallus.UCSC.galGal6 is not installed/found!
```

如果BSgenome.Ggallus.UCSC.galGal6的包在R默认的路径下，再load一次`library(BSgenome.Ggallus.UCSC.galGal6)`，可以解决问题。([issues2144](https://github.com/GreenleafLab/ArchR/issues/2144))

如果包BSgenome.Ggallus.UCSC.galGal6在用户指定路径下，由于.requirePackage(genome)只会在默认路径下搜索，所以需要将用户路径加入到R的库搜索路径：

```
> x = "BSgenome.Ggallus.UCSC.galGal6"
> x %in% rownames(installed.packages())
[1] FALSE
> .libPaths()
[1] "/SGRNJ01/Public/Software/conda_env/ArchR/lib/R/library"
> .libPaths(c("/Personal/wangxiang/R_lib/", .libPaths()))
> .libPaths()
[1] "/Personal/wangxiang/R_lib"
[2] "/SGRNJ01/Public/Software/conda_env/ArchR/lib/R/library"
> x %in% rownames(installed.packages())
[1] TRUE
```

**流程优化建议：非内置物种的BSgenome包放在统一路径下，代码将该路径加入到默认库路径。**

</br>

# `addReproduciblePeakSet`报错问题解决

addReproduciblePeakSet的报错，主要与其调用的MACS2软件相关。([issues239](https://github.com/GreenleafLab/ArchR/issues/239#issuecomment-656092828))

在排除MACS2环境问题后，报错主要与基因组大小设置相关，非内置物种，需要设置`genomeSize`参数。

目前遇到过两种报错，一种是未设置`genomeSize`参数时，

![image](https://github.com/user-attachments/assets/dda3b948-cddd-438b-a771-9d37b6b4c9cb)

另一种是参数设置不合理时，

![image](https://github.com/user-attachments/assets/d5d3a22b-59ea-40d7-87ef-94aa1a7f3ceb)

查看使用的物种基因组大小

```
> library(BSgenome.Hsapiens.UCSC.hg38)
> seq_info <- seqinfo(Hsapiens)
> total_size <- sum(seqlengths(seq_info))
> total_size
[1] 3257347282
```

但是软件默认的基因组大小要比这个值小

![image](https://github.com/user-attachments/assets/12dd52fe-3e7b-4270-a9cb-afd1a0a8d42e)

关于合理设置genomeSize，在[MACS2的使用](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)中进行了讨论，
需要设置的是 *Effective genome Size*，即是可比对区、低重复区的基因组大小，对于Hard mask的基因组来说，就是非N区的大小。
关于Effective genome Size的定义和常见物种的大小，参见[deepTools](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html)。

```
projHeme4 <- addReproduciblePeakSet(
        ArchRProj = projHeme4,
        groupBy = "Clusters2",
        pathToMacs2 = pathToMacs2,
        minCells = 50, 
        genomeSize = 1.04e+09)
```

**流程优化建议：可以根据deepTools的两种方案计算后参数传入；也可以简单点设置为基因组大小的80%，特别复杂物种（如高重复，gap较多的）可适当降低，如70%。**

</br>

# H5文件读取报错问题解决

在使用`getMarkerFeatures`时报错，

![image](https://github.com/user-attachments/assets/7fcf292b-b6b8-43e7-b5e5-71da9ff98255)

查看代码，多线程会有同时打开或关闭h5文件的风险

https://github.com/GreenleafLab/ArchR/blob/6feec354ad6c8052ddbc4626a2ca2d858ed465bf/R/ArrowUtils.R#L197

设置单线程后，问题解决（其他命令如果遇到h5的问题也可以这样解决）

```
 markersPeaks <- getMarkerFeatures(
        ArchRProj = projHeme4,
        useMatrix = "PeakMatrix",
        groupBy = "Clusters2",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = testmethod,
        threads = 1)
```

**流程优化建议：这个软件运行速度很快，为避免中断，建议全部设为1。**
