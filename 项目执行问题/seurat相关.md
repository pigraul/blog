# 所有细胞表达量为0时FeaturePlot的bug

- 问题描述：如果某个gene在所有细胞中表达量为0，应该所有点的颜色为灰色（假设色度设置：最小值->灰色，最大值->红色），但是有的时候会显示全部为红色，尤其是所有gene在一起展示的时候，会出现歧义，如图CA6基因
 ![10410b49980cac0afab2991c325499f9](https://github.com/user-attachments/assets/cbf42a20-1d0a-4426-8207-14260cafd1d9)

- 解决方法，添加`keep.scale = "all"`参数

  ```
  p <- FeaturePlot(PRO, features=marker,cols = c("lightgrey","red"),
        min.cutoff=0,max.cutoff="q90",pt.size =0.1,
        reduction="umap",slot = "data",raster=FALSE,repel =TRUE,ncol=3, keep.scale = "all")
  ```
  参考 https://www.biostars.org/p/9467039/
