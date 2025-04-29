# 图形释义

边的颜色/权重、节点的颜色/大小/形状的解释：在所有的可视化图表中，边的颜色与作为发送方的源相对应，并且边的权重与相互作用的强度成正比。更粗的边线表示更强的信号。

在层次图和圆形图中，圆形的大小与每个细胞组中的细胞数量成正比。在层次图中，实心圆和空心圆分别代表源和目标。

在弦图中，内部较细的条形的颜色表示接收来自相应外部条形信号的目标。内部条形的大小与目标接收到的信号强度成正比。这样的内部条形有助于解读复杂的弦图。**需要注意的是，对于某些细胞组，存在一些没有任何弦连接的内部条形，请忽略它们，因为这是一个尚未被 circlize 包解决的问题。**

Explnations of edge color/weight, node color/size/shape: In all visualization plots, edge colors are consistent with the sources as sender, and edge weights are proportional to the interaction strength. Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group. In the hierarchy plot, solid and open circles represent source and target, respectively. In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note that there exist some inner bars without any chord for some cell groups, please just igore it because this is an issue that has not been addressed by circlize package.


</br>

# cellcat不同版本相关问题

## cellchat的v1和v2不兼容
-  v1无法读取v2生成的rds
-  v1版本中有些slot在v2中不存在或者变更

##  命令使用有所不同
主要是在一些作图中，`from / to`变为`sources.use / targets.use`

比如`netVisual_bubble`，在v1中为`netVisual_bubble(cellchat, from = ligands, to = recepters)`，在v2中为`netVisual_bubble(cellchat, sources.use = ligands, targets.use = recepters)`

</br>


# chord弦图相关问题
默认情况下，每个细胞类型的弦的宽度表示信号强度，宽度自适应，有的老师想使每个细胞类型的弦的宽度一致

我们默认画弦图使用`netVisual_aggregate(cellchat, signaling = "TGFb", layout = "chord"`，这个命令无法修改弦的宽度；cellchat提供了另外一个画弦图的命令`netVisual_chord_cell(cellchat, signaling = "TGFb", scale = TRUE)`，将`scale = TRUE`则可以画相同宽度，`scale = FALSE`结果与`netVisual_aggregate`一致

注意一旦scale了，通过图形，只能比较细胞内部的强度，而无法比较细胞间的强度。
![image](https://github.com/user-attachments/assets/0b00c180-ac3c-4f44-97e5-36962a2e2b4b)

![image](https://github.com/user-attachments/assets/818d788d-a32f-48cc-87ab-5a8b67e3ef6a)

- 需要注意的一点，如果存在孤立的细胞类型且需要显示该类型（默认参数），则程序会自动scale，无论是否设置；如果将`remove.isolate = TRUE`，那么宽度会随着信号强度改变
  https://github.com/jinworks/CellChat/blob/346fb615aaefa3c1a2830ec2fe3a933c5f723c83/R/visualization.R#L2722

</br>

# netVisual_heatmap()使用问题

- 统一尺度问题：行/列注释的barplot的y坐标，heatmap的颜色

并行展示多个样本的heatmap，用作者提供的方案，虽然legend只展示了一个，但是每个图片还是自己的尺度，并没有统一（We did not control the color scale in different plots.）

https://github.com/jinworks/CellChat/issues/20 &nbsp;&nbsp;&nbsp;&nbsp; https://github.com/sqjin/CellChat/issues/414 &nbsp;&nbsp;&nbsp;&nbsp; https://github.com/sqjin/CellChat/issues/361

- 解决方案，自己修改netVisual_heatmap代码

1. 统一行/列barplot的y轴

修改 https://github.com/jinworks/CellChat/blob/main/R/visualization.R#L1971 部分（ha2同），在`gp = gpar(fill = color.use.row, col=color.use.row)`中添加ylim，如`gp = gpar(fill = color.use.row, col=color.use.row),  ylim = c(0, 0.03))`

2. 统一legend色度

只在`heatmap_legend_param`中设置`at = colorbar.break`不能解决问题，虽然`colorbar.break`设置了范围0 - 0.015，但是颜色显示还会自动调整为最大值显示最深颜色，而不是中间色（假设最大值为0.01），`color.heatmap.use`自动映射颜色，所以需要将`color.heatmap.use`与固定范围结合，即在`color.heatmap.use`定义完后加一步固定范围的语句，如在 https://github.com/jinworks/CellChat/blob/main/R/visualization.R#L1978 完成后，加入

```
  color.heatmap.use <- circlize::colorRamp2(
  seq(0, 0.015, length.out = 100),  # 固定范围 [0, 0.015]
  color.heatmap.use                 # 使用你的颜色向量
```

修改前
![image](https://github.com/user-attachments/assets/28aff27d-9cc8-41b3-8002-abff228aa25a)

修改后
![image](https://github.com/user-attachments/assets/1a6c69bf-15c8-49c7-a16b-335e0682f994)



</br>

# netAnalysis_signalingRole_scatter统一坐标

与heatmap类似，修改代码 https://github.com/jinworks/CellChat/blob/346fb615aaefa3c1a2830ec2fe3a933c5f723c83/R/analysis.R#L2407 （或L2405，根据参数需要）
```
    gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size) +
      scale_x_continuous(limits = c(0,25)) +  # 设置 x 轴范围
      scale_y_continuous(limits = c(0,15))   # 设置 y 轴范围
```

# 提取特定pathway的基因

cellchat自带画基因表达量的命令`plotGeneExpression`，用来画指定通路或者L-R pairs基因的表达量的“violin”, “dot”, “bar”。但是如果对表达量作图有个性化的需求，可以用`extractEnrichedLR`来提取出相应的基因列表，再自己画图，注意加参数`geneLR.return = TRUE`（如果不加，只会返回$pairLR），例如`extractEnrichedLR(cellchat_bi, signaling = "TGFb" ,geneLR.return = TRUE)`， 返回结果为
```
$pairLR
    A data.frame: 6 × 1 interaction_name
    <chr>
    TGFB1_TGFBR1_TGFBR2
    TGFB2_TGFBR1_TGFBR2
    TGFB3_TGFBR1_TGFBR2
    TGFB1_ACVR1_TGFBR1
    TGFB2_ACVR1_TGFBR1
    TGFB3_ACVR1_TGFBR1
$geneLR
        'Tgfb1''Tgfb2''Tgfb3''Tgfbr1''Acvr1''Tgfbr2'
```
