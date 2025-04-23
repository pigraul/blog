# 图形释义

边的颜色/权重、节点的颜色/大小/形状的解释：在所有的可视化图表中，边的颜色与作为发送方的源相对应，并且边的权重与相互作用的强度成正比。更粗的边线表示更强的信号。

在层次图和圆形图中，圆形的大小与每个细胞组中的细胞数量成正比。在层次图中，实心圆和空心圆分别代表源和目标。

在弦图中，内部较细的条形的颜色表示接收来自相应外部条形信号的目标。内部条形的大小与目标接收到的信号强度成正比。（下图红色框标记的）这样的内部条形有助于解读复杂的弦图。**需要注意的是，对于某些细胞组，存在一些没有任何弦连接的内部条形，请忽略它们，因为这是一个尚未被 circlize 包解决的问题。**

Explnations of edge color/weight, node color/size/shape: In all visualization plots, edge colors are consistent with the sources as sender, and edge weights are proportional to the interaction strength. Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group. In the hierarchy plot, solid and open circles represent source and target, respectively. In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note that there exist some inner bars without any chord for some cell groups, please just igore it because this is an issue that has not been addressed by circlize package.

![image](https://github.com/user-attachments/assets/508e4eb1-e4b6-4cab-b5e5-f9dfee31cbae)

</br>

# cellcat不同版本相关问题

## cellchat的v1和v2不兼容
-  v1无法读取v2生成的rds
-  v1版本中有些slot在v2中不存在或者变更

##  命令使用有所不同
主要是在一些作图中，`from / to`变为`sources.use / targets.use`

比如`netVisual_bubble`，在v1中为`netVisual_bubble(cellchat, from = ligands, to = recepters)`，在v2中为`netVisual_bubble(cellchat, sources.use = ligands, targets.use = recepters)`

</br>

# netVisual_heatmap()使用问题

统一尺度问题

并行展示多个样本的heatmap，虽然legend只展示了一个，但是每个图片还是自己的scale尺度

https://github.com/jinworks/CellChat/issues/20

https://github.com/sqjin/CellChat/issues/414
