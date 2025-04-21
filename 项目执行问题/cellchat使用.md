# cellcat不同版本相关问题

## cellchat的v1和v2不兼容
-  v1无法读取v2生成的rds
-  v1版本中有些slot在v2中不存在或者变更

##  命令使用有所不同
主要是在一些作图中，`from / to`变为`sources.use / targets.use`

比如`netVisual_bubble`，在v1中为`netVisual_bubble(cellchat, from = ligands, to = recepters)`，在v2中为`netVisual_bubble(cellchat, sources.use = ligands, targets.use = recepters)`


# netVisual_heatmap()使用问题

统一尺度问题

并行展示多个样本的heatmap，虽然legend只展示了一个，但是每个图片还是自己的scale尺度

https://github.com/jinworks/CellChat/issues/20

https://github.com/sqjin/CellChat/issues/414
