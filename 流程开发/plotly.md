# 分面画图设置

- 自定义legend range，用于几个分面数据范围不一致，统一color scale颜色
```
fig = px.scatter(vol_df, x="-log10(p_val)", y="avg_logFC", color="gene_tor",facet_col="cluster",
                color_continuous_scale=px.colors.sequential.PuRd,
	     range_color=(0,0.5) ## 可以自定义legend range的范围
)
```

- 分面title设置
```
for i, facet in enumerate(fig.layout.annotations):  ## 分面title部分
    facet.textangle = -45  # 设置标题角度
    facet.font.size = 12    # 设置字体大小

fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
```

- 分面坐标轴相关
```
fig.update_xaxes(matches=None)   ## x轴不共享坐标范围

fig.update_xaxes(title=None)  ## 不显示x轴title

for axis in fig.select_xaxes():  ## 不显示x轴坐标刻度，如果不for,只会修改第一个
    axis.update(showticklabels=False)
```


# histogram
- histogram作图时，如果bin的范围总是跨越最大或者最小值，update_trace来调整 (https://community.plotly.com/t/histogram-bin-size-with-plotly-express/38927/5)

![image](https://github.com/user-attachments/assets/d0c23a83-9542-48fa-8308-8929633988df)


# Number Formatting

https://dash.plotly.com/datatable/data-formatting



# 在图像外加上border
https://plotly.com/python/axes/

对x/y轴画上镜像，显示line，且要给line颜色
```
fig.update_xaxes(mirror=True, showline=True,linecolor='grey')
fig.update_yaxes(mirror=True, showline=True,linecolor='grey')
```



# violin核密度函数设置
`fig.update_traces(spanmode=<VALUE>)`

**spanmode参数** 的功能解释，核心是定义 “密度函数计算的数据范围” 如何确定，分三种模式：

- spanmode="soft"（默认）

密度函数的计算范围会超出数据本身的极值：从 “样本最小值 - 2 个带宽” 到 “样本最大值 + 2 个带宽”。

（注：“带宽” 是计算核密度估计的参数，默认由 Silverman 法则自动确定，此模式能避免密度函数在数据边缘突然截断，让小提琴图两端更平滑。）

- spanmode="hard"

密度函数的计算范围严格限定在数据自身的极值内：从 “样本最小值” 到 “样本最大值”。

此模式下，密度函数会在数据边缘直接终止，小提琴图两端会呈现 “截平” 的效果，与ggplot2的默认的小提琴图行为一致。

- spanmode="manual"

自定义密度函数的计算范围，需配合span参数使用（span需传入一个列表，如[min_val, max_val]，明确指定范围的上下限），适用于需要精准控制可视化范围的场景。