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
