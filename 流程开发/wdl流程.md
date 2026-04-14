# 参考

[quickstart](https://docs.openwdl.org/getting-started/quickstart.html)

[vsmalladi](https://vsmalladi.github.io/openwdl.github.io/getting-started/)

[WDL学习笔记](https://pzweuj.github.io/posts/WDL)

[WDL Example](https://academy.dnanexus.com/buildingworkflows/wdl/wdl_hello)


# wdl流程注释

`#`注释


# wdl布尔型值

`true` or `false` 均为小写，作为变量传递引用例如

`~{if plot_flag then "--all_type_plot" else ""}`


# 可选输入判断

`defined(input_opt)` 作为变量传递引用例如

`~{if defined(input_opt) then input_opt else ""}`




