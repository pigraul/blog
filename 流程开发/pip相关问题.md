某个目录（假设`b`）没有被pip install安装
| 问题原因 | 解决方法 |
|---------|---------|
| `b`缺少 `__init__.py` | 在`b`下创建`__init__.py` |
| 未在`setup.py`中声明 | 使用`find_packages(include=["my_package*"])` |
| `b`是非 Python 文件     | 添加`MANIFEST.in`并设置`include_package_data=True`     |
| 安装缓存问题     | 先卸载旧包 (`pip uninstall my_test`)，再重新安装     |
	
