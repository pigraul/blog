## csc_matrix统计非0元素个数
要计算 csc_matrix 每一行或者每一列大于0的元素个数，你可以使用 getnnz(axis) 方法，其中 axis 参数指定计算的方向。对于稀疏矩阵，getnnz(axis=0) 表示计算每列大于0的元素个数，而 getnnz(axis=1) 表示计算每行大于0的元素个数。
```
from scipy.sparse import csc_matrix

data = [[0, 1, 0],
        [0, 0, 1],
        [1, 0, 0]]
matrix = csc_matrix(data)

# 计算每列大于0的元素个数
non_zero_per_column = matrix.getnnz(axis=0)

# 计算每行大于0的元素个数
non_zero_per_row = matrix.getnnz(axis=1)
```
