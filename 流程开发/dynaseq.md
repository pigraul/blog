# conversion提速

- 问题1：read.get_aligned_pairs会遍历read的每个位点，如果该read不存在突变，会浪费很多时间，且大多数reads都不存在突变

  解决方案：`nM` number of mismatches. 通过判断 nM>0 来确定是否需要进一步read.get_aligned_pairs

- 问题2：此前按照细胞拆分bam，部分细胞数据多导致运行时间过长
  
  解决方案：使用picard根据reads数拆分bam，确保每个bam的运行时间差不多（2-3分钟）；replacement步骤使用同一方法提速（先暂时保留临时拆分的文件，replacement步骤统一删除）
