## 切分bam，并行计算，减少时间
 (https://nf-co.re/subworkflows/bam_split_by_region)
 
切分等级考虑（目前选择4）：

	1. 按chr：单个文件大，运行时间过长
 
	2. 按gene：文件过多，io问题，后续合并文件数目限制问题
 
	3. 按block：将gene分成多个block切分；需要设置block原则，份数或者大小，不同基因组考虑较多；index时设置
 
	4. 按reads number且，目前选择1M reads一个文件（参数可调整，大文件调大，时间换空间），大多在300-500个文件


## Bulk slamseq流程
(https://github.com/nf-core/slamseq)

</br>


## 只想在自己的仓库中保持原始仓库的更新，而不做任何修改或合并请求（Pull Request），可以通过以下步骤实现：
添加远程原始仓库</br>
首先，将原始仓库添加为远程仓库：

### 在本地克隆你的仓库 AB
```
git clone https://github.com/your-username/AB.git
cd AB
```

### 添加原始仓库 AA 作为一个远程仓库
`git remote add upstream https://github.com/original-repo/AA.git`

同步更新
每当原始仓库 AA 有更新时，你需要同步这些更新到你的仓库 AB。

### 获取原始仓库的最新修改
`git fetch upstream`

### 切换到你的主分支（假设是 master，如果不是，请替换为你的主分支名称）
`git checkout master`

### 将原始仓库的更新合并到你的主分支
`git merge upstream/master`

这一步如果报错：
fatal: refusing to merge unrelated histories
这个错误通常发生在尝试合并两个不相关的 Git 历史（unrelated histories）时。这种情况通常发生在以下几种情况下：
1. 两个仓库的历史没有共同的祖先：这意味着 Git 无法找到两个仓库历史中的共同提交点，因此默认情况下拒绝合并。
2. 新的仓库：如果你正在尝试将一个全新的仓库合并到现有的仓库中，它们几乎肯定是没有共同历史的。
解决方法通常取决于你的具体情况：
- 解决方法一：允许合并不相关历史
如果你确定这两个仓库确实应该合并，你可以在执行 git merge 时加上 --allow-unrelated-histories 参数：

`git merge --allow-unrelated-histories <branch_name>`

这样做将允许 Git 合并两个不相关的历史。合并后，你可能需要解决任何冲突，并确认合并的结果是你期望的。

- 解决方法二：考虑使用 git pull 而不是 git merge
有时候你可能想要执行 git pull 而不是直接的 git merge。git pull 是 git fetch 后跟 git merge 的缩写。你可以尝试：

`git pull origin <branch_name> --allow-unrelated-histories`

这将先获取远程分支的内容，然后尝试合并。

- 解决方法三：考虑使用 git rebase 代替合并
另一个选择是使用 git rebase 将另一个分支的提交在你的分支上重放。这种方法通常用于使提交历史更加清晰，但要注意，git rebase 也可能会产生冲突需要解决。

`git rebase origin/<branch_name>`

在使用 git rebase 时也可能需要解决冲突，特别是在有共同修改的情况下。

注意事项

• 在执行任何这些操作之前，确保你理解你正在做的改变以及它们对你的项目的影响。

• 如果你不确定两个仓库是否应该合并，最好先备份或者在本地进行试验。

通过上述方法之一，你应该能够解决 "refusing to merge unrelated histories" 错误。






</br>

## Starsolo Statistics Explanation
https://github.com/alexdobin/STAR/issues/1887


</br>

## 怎么使用命令行的输出作为一个变量在script里使用

(https://stackoverflow.com/questions/66568781/how-to-call-a-variable-created-in-the-script-in-nextflow)


</br>

## 删除临时文件
`process.scratch = true`

(https://github.com/nextflow-io/nextflow/issues/165)





