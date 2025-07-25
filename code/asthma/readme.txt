lcczscore_pip_simplify.ipynb
	## 1. 简化计算单细胞lcczscore、
      		1.1 pip>1时，正常随机 ，只记录单细胞中基因数量有变化的点  【高扰动】
     		1.2 pip = 0时的lcczscore和weighted-lcczscore，expressmean-lcczscore 【低扰动】
			* lcczscore:  最大连通分支上节点数量的显著性
			* weighted-lcczscore: 最大连通分支上表达值的显著性
			* expressmean-lcczscore: 最大连通分支上平均每个基因表达的显著性
	## 2. 根据 lcczscore pip 计算得到的结果提取每个细胞的6个点
		对于每个细胞
			* 统计pip>1 的lcczscore最大值，以及该位置对应的weightedlcczscore，expresszscore
			* 统计pip=0处(所有基因) 的lcczscore，以及该位置对应的weightedlcczscore，expresszscore


lcczscore pip-simplify-CoreandPeri.ipynb
	不画每个细胞的lcczsore变化曲线，只统计性状核心、外围基因在单个细胞内的表达情况
