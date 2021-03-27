library(piecewiseSEM)

#数据集，详情 ?keeley
data(keeley)
head(keeley)

#分段 SEM，详情 ?psem
#这里，对于每个独立的响应方程，直接使用简单线性回归，确定响应关系
#其它情况，如果已知变量间的某种非线性关系，可以使用非线性模型
keeley_psem <- psem(
	lm(firesev ~ age, data = keeley),
	lm(cover ~ firesev, data = keeley),
	data = keeley)

summary(keeley_psem)
plot(keeley_psem)

#将 age 和 cover 的因果关系考虑在内，创建新模型
keeley_psem2 <- psem(
  lm(cover ~ firesev + age, data = keeley),
  lm(firesev ~ age, data = keeley),
  data = keeley
)

#比较新（含 age 和 cover 的关系）旧（不含 age 和 cover 的关系）模型
#全局模型的卡方检验
anova(keeley_psem, keeley_psem2)

#参数估计值（回归系数）
coefs(keeley_psem, standardize = 'none')
coefs(keeley_psem, standardize = 'scale')	#可显示标准化后的

#定向分离测试
dSep(keeley_psem, .progressBar = FALSE)

#Fisher’s C 统计量，p 值等同于上述定向分离测试
fisherC(keeley_psem)

#AIC 值，用于比较多个模型
AIC(keeley_psem, keeley_psem2)

#一个带广义线性模型的 SEM 示例
set.seed(1)

dat <- data.frame(
	x = runif(100),
	y1 = runif(100),
	y2 = rpois(100, 1),
	y3 = runif(100)
)

modelList <- psem(
	lm(y1 ~ x, data = dat),
	glm(y2 ~ x, family = 'poisson', data = dat),
	lm(y3 ~ y1 + y2, data = dat),
	data = dat
)

summary(modelList, conserve = TRUE)
