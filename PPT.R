#Program Permutation Test
#Author: Feixiao_L
#------------------------------------------------------------
#置换检验
#1. 抽样检测法
#一次试验样本容量估计（使用CLT简化计算）
sample_size = function(N, p, z, delta){
	n = z * z * p * (1 - p) / (delta * delta)
	return(n)
}
#一次实验给出估计结果
single_exp = function(X, Y, size){
	nx = length(X)
	ny = length(Y)
	A = c(X, Y)
	S = sum(A)
	D = mean(X) - mean(Y)
	result = numeric(size)
	for (i in 1:size){
		result[i] = sum(sample(A, nx))
	}
	outa = mean(as.numeric((result/nx - (S - result)/ny) >= D))
	outb = mean(as.numeric((result/nx - (S - result)/ny) <= D))
	return(min(outa, outb))
}
#检验主体函数
#计算概率
#输入：样本A, 样本B, 迭代次数（一次试验的样本容量）默认为20000，严格不超过15w
#输出：单侧检验的概率结果，只返回概率度小于等于0.5的那边
#计算模式：
#a. 简单一次实验（不控制误差）
#b. 自动调整误差一次试验（控制误差，默认95%控制在p的下一个数量级但最大不超过0.001）
RPT = function(X, Y, size = 20000){
	nx = length(X)
	ny = length(Y)
	A = c(X, Y)
	D = sum(X)
	result = 0
	for (i in 1:size){
		if (sum(sample(A, nx)) <= D){
			result = result + 1
		}
	}
	return(c(result, size, result/size, 1.95 * sqrt(result/size * (1 - result/size) / size)))
}

#------------------------------------------------------------
#2. 枚举法
#利用'combn'函数枚举出所有组合数值
give_combn_diff = function(X, Y){
	nx = length(X)
	ny = length(Y)
	A = c(X, Y)
	S = sum(A)
	mean_diff = function(vec){
		return(mean(vec) - (S - sum(vec)) / ny)
	}
	return(combn(A, nx, mean_diff))
}
#跳跃法减少计算次数
EPT = function(Sa, Sb, front = 'left', include = TRUE){
	#初始化
	A = sum(Sa)
	S = sort(c(Sa, Sb))
	N = length(S)
	R = length(Sa)
	C = c(1:R)
	C_stop = c((N - R + 1):N)
	C_lim = C_stop
	calc = 0
	#定义子函数
	partsum = function(vec, index){s = 0; for (i in index){s = s + vec[i]}; return(s)}
	judge = function(index1, index2, len){out = FALSE; for(i in len:1){if (index1[i] < index2[i]){out = TRUE; break}}; return(out)}
	j = 0
	#主循环
	while (TRUE){
		j = j+1
		if (judge(C, C_stop, R)){
			if(judge(C, C_lim, R)){
				if (partsum(S, C) > A){
					if (C[R] < C_lim[R]){C_lim = C; C_stop = c((C[R] - R + 1):C[R])}
					for (i in (R - 1):1){if (C[i] < N - R + i){break}}
					C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
				}else{
					calc = calc + 1
					for (i in R:1){if (C[i] < N - R + i){break}}
					C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
				}
			}else{
				for (i in R:1){if (C[i] < C_stop[i]){break}}
				C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
			}
		}else{
			break
		}
	}
	if (partsum(S, C) <= A){calc = calc + 1}
	return(c(calc, j))
}

EEPT = function(Sa, Sb, front = 'left', include = TRUE){
	#初始化
	A = sum(Sa)
	S = sort(c(Sa, Sb))
	N = length(S)
	R = length(Sa)
	C = c(1:R)
	C_stop = c((N - R + 1):N)
	calc = 0
	#定义子函数
	partsum = function(vec, index){s = 0; for (i in index){s = s + vec[i]}; return(s)}
	judge = function(index1, index2, len){out = FALSE; for(i in len:1){if (index1[i] < index2[i]){out = TRUE; break}}; return(out)}
	j = 0
	#主循环
	while (TRUE){
		j = j+1
		if (judge(C, C_stop, R)){
			if (partsum(S, C) <= A){calc = calc + 1}
			for (i in R:1){if (C[i] < C_stop[i]){break}}
			C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
		}else{
			break
		}
	}
	if (partsum(S, C) <= A){calc = calc + 1}
	return(c(calc, j))
}

EEEPT = function(Sa, Sb, front = 'left', include = TRUE){
	#初始化
	A = sum(Sa)
	S = sort(c(Sa, Sb))
	N = length(S)
	R = length(Sa)
	C = c(1:R)
	C_stop = c((N - R + 1):N)
	C_lim = C_stop
	calc = 0
	#定义子函数
	partsum = function(vec, index){s = 0; for (i in index){s = s + vec[i]}; return(s)}
	judge = function(index1, index2, len){out = FALSE; for(i in len:1){if (index1[i] < index2[i]){out = TRUE; break}}; return(out)}
	j = 0
	#主循环
	while (TRUE){
		j = j+1
		if (judge(C, C_stop, R)){
			if(judge(C, C_lim, R)){
				if (partsum(S, C) > A){
					if (C[R] < C_lim[R]){C_lim = C; C_stop = c((C[R] - R + 1):C[R])}
					for (i in (R - 1):1){if (C[i] < C[R] - R + i){break}}
					C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
				}else{
					calc = calc + 1
					for (i in R:1){if (C[i] < N - R + i){break}}
					C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
				}
			}else{
				for (i in R:1){if (C[i] < C_stop[i]){break}}
				C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
			}
		}else{
			break
		}
	}
	if (partsum(S, C) <= A){calc = calc + 1}
	return(c(calc, j))
}


#C_stop == c((N-R+1):N) | C == c((C_lim[R]-R+1):C_lim[R])
#------------------------------------------------------------
#测试样例生成
set.seed(1234)
S = 20000
N = 10
X = round(rnorm(N, 10, 5), 4) #* 10000
Y = round(rnorm(N, 15, 5), 4) #* 10000
#测试用公用代码
t.test(X, Y)
sum(as.numeric(combn(c(X, Y), N, sum) <= sum(X)))
#EPT(X, Y)
EEPT(X, Y)
EEEPT(X, Y)
RPT(X, Y, S)
#system.time()
#system.time(EPT(X, Y))
system.time(EEPT(X, Y))
system.time(EEEPT(X, Y))
system.time(RPT(X, Y, S))
