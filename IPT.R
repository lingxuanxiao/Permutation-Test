#Program Permutation Test (improved ed.)
#Author: Feixiao_L
#======================================================================================
#有关置换检验的若干脚本（增进版）
#======================================================================================
#检验算法：
#1. Randomized permutation test
RPT = function(Sa, Sb, size = 20000){
	A = sum(Sa)
	S = c(Sa, Sb)
	N = length(Sa)
	result = numeric(size)
	for (i in 1:size){result[i] = sum(sample(S, N))}
	return(mean(as.numeric(result <= A)))
}
#2. Exact permutation test
EPT = function(Sa, Sb){
	A = sum(Sa)
	S = sort(c(Sa, Sb))
	N = length(S)
	R = length(Sa)
	C = c(1:R)
	calc = 0
	j = 0
	while (TRUE){
		j = j+1
		if (C[1] < N-R+1){
			if (sum(S[C]) <= A){calc = calc + 1}
			for (i in R:1){if (C[i] < N-R+i){break}}
			C = c(C[1:i-1], c((C[i]+1):(C[i]+1+R-i)))
		}else{
			break
		}
	}
	if (sum(S[C]) <= A){calc = calc + 1}
	return(c(calc, j))
}
#3. Improved exact permutation test
IEPT = function(Sa, Sb){
	A = sum(Sa)
	S = sort(c(Sa, Sb))
	N = length(S)
	R = length(Sa)
	C = c(1:R)
	calc = 0
	j = 0
	Sup = N
	C_lim = c((N - R + 1):N)
	judge = function(index1, index2){out = FALSE; for(i in R:1){if (index1[i] < index2[i]){out = TRUE; break}}; return(out)}
	while (TRUE){
		j = j+1
		if (C[1] < Sup-R+1){
			if (judge(C, C_lim)){
				if (sum(S[C]) <= A){
					calc = calc + 1
					for (i in R:1){if (C[i] < N-R+i){break}}
					C = c(C[1:i-1], c((C[i]+1):(C[i]+1+R-i)))
				}else{
					if (C[R] < C_lim[R]){C_lim = C; Sup = C[R]}
					for (i in R:1){if (C[i] < C[R]-R+i){break}}
					C = c(C[1:i-1], c((C[i]+1):(C[i]+R-i+1)))
				}
			}else{
				for (i in R:1){if (C[i] < Sup-R+i){break}}
				C = c(C[1:i-1], c((C[i]+1):(C[i]+1+R-i)))
			}
		}else{
			break
		}
	}
	if (sum(S[C]) <= A){calc = calc + 1}
	return(c(calc, j))
}
#======================================================================================
#测试样例
set.seed(1234)
S = 20000
N = 10
X = round(rnorm(N, 10, 5), 4) #* 10000
Y = round(rnorm(N, 15, 5), 4) #* 10000
sum(as.numeric(combn(c(X, Y), N, sum) <= sum(X)))
round(RPT(X, Y, S) * choose(2*N, N))
EPT(X, Y)
IEPT(X, Y)
system.time(sum(as.numeric(combn(c(X, Y), N, sum) <= sum(X))))
system.time(round(RPT(X, Y, S) * choose(2*N, N)))
system.time(EPT(X, Y))
system.time(IEPT(X, Y))
#======================================================================================
#数据实验
#实验1：不同信噪比下给出置信度的比较
#1. 测试集生成
set.seed(12345)
N = 10
dB = -4:10			#dB是表示功率量之比的一种单位，等于功率强度之比的常用对数的10倍
sigma = 1
mu = 10^(dB/20) * sigma
mu = c(0, mu)
TestSet = list()
for (i in 1:16){TestSet[[i]] = round(rnorm(N * 1000, mu[i], sigma), 4)}
#2. 进行测试
set.seed(54321)
Tresult = matrix(0, 1000, 15)
Eresult = matrix(0, 1000, 15)
Rresult_1 = matrix(0, 1000, 15)
Rresult_2 = matrix(0, 1000, 15)
Rresult_3 = matrix(0, 1000, 15)
Iresult = matrix(0, 1000, 15)
for (i in 1:1000){
	X = TestSet[[1]][1:10 + (i-1)*10]
	for (j in 2:16){
		Y = TestSet[[j]][1:10 + (i-1)*10]
		#Tresult[i, (j-1)] = t.test(X, Y, 'less')$p.value
		#Eresult[i, (j-1)] = EPT(X, Y)[1]/choose(20, 10)
		#Rresult_1[i, (j-1)] = RPT(X, Y, 20000)
		#Rresult_2[i, (j-1)] = RPT(X, Y, 50000)
		#Rresult_3[i, (j-1)] = RPT(X, Y, 100000)
		#Iresult[i, (j-1)] = IEPT(X, Y)[1]/choose(20, 10)
	}
}
write.csv(Tresult, file = "C:/Users/Feixiao_L/Desktop/PTResult/T.csv", row.names = F, quote = F)
write.csv(Eresult, file = "C:/Users/Feixiao_L/Desktop/PTResult/E.csv", row.names = F, quote = F)
write.csv(Rresult_1, file = "C:/Users/Feixiao_L/Desktop/PTResult/R1.csv", row.names = F, quote = F)
write.csv(Rresult_2, file = "C:/Users/Feixiao_L/Desktop/PTResult/R2.csv", row.names = F, quote = F)
write.csv(Rresult_3, file = "C:/Users/Feixiao_L/Desktop/PTResult/R3.csv", row.names = F, quote = F)
write.csv(Iresult, file = "C:/Users/Feixiao_L/Desktop/PTResult/I.csv", row.names = F, quote = F)
#实验2：不同信噪比下计算时间的比较
#实验3：不同样本容量下计算时间的比较
