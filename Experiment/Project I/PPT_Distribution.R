#PPT Distribution(Rough)
#Edition: losing accuracy & reducing work time
#Program of permutation test(PPT)
#Quasi random number vs. Pseudo random number
#Author: Feixiao_L
#------------------------------------------------------------
#1. Quasi random number
#Halton
#Halton sequence base on p
#input: p is a prime number, n is the length of sequence
#output: quasi random number in [0, 1]
Halton = function(n, p = 2){
	out = numeric(n)
	for (j in 1:n){
		i = j
		k = 1
		b = 0
		#i is the origin number
		#base conversion
		while (i>0){
			b = i %% p * p^(-k) + b
			i = i %/% p
			k = k + 1
		}
		out[j] = b
	}
	return(out)
}
#------------------------------------------------------------
#2. Producing Combination
#Based on 'combn'
#Give the combn different result
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
#------------------------------------------------------------
#3. Calculate the similarity
#Give the prop table of the result, based on 'table'
give_table = function(X, bin = NULL){
	X = table(X)
	if (is.null(bin)){
		out = X / sum(X)
	}else{
		index = numeric(length(bin))
		for (i in 1:length(as.numeric(names(X)))){
			index[which(bin == as.numeric(names(X))[i])] = as.numeric(X)[i]
		}
		out = index / sum(index)
		names(out) = bin
	}
	return(out)
}
#Give the distance of 2 vectors
give_distance = function(X, Y, p = 2){
	out = (sum(abs(X - Y) ^ p)) ^ (1 / p)
	return(out)
}
#Give the angle similarity of 2 vectors
give_sim = function(X, Y, standrad = FALSE){
	if (standrad){
		return(1 - abs(cor(X, Y, method = 'pearson')))
	}else{
		return(as.numeric(1 - crossprod(X, Y)/sqrt(crossprod(X, X)*crossprod(Y, Y))))
	}
}
#Give the Kullback–Leibler divergence of 2 vectors
give_KL = function(X, Y){
#calculate the KL(X||Y)
#lim(x->0)(x*log(x)) = 0
	n = length(X)
	out = 0
	for (i in 1:n){
		if (X[i] != 0){
			out = out + X[i] * log(X[i] / Y[i])
		}
	}
	return(as.numeric(out))
}
#------------------------------------------------------------
#Test Sample
#Const
digit = 1
ilter = 10000
#Give the sample
set.seed(1234)
X = rnorm(10, 10, 5)
Y = rnorm(10, 15, 5)
A = give_combn_diff(X, Y)
A = round(A, digit)
SI = as.numeric(give_table(A))
BIN = as.numeric(names(give_table(A)))
#Random List
set.seed(1234)
PRand_index = round(runif(ilter) * choose(20, 10))
PRand = numeric(ilter)
for (i in 1:ilter){
	PRand[i] = A[PRand_index[i]]
}
QRand_index = round(Halton(ilter) * choose(20, 10))
QRand = numeric(ilter)
for (i in 1:ilter){
	QRand[i] = A[QRand_index[i]]
}
set.seed(1234)
SRand = sample(A, ilter)
#initialize the result
Manhattan = matrix(numeric(4 * ilter), nrow = 4)
Euclid = matrix(numeric(4 * ilter), nrow = 4)
Angle = matrix(numeric(4 * ilter), nrow = 4)
Corelation = matrix(numeric(4 * ilter), nrow = 4)
KL = matrix(numeric(4 * ilter), nrow = 4)
#Calculate
for (i in 1:ilter){
	#Order
	O = give_table(A[1:i], bin = BIN)
	Manhattan[1, i] = give_distance(O, SI, p = 1)
	Euclid[1, i] = give_distance(O, SI)
	Angle[1, i] = give_sim(O, SI)
	Corelation[1, i] = give_sim(O, SI, standrad = TRUE)
	KL[1, i] = give_KL(O, SI)
	#Pseudo Random
	P = give_table(PRand[1:i], bin = BIN)
	Manhattan[2, i] = give_distance(P, SI, p = 1)
	Euclid[2, i] = give_distance(P, SI)
	Angle[2, i] = give_sim(P, SI)
	Corelation[2, i] = give_sim(P, SI, standrad = TRUE)
	KL[2, i] = give_KL(P, SI)
	#Quasi Random
	Q = give_table(QRand[1:i], bin = BIN)
	Manhattan[3, i] = give_distance(Q, SI, p = 1)
	Euclid[3, i] = give_distance(Q, SI)
	Angle[3, i] = give_sim(Q, SI)
	Corelation[3, i] = give_sim(Q, SI, standrad = TRUE)
	KL[3, i] = give_KL(Q, SI)
	#Sample
	S = give_table(SRand[1:i], bin = BIN)
	Manhattan[4, i] = give_distance(S, SI, p = 1)
	Euclid[4, i] = give_distance(S, SI)
	Angle[4, i] = give_sim(S, SI)
	Corelation[4, i] = give_sim(S, SI, standrad = TRUE)
	KL[4, i] = give_KL(S, SI)
}
#Draw
pdf('C:/Users/Feixiao_L/Desktop/temp.pdf')
color = c('dark red', 'dark blue', 'forestgreen', 'yellow')
simul_time = c(1:ilter)
#原始图
plot(simul_time, numeric(ilter), 'l', col = 'black', ylim = c(0, 1))
for (i in c(1:4)){
	lines(simul_time, Manhattan[i, ], 'l', col = color[i])
}
plot(simul_time, numeric(ilter), 'l', col = 'black', ylim = c(0, 1))
for (i in c(1:4)){
	lines(simul_time, Euclid[i, ], 'l', col = color[i])
}
plot(simul_time, numeric(ilter), 'l', col = 'black', ylim = c(0, 1))
for (i in c(1:4)){
	lines(simul_time, Angle[i, ], 'l', col = color[i])
}
plot(simul_time, numeric(ilter), 'l', col = 'black', ylim = c(0, 1))
for (i in c(1:4)){
	lines(simul_time, Corelation[i, ], 'l', col = color[i])
}
plot(simul_time, numeric(ilter), 'l', col = 'black', ylim = c(0, 1))
for (i in c(1:4)){
	lines(simul_time, KL[i, ], 'l', col = color[i])
}
#对数作图
#plot(NA, type = "n", log = "xy", ylim = c(1e-10,1), xlim = c(1e-10,1))


dev.off()









