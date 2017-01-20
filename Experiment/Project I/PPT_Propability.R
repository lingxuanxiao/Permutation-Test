#PPT Propability(Rough)
#Program of permutation test(PPT)
#Quasi random number vs. Pseudo random number
#Author: Feixiao_L
#------------------------------------------------------------
#1. Quasi random number
#Halton
#Halton sequence base on p
#input: p is a prime number, n is the length of sequence
#output: quasi random number in [0, 1]
Halton = function(n, move = 0, p = 2){
	out = numeric(n)
	for (j in 1:n){
		i = j + move
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
#Test Sample
#Give the sample
set.seed(1234)
X = rnorm(10, 10, 5)
Y = rnorm(10, 15, 5)
D = mean(X) - mean(Y)
A = give_combn_diff(X, Y)
ilter = 100000
#Monte Carlo
#Pseudo Random List
set.seed(1234)
PRand_index = round(runif(ilter) * choose(20, 10))
PRand = numeric(ilter)
for (i in 1:ilter){
	PRand[i] = A[PRand_index[i]]
}
#Quasi Random List
QRand_index = round(Halton(ilter, 3) * choose(20, 10))
QRand = numeric(ilter)
for (i in 1:ilter){
	QRand[i] = A[QRand_index[i]]
}
#List by sample
set.seed(1234)
SRand = sample(A, ilter)
#Calculate
#Calculate the real prop: P(X>=D)
Prop = mean(as.numeric(A >= D))
#Calculate the simulate prop: P(X>=D)
O = as.numeric(A[1:ilter] >= D)
P = as.numeric(PRand >= D)
Q = as.numeric(QRand >= D)
S = as.numeric(SRand >= D)
prop_result = matrix(0, 4, ilter)
row.names(prop_result) = c('Order', 'Pseudo', 'Quasi', 'Sample')
for (i in 1:ilter){
	prop_result[1, i] = abs(mean(O[1:i]) - Prop)
	prop_result[2, i] = abs(mean(P[1:i]) - Prop)
	prop_result[3, i] = abs(mean(Q[1:i]) - Prop)
	prop_result[4, i] = abs(mean(S[1:i]) - Prop)
}
#Draw
pdf('C:/Users/Feixiao_L/Desktop/temp.pdf')
color = c('dark red', 'dark blue', 'forestgreen', 'yellow')
simul_time = c(1:ilter)
plot(NA, type = "n", log = "xy", xlim = c(1, ilter), ylim = c(1e-10, 0.01))
for (i in c(1:4)){
	lines(simul_time, prop_result[i, ], type = 'l', col = color[i])
}
dev.off()









