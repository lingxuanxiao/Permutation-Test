#PPT Propability(Adv Ed.)
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
#Give the test sample
set.seed(12345)
X = rnorm(10, 10, 5)
Y = rnorm(10, 15, 5)
Diff = mean(X) - mean(Y)
All = as.numeric(give_combn_diff(X, Y) >= Diff)
#------------------------------------------------------------
#Const
Prop = mean(All)
ilter = 100
size = 2 ^ c(7:16) - 1
S = choose(20, 10)
#------------------------------------------------------------
#Monte Carlo
#initialization
set.seed(54321)
p_pseudo_result = matrix(0, 10, 100)
p_quasiR_result = matrix(0, 10, 100)
p_quasiF_result = numeric(10)
p_sample_result = matrix(0, 10, 100)
#simulate
for (i in 1:10){
  for (j in 1:ilter){
    index1 = ceiling(runif(size[i]) * S)
    index2 = ceiling(Halton(size[i], round(runif(1) * 128)) * S)
    Pseudo = numeric(size[i])
    QuasiR = numeric(size[i])
    for (k in 1:size[i]){
      Pseudo[k] = All[index1[k]]
      QuasiR[k] = All[index2[k]]
    }
    p_pseudo_result[i, j] = mean(Pseudo)
    p_quasiR_result[i, j] = mean(QuasiR)
    p_sample_result[i, j] = mean(sample(All, size[i]))
  }
  index =  ceiling(Halton(size[i]) * S)
  QuasiF = numeric(size[i])
  for (k in 1:size[i]){
    QuasiF[k] = All[index[k]]
  }
  p_quasiF_result[i] = mean(QuasiF)
}
#------------------------------------------------------------
#OutFile
write.csv(t(p_pseudo_result), 'C:/Users/Feixiao_L/Desktop/temp1.csv', row.names = FALSE, quote = FALSE)
write.csv(t(p_quasiR_result), 'C:/Users/Feixiao_L/Desktop/temp2.csv', row.names = FALSE, quote = FALSE)
write.csv(t(p_sample_result), 'C:/Users/Feixiao_L/Desktop/temp3.csv', row.names = FALSE, quote = FALSE)
write.csv(t(p_quasiF_result), 'C:/Users/Feixiao_L/Desktop/temp4.csv', row.names = FALSE, quote = FALSE)
#------------------------------------------------------------
#Mean

#------------------------------------------------------------
#Draw
pdf('C:/Users/Feixiao_L/Desktop/temp.pdf')
color = c('dark red', 'dark blue', 'forestgreen', 'yellow')
simul_time = c(1:ilter)
plot(NA, type = "n", log = "xy", xlim = c(1, ilter), ylim = c(1e-10, 0.01))
for (i in c(1:4)){
	lines(simul_time, prop_result[i, ], type = 'l', col = color[i])
}
dev.off()