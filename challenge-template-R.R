install.packages("data.table", repos="cloud.r-project.org")
library(data.table)

X.train <- as.matrix(data.frame(fread("data/ancestry_train.data")))
y.train <- as.matrix(data.frame(fread("data/ancestry_train.solution")))

X.test <- as.matrix(data.frame(fread("data/ancestry_test.data")))

#cost function
cost = function(y.true, y.hat){
	mse = -log10(mean((y.true-y.hat)**2)+1e-5)
	return(mse)
}
	
#baseline
base = matrix(1/3, nrow(y.train), ncol(y.train))
cost(y.train, base) 

y.test = matrix(1/3, nrow(X.test), ncol(y.train))

#ridge regression
beta = function(genotrain, phenotrain, lambda) {
	SNP.count = ncol(genotrain)
	imatrix = diag(SNP.count) * lambda
	tgenotrain = t(genotrain)

	bhats = solve((tgenotrain %*% genotrain) + imatrix) %*% (tgenotrain %*% phenotrain)
	return(bhats)
}


#bonferroni procedure
data.train <- cbind(y.train[,1], X.train)
model <- lm(y ~ ., data.frame(data.train))
pvalues <- summary(model)$coefficients[,4]
filterSNPS <- c()
for (i in 2:length(pvalues)) {
	if (pvalues[i] <= 0.05/ncol(X.train)) {
		filterSNPS = append(filterSNPS,i-1)
	}
}

#train on selected SNPs
bfX.train <- cbind(numeric(2000), X.train[,7], X.train[,23], X.train[,55], X.train[,60], X.train[,132], X.train[,178])
bfX.test <- cbind(numeric(3000), X.test[,7], X.test[,23], X.test[,55], X.test[,60], X.test[,132], X.test[,178])

result <- bfX.test %*% beta(bfX.train, y.train, 0.001)

#save and zip the file
fwrite(as.data.frame(result), file = "predictions.csv",  
			 sep = " ", quote=FALSE, row.names = F, col.names = F)

system("zip -r predictions.zip predictions.csv")
