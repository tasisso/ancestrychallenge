install.packages("data.table", repos="cloud.r-project.org")
library(data.table)

X.train <- as.matrix(data.frame(fread("data/ancestry_train.data")))
y.train <- as.matrix(data.frame(fread("data/ancestry_train.solution")))

X.test <- as.matrix(data.frame(fread("data/ancestry_test.data")))

#Augment training and testing genomatrices
cbind(a = 0, df)
X.train <- cbind(V0=1, X.train)
X.test <- cbind(V0=1, X.test)

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

	#matrix operations
	bhats = solve((tgenotrain %*% genotrain) + imatrix) %*% (tgenotrain %*% phenotrain)
	return(bhats)
}

#predict y
result <- X.test %*% beta(X.train, y.train, 0.001)

#save and zip the file
fwrite(as.data.frame(result), file = "predictions.csv",  
			 sep = " ", quote=FALSE, row.names = F, col.names = F)

system("zip -r predictions.zip predictions.csv")
