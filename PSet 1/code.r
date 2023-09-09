a <- read.csv("gamma-ray.csv")
intervals_total <- sum((a[, 1]))
instance_total <- sum((a[, 2]))
lambda_0 <- 0.003880851

total_time <- c(a[, 1])
instances <- c(a[, 2])
lambdas <- array(c(instances / total_time, instances, total_time), dim = c(100, 3))
r_t <- sum(instances)
likelihood_1 <- function(t_i, g_i, lambda) {
    exp_coeff <- 0
    lambda_coeff <- 0
    numerator_prod <- 1
    denom_prod <- 1
    for (x in 1:100) {
        item1 <- t_i[x]^g_i[x]
        numerator_prod <- numerator_prod * item1
        denom_prod <- denom_prod * factorial(g_i[x])
        exp_coeff <- exp_coeff +  t_i[x]
        lambda_coeff <- lambda_coeff +  g_i[x]

    }
    exp_coeff <- -exp_coeff*lambda
    result <- exp(exp_coeff) * lambda^lambda_coeff * numerator_prod / denom_prod
    print(result)
}

likelihood_2 <- function(data) {
    # print(data)c(lambda[, 3]), c(lambda[, 2]),
    t_i <- c(data[, 3])
    g_i <- c(data[, 2])
    lams <- c(data[, 1])

    exp_coeff <- 0
    lambda_prod <- 1
    numerator_prod <- 1
    denom_prod <- 1
    for (x in 1:100) {
        
        numerator_prod <- numerator_prod * t_i[x]^g_i[x]
        denom_prod <- denom_prod * factorial(g_i[x])
        lambda_prod <-lambda_prod * lams[x]^g_i[x]

        exp_coeff <- exp_coeff +  lams[x] * t_i[x]

    }
    print(numerator_prod)
    exp_coeff <- -exp_coeff
    result <- exp(exp_coeff) * lambda_prod * numerator_prod / denom_prod
}
like1 <- likelihood_1(c(lambdas[, 3]), c(lambdas[, 2]),lambda_0 )
like2 <- likelihood_2(lambdas)
print(-2*log(like1/like2))
curve(dchisq(x, df = 99), from = 60, to = 140)
abline(v = 123)
points(104.3979,0) 
print(qchisq(.05, 99, lower.tail=FALSE))
print(pchisq(104.3979,99,lower.tail=FALSE))

data <- read.csv(file="golub_data/golub.csv",head=TRUE,sep=",")
p_vector = c()
 for (i in 1:3051) {
 x <- as.numeric(data[i:i,2:28])
 y <- as.numeric(data[i:i,29:39])
 res <- wilcox.test(x, y, alternative = "two.sided")
 p <- res$p.value
 p_vector <- c(p_vector,p)
 }

p_vector_BH <- p.adjust(p_vector,"BH")
print(sum(p_vector_BH<0.05))
print(sum(p_vector<0.05))

 library(MASS)
 X_data <- read.csv(file="syn_X.csv",header=FALSE,sep=",")
 y_data <- read.csv(file="syn_y.csv",header=FALSE,sep=",")
 X <- cbind(as.matrix(X_data),matrix(1,nrow=100))
y <- as.matrix(y_data)
 B_ideal <- ginv(t(X) %*% X) %*% t(X) %*% y
 print(B_ideal)


M_data <- read.csv(file="mortality.csv",header=TRUE,sep=",")
   
  X = as.matrix(M_data[,3:16])
  y = as.matrix(M_data[,2])
  X <- cbind(as.matrix(X),matrix(1,nrow=59))
  for (i in c(6,7,9,12,13,14)){
 X[,i] <- log(X[,i])
 }
  for (i in 1:14){
 X[,i] <- scale(X[,i])
 }
  y[,1] <- scale(y[,1])
  library(MASS)
  B_ideal <- ginv(t(X) %*% X) %*% t(X) %*% y
  # Initialize some things for the algorithm
  B <- matrix(0,nrow=15)
  MSE <- c(sum((y - X %*% B)^2)/length(y))
  d <- c(sqrt(sum((B - B_ideal)^2)))
  not_close <- TRUE
  epsilon <- 0.001
  i <- 0
  #Find the optimal L, and in turn, alpha
  L <- max(eigen(t(X) %*% X)$values)
  alpha <- 1/L
  #Actual algorithm; continues until we are changing by less then epsilon.
  while (not_close){
 i <- i + 1
 B_new <- B + 2*alpha*t(X) %*% (y - X %*% B)
 if (abs(sum((y - X %*% B_new)^2) - sum((y - X %*% B)^2)) < epsilon){
 not_close = FALSE
 }
 B <- B_new
 MSE <- c(MSE, sum((y - X %*% B)^2)/length(y))
 d <- c(d, sqrt(sum((B - B_ideal)^2)))
 }
# plot(qqline(M_data$Rain))
xvals = 0:i
#Plot MSE


#Plot distance between B and B_ideal


lr = lm(Mortality ~ JanTemp + JulyTemp +
            + RelHum + Rain + Educ + log(Dens) +
            + log(NonWhite) + WhiteCollar + log(Pop) +
            + House + Income + log(HC) + log(NOx) +
            + log(SO2), M_data)
 print(summary(lr))
 
 
 