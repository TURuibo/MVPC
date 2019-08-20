# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

n = 100000
a =  rnorm(n, 0, 1)
b =  a + rnorm(n, 0, 1)
b2 = a + rnorm(n, 0, 1)
c =  b + b2 + rnorm(n, 0, 1)

e =  rnorm(n, 0, 1)

data_t = data.frame(a)
data_t[,2]=b
data_t[,3]=b2
data_t[,4]=c

cor(data_t)
suffStat_t  = list(C=cor(data_t), n=n)
gaussCItest(1, 4, c(2,3), suffStat_t)

# Generating missing data
w  = a + c + e 
r = w >  0
data_t[,5]=w

data_t2 = data_t
data_t2[!r,1] = NA

suffStat_t2 = list(data=data_t2)
gaussCItest_td(1, 4, c(2,3),suffStat_t2)

## Permutation-based correction methods 

prt = list()
m = c(1)
prt[[1]]=c(5)

prt_m<-data.frame(m=m)
prt_m[['prt']]<-prt

suffStat = list(data = data_t2, prt_m = prt_m)

PermCCItest(1, 4, c(2,3), suffStat)



# Step 1: Learning generaive model for {X, Y, S} to impute X, Y, and S
data <- test_wise_deletion(c(1, 2, 3), data_t2)
xnam <- paste0("data[,", 3,"]")
fmla <- as.formula(paste("data[,1] ~ 0 + ", paste(xnam, collapse= "+")))
fit <- lm(formula = fmla, data = data)
res <- residuals(fit) # residuals

fmla2 <- as.formula(paste("data[,2] ~ 0 + ", paste(xnam, collapse= "+")))
fit2 <- lm(formula = fmla2, data = data)
res2 <- residuals(fit2) # residuals

# Step 2: Shuffle the source "W" -- parents of the missingness indicators
w_p = sample(w)
w_p = w_p[1:length(res)]
data[,3] = w_p
# Step 3: Generate the virtual data follows the full data distribution P(X, Y, S)

vir = predict(fit,list(data)) + res
vir2 = predict(fit2,list(data)) + res2

data_t3 = data.frame(vir)
data_t3[,2]=vir2

suffStat_t3 = list(data=data_t3)
gaussCItest_td(1, 2, c(), suffStat_t3)
