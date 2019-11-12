proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

td_pc = list()
mvpc_drw = list()
mvpc_permc = list()
ref = list()
pc = list()

std_td_pc = list()
std_mvpc_drw = list()
std_mvpc_permc = list()
std_ref = list()
std_pc = list()

r_td_pc = list()
r_mvpc_drw = list()
r_mvpc_permc = list()
r_ref = list()
r_pc = list()

## ********* Synthethic Binary Data Generation ********* 
sz  <-600000

count = 1
shd_mvpc_permc = c()
shd_mvpc_drw = c()

shd_ref = c()
shd_pc = c()
shd_td_pc = c()

rp_mvpc_drw = list()
rp_mvpc_permc = list()

rp_ref = list()
rp_pc = list()
rp_td_pc = list()


i_g = 13
print(paste('graph=', i_g,'number of samples: ', sz))
data.name = paste('data',i_g,sep='')
graph.name = paste('graph',i_g,sep='')

data.file <- paste(data_path,'/',data.name,'.txt',sep="")
graph.file<- paste(data_path,'/',graph.name,'.txt',sep="")

data = load_bin_data(data.file)
data = data[1:sz,]
DAG = load_bin_graph(graph.file)
CPDAG = dag2cpdag(DAG)

# Detect colliders and  Colliders' parents
cldr <- detect_colliders(DAG)
cldr_prt <- detect_colliders_prt(DAG, cldr)
# Choose missingness inidcator and their parents


p_m <- create_mnar_ind(cldr,cldr_prt,
                       num_var=20, 
                       num_extra_e=5, 
                       num_m = 10) 
ms = p_m$ms
prt_ms = p_m$prt_ms

# Generate missing values 
mask = data != data

for(i in 1:length(ms)){
  nsample = nrow(data)
  pr1 <- plogis( 2*data[,prt_ms[i]]-1)
  r <- rbinom(nsample, 1, prob = pr1)==1
  mask[,ms[i]] = r
}
data_m = data
data_m[mask] = NA


prt_m<-data.frame(m=ms)

prt = list()
for(i in 1:length(prt_ms)){
  prt[[i]] = c(prt_ms[i])  
}

prt_m[['prt']]<-prt

# # ********* Correction *********
suffStat = list(data=data_m,adaptDF=FALSE)
res_mvpc_permc<-mvpc(suffStat, binCItest_td, binCItest.permc,prt_m, alpha=0.05, p=20)
shd(res_mvpc_permc,CPDAG)

## ********* Correction *********
suffStat = list(data=data_m,adaptDF=FALSE)
res_mvpc_drw<-mvpc(suffStat, binCItest_td, binCItest.drw, prt_m, alpha=0.05, p=20)
shd(res_mvpc_drw,CPDAG)

suffStat = list(data=data_m,prt_m=prt_m,adaptDF=FALSE)

# binCItest.permc(7,15,c(),suffStat)
# binCItest_td(7,15,c(),suffStat)
binCItest.drw(7,15,c(),suffStat)
binCItest.drw(7,19,c(),suffStat)
prt_m


x_ <- data[,7]
# py <- plogis( x_ ); y_ <- rbinom(n, 1, prob = py) # {0,1}
y_ <- data[,15]
dat <- cbind(x_,y_)

binCItest(1,2,c(), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.

# Add missingness mechanism
w_ = data[,19]
# w2_ <- rbinom(n, 1, pr=1/2)
# pw2 <- plogis( 2*x_-0.5); w2_ <- rbinom(n, 1, prob = pw2) # {0,1}

pr1 <- plogis( 2*w_-1); yr_ <- rbinom(n, 1, prob = pr1)
pr <- plogis( 2*y_-1); wr_ <- rbinom(n, 1, prob = pr)

x = x_
y = y_
w = w_
# w2 = w2_

y[yr_==1] = NA
w[wr_==1] = NA

dat = cbind(x,y,w)
prt = list()
prt[[1]] = c(3)
prt[[2]] = c(2)

prt_m<-data.frame(m=c(2,3))
prt_m[['prt']]<-prt

suffStat = list(data = dat, prt_m=prt_m, adaptDF = FALSE)
binCItest_td(1, 2, c(), suffStat)
binCItest.drw(1, 2, c(), suffStat)










dat = data_m[,c(7,19, 15)]
r2 = is.na(dat[,2]) * 1 
r3 = is.na(dat[,3]) * 1 
# P(R2 = 0 , R3 = 0)/ P(R2 = 0 | x3, R3=0) P(R3 = 0 | x2, R2 = 0 )

pr0 = sum((r2==0) & (r3==0))/nrow(dat)

pr3x20r20 = sum(r3[r2==0] == 0 & (dat[r2==0,2]==0) ) / sum(dat[r2==0,2]==0)
pr3x21r20 = sum(r3[r2==0] == 0 & (dat[r2==0,2]==1) ) / sum(dat[r2==0,2]==1)

pr20x30r30 = sum(r2[r3==0]==0 & dat[r3==0,3]==0)/sum(dat[r3==0,3]==0)
pr20x31r30 = sum(r2[r3==0]==0 & dat[r3==0,3]==1)/sum(dat[r3==0,3]==1)

px30x20 = pr0/(pr3x20r20*pr20x30r30)
px30x21 = pr0/(pr3x21r20*pr20x30r30)
px31x20 = pr0/(pr3x20r20*pr20x31r30)
px31x21 = pr0/(pr3x21r20*pr20x31r30)

na.w = is.na(dat[,2]) | is.na(dat[,3])
a = (dat[,2] == 0) & (dat[,3] == 0 )
a[na.w] = FALSE
b = (dat[,2] == 0) & (dat[,3] == 1 )
b[na.w] = FALSE
c = (dat[,2] == 1) & (dat[,3] == 0 )
c[na.w] = FALSE
d = (dat[,2] == 1) & (dat[,3] == 1 )
d[na.w] = FALSE

weights = rep(1, nrow(dat))
weights[a] = px30x20
weights[b] = px31x20
weights[c] = px30x21
weights[d] = px31x21
weights[na.w] = NA

data = cbind(dat[,c(1,2)],weights)
data = test_wise_deletion(c(1,2,3), data)
binCItest_w(1, 2, c(), data[,c(3)], list(dm = data[,c(1,2)], adaptDF = FALSE))




dat = data_m[,c(7,15, 19)]
r2 = is.na(dat[,2]) * 1 
r3 = is.na(dat[,3]) * 1 
# P(R2 = 0 , R3 = 0)/ P(R2 = 0 | x3, R3=0) P(R3 = 0 | x2, R2 = 0 )

pr0 = sum((r2==0) & (r3==0))/nrow(dat)

pr3x20r20 = sum(r3[r2==0] == 0 & (dat[r2==0,2]==0) ) / sum(dat[r2==0,2]==0)
pr3x21r20 = sum(r3[r2==0] == 0 & (dat[r2==0,2]==1) ) / sum(dat[r2==0,2]==1)

pr20x30r30 = sum(r2[r3==0]==0 & dat[r3==0,3]==0)/sum(dat[r3==0,3]==0)
pr20x31r30 = sum(r2[r3==0]==0 & dat[r3==0,3]==1)/sum(dat[r3==0,3]==1)

px30x20 = pr0/(pr3x20r20*pr20x30r30)
px30x21 = pr0/(pr3x21r20*pr20x30r30)
px31x20 = pr0/(pr3x20r20*pr20x31r30)
px31x21 = pr0/(pr3x21r20*pr20x31r30)

na.w = is.na(dat[,2]) | is.na(dat[,3])
a = (dat[,2] == 0) & (dat[,3] == 0 )
a[na.w] = FALSE
b = (dat[,2] == 0) & (dat[,3] == 1 )
b[na.w] = FALSE
c = (dat[,2] == 1) & (dat[,3] == 0 )
c[na.w] = FALSE
d = (dat[,2] == 1) & (dat[,3] == 1 )
d[na.w] = FALSE

weights = rep(1, nrow(dat))
weights[a] = px30x20
weights[b] = px31x20
weights[c] = px30x21
weights[d] = px31x21
weights[na.w] = NA

data = cbind(dat[,c(1,2)],weights)
data = test_wise_deletion(c(1,2,3), data)
binCItest_w(1, 2, c(), data[,c(3)], list(dm = data[,c(1,2)], adaptDF = FALSE))


# data_m2 = data_m2[!r, ]
# suffStat = list(data=data_m2,prt_m=prt_m,adaptDF=FALSE)
# binCItest.drw(7,15,c(),suffStat)
# #
