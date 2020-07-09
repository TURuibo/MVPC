#****************** Evaluation ******************
comp_com_td_mvpc<- function(exp="mar_cor",
                            n=1000, 
                            times=10,
                            p=20,
                            data_path="/Users/ruibo/Desktop/mvpc/mvpc-xyz/causality/R-proj/data"
){
  shd_comp= rep(0,times)
  shd_mvpc = rep(0,times)
  shd_test = rep(0,times)
  
  for(i in 1:times){
    print(i)
    # i=10
    load(paste(data_path,"/syn/",exp,n,"_",i,".rda",sep=""))
    
    data_mcar = data_all$data_mcar
    data_c = data_all$data_c
    suffStat = list(data_m = data_all$suffStat$data)
    myCPDAG=data_all$cpdag
    myDAG=data_all$dag
    
    res_mvpc<-mvpc(suffStat, gaussCItest_td,PermCCItest, alpha=0.01, p=p)
    res_pc<-pc(suffStat, gaussCItest_td, alpha=0.01, p=p)
    
    suffStat  = list(C = cor(data_c),n=length(data_c[,1]))
    res_ref<-pc(suffStat, gaussCItest, alpha=0.01, p=p)
    
    shd_comp[i] = shd(res_ref, myCPDAG)
    shd_mvpc[i] = shd(res_mvpc, myCPDAG)
    shd_test[i] = shd(res_pc, myCPDAG)  
  }
  return(list(shd_comp=shd_comp,shd_mvpc=shd_mvpc,shd_test=shd_test))
}

eva.detection<-function(prt1,prt2){
  count = 0
  for(i in 1:length(prt1)){
    if(length(prt1[i])!=length(prt2[[i]])){count = count + 1 }
    else if(prt1[i]!=prt2[[i]]){count = count + 1 }
  }
  return(count)
}

test_adj<-function(myCPDAG,res){
  ## Ajacency
  
  true_skel<-get_edge_pairs(myCPDAG)
  our_skel<-get_edge_pairs(res)
  
  tran_true_skel<-get_trans_edge_pairs(myCPDAG)
  tran_our_skel<-get_trans_edge_pairs(res)
  
  true_pair <- complex(real = true_skel[,1], imaginary = true_skel[,2])
  ours_pair <- complex(real = our_skel[,1], imaginary = our_skel[,2])
  
  tran_true_pair <- complex(real = tran_true_skel[,1], imaginary = tran_true_skel[,2])
  tran_ours_pair <- complex(real = tran_our_skel[,1], imaginary = tran_our_skel[,2])
  
  # all pairs:
  all_pair<-union(true_pair,tran_true_pair)
  all_pair_us <-union(ours_pair,tran_ours_pair)
  
  num_our<-length(all_pair_us)
  num_true<-length(all_pair)
  num_us_correct<-length(intersect(all_pair_us,all_pair))
  recall<-num_us_correct/num_true
  precision<-num_us_correct/num_our
  return(list(recall=recall, precision=precision))
}

get_edge_pairs<-function(G){
  g<-as(G,'matrix')
  g[lower.tri(g)]<-0
  g<-g==1
  which(g, arr.ind = TRUE)
}

get_trans_edge_pairs<-function(G){
  g<-t(as(G,'matrix'))
  g[lower.tri(g)]<-0
  g<-g==1
  which(g, arr.ind = TRUE)
}

compute_rp<-function(rp_){
  recall = c()
  precision = c()
  count = 1
  for(i in 1:length(rp_)){
    if(i == 6 | i ==12 | i==18 | i==24){next}
    recall[count] = rp_[[i]]$recall
    precision[count] = rp_[[i]]$precision
    count = count + 1
  }
  return(c(mean(recall), sd(recall), mean(precision),sd(precision)))
}

compute_f1<-function(rp_){
  recall = c()
  precision = c()
  f1 = c()
  count = 1
  for(i in 1:length(rp_)){
    if(i == 6 | i ==12 | i==18 | i==24){next}
    recall[count] = rp_[[i]]$recall
    precision[count] = rp_[[i]]$precision
    f1[count] = 2 * ( recall[count] * precision[count])/( recall[count] + precision[count])
    count = count + 1
  }
  return(c(mean(f1), sd(f1)))
}