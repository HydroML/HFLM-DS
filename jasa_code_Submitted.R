rm(list = ls())

#import functions
library(glmnet)
library(Matrix)
library(ggplot2)
library(lubridate)
library(caTools)
library(forecast)
library(zoo)
library(ranger)
library(missRanger)
library(GPfit)
library(lattice)
library(RColorBrewer)
library(segmented)

#get order and location of all coefficient function nodes
get_index_mat <- function(M = 4,max_lag=3) {
  K <- (M+1)*max_lag # total number of nodes
  what<-expand.grid(0:(max_lag-1),0:M)
  data.frame("Node" = 1:K, "Y_coord" = what$Var1, "X_coord" = what$Var2)
}


#build smoothness regularization matrices
get_D <- function(M,max_lag) {
  ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  K <- nrow(ind_mat)
  num_pairs <- (M+1)*(max_lag-1)
  
  # First we build the D_h matrix
  D_v <- Matrix(0, nrow = num_pairs, ncol = K)
  ind_diff_v <- which(diff(ind_mat[,3]) == 0)
  
  D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v])] <- -1
  D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v]+1)] <- 1
  D_v<-rbind(D_v,sparseMatrix(i=1:(M+1),j=which(ind_mat$Y_coord == max(ind_mat$Y_coord)), x=rep(1,(M+1)),dims = c(M+1,K) )) # add vanishing boundary to top
  
  # Secondly we have the D_v matrix
  num_pairs <- (M)*(max_lag)
  D_h <- Matrix(0, nrow = num_pairs, ncol = K)
  ind_mat2 <- ind_mat[order(ind_mat[,2]),]
  ind_diff_h<- which(diff(ind_mat2[,2]) == 0)
  D_h[cbind(1:num_pairs, sort(ind_mat2[ind_diff_h,1]))] <- -1
  D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h+1,1]))] <- 1
  D_h<-rbind(D_h,sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == min(ind_mat$X_coord)), x=rep(1,max_lag),dims = c(max_lag,K) ) + sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == max(ind_mat$X_coord)), x=rep(-1,max_lag),dims = c(max_lag,K) )) # add periodicity condition
  
  return(list(D_h,D_v))
}


get_grp_normsH <- function(bvec, M,max_lag) {
  grps <- as.factor(rep(0:M, rep(max_lag,M+1)))
  norms <- tapply(bvec, grps, FUN = function(x){
    cumsum(as.numeric(x[length(x):1])^2)[length(x):1]
  })
  
  as.numeric(do.call("c", norms))
}


#plotting beta(s,t)
plot_bvec_better<-function(bvec,M,max_lag){
  index_dat<-get_index_mat(M,max_lag)
  delta<-compute_delta(bvec,M,max_lag = max_lag)
  
  mat<-matrix(0,ncol = M+1,nrow = max_lag)
  mat[cbind(index_dat$Y_coord+1,index_dat$X_coord+1)]<-bvec
  mat<-matrix(mat[1:min(ceiling(max(delta)/10)*10,max_lag),],ncol = M+1)
  range_p<-max(abs(mat[mat>0]),2e-15)
  range_n<-max(abs(mat[mat<0]),2e-15)
  max_range<-max(abs(mat))
  tot_range<-range_p+ range_n
  totalcols<-35
  rangepercol<-tot_range/totalcols
  numblues<-ceiling(range_p/rangepercol)
  numreds<-totalcols- numblues +1
  
  coul <- colorRampPalette(brewer.pal(8, "Blues"))(max(numblues+2,ceiling(totalcols/2)))
  coul2 <- colorRampPalette(brewer.pal(8, "Reds"))(max(numreds+2,ceiling(totalcols/2)))
  colorvec<-c(coul2[(numreds-1+2):2],"#FFFFFF",coul[3:(numblues-1+3)])
  
  rownames(mat)<-0:(nrow(mat)-1)
  mat<-setNames(reshape2::melt(t(mat)),c("x","y","z"))
  
  levelplot(z~x*y,mat,col.regions = colorvec,at=sort(c(seq(from=-1e-15,to=-range_n,length.out=numreds+1 ),seq(from=1e-15,to=range_p,length.out=numblues+1))),
            ylab=list("Lag (Days)",cex=1.4),xlim = c(0.5,M+1+0.5),aspect = "xy",xlab=list("Day of the year" ,cex=1.4),colorkey=list(axis.text=list(cex=1.2)),
            scales=list(x=list(cex=1.2),y=list(cex=1.2)))
}


# form matrix of input var
get_x<-function(x_ts,M,max_lag,first_day){
  K=(M+1)*max_lag
  nrowX<-(length(x_ts)-max_lag+1)
  i_vec<-rep(0, max_lag*nrowX)
  j_vec<-rep(0, max_lag*nrowX)
  x_vec<-rep(0, max_lag*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+max_lag)]<-rep(i,max_lag)
    j_vec[(indexCounter+1):(indexCounter+max_lag)]<-(((indexCounter + (first_day-1)*max_lag) %% K):((indexCounter+max_lag-1 + (first_day-1)*max_lag) %% K))+1
    x_vec[(indexCounter+1):(indexCounter+max_lag)]<-x_ts[(i+max_lag-1):(i)]
    indexCounter<-indexCounter+max_lag
  }
  
  sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, K))
}

compute_delta<-function(bvec,M,max_lag){
  deltavec<-rep(0,M+1)
  
  index_dat<-get_index_mat(M,max_lag)
  
  for(i in 0:M){
    curobs<-bvec[index_dat$X_coord==i]
    curobs<-curobs[length(curobs):1]
    deltavec[i+1]<-max_lag-max(which(cumsum(curobs)==0),0)
  }
  deltavec
}


#optimize w_h and w_v hyperparameters with bayesian optimization
optimize_hp<-function(hp_range,XTX,hTh,vTv,XTy,Xval,ytrain,yval,n_hp,To_Zero){
  n_init<-30 #get starting points
  R2<-rep(0,round(n_init))
  
  max_w_h<-hp_range[[1]][2]
  max_w_v<-hp_range[[2]][2]
  
  min_w_h<-hp_range[[1]][1]
  min_w_v<-hp_range[[2]][1]
  
  w_h<-runif(length(R2),min=min_w_h,max=max_w_h)
  w_v<-runif(length(R2),min=min_w_v,max=max_w_v)
  
  
  for(h in 1:length(R2)){
    pred_beta<-rep(0,ncol(Xval))
    pred_beta[!To_Zero] <- Matrix::solve(XTX+ exp(w_h[h])*hTh + exp(w_v[h])*vTv, XTy, sparse = TRUE)
    Y_hat_val=Xval %*% pred_beta
    R2[h]<-1-sum((yval- Y_hat_val)^2)/sum((mean(yval)-yval)^2)
    R2[h]<-max(R2[h],-1)
    print(paste0("This is R2[h]: ",R2[h], " , and h is: ",h," and w_h is: ",w_h[h], " and w_v is: ", w_v[h]))
  }
  
  while(h<n_hp){ #start gaussian process
    
    w_h_scaled<-(w_h- min_w_h)/(max_w_h- min_w_h)
    w_v_scaled<-(w_v- min_w_v)/(max_w_v- min_w_v)

    mod<-GP_fit(X=data.frame(w_h=w_h_scaled, w_v=w_v_scaled),Y=R2)
    
    w_h_long<-runif(1000,min=min_w_h,max=max_w_h)
    w_v_long<-runif(1000,min=min_w_v,max=max_w_v)
    
    w_h_long_scaled<-(w_h_long- min_w_h)/(max_w_h- min_w_h)
    w_v_long_scaled<-(w_v_long- min_w_v)/(max_w_v- min_w_v)
    
    pred <- predict.GP(mod, xnew = data.frame(w_h=w_h_long_scaled,w_v=w_v_long_scaled))
    mu <- pred$Y_hat
    sigma <- sqrt(pred$MSE)
    Z <- (mu - max(pred$Y_hat))/sigma
    expected_imp <- sigma*(Z  * pnorm(Z) + dnorm(Z))
    expected_imp[is.na(expected_imp)]<-0
    
    nextone<-which.max(expected_imp)
    w_h<-c(w_h,w_h_long[nextone])
    w_v<-c(w_v,w_v_long[nextone])
    
    h=h+1
    pred_beta<-rep(0,ncol(Xval))
    pred_beta[!To_Zero] <- Matrix::solve(XTX+ exp(w_h[h])*hTh + exp(w_v[h])*vTv, XTy, sparse = TRUE)
    Y_hat_val=Xval %*% pred_beta
    curR2<- 1-sum((yval- Y_hat_val)^2)/sum((mean(yval)-yval)^2)
    curR2<-max(curR2,-1)
    R2<-c(R2,curR2)
    print(paste0("Gaussian process: This is R2[h]: ",R2[h], " , and h is: ",h," and w_h is: ",w_h[h], " and w_v is: ", w_v[h]))
  }
  list(w_h,w_v,R2)
}

#find elbow point of monotone decreasing plot
get_elbow<-function(qVec,R2Vec){
  qVecScaled<-(qVec-min(qVec))/(max(qVec)-min(qVec))
  R2VecScaled<- (R2Vec-min(R2Vec))/(max(R2Vec)-min(R2Vec))
  kneedledist<-R2VecScaled - (1-qVecScaled) 
  kneedledist
}


compute_beta_R2<-function(est,true){
  1-sum((est-true)^2)/sum((true-mean(true))^2)
}


optimize_q<-function(grpnorms,Nq=500,Xmat_train,Xmat_val,y_train,y_val,hMat,vMat,w_h_best,w_v_best){
  q_ngn<-seq(from=0,to=1-365/dim(Xmat_train)[2],length.out=100) #important to select this carefully
  R2<-rep(0,length(q_ngn))
  for(h in 1:length(R2)){
    To_Zero<-grpnorms<quantile(grpnorms,q_ngn[h])
    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(hMat[,!To_Zero])
    vProd_sparse<-crossprod(vMat[,!To_Zero])
    
    pred_beta<-rep(0,length(pred_beta))
    pred_beta[!To_Zero]<- Matrix::solve(XXProd_sparse+ w_h_best*hProd_sparse + w_v_best*vProd_sparse, XyProd_sparse, sparse = TRUE)
    Y_hat_val=Xmat_val %*% pred_beta
    R2[h]<-1-sum((y_val- Y_hat_val)^2)/sum((y_val-mean(y_val))^2)
    print(paste0("This is sparse R2[h]: ",R2[h], " , and h is: ",h))
  }
  
  whichbestq<-which.max(get_elbow(q_ngn,R2))
  new_q<-q_ngn[q_ngn<= q_ngn[whichbestq]]
  new_R2<-R2[q_ngn<= q_ngn[whichbestq]]
  
  whichbestq2<-which.max(get_elbow(new_q,new_R2))
  finalq2<- new_q[whichbestq2]
  finalq1<- q_ngn[whichbestq]
  searchRange<-c(max(0,finalq2-0.05),min(1,finalq1+0.05))
  
  q_ngn2<-runif(Nq- length(R2),min=searchRange[1],max=searchRange[2]) #important to select this carefully
  R22<-rep(0,length(q_ngn2))
  for(h in 1:length(R22)){
    To_Zero<-grpnorms<quantile(grpnorms,q_ngn2[h])

    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(hMat[,!To_Zero])
    vProd_sparse<-crossprod(vMat[,!To_Zero])
    
    pred_beta<-rep(0,length(pred_beta))
    pred_beta[!To_Zero]<- Matrix::solve(XXProd_sparse+ w_h_best*hProd_sparse + w_v_best*vProd_sparse, XyProd_sparse, sparse = TRUE)
    Y_hat_val=Xmat_val %*% pred_beta
    R22[h]<-1-sum((y_val- Y_hat_val)^2)/sum((y_val-mean(y_val))^2)
    if(h%%50==0) print(h/length(R22))
  }
  q_all<-c(q_ngn,q_ngn2)
  R2_all<-c(R2,R22)
  
  return(list(q_all,R2_all))
}

get_x_kir<-function(x_ts,y_ts,M,max_lag,first_day,ars){
  differ<-length(x_ts)- length(y_ts)
  max_lag=max_lag+ars
  first_day= (first_day + ars) %% (M+1)
  nrowX<- length(x_ts)-max(max_lag-1,ars)
  
  K=(M+1)*(max_lag)
  i_vec<-rep(0, max_lag*nrowX)
  j_vec<-rep(0, max_lag*nrowX)
  x_vec<-rep(0, max_lag*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+max_lag)]<-rep(i,max_lag)
    j_vec[(indexCounter+1):(indexCounter+max_lag)]<-(((indexCounter + first_day*max_lag) %% K):((indexCounter+max_lag-1 + first_day*max_lag) %% K))+1
    x_vec[(indexCounter+1):(indexCounter+max_lag)]<-x_ts[(i+max_lag-1):i]
    indexCounter<-indexCounter+max_lag
  }
  
  AllX=sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, K))
  
  
  i_vec<-rep(0, ars*nrowX)
  j_vec<-rep(0, ars*nrowX)
  x_vec<-rep(0, ars*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+ars)]<-rep(i,ars)
    j_vec[(indexCounter+1):(indexCounter+ars)]<-1:ars
    x_vec[(indexCounter+1):(indexCounter+ars)]<-y_ts[(i+ars-1):(i)]
    indexCounter<-indexCounter+ars
  }
  
  AllY=sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, ars))
  
  Allboth<-cbind2(AllX,AllY)
  Allboth
}


get_truebeta<-function(pred_beta,indexmat,ars){
  y_cords<-unique(indexmat$Y_coord)
  big_indexmat<-get_index_mat(M=max(indexmat$X_coord),max_lag = max(y_cords+1)+ars)
  big_indexmat$val<-pred_beta[1:(length(pred_beta)-ars)]
  indexmat<-merge(indexmat,big_indexmat,by=c("X_coord","Y_coord"),all.x = T)
  indexmat<-indexmat[order(indexmat$Node.x),]
  arcs<-pred_beta[(length(pred_beta)-ars+1):length(pred_beta)]
  real_beta<-indexmat$val
  
  for(y_cord in y_cords[2:length(y_cords)]){
    aligned_beta<-real_beta[which(indexmat$Y_coord==y_cord-1)]
    aligned_beta<-aligned_beta[c(length(aligned_beta),1:(length(aligned_beta)-1))]
    real_beta[indexmat$Y_coord==y_cord]<-real_beta[indexmat$Y_coord==y_cord]+ arcs[1]* aligned_beta
    if(y_cord>1 & length(arcs)>1){
      aligned_beta<-real_beta[which(indexmat$Y_coord==y_cord-2)]
      aligned_beta<-aligned_beta[c(length(aligned_beta)-1,length(aligned_beta),1:(length(aligned_beta)-2))]
      real_beta[indexmat$Y_coord==y_cord]<-real_beta[indexmat$Y_coord==y_cord]+ arcs[2]* aligned_beta
    }
  }
  real_beta
}


#build smoothness regularization matrices
get_D_kir <- function(M,max_lag,est_arcoefs) {
  ars=length(est_arcoefs)
  max_lag=max_lag+ars
  ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  K <- nrow(ind_mat)
  num_pairs <- (M+1)*(max_lag-1)
  
  D_v<-bandSparse(K,k=c(0,1),diagonals =list(rep(-1,K), rep(c(rep(1,max_lag-1),0),M+1)[-K]))
  #D_v<-bandSparse(K,k=c(-1,0,1),diagonals =list(rep(est_arcoefs[1]^2+est_arcoefs[2]-est_arcoefs[1],K-1),rep(est_arcoefs[1]-1,K), rep(c(rep(1,max_lag-1),0),M+1)[-K]))
  extracols<-Matrix(0,nrow = nrow(D_v),ncol = ars)
  D_v<-cbind2(D_v,extracols)
  extrarows<-Matrix(0,nrow = ars,ncol = ncol(D_v))
  D_v<-rbind2(D_v,extrarows)
  
  
  #D_v[1,1]<- est_arcoefs[1]-1
  #D_v[2,1]<- est_arcoefs[1]^2+ est_arcoefs[2]- est_arcoefs[1]
  #D_v[2,2]<- est_arcoefs[1]-1
  
  # max_lag=max_lag+ars
  # ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  # K <- nrow(ind_mat)
  # Secondly we have the D_v matrix
  if(M>0){
    num_pairs <- (M)*(max_lag)
    D_h <- Matrix(0, nrow = num_pairs, ncol = K)
    ind_mat2 <- ind_mat[order(ind_mat[,2]),]
    ind_diff_h<- which(diff(ind_mat2[,2]) == 0)
    D_h[cbind(1:num_pairs, sort(ind_mat2[ind_diff_h,1]))] <- -1
    D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h+1,1]))] <- 1
    D_h<-rbind(D_h,sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == min(ind_mat$X_coord)), x=rep(1,max_lag),dims = c(max_lag,K) ) + sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == max(ind_mat$X_coord)), x=rep(-1,max_lag),dims = c(max_lag,K) )) # add periodicity condition
  } else{
    D_h= sparseMatrix(i=1,j=1,x=0,dims = c(1,K))
  }
  extracols<-Matrix(0,nrow = nrow(D_h),ncol = ars)
  D_h<-cbind2(D_h,extracols)
  extrarows<-Matrix(0,nrow = ars,ncol = ncol(D_h))
  D_h<-rbind2(D_h,extrarows)
  
  return(list(D_h,D_v))
  
}

Zero_out_pred_beta<- function(b,d){
  newb<-b
  M<-length(d)-1
  max_lag<-length(b)/(M+1)
  ind_mat<-get_index_mat(M,max_lag)
  for(m in 1:(M+1)){
    newb[ind_mat$Y_coord>(d[m]-1) & ind_mat$X_coord==(m-1)]<-0
  }
  newb
}


###########################################################################################################
######################################   process data   ###################################################
############################################################################################################
setwd("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_aoas/")

#do 1344 (British Columbia, Canada), 1759 (Florida, United States)
gridcode<-1759

simu_x<-readRDS(file=paste0("grid",gridcode,"_rain.rds"))
Y<-readRDS(file = paste0("grid",gridcode,"_streamflow.rds"))

M = ncol(simu_x)-1
max_lag<-150
K = (M+1)*max_lag
n_hp<-65
ars=2

#regularization
RegMat<-get_D(M,max_lag)
hProd<-crossprod(RegMat[[1]])
vProd<-crossprod(RegMat[[2]])

x<-as.numeric(t(simu_x))
y<-as.numeric(t(Y))
y<-y[max_lag:length(y)]

trainfrac<-0.6
valfrac<-0.2
testfrac<-0.2


y_train<-y[1:round(length(y)*trainfrac)]
y_val<-y[(round(length(y)*trainfrac)+1):round(length(y)*(trainfrac+valfrac))]
y_test<-y[(round(length(y)*(trainfrac+valfrac))+1):length(y)]
y_full<-y[1:round(length(y)*(trainfrac+valfrac))]

x_train<-x[1:(max_lag+length(y_train)-1)]
x_val<-x[(length(y_train)+1):(max_lag+length(y_train)+length(y_val)-1)]
x_test<-x[(length(y_train)+length(y_val)+1):length(x)]
x_full<-x[1:(max_lag+length(y_train)+length(y_val)-1)]

first_day<-max_lag
first_day_train<-max_lag
first_day_val<-(length(y_train)+max_lag) %% 365
first_day_test<-(length(y_train)+length(y_val)+max_lag) %% 365
first_day_full<-max_lag

Xmat_train<-get_x(x_train,M,max_lag,first_day_train)
Xmat_val<-get_x(x_val,M,max_lag,first_day = first_day_val)
Xmat_test<-get_x(x_test,M,max_lag,first_day = first_day_test)
Xmat_full<-get_x(x_full,M,max_lag,first_day = first_day_full)
Xmat<-get_x(x,M,max_lag,first_day = first_day)

XXProd_train<-crossprod(Xmat_train)
XXProd_val<-crossprod(Xmat_val)
XXProd_test<-crossprod(Xmat_test)
XXProd_full<-crossprod(Xmat_full)
XXProd<-crossprod(Xmat)

XyProd_train<-crossprod(Xmat_train,matrix(y_train,ncol = 1))
XyProd_val<-crossprod(Xmat_val,matrix(y_val,ncol = 1))
XyProd_full<-crossprod(Xmat_full,matrix(y_full,ncol = 1))
XyProd_test<-crossprod(Xmat_test,matrix(y_test,ncol = 1))
XyProd<-crossprod(Xmat,matrix(y,ncol = 1))



############################################################################################################################################3
################################################    test on real streamflow data    #########################################################
#############################################################################################################################################

To_Zero<-rep(F,K)
result<-optimize_hp(hp_range=list(c(10,20),c(-5,15)),XTX=XXProd_train,hTh=hProd,vTv=vProd,XTy=XyProd_train,Xval=Xmat_val
                    ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)


best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
plot(result[[1]],result[[3]])
plot(result[[2]],result[[3]])


pred_beta <- as.numeric(Matrix::solve(XXProd_full+ w_h[best]*hProd + w_v[best]*vProd, XyProd_full, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

#detect where low betas are
grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
plot_bvec_better(grpnorms,M,max_lag)


#optimize quantile for thresholding
q_result<-optimize_q(grpnorms = grpnorms,Nq=500,Xmat_train = Xmat_full,Xmat_val = Xmat_full,
                     y_train = y_full,y_val = y_full,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
q_ngn<-q_result[[1]]
R2<-q_result[[2]]
plot(q_ngn,R2)
distances1<-get_elbow(q_ngn,R2)
whichbestq<-which.max(distances1)
points(q_ngn[whichbestq],R2[whichbestq],col="red")
finalq<-q_ngn[whichbestq]
print(finalq)
To_Zero<-grpnorms<quantile(grpnorms,finalq)


big_indexmat<-get_index_mat(M=M,max_lag = ars+ max_lag)
small_indexmat<-get_index_mat(M=M,max_lag = max_lag)
To_Zero_big<-rep(F,nrow(big_indexmat))
for(m in 1:(M+1)){
  tot_in_col<-sum(To_Zero[small_indexmat$X_coord==(m-1)])
  To_Zero_big[big_indexmat$Y_coord> max(big_indexmat$Y_coord)-tot_in_col & big_indexmat$X_coord==(m-1)]<- T
}
To_Zero<-c(To_Zero_big,rep(F,ars))

Xmat<-get_x_kir(x,y,M,max_lag,first_day,ars = ars)
Xmat_train<- Xmat[1:(nrow(Xmat_train)),]
Xmat_val<-Xmat[(nrow(Xmat_train)+1):nrow(Xmat_full),]
Xmat_full<-rbind2(Xmat_train,Xmat_val)

y_train<-y[(ars+1):(nrow(Xmat_train)+ars)]
y_val<-y[(nrow(Xmat_train)+ars+1):(nrow(Xmat_full)+ars)]
y_full<-c(y_train,y_val)
y<-y[(ars+1):length(y)]

RegMat<-get_D_kir(M,max_lag,est_arcoefs = rep(0,ars))
XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

result<-optimize_hp(hp_range=list(c(10,20),c(-5,15)),XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_val
                    ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)

best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
#refit model and update beta
XyProd_sparse<-crossprod(Xmat_full[,!To_Zero],matrix(y_full,ncol = 1))
XXProd_sparse<-crossprod(Xmat_full[,!To_Zero])


pred_beta<-rep(0,ncol(Xmat))
pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*hProd_sparse + w_v[best]*vProd_sparse, XyProd_sparse, sparse = TRUE))
plot_bvec_better(pred_beta[1:(length(pred_beta)-ars)],M,max_lag+ars)
est_delta<-compute_delta(pred_beta[1:(length(pred_beta)-ars)],M,max_lag = max_lag+ars)-ars

indexmat_new<-get_index_mat(M,max_lag)
pred_beta<-get_truebeta(pred_beta,indexmat_new,ars=ars)
plot_bvec_better(pred_beta,M,max_lag)
pred_beta<-Zero_out_pred_beta(b=pred_beta,d=est_delta)

Y_hat_test=Xmat_test %*% pred_beta
1-sum((y_test- Y_hat_test)^2)/sum((y_test- mean(y_full))^2)
plot_bvec_better(pred_beta,M,max_lag)



############################################################################################################################################3
################################################    final prediction    #########################################################
#############################################################################################################################################


set.seed(111)
To_Zero<-rep(F,K)
result<-optimize_hp(hp_range=list(c(10,20),c(-5,15)),XTX=XXProd_full,hTh=hProd,vTv=vProd,XTy=XyProd_full,Xval=Xmat_test
                    ,yval=y_test,ytrain = y_full,n_hp=n_hp,To_Zero = To_Zero)


best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
plot(result[[1]],result[[3]])
plot(result[[2]],result[[3]])

pred_beta <- as.numeric(Matrix::solve(XXProd+ w_h[best]*hProd + w_v[best]*vProd, XyProd, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

#detect where low betas are
grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
plot_bvec_better(grpnorms,M,max_lag)


#optimize quantile for thresholding
q_result<-optimize_q(grpnorms = grpnorms,Nq=500,Xmat_train = Xmat,Xmat_val = Xmat,
                     y_train = y,y_val = y,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
q_ngn<-q_result[[1]]
R2<-q_result[[2]]
plot(q_ngn,R2)
distances1<-get_elbow(q_ngn,R2)
whichbestq<-which.max(distances1)
points(q_ngn[whichbestq],R2[whichbestq],col="red")
finalq<-q_ngn[whichbestq]
print(finalq)
To_Zero<-grpnorms<quantile(grpnorms,finalq)



big_indexmat<-get_index_mat(M=M,max_lag = ars+ max_lag)
small_indexmat<-get_index_mat(M=M,max_lag = max_lag)
To_Zero_big<-rep(F,nrow(big_indexmat))
for(m in 1:(M+1)){
  tot_in_col<-sum(To_Zero[small_indexmat$X_coord==(m-1)])
  To_Zero_big[big_indexmat$Y_coord> max(big_indexmat$Y_coord)-tot_in_col & big_indexmat$X_coord==(m-1)]<- T
}
To_Zero<-c(To_Zero_big,rep(F,ars))

Xmat<-get_x_kir(x,y,M,max_lag,first_day,ars = ars)
Xmat_train<- Xmat[1:(nrow(Xmat_train)),]
Xmat_val<-Xmat[(nrow(Xmat_train)+1):nrow(Xmat_full),]
Xmat_test<-Xmat[(nrow(Xmat_full)+1):nrow(Xmat),]
Xmat_full<-rbind2(Xmat_train,Xmat_val)

y_train<-y[(ars+1):(nrow(Xmat_train)+ars)]
y_val<-y[(nrow(Xmat_train)+ars+1):(nrow(Xmat_full)+ars)]
y_test<-y[(nrow(Xmat_full)+ars+1):length(y)]
y_full<-c(y_train,y_val)
y<-y[(ars+1):length(y)]

RegMat<-get_D_kir(M,max_lag,est_arcoefs = rep(0,ars))
XyProd_sparse<-crossprod(Xmat_full[,!To_Zero],matrix(y_full,ncol = 1))
XXProd_sparse<-crossprod(Xmat_full[,!To_Zero])
hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

result<-optimize_hp(hp_range=list(c(10,20),c(-5,15)),XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_test
                    ,yval=y_test,ytrain = y_full,n_hp=n_hp,To_Zero = To_Zero)

best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
#refit model and update beta
XyProd_sparse<-crossprod(Xmat[,!To_Zero],matrix(y,ncol = 1))
XXProd_sparse<-crossprod(Xmat[,!To_Zero])


pred_beta<-rep(0,ncol(Xmat))
pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*hProd_sparse + w_v[best]*vProd_sparse, XyProd_sparse, sparse = TRUE))
plot_bvec_better(pred_beta[1:(length(pred_beta)-ars)],M,max_lag+ars)
est_delta<-compute_delta(pred_beta[1:(length(pred_beta)-ars)],M,max_lag = max_lag+ars)-ars

indexmat_new<-get_index_mat(M,max_lag)
pred_beta<-get_truebeta(pred_beta,indexmat_new,ars=ars)
plot_bvec_better(pred_beta,M,max_lag)
pred_beta<-Zero_out_pred_beta(b=pred_beta,d=est_delta)
plot_bvec_better(pred_beta,M,max_lag)

jpeg(paste0("groundtruth",gridcode,"_2024.jpeg"),width = 8, height = 6.5, units = 'in',res = 500)
plot_bvec_better(pred_beta,M,max_lag)
dev.off()

write.csv(data.frame(pred_beta=pred_beta),paste0("groundtruth",gridcode,"_2024.csv"),row.names = F)



############################################################################################################################################3
################################################    run simulation study       #########################################################
#############################################################################################################################################
setwd("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_aoas/")
gridcode<-1344
desiredR2<-0.4


simu_x<-readRDS(file=paste0("grid",gridcode,"_rain.rds"))
x<-as.numeric(t(simu_x))

M = ncol(simu_x)-1
max_lag<-150
K = (M+1)*max_lag
trainfrac<-0.8
valfrac<-0.2

bvec<-read.csv(paste0("groundtruth",gridcode,"_2024.csv"))
bvec<-bvec$pred_beta
true_delta<-compute_delta(bvec,M,max_lag)

n_hp<-65
nsim<-30
beta_R2<-rep(0,nsim)
delta_bias<-rep(0,nsim)
delta_cor<-rep(0,nsim)
arcfs<-c(1.5,-0.52)
ars=2
hp_range<-list(c(8,24),c(-5,15))
for(i in 1:nsim){
  RegMat<-get_D(M,max_lag)
  hProd<-crossprod(RegMat[[1]])
  vProd<-crossprod(RegMat[[2]])
  
  first_day<- max_lag %% (M+1)
  Xmat<-get_x(x,M,max_lag,first_day)
  Xmat_train<-Xmat[1:round(nrow(Xmat)*trainfrac),]
  Xmat_val<-Xmat[(1+round(nrow(Xmat)*trainfrac)):nrow(Xmat),]
  
  y=as.numeric(Xmat %*% bvec)
  noiseVar<-(var(y)-desiredR2*var(y))/desiredR2
  set.seed(i)
  noise<-as.numeric(arima.sim(list(order=c(length(arcfs),0,0),ar=arcfs),n=length(y),sd=1))
  noise<- (noise - mean(noise))* sqrt(noiseVar)/sd(noise)
  y=y+noise
  y_train<-y[1:nrow(Xmat_train)]
  y_val<-y[(nrow(Xmat_train)+1):length(y)]
  
  XyProd_train<-crossprod(Xmat_train,matrix(y_train,ncol = 1))
  XyProd_val<-crossprod(Xmat_val,matrix(y_val,ncol = 1))
  XyProd<-crossprod(Xmat,matrix(y,ncol = 1))
  XXProd_train<-crossprod(Xmat_train)
  XXProd_val<-crossprod(Xmat_val)
  XXProd<-crossprod(Xmat)
  
  
  To_Zero<-rep(F,K)
  result<-optimize_hp(hp_range=hp_range,XTX=XXProd_train,hTh=hProd,vTv=vProd,XTy=XyProd_train,Xval=Xmat_val
                      ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)
  
  best<-which.max(result[[3]])
  bestR2<-max(result[[3]])
  w_h<-exp(result[[1]])
  w_v<-exp(result[[2]])
  pred_beta <- as.numeric(Matrix::solve(XXProd+ w_h[best]*hProd + w_v[best]*vProd, XyProd, sparse = TRUE))
  #plot_bvec_better(pred_beta,M,max_lag)
  print(paste0("non sparese w_h and w_v are: ",round(w_h[best],3)," and ",round(w_v[best],3)))
  
  #detect where low betas are
  grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
  #plot_bvec_better(grpnorms,M,max_lag)
  
  q_result<-optimize_q(grpnorms = grpnorms,Nq=500,Xmat_train = Xmat,Xmat_val = Xmat,
                       y_train = y,y_val = y,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
  q_ngn<-q_result[[1]]
  R2<-q_result[[2]]
  plot(q_ngn,R2)
  distances1<-get_elbow(q_ngn,R2)
  whichbestq<-which.max(distances1)
  points(q_ngn[whichbestq],R2[whichbestq],col="red")
  finalq<-q_ngn[whichbestq]
  print(finalq)
  To_Zero<-grpnorms<quantile(grpnorms,finalq)
  
  #redo smoothness
  big_indexmat<-get_index_mat(M=M,max_lag = ars+ max_lag)
  small_indexmat<-get_index_mat(M=M,max_lag = max_lag)
  To_Zero_big<-rep(F,nrow(big_indexmat))
  for(m in 1:(M+1)){
    tot_in_col<-sum(To_Zero[small_indexmat$X_coord==(m-1)])
    To_Zero_big[big_indexmat$Y_coord> max(big_indexmat$Y_coord)-tot_in_col & big_indexmat$X_coord==(m-1)]<- T
  }
  To_Zero<-c(To_Zero_big,rep(F,ars))
  
  
  Xmat<-get_x_kir(x,y,M,max_lag,first_day,ars = ars)
  Xmat_train<- Xmat[1:nrow(Xmat_train),]
  Xmat_val<-Xmat[(nrow(Xmat_train)+1):nrow(Xmat),]
  
  y_train<-y[(ars+1):(nrow(Xmat_train)+ars)]
  y_val<-y[(nrow(Xmat_train)+ars+1):length(y)]
  y<-y[(ars+1):length(y)]
  
  RegMat<-get_D_kir(M,max_lag,est_arcoefs = arcfs[1:ars])
  XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
  XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
  hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
  vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])
  
  result<-optimize_hp(hp_range=hp_range,XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_val
                      ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)
  
  best<-which.max(result[[3]])
  w_h<-exp(result[[1]])
  w_v<-exp(result[[2]])
  #refit model and update beta
  XyProd_sparse<-crossprod(Xmat[,!To_Zero],matrix(y,ncol = 1))
  XXProd_sparse<-crossprod(Xmat[,!To_Zero])
  hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
  vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])
  
  
  pred_beta<-rep(0,ncol(Xmat))
  pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*hProd_sparse + w_v[best]*vProd_sparse, XyProd_sparse, sparse = TRUE))
  plot_bvec_better(pred_beta[1:(length(pred_beta)-ars)],M,max_lag+ars)
  est_delta<-compute_delta(pred_beta[1:(length(pred_beta)-ars)],M,max_lag = max_lag+ars)-ars
  
  indexmat_new<-get_index_mat(M,max_lag)
  pred_beta<-get_truebeta(pred_beta,indexmat_new,ars=ars)
  plot_bvec_better(pred_beta,M,max_lag)
  pred_beta<-Zero_out_pred_beta(b=pred_beta,d=est_delta)
  plot_bvec_better(pred_beta,M,max_lag)
  
  beta_R2[i]<-compute_beta_R2(pred_beta,bvec)
  delta_bias[i]<-mean(est_delta-true_delta)
  delta_cor[i]<-cor(est_delta,true_delta)
  
  write.csv(data.frame(bR2=beta_R2,db=delta_bias,dc=delta_cor),
            paste0("SimStudy",gridcode,"_R2_",desiredR2,arcfs[1],arcfs[2],"_endkir.csv"),row.names = F)
  print(i)
}


SimStudy1344_R2_0.4 <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_aoas/SimStudy1344_R2_0.4.csv")
mean(SimStudy1344_R2_0.4$bR2)
sd(SimStudy1344_R2_0.4$bR2)
mean(SimStudy1344_R2_0.4$db)
sd(SimStudy1344_R2_0.4$db)
mean(SimStudy1344_R2_0.4$dc)
sd(SimStudy1344_R2_0.4$dc)

SimStudy1344_R2_0.8 <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_aoas/SimStudy1344_R2_0.8.csv")
mean(SimStudy1344_R2_0.8$bR2)
sd(SimStudy1344_R2_0.8$bR2)
mean(SimStudy1344_R2_0.8$db)
sd(SimStudy1344_R2_0.8$db)
mean(SimStudy1344_R2_0.8$dc)
sd(SimStudy1344_R2_0.8$dc)


SimStudy1759_R2_0.4 <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_aoas/SimStudy1759_R2_0.4.csv")
mean(SimStudy1759_R2_0.4$bR2)
sd(SimStudy1759_R2_0.4$bR2)
mean(SimStudy1759_R2_0.4$db)
sd(SimStudy1759_R2_0.4$db)
mean(SimStudy1759_R2_0.4$dc)
sd(SimStudy1759_R2_0.4$dc)

SimStudy1759_R2_0.8 <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_aoas/SimStudy1759_R2_0.8.csv")
mean(SimStudy1759_R2_0.8$bR2)
sd(SimStudy1759_R2_0.8$bR2)
mean(SimStudy1759_R2_0.8$db)
sd(SimStudy1759_R2_0.8$db)
mean(SimStudy1759_R2_0.8$dc)
sd(SimStudy1759_R2_0.8$dc)



####################################################################################3333
# plot average and year streamflow/precip
###########################################################################################
rm(list = ls())
setwd("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/Code_and_data_jasa/")
library(lubridate)
gridcode<-1344
cur_data <- read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/WaterInputOutputData_with_HBV_sim_streamflow/WaterInputOutputData_with_HBV_sim_streamflow/water_input_output_data_with_HBV_sim_streamflow_",gridcode,".csv"))
cur_data$Observed_Streamflow_mm.d[is.na(cur_data$Observed_Streamflow_mm.d)]<-cur_data$Simulated_Streamflow_mm.d[is.na(cur_data$Observed_Streamflow_mm.d)]
datevec<-seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day")


simu_x<-cur_data$Forcing_Precipitation_mm.d
simu_x<-simu_x[!(month(datevec)==2 & day(datevec)==29)]
simu_x<-matrix(simu_x,ncol = 365,byrow = T)
Y<-cur_data$Observed_Streamflow_mm.d
Y<-Y[!(month(datevec)==2 & day(datevec)==29)]
Y<-matrix(Y,ncol = 365,byrow = T)
Y[Y<0]<-0


jpeg(paste0("AvgRainStream",gridcode,".jpeg"),width = 6, height = 4, units = 'in',res = 500)
par(mar = c(5, 5, 2, 4.5))
plot(colMeans(simu_x),type="l",xlab = "Day of the year",ylab="Average rainfall (mm)",lwd=1.5)
par(new = TRUE)
plot(colMeans(Y),type="l",col="steelblue2",ylab = "",xlab = "",axes = FALSE, bty = "n",lwd=1.6)
axis(side=4, at = pretty(range(colMeans(Y))),col="steelblue2",col.axis="steelblue2")
mtext("Average streamflow (mm)", side=4, line=3,col = "steelblue2")
dev.off()


jpeg(paste0("OneYearRainStream",gridcode,".jpeg"),width = 6, height = 4, units = 'in',res = 500)
par(mar = c(5, 5, 2, 4.5))
plot(simu_x[20,],type="l",xlab = "Day of the year",ylab="Average rainfall (mm)",lwd=1.5)
par(new = TRUE)
plot(Y[20,],type="l",col="steelblue2",ylab = "",xlab = "",axes = FALSE, bty = "n",lwd=1.6)
axis(side=4, at = pretty(range(Y[20,])),col="steelblue2",col.axis="steelblue2")
mtext("Average streamflow (mm)", side=4, line=3,col = "steelblue2")
dev.off()





gridcode<-1759
cur_data <- read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/WaterInputOutputData_with_HBV_sim_streamflow/WaterInputOutputData_with_HBV_sim_streamflow/water_input_output_data_with_HBV_sim_streamflow_",gridcode,".csv"))
cur_data$Observed_Streamflow_mm.d[is.na(cur_data$Observed_Streamflow_mm.d)]<-cur_data$Simulated_Streamflow_mm.d[is.na(cur_data$Observed_Streamflow_mm.d)]
datevec<-seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day")


simu_x<-cur_data$Forcing_Precipitation_mm.d
simu_x<-simu_x[!(month(datevec)==2 & day(datevec)==29)]
simu_x<-matrix(simu_x,ncol = 365,byrow = T)
Y<-cur_data$Observed_Streamflow_mm.d
Y<-Y[!(month(datevec)==2 & day(datevec)==29)]
Y<-matrix(Y,ncol = 365,byrow = T)
Y[Y<0]<-0


jpeg(paste0("AvgRainStream",gridcode,".jpeg"),width = 6, height = 4, units = 'in',res = 500)
par(mar = c(5, 5, 2, 4.5))
plot(colMeans(simu_x),type="l",xlab = "Day of the year",ylab="Average rainfall (mm)",lwd=1.5)
par(new = TRUE)
plot(colMeans(Y),type="l",col="steelblue2",ylab = "",xlab = "",axes = FALSE, bty = "n",lwd=1.6)
axis(side=4, at = pretty(range(colMeans(Y))),col="steelblue2",col.axis="steelblue2")
mtext("Average streamflow (mm)", side=4, line=3,col = "steelblue2")
dev.off()


jpeg(paste0("OneYearRainStream",gridcode,".jpeg"),width = 6, height = 4, units = 'in',res = 500)
par(mar = c(5, 5, 2, 4.5))
plot(simu_x[20,],type="l",xlab = "Day of the year",ylab="Average rainfall (mm)",lwd=1.5)
par(new = TRUE)
plot(Y[20,],type="l",col="steelblue2",ylab = "",xlab = "",axes = FALSE, bty = "n",lwd=1.6)
axis(side=4, at = pretty(range(Y[20,])),col="steelblue2",col.axis="steelblue2")
mtext("Average streamflow (mm)", side=4, line=3,col = "steelblue2")
dev.off()



