clear all
clc
% num_sp = 1000;
% X=randn(num_sp,1);
% Y=X+0.5*randn(num_sp,1);
% Z=Y+0.5*randn(num_sp,1);
% 
% [p_val stat]=indtest_new(X,Z,[],[]);
% p_val % X and Z should be dependent
% 
% [p_val stat]=indtest_new(X,Z,Y,[]);
% p_val
% 
% ind  = Z > 1;
% [p_val stat]=indtest_new(X(ind),Z(ind),Y(ind),[]);
% p_val

X= randn(100,1);
X_target = randn(100,1);

width = 0.4 * std(X);
[beta_cs EXITFLAG_cs] = betaKMM_improved(X, X_target, width, 0, 0);