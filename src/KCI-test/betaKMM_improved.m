% calculate betas solving min. problem (12) in paper
function [beta EXITFLAG] = betaKMM_improved(X, Xtst, sigma, doingRealTraining, regression);

nsamples = size(X,1);  % number of train samples
ntestsamples = size(Xtst,1);  % number of test samples
Centering_Kx = 0;

% variables: (here in the program / in (12) in the paper)
% H is K
% f is kappa
%

% minimize...
% 'calculating H=K...'
H = rbf_dot(X,X,sigma,0);
H=(H+H')/2; %make the matrix symmetric (it isn't symmetric before because of bad precision)

% 'calculating f=kappa...'
R3 = rbf_dot(X,Xtst,sigma,0);
f=(R3*ones(ntestsamples, 1));
f=-nsamples/ntestsamples*f;

% did the same, but slowlier:
% f=-nsamples/ntestsamples*ones(nsamples,1);
% for i=1:nsamples
%     fi=0;
%     for j=1:ntestsamples
%         fi = fi + rbf_dot(X(i,:),Xtst(j,:),sigma);
%     end
%     f(i,1) = f(i,1)*fi;
% end
%
% do they really the same?
%'different f?'
%[f1 f]

% subject to...
% abs(sum(beta_i) - m) <= m*eps
% which is equivalent to A*beta <= b where A=[1,...1;-1,...,-1] and b=[m*(eps+1);m*(eps-1)]
eps = (sqrt(nsamples)-1)/sqrt(nsamples);
%eps=1000/sqrt(nsamples);
% A=ones(1,nsamples);
% A(2,:)=-ones(1,nsamples);
% b=[nsamples*(eps+1); nsamples*(eps-1)];

Aeq = [];
beq = [];
% 0 <= beta_i <= 1000 for all i
LB = zeros(nsamples,1);
UB = ones(nsamples,1).*100; % 1000

% reparameterizatio of beta to make it smooth
if nsamples < 200
    width = 0.8;
elseif nsamples < 1200
    width = 0.5;
    %        width = 0.8;
else
    width = 0.3; % 0.3
end
lambda2 = 0.1;
Thresh = 1E-5;

Kx = rbf_dot(X,X,width,0);
if Centering_Kx
    Hc1 = eye(nsamples) - ones(nsamples,nsamples)/nsamples;
    %     KY_x = Hc1 * KY_x * Hc1; %%%??? or
    
    Kx = Hc1 * Kx * Hc1;
end

% will be used in the iterations
pdinv_Kx = pdinv(Kx + lambda2*eye(nsamples));

% [UU1,SS1,VV1] = svd(Kx *pdinv_Kx);
% eig1 = diag(SS1);
% II1 = find(eig1 > max(eig1) * Thresh);
% Q = Kx*pdinv_Kx * VV1(:,II1);
Q = Kx * pdinv_Kx;
alpha0 = pdinv(Q'*Q) * Q'* ones(nsamples,1);


Rc = Q;

J = Rc' * H * Rc; J = (J+J')/2;
Thresh_beta = 1E-3;
UB_beta = 10;
B = 10; lb_sum_beta = 1-B/sqrt(nsamples)/4; ub_sum_beta =  1+B/sqrt(nsamples)/4;


A = [-Rc; Rc; ones(1,nsamples)*Rc; -ones(1,nsamples)*Rc];
b = [zeros(nsamples, 1)-Thresh_beta; ones(nsamples,1)*UB_beta; ub_sum_beta*nsamples; -lb_sum_beta*nsamples];
options = optimset('MaxIter', 1000); % 5000

% X=QUADPROG(H,f,A,b,Aeq,beq,LB,UB) attempts to solve the quadratic programming problem:
%              min 0.5*x'*H*x + f'*x
% subject to:  A*x <= b
%              Aeq*x = beq
%              LB <= x <= UB

% 'solving quadprog for betas...'
[alpha,FVAL,EXITFLAG] = quadprog(J, Rc' * f, A,b,Aeq,beq,-1E4*ones(size(Rc,2),1), 1E4*ones(size(Rc,2),1),alpha0, options);
% [beta,FVAL,EXITFLAG] = quadprog(H,f,A,b,Aeq,beq,LB,UB);
% EXITFLAG
if ((EXITFLAG==0 ) && (doingRealTraining==1))
    %[beta,FVAL,EXITFLAG] = quadprog(H,f,A,b,Aeq,beq,LB,UB,beta,optimset('MaxIter',1e4));
    EXITFLAG
end

beta = Rc * alpha;

if (regression==0)
    % guarantee that all beta greater than 0
    threshold=0.01*abs(median(beta));
    beta(find(beta<threshold)) = threshold;
%     sprintf('number of beta <~ %f: %d (0 is good)', threshold, length(find(beta<threshold)))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXITFLAG:
%       1  QUADPROG converged with a solution X.
%       3  Change in objective function value smaller than the specified tolerance.
%       4  Local minimizer found.
%       0  Maximum number of iterations exceeded.
%      -2  No feasible point found.
%      -3  Problem is unbounded.
%      -4  Current search direction is not a direction of descent; no further
%           progress can be made.
%      -7  Magnitude of search direction became too small; no further progress can
%           be made.
