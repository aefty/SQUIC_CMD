clc;
clear all;
rng(2);

p=10;
n=1000;

density=.2;
r=abs(rand(p,1)+1);
iCov= sprandsym(p,density,r);

mu=rand(p,1);
Y=mvnrnd(mu,inv(iCov),n);

lambda_scalar =100;
lambda_matrix=-lambda_scalar* double(boolean(iCov));

%lambda_full_QUIC = full(lambda_matrix);
%lambda_full_QUIC(lambda_full_QUIC==0)=lambda_scalar;

max_iter=200;
drop_tol=1e-6;
term_tol=1e-6;

verbose_in=0;

[X,W,squic_info]=SQUIC(Y,lambda_scalar,lambda_matrix,max_iter,drop_tol,term_tol,verbose_in);
[X_QUIC,W_QUIC]=QUICi('default',cov(Y),lambda_scalar+full(lambda_matrix),term_tol, verbose_in, max_iter,eye(p),eye(p));

disp('X')
full(X)
disp('X_QUIC')
full(X_QUIC)
disp('iCov')
full(iCov)