clc;
clear all;
close all
rng(2);

time_squic=[];
time_squic_core=[];
time_quic=[];

err_quic=[];
err_squic=[];

n=100;
p_set=100:100:1000;

for p=p_set;
    
    p
    
    density=5/p;
    r=abs(rand(p,1)+1);
    iCov= sprandsym(p,density,r);
    
    mu=rand(p,1);
    Y=mvnrnd(mu,inv(iCov),n);
    
    % L = lambda_matrix + lambda_scalar*(1-pattern(lambda_matrix))
    lambda_scalar = 5;
    lambda_matrix = lambda_scalar/1e12 * double(boolean(iCov));
    
    lambda_QUIC                 = lambda_matrix;
    lambda_QUIC(lambda_QUIC==0) = lambda_scalar;
    
    %lambda_full_QUIC = full(lambda_matrix);
    %lambda_full_QUIC(lambda_full_QUIC==0)=lambda_scalar;
    
    max_iter=100;
    drop_tol=1e-4;
    term_tol=1e-4;
    
    verbose_in=0;
    
    disp('SQUIC')
    tic
    [X,W,squic_info]=SQUIC(Y,lambda_scalar,lambda_matrix,max_iter,drop_tol,term_tol,verbose_in);
    time_squic(end+1)=toc
    
    time_squic_core(end+1)=squic_info.time_total
    
    disp('QUIC')
    tic
    [X_QUIC,W_QUIC]=QUICi('default',cov(Y),lambda_QUIC,term_tol, verbose_in, max_iter,eye(p),eye(p));
    time_quic(end+1)=toc
    
    err_squic(end+1)=full(sum(sum(abs(X-iCov))));
    err_quic(end+1)=sum(sum(abs(X_QUIC-iCov)));
    
end

figure
subplot(1,2,1)
plot(p_set,err_squic);hold on;
plot(p_set,err_quic);hold on;
legend('squic','quic')

subplot(1,2,2)
plot(p_set,time_squic);hold on;
plot(p_set,time_squic_core);hold on;
plot(p_set,time_quic);hold on;
legend('squic','squic core','quic')