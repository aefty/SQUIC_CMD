clc;
clear all;
close all

addpath QUIC


DATA={};
DATA.p_set=[];

DATA.time_squic=[];
DATA.time_squic_core=[];
DATA.time_quic=[];

DATA.err_l2_squic=[];
DATA.err_l2_quic=[];

DATA.err_l0_squic=[];
DATA.err_l0_quic=[];

n=100;
DATA.p_set=[20,40,80,160,320,640,1280];

for i=1:length(DATA.p_set);

    p=DATA.p_set(i);
    
    for j=1:10

        rng(p+j);

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
        DATA.time_squic(j,i)=toc
        
        DATA.time_squic_core(j,i)=squic_info.time_total
        
        disp('QUIC')
        tic
        [X_QUIC,W_QUIC]=QUIC('default',cov(Y),lambda_QUIC,term_tol, verbose_in, max_iter,eye(p),eye(p));
        DATA.time_quic(j,i)=toc
        
        DATA.err_l2_squic(j,i)=full(sum(sum(abs(X-iCov))));
        DATA.err_l2_quic(j,i)=sum(sum(abs(X_QUIC-iCov)));

        DATA.err_l0_squic(j,i)=full(sum(sum(abs(boolean(X)-boolean(iCov)))));
        DATA.err_l0_quic(j,i)=sum(sum(abs(boolean(X_QUIC)-boolean(iCov))));
    end  
end

save("results","DATA")

