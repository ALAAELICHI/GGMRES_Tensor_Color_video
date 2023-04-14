%%
clear all;
clc;
EXACTSOL=double(imread('enamel.tif'));
EXACTSOL=EXACTSOL(1:256,1:256);
sigma=2.5;
t=6;
[n,n]=size(EXACTSOL);
for i=1:n
    for j=1:n
        if abs(i-j)<= t;
            H1(i,j)=(1/(2*sqrt(2*pi)*sigma))*exp(-(i-j)^2/(2*sigma^2));
        else 
           H1(i,j)=0;
        end
    end
end
H1=sparse(H1);
H2=H1;
B_exact=H1*EXACTSOL*H2'; 
E = randn (size(B_exact));
E = E / norm(E,'fro') ; 
B = B_exact+0.001*norm(B_exact,'fro')*E ;
%*******************************************************************

s=norm(0.001*norm(B_exact,'fro')*E,'fro');

eta=1.1;
mu0=10;

tic
[X,mu,err_r,l]=GGKB(B,EXACTSOL,n,H1,H2,mu0,eta,s);
time=toc
GGKB_steps=l
err_r
mu
mesh(X)
%*******************************************************************

