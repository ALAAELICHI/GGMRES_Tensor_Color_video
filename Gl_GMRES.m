function [X,alpha]=Gl_GMRES(A,BB,B,X,restrt,max_it)
R=B;
m=restrt;
k=0;
alpha=zeros(1,max_it);
while(k<max_it)
    k=k+1;
    [H,K]=Global_Arnoldi(A,BB,R,m);
    betta=tnorm(R);
    [U,S,V] = svd(H);
    g=diag(S);
    f=U'*eye(m+1,1);
    alpha(k) = gcv_tik(f,g,betta);
    M=H'*H+alpha(k)^2*eye(m);           
    b=betta*H'*eye(m+1,1);
    y=M\b;
    for jj=1:m
        X=X+K(:,:,:,jj)*y(jj);
    end
    R=B-operator(A,BB,X);
end