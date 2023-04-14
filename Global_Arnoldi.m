function [H,K]=Global_Arnoldi(A,B,R0,m)
[M,N,k]=size(double(R0));
K=zeros(M,N,k,m+1);
H=zeros(m+1,m);
V1=R0/tnorm(R0);
K(:,:,:,1)=V1;
for i=1:m
    V=operator(A,B,V1);
    for j=1:i
        H(j,i)=V(:)'*V1(:);
        V=V-H(j,i)*K(:,:,:,j);
    end
    H(i+1,i)=tnorm(V);
    V1=V/H(i+1,i);
    K(:,:,:,i+1)=V1;
end
end