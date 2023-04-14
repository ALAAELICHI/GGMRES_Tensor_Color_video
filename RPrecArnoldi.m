function [H, W] = RPrecArnoldi(AA,BB, H, V, psf)
% 
% One step of the Arnoldi algorithm applied to the right-preconditioned
% matrix A*M
%
% Inputs:
%       - A:        the system matrix, that can also be a psfMatrix
%       - M:        the preconditioner
%       - H:        the Hessenber matrix so far computed
%       - V:        the orthonormal basis for the Krylov subspace so far computed
%       - psf:      if the matrix A is a psfMatrix (and therefore we should work with 2D arrays): psf=1; otherwise, psf=0;
%
% Outputs:
%       - H:        updated Hessenberg matrix
%       - W:        new basis vector for the updated Krylov subspace
%
[p,q,ell]=size(AA);
k = size(H,2)+1;
h=zeros(k+1,1);
m=size(V,1);
if psf    
    proj=zeros(m);
    temp1=M.*cell2mat(V(k));
    W=A*temp1;
    for l=1:k
        temp2=cell2mat(V(l));
        h(l)=sum(sum(temp2.*W));
        proj=proj+h(l)*temp2;
    end
    W=W-proj;
    h(k+1)=norm(W,'fro');
else
    proj=zeros(m,1);
    %W=M*V(:,k);
    W=vec(operator(AA,BB,reshape(V(:,k),[p,q,ell])));%A*V(:,k);
    for l=1:k
        h(l)=V(:,l)'*W;
        proj=proj+h(l)*V(:,l);
    end
    W=W-proj;
    h(k+1)=norm(W);
end
if h(k+1)~=0
    W=W/h(k+1);
else
    display('Breakdown of Arnoldi')
end
if k==1
    H=[H h];
else
    H = [H; zeros(1,k-1)];
    H=[H h];
end
