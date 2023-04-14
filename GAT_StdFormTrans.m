function [xapprox,Res,Errel,Lambda,stopAT]=GAT_StdFormTrans(AA,BB,b,lambda,m,s,sigma,eta,x0,x)
%
% Arnoldi-Tikhonov method applied to the problem transformed into standard form. 
% Parameter selection strategy: discrepancy principle; at each iteration:
% secant-update method.
%
% Inputs:
%       - A:        the system matrix, that can also be a psfMatrix
%       - b:        noisy right hand side that can both be a 1D and 2D array
%       - M:        the "inverse" of the regularization matrix, that acts as a preconditioner 
%       - lambda:   initial value for the regularization parameter
%       - m:        maximum number of iterations
%       - s:        decide whether to stop (s='y') or not (s='n') once the discrepancy principle is satisfied
%       - sigma:    noise level
%       - eta:      safety parameter to define the discrepancy principle; usually eta=1.01
%       - x0:       initial guess for the solution (default: x0=0)
%       - x:        exact solution
%
% Outputs:
%       - xapprox:  sequence of all m+1 the intermediate solution. The first one is the initial guess, 
%                   the one computed when the discrepancy principle is satisfied is at position stopFT+1.
%                   If a 1D problem is considered, xapprox is an array of size (length(xex))x(m+1);
%                   if a 2D problem is considered, xapprox is a cell array of size (m+1)x1.
%       - Res:      sequence of the quantities norm(b-A*xapprox(j+1),2)/norm(b,2)
%       - Errel:    sequence of the quantities norm(xapprox(j+1)-xex,2)/norm(xex)
%       - Lambda:   sequence of the regularization parameters employed at each iteration
%       - stopAT:   iteration at which the discrepancy principle is fulfilled;
%
[p,q,ell]=size(AA);
Res=zeros(m,1);
Errel=zeros(m,1);
Lambda=zeros(m,1);
stopAT=m;
r0=b-vec(operator(AA,BB,reshape(x0,[p,q,ell])));%b-A*x0;
% if size(b,2)>1
%     psf=1;    
%     n=size(b,2);
%     xapprox=cell(1,m+1);
%     n_r=norm(r0,'fro');
%     n_x=norm(x,'fro');
%     n_b=norm(b,'fro');
%     xapprox{1}=x0;
%     V{1}=r0/n_r;
% else
    psf=0;
    N=length(b);
    xapprox=zeros(N,m+1);
    n_r=norm(r0);
    n_x=norm(x);
    n_b=norm(b);
    xapprox(:,1)=x0;
    V(:,1)=r0/n_r;
%end
Lambda(1)=lambda;
H=[];
for k=1:m
    [H, w] = RPrecArnoldi(AA,BB,H,V,psf);
    if psf
        V{k+1}=w;
    else
        V=[V, w];
    end
    d=[n_r;zeros(k,1)];
    yGMRES=H\d;
    n_rGMRES=norm(H*yGMRES-d);
    Areg=[H; sqrt(lambda)*eye(k)]; % Areg=[H; sqrt(lambda)*M*v(:,1:k)]; % Areg=[H; sqrt(lambda)*v(:,1:k)]; % TEST
    breg=[d; zeros(k,1)]; % breg=[nr;zeros(k,1);zeros(N,1)]; % 
    yk=Areg\breg;
    if psf
        xk=zeros(n);
        for i=1:k
            xk=xk+yk(i)*cell2mat(V(i));
        end
        xk=x0+M.*xk;
        xapprox{k+1}=xk;
    else
        xk=V(:,1:k)*yk; % xj=x0+L*Z*yj; % TEST % seems correst like this
        xk=x0+xk;
        xapprox(:,k+1)=xk;
    end
    res=H*yk-d;  
    nresb=norm(res)/n_b;
    Res(k)=nresb;
    if psf
        Errel(k)=norm(xk-x,'fro')/n_x;
    else
        Errel(k)=norm(xk-x)/n_x; % Errel(j)=norm(xjS-xS)/norm(xS); % 
    end
    if nresb<sigma*eta && stopAT==m
        stopAT=k;
        if s=='y'
            xapprox=xapprox(:,1:k+1); Res=Res(1:k); 
            Errel=Errel(1:k); Lambda=Lambda(1:k);
            return
        end
    end
    mul=(eta*n_b*sigma-n_rGMRES)/(norm(res)-n_rGMRES); 
    mul=abs(mul);
    lambda=mul*lambda;
    Lambda(k+1)=lambda;
end
