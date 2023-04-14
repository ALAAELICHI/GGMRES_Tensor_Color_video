function [xm,lambda]=restArnoldi(AA,BB,b,m,x0)
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
r0=b-vec(operator(AA,BB,reshape(x0,[p,q,ell])));%b-A*x0;
psf=0;
n_r=norm(r0);
V(:,1)=r0/n_r;
H=[];
for k=1:m
    [H, w] = RPrecArnoldi(AA,BB,H,V,psf);
    if psf
        V{k+1}=w;
    else
        V=[V, w];
    end

end
d=[n_r;zeros(m,1)];
[Uh,Sh,Vh] = svd(H,0);
s=diag(Sh);
%bhat=Uh'*d;
%lambda = gcv_tikhonov(s, bhat);
% size(Uh)
 f=Uh'*eye(m+1,1);
  lambda = gcv_tik(f,s,n_r);
%lambda=0.0001;
%[lambda,G,reg_param] = gcv(Uh,s,d);
Areg=[H; sqrt(lambda)*eye(m)]; % Areg=[H; sqrt(lambda)*M*v(:,1:k)]; % Areg=[H; sqrt(lambda)*v(:,1:k)]; % TEST
breg=[d; zeros(m,1)]; % breg=[nr;zeros(k,1);zeros(N,1)]; %
ym=Areg\breg;
xm=V(:,1:m)*ym; % xj=x0+L*Z*yj; % TEST % seems correst like this
xm=x0+xm;