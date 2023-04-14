function [xapprox,Lambda] = RestGmres(AA,BB,b,ne,ni)
%
% Algorithm 2 along with some variants.
%
% Inputs:
%       - A:        the system matrix, that can also be a psfMatrix
%       - b:        noisy right hand side that can both be a 1D and 2D array
%       - sigma:    noise level
%       - eta:      safety parameter to define the discrepancy principle; usually eta=1.01
%       - optn:     specifies whcih variant has to be applied at each restarts
%                 - optn=1 : set to zero the negative entries of the last solution
%                 - optn=2 : take the absolute value of the last solution
%                 - optn=3 : take the new initial guess always equal to 0
%                 - optn=4 : don't modify the last solution
%                 - optn=5: project to the nonnegative subspace, with a
%                   given volume
%       - ne:       maximum number of outer iterations
%       - ni:       maximum number of inner iterations
%       - x:        exact solution
%
% Outputs:
%       - xapprox:  sequence of all the m+1 intermediate solution. The first one is the initial guess, 
%                   the one computed when the discrepancy principle is satisfied is at position stopFT+1.
%                   If a 1D problem is considered, xapprox is an array of size (length(xex))x(m+1);
%                   if a 2D problem is considered, xapprox is a cell array of size (m+1)x1.
%       - Errel:    sequence of the quantities norm(xapprox(j+1)-xex,2)/norm(xex)
%       - Discr:    sequence of the quantities norm(b-A*xapprox(j+1),2)/norm(b,2)
%       - Lambda:   sequence of the regularization parameters employed at each iteration
%       - volume:   volume of the images (sum of the pixels), at each
%                   iteration
%       - TotIt:    total number of iterations performed (obtained adding the number of iterations performed ad each restart)
%       - stopReSt: vector containing the number of iterations performed at each restart
%
% Errel=zeros(ne*ni,1);
% Discr=zeros(ne*ni,1);
 Lambda=zeros(ne,1);
% stopReSt=zeros(ne,1);
% volume=zeros(ne,1);
%b(b<0)=0; % first attempt to force nonnegativity
% TotIt=0;
% lambda=1;

N=length(b);
xapprox=zeros(N,ne);
x0=zeros(N,1);

for extit=1:ne
    display(extit)
    %[xapprox_temp,Res_temp,Errel_temp,Lambda_temp,stopReSt_temp]=GAT_StdFormTrans(AA,BB,b,lambda,ni,s,sigma,eta,x0);
    [xm,lambda]=restArnoldi(AA,BB,b,ni,x0);
    Lambda(extit)=lambda;
    x0=xm;
    xapprox(:,extit)=xm;
end
