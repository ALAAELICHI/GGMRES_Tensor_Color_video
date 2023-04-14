function [xapprox,Errel,Discr,Lambda,volume,TotIt,stopReSt] = ReStart(AA,BB,b,sigma,eta,ne,ni,x)
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
Errel=zeros(ne*ni,1);
Discr=zeros(ne*ni,1);
Lambda=zeros(ne*ni,1);
stopReSt=zeros(ne,1);
volume=zeros(ne,1);
%b(b<0)=0; % first attempt to force nonnegativity
v=sum(sum(b));
TotIt=0;
lambda=1;
% if size(b,2)>1
%     psf=1;
%     n=size(b,2);
%     xapprox=cell(ne*ni,1);
%     x0=zeros(n);
%     M=ones(n);
% else
%     psf=0;
    N=length(b);
    xapprox=zeros(N,ne*ni);
    x0=zeros(N,1);
   % M=speye(N);
%end
for extit=1:ne
    display(extit)
    s='y'; % stop as son as the discrepancy is satisfied
    [xapprox_temp,Res_temp,Errel_temp,Lambda_temp,stopReSt_temp]=GAT_StdFormTrans(AA,BB,b,lambda,ni,s,sigma,eta,x0,x);
    EfIt=stopReSt_temp;
%     if psf
%         xapprox((TotIt+1):(TotIt+EfIt))=xapprox_temp(2:(EfIt+1));
%         x0=cell2mat(xapprox_temp(EfIt+1));
%         switch optn
%             case 1
%                 x0(x0<0)=0;
%                 M=sqrt(x0);
%             case 2
%                 x0=abs(x0);
%                 M=sqrt(x0);
%             case 3
%                 x0(x0<0)=0;
%                 M=sqrt(x0);
%                 x0=zeros(n);
%             case 4
%                 M=sqrt(abs(x0));
%         end
%     else
%         xapprox(:,(TotIt+1):(TotIt+EfIt))=xapprox_temp(:,2:(EfIt+1));
%         x0=xapprox_temp(:,EfIt+1);
%         switch optn
%             case 1
%                 x0(x0<0)=0;
%                 M=spdiags(sqrt(x0),0:0,N,N);
%             case 2
%                 x0=abs(x0);
%                 M=spdiags(sqrt(x0),0:0,N,N);
%             case 3
%                 x0(x0<0)=0;
%                 M=spdiags(sqrt(x0),0:0,N,N);
%                 x0=zeros(N,1);
%             case 4
%                 M=spdiags(sqrt(abs(x0)),0:0,N,N);
%         end
%     end    
    EfIt=stopReSt_temp;
    Lambda((TotIt+1):(TotIt+EfIt))=Lambda_temp(1:EfIt);
    Discr((TotIt+1):(TotIt+EfIt))=Res_temp(1:EfIt);
    Errel((TotIt+1):(TotIt+EfIt))=Errel_temp(1:EfIt);
    stopReSt(extit)=stopReSt_temp;
    volume(extit)=sum(sum(x0));
    TotIt=TotIt+EfIt;
    lambda=Lambda_temp(EfIt);    
end
% if psf
%     xapprox=xapprox(1:TotIt);
% else
    xapprox=xapprox(:,1:TotIt);
%end
Lambda=Lambda(1:TotIt);
Discr=Discr(1:TotIt);
Errel=Errel(1:TotIt);
