function [A,At]=tensorGauss(N,band,sigma)
%tensorBlur function provides a N*N*N*N blur tensor
%
%        A=tensorBlur(P,center,BC)
%        A=tensorBlur(P,center)
%
%  Input:
%        P  Array containing the point spread function.
%   center  [row, col] = indices of center of PSF, P.
%       BC  String indicating boundary condition.
%             ('zero', 'reflexive', or 'periodic')
%           Default is 'zero'.
%        N  The fourth-order bluriing tensor size
%
%  Output:
%   A  fourth-order blurring tenor 


% Check inputs and set default parameters.
%
% Initialization.
if (nargin < 2), band = 3; end
band = min(band,N);
if (nargin < 3), sigma = 0.7; end

% Construct the matrix as a Kronecker product.
z = [exp(-([0:band-1].^2)/(2*sigma^2)),zeros(1,N-band)];
A = round(toeplitz(z),4);
%A = (1/sqrt(2*pi*sigma^2))*A;
B = (1/(2*pi*sigma^2))*kron(sparse(A),sparse(A));
A=tensor(B);
%A=ndsparse(A,[N,N,N,N]);
A=spreshape(A,[N,N,N,N]);
B=tensor(B');
At=spreshape(B,[N,N,N,N]);
%A=reshape(A,[N,N,N,N]);
%A=full(ndsparse(A,[N,N,N,N]));
