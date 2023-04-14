clear all
X1 = im2double(imread('peppers.png'));
X2 = im2double(imread('papav256.pbm'));
 subplot(121); imshow(double(X1),[]); 
 subplot(122);imshow(double(X2),[]);
 figure
 A=blur(512,6,4);
AA(:,:,1)=0.70*A; AA(:,:,2)=0.15*A; AA(:,:,3)=0.15*A;
%-------------------------------------------------------------------------
B_exact= tprod(AA,X1); 
E = randn(size(B_exact));
E = E/tnorm(E) ; 
nu=0.001*tnorm(B_exact)*E;
B1 = B_exact+nu;
 A=blur(256,6,4);
 AAA(:,:,1)=0.70*A; AAA(:,:,2)=0.15*A; AAA(:,:,3)=0.15*A;
B_exact= tprod(AAA,X2); 
E = randn(size(B_exact));
E = E/tnorm(E) ; 
nu=0.001*tnorm(B_exact)*E;
B2 = B_exact+nu;
subplot(121); imshow(double(B1),[]); 
 subplot(122);imshow(double(B2),[]);