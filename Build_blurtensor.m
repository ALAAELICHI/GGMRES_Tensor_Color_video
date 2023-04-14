addpath(pwd)
N=256;
P=psfGauss(5,2);
A=sptensor(zeros(3,3,5,5));
for i=1:N
    for j=1:N
    A(i,j,i,j)=P(5,5);A(i,j,i,j+1)=P(5,4);A(i,j,i,j+2)=P(5,3);A(i,j,i,j+3)=P(5,2);A(i,j,i,j+4)=P(5,1);
    A(i,j,i+1,j)=P(4,5);A(i,j,i+1,j+1)=P(4,4);A(i,j,i+1,j+2)=P(4,3);A(i,j,i+1,j+3)=P(4,2);A(i,j,i+1,j+4)=P(4,1);
    A(i,j,i+2,j)=P(3,5);A(i,j,i+2,j+1)=P(3,4);A(i,j,i+2,j+2)=P(3,3);A(i,j,i+2,j+3)=P(3,2);A(i,j,i+2,j+4)=P(3,1);
    A(i,j,i+3,j)=P(2,5);A(i,j,i+3,j+1)=P(2,4);A(i,j,i+3,j+2)=P(2,3);A(i,j,i+3,j+3)=P(2,2);A(i,j,i+3,j+4)=P(2,1);
    A(i,j,i+4,j)=P(1,5);A(i,j,i+4,j+1)=P(1,4);A(i,j,i+4,j+2)=P(1,3);A(i,j,i+4,j+3)=P(1,2);A(i,j,i+4,j+4)=P(1,1);
    end
end
% 
% 
save blurtensor A
  