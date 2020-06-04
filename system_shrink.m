function [Ar] = system_shrink (N,n,Af)

%Coarse_Reference  Generates an expanded coarse coefficient matrix and solution vector
%[Ar] = system_shrink (N,n,Af);
%
% Input:
%   N         The domain of fine image that consists of 
%             N^2 cells.   
%   n         The domain of coarse image that consists of 
%             n^2 cells.
%   Af        Fine coefficient matrix with N^2 columns         
%
% Output:
%   Ar        Coefficient matrix with n^2 columns that is built from A
%             Fine with N^2 columns.



% Reference: 
%
% Code written by: Aydin Torkabadi
% 
% Copyright 2019 

%Shrinking the fine Matrix

nM = size(Af,1);
k = N/n;
[j,i,J,I,counter] = deal(0);
Ar = zeros(nM,n*n);
indx = zeros(N,N);
nn = n^2;
     
for mm = 1:nM
    for j=1:n
        for i=1:n
        indx(i,j) = 1+(j-1)*k*N + (i-1) * k;
        counter = 0;
            for J = 1:k
                for I =1:k
                counter = Af(mm,indx(i,j)+(J-1)*N+(I-1));
                Ar(mm,((j-1)*n+i))= counter + Ar(mm,((j-1)*n+i));
            %                                  %jumper %SntSv     
                end
            end
        end
    end
end