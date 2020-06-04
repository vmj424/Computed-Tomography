function [Ar] = system_expand (N,n,Ac)

%Coarse_Reference  Generates an expanded coarse coefficient matrix and solution vector
%[Ar] = system_expand (N,n,Ac);

%
% Input:
%   N         The domain of fine image that consists of 
%             N^2 cells.   
%   n         The domain of coarse image that consists of 
%             n^2 cells.
%   Ac        Coarse coefficient matrix with n^2 columns 
         
%
% Output:
%   Ar        Coefficient matrix with N^2 columns that is built from A
%             Coarse with n^2 columns.

% Reference: 
%
% Code written by: Aydin M. Torkabadi, Esam Hussein
% 
% Copyright 2019 Aydin Torkabadi, Esam Hussein

%Expansion of Coarse Matrix
nM = size(Ac,1);
k = N/n;
[j,i,J,I] = deal(0);
Ar = zeros(nM,N*N);
index = zeros(n,n);
     
    for j=1:n
        for i=1:n
            index(i,j) = 1+(j-1)*k*N + (i-1) * k;
            for J = 1:k
                for I =1:k
                    Ar(:,index(i,j)+(J-1)*N+(I-1))=(1/k^2).*Ac(:,((j-1)*n+i));
                    %                      %SntSv     
                end
            end
        end
    end

% % Create rhs
% Yr = Ar * Xf;