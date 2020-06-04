clear
clc
close all 

%% parameter set
p = 185;
Maxtheta = 179;
    theta = 0:2.4:Maxtheta;
    nA = numel(theta); % Number of angles
    nM = p*nA; % Number of measurments
    gamma = 0.01;
    % Output variables and related parameters
    n = 64; %  
    k = 2;
    N = n*2; % Fine
%     d_f = round(sqrt(2)*N);
    d_c = round(sqrt(2)*n); % first to last beam distance in coarse 
    d_f = d_c*2;            % first to last beam distance in fine 
    
    DoD1 = p*nA/(n*n);
    DoD2 = p*nA/(N*N);

%% System Matrices configuration

% Construction of Fine matrix 
[Af,yf,xf,thetaf,pf,df] = paralleltomo(N,theta,p,d_f);  
  
% Construction of coarse matrix 
[Ac,yc,xc,theta,pc,d] = paralleltomo(n,theta,p,d_c);

%% noise 
rng(1);
aaa= 0;
bbb= 0; %choose accordingly 0.02, 0.05 0.1
ccc = aaa+randn(size(yf))*bbb;
yf = yf+ccc;

yc = yc+ccc;
%% load stored A_collapsed matrix 
load('CT_128by128.mat')
%%
% % Construction of adjusted coarse A_collapsed
% [A_collapsed] = system_shrink(N,n,Af);
% 
% % Construction of adjusted fine A_expanded
% [A_expanded] = system_expand(N,n,Ac);

%% Image adjustments

% Enlarge coarse to obtain adjusted fine 
x_enlarged = kron(reshape(xc,n,n),ones(2)); % Fine image from actual coarse phantom
x_enlarged = x_enlarged(:);
% or
%     Af_trans = transpose(Af);
%     ATAf = (Af_trans*Af);
% xf_enlarged2 = ATAf\Af_trans*Ac_collaps*xc; % Coarse image from actual fine phantom
%% reconstruction

% Track 1 *****************************************************************

% Reconstruct fine
Xf_t1 = lsqr(Af,A_collapsed*xc); 

% Reconstruct coarse 
    yt1 = Af*x_enlarged;
    yt1 = yt1+ccc;
    F = transpose(A_collapsed);
    G = (F*A_collapsed)\F; %(A^T*A)^-1*A^T
Xc_t1 = G*yt1;
%%
figure(01)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(Xc_t1,n,n)), xlabel('Coarse'),set(gca,'FontSize',14)
colormap gray

figure(02)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(Xc_t1,100),title('Coarse'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
%% 
% sinogram = reshape(yf,[p,nA]);
% fbp = iradon(sinogram,theta);
% U = fbp(1:N,1:N);
% U = U(:);

%% 



%%
RRMSE = (norm(xf-U(:))/norm(xf))*100;

psnr_U = psnr(U(:),xf);
ssim_U = ssim(U(:),xf);


figure(11)
subplot(1,2,2),imagesc(reshape(U,N,N)), xlabel('TV-regularization');
title(['RRMSE (%) = ', num2str(RRMSE)])%, xlabel('10643 elements in A with 29 projections')
subplot(1,2,1),imagesc(reshape(xf,N,N)),xlabel('Origional');
colormap gray

 
figure(12)
subplot(1,2,1), histogram(xf,100), xlabel('Original')
subplot(1,2,2), histogram(U(:),100), xlabel('Reconstruction')
%%                                           Enlargement

% Enlarge coarse to obtain adjusted fine 
x = kron(reshape(Xc_t1,n,n),ones(2)); % Fine image from actual coarse phantom
x = x(:);

figure(222111)
subplot(1,4,1),imagesc(reshape(Xc_t1,n,n)), xlabel('Coarse')
subplot(1,4,2),imagesc(reshape(xc,n,n)), xlabel('Original Coarse')
subplot(1,4,3), imagesc(reshape(x,N,N)), xlabel('Extended Coarse')
subplot(1,4,4), imagesc(reshape(xf,N,N)), xlabel('Original fine')
colormap gray
%%                                           FUSION
x1 = U(:);
% create a matrix for fusion 
B = [x,x1];

W1 = 0.75*(B(:,1)) + 0.25*(B(:,2));
W2 = 0.5*(B(:,1)) + 0.5*(B(:,2));
W3 = 0.25*(B(:,1)) + 0.75*(B(:,2));
W4 = 0.05*(B(:,1)) + 0.95*(B(:,2));
W5 = 0.01*(B(:,1)) + 0.99*(B(:,2));


A1 = (norm(xf-W1)/norm(xf))*100;
A2 = (norm(xf-W2)/norm(xf))*100;
A3 = (norm(xf-W3)/norm(xf))*100;
A4 = (norm(xf-W4)/norm(xf))*100;
A5 = (norm(xf-W5)/norm(xf))*100;
A6 = (norm(xf-U(:))/norm(xf))*100;
Accuracy_Error = [A1,A2,A3,A4,A5,A6]';

B1 = ssim(W1,xf);
B2 = ssim(W2,xf);
B3 = ssim(W3,xf);
B4 = ssim(W4,xf);
B5 = ssim(W5,xf);
B6 = ssim(U(:),xf);

SSIM = [B1,B2,B3,B4,B5,B6]';

C1 = psnr(W1,xf);
C2 = psnr(W2,xf);
C3 = psnr(W3,xf);
C4 = psnr(W4,xf);
C5 = psnr(W5,xf);
C6 = psnr(U(:),xf);

PSNR = [C1,C2,C3,C4,C5,C6]';

Weight = [{'25'}, {'50'}, {'75'}, {'95'},{'99'}, {'100'}]';

T = table(Weight, Accuracy_Error, SSIM, PSNR);
T(1:6,:)
%%
figure(235)
subplot(1,5,1), imshow(reshape(W1,N,N)), xlabel('Fine weight 25%'), title(['Accuracy Error = ' , num2str(A1)])
subplot(1,5,2), imshow(reshape(W2,N,N)), xlabel('Fine weight 50%'), title(['Accuracy Error = ' , num2str(A2)])
subplot(1,5,3), imshow(reshape(W3,N,N)), xlabel('Fine weight 75%'), title(['Accuracy Error = ' , num2str(A3)])
subplot(1,5,4), imshow(reshape(W4,N,N)), xlabel('Fine weight 95%'), title(['Accuracy Error = ' , num2str(A4)])
subplot(1,5,5), imshow(reshape(W5,N,N)), xlabel('Fine weight 99%'), title(['Accuracy Error = ' , num2str(A5)])
%subplot(1,6,6), imshow(reshape(U(:),N,N)), xlabel('Fine weight 100%'), title(['Accuracy Error = ' , num2str(A6)])