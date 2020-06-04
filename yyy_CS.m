
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
    n = 64; % Coarse
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
    F = transpose(A_collapsed);
    G = (F*A_collapsed)\F; %(A^T*A)^-1*A^T
Xc_t1 = G*yt1;
%%
figure(211)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(Xc_t1,n,n)), xlabel('Coarse'),set(gca,'FontSize',14)
colormap gray

figure(221)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(Xc_t1,100),title('Coarse'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
%% 
psi = real(fft2(eye(N^2)));
load('s1_no_noise.mat');
U = psi*s1;

%%


RRMSE = (norm(xf-U(:))/norm(xf))*100;
TV_y = Af*U(:);

M_RRMSE = (norm(TV_y-yf)/norm(yf))*100;

figure(14)
subplot(1,2,2),imagesc(reshape(U,N,N)), xlabel('TV-regularization');
title(['RRMSE (%) = ', num2str(RRMSE)])%, xlabel('10643 elements in A with 29 projections')
subplot(1,2,1),imagesc(reshape(xf,N,N)),xlabel('Origional');
colormap gray

 

%%                                           Enlargement

% Enlarge coarse to obtain adjusted fine 
x = kron(reshape(Xc_t1,n,n),ones(2)); % Fine image from actual coarse phantom
x = x(:);
y_coarse = Af*x;

figure(222111)
subplot(1,4,1),imagesc(reshape(Xc_t1,n,n)), xlabel('Coarse')
subplot(1,4,2),imagesc(reshape(xc,n,n)), xlabel('Original Coarse')
subplot(1,4,3), imagesc(reshape(x,N,N)), xlabel('Extended Coarse')
subplot(1,4,4), imagesc(reshape(xf,N,N)), xlabel('Original fine')
colormap gray
 %%                                 DATA VISUALATION
 
mu = sort(U(:));
lu = sort(x(:));
figure(1126)
ylim([-0.2 1.2])
stairs(lu);
hold on
stairs(mu);
legend('Enlarged Coarse', 'TV')

figure(1127)
subplot(1,2,1), histogram(U(:),100), xlabel('TV-regularization')
subplot(1,2,2), histogram(x,100), xlabel('Enlarged coarse')

%%                                           FUSION
x1 = U(:);
% create a matrix for fusion 
B = [x,x1];
    
% different fusion methods 
Mean = mean(B,2);
Weighted_average = 0.01*(B(:,1)) + 0.99*(B(:,2));
Product = real(sqrt(dot(B(:,1),B(:,2),2)));
Min = min(B,[],2);
Median = median(B,2); 
Max = max(B,[],2);

% fuzzy Fusion 
fis = readfis('Triangle_FUSION2');
input =  [x(:),x1(:)];
fuzzy1 = evalfis(input,fis);

y_fuzzy1 = Af*fuzzy1;

fis2 = readfis('SingleFIS2');
input =  [x1(:)];
fuzzy2 = evalfis(input,fis2);

y_fuzzy2 = Af*fuzzy2;


figure(31)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(fuzzy2,N,N)), xlabel('Fuzzy2'),set(gca,'FontSize',14)
colormap gray

figure(32)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(fuzzy2,100),title('Fuzzy2'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)

%%
                                            % Accuracy Error 


                                            
A1 = (norm(xf-U(:))/norm(xf))*100;
A2 = (norm(xf-x)/norm(xf))*100;
A3 = (norm(xf-Mean)/norm(xf))*100;
A4 = (norm(xf-Weighted_average)/norm(xf))*100;
A5 = (norm(xf-Product)/norm(xf))*100;
A6 = (norm(xf-Min)/norm(xf))*100;
% A7 = (norm(xf-Median)/norm(xf))*100;
A8 = (norm(xf-Max)/norm(xf))*100;
A9 = (norm(xf-fuzzy1)/norm(xf))*100;
A10 = (norm(xf-fuzzy2)/norm(xf))*100;
Accuracy_Error = [A1,A2,A3,A4,A5,A6,A8,A9,A10]';


                                           


%

% figure(21)
% subplot(2,5,1), imshow((reshape(xf,[N,N]))), xlabel('Actual'),set(gca,'FontSize',14)
% %title(['DoD =' , num2str(DoD2)])
% subplot(2,5,2), imshow((reshape(U(:),[N,N]))), xlabel('Compressed Sensing'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A1)])
% ylabel(['Measurement =', num2str(M1)])
% subplot(2,5,3), imshow((reshape(Xc_t1,[n,n]))), xlabel('Enlarged Coarse'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A2)])
% ylabel(['Measurement = ', num2str(M2)])
% subplot(2,5,4), imshow((reshape(Mean,[N,N]))), xlabel('Mean'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A3)])
% ylabel(['Measurement = ', num2str(M3)])
% subplot(2,5,5), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A4)])
% ylabel(['Measurement = ', num2str(M4)])
% subplot(2,5,6), imshow((reshape(Product,[N,N]))), xlabel('Product'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A5)])
% ylabel(['Measurement = ', num2str(M5)])
% subplot(2,5,7), imshow((reshape(Min,[N,N]))), xlabel('Min'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A6)])
% ylabel(['Measurement = ', num2str(M6)])
% subplot(2,5,8), imshow((reshape(Max,[N,N]))), xlabel('Max'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A8)])
% ylabel(['Measurement = ', num2str(M8)])
% subplot(2,5,9), imshow((reshape(fuzzy1,[N,N]))), xlabel('Fuzzy 1'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A9)])
% ylabel(['Measurement = ', num2str(M9)])
% subplot(2,5,10), imshow((reshape(fuzzy2,[N,N]))), xlabel('Image Refinement'),set(gca,'FontSize',14)
% title(['Accuracy Error = ' , num2str(A10)])
% ylabel(['Measurement = ', num2str(M10)])
% colormap gray


figure(22)
subplot(2,5,1), imshow((reshape(xf,[N,N]))), xlabel('Actual'),set(gca,'FontSize',16)
%title(['DoD =' , num2str(DoD2)])
subplot(2,5,2), imshow((reshape(U(:),[N,N]))), xlabel('Compressed Sensing'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A1)])
% ylabel(['Measurement Error =', num2str(M1)])
subplot(2,5,3), imshow((reshape(Xc_t1,[n,n]))), xlabel('Enlarged Coarse'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A2)])
% ylabel(['Measurement Error = ', num2str(M2)])
subplot(2,5,4), imshow((reshape(Mean,[N,N]))), xlabel('Arithmetic Mean'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A3)])
% ylabel(['Measurement Error = ', num2str(M3)])
subplot(2,5,5), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A4)])
% ylabel(['Measurement Error = ', num2str(M4)])
subplot(2,5,6), imshow((reshape(Product,[N,N]))), xlabel('Geometric Mean'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A5)])
% ylabel(['Measurement Error = ', num2str(M5)])
subplot(2,5,7), imshow((reshape(Min,[N,N]))), xlabel('Minimum'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A6)])
% ylabel(['Measurement Error = ', num2str(M6)])
subplot(2,5,8), imshow((reshape(Max,[N,N]))), xlabel('Maximum'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A8)])
% ylabel(['Measurement Error = ', num2str(M8)])
subplot(2,5,9), imshow((reshape(fuzzy1,[N,N]))), xlabel('Fuzzy'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A9)])
% ylabel(['Measurement Error = ', num2str(M9)])
subplot(2,5,10), imshow((reshape(fuzzy2,[N,N]))), xlabel('Image Refinement'),set(gca,'FontSize',16)
% title(['Accuracy Error = ' , num2str(A10)])
% ylabel(['Measurement Error = ', num2str(M10)])
colormap gray


%

sss1 = ssim(U(:),xf);% TV
sss2 = ssim(x(:),xf);% coarse extended
sss3 = ssim(Mean,xf);
sss4 = ssim(Weighted_average,xf);
sss5 = ssim(Product,xf);
sss6 = ssim(Min,xf);
%sss7 = ssim(Median,xf);
sss8 = ssim(Max,xf);
sss9 = ssim(fuzzy1,xf);
sss10 = ssim(fuzzy2,xf);
SSIM = [sss1,sss2,sss3,sss4,sss5,sss6,sss8,sss9,sss10]';


rrr1 = psnr(U(:),xf);% TV
rrr2 = psnr(x(:),xf);% coarse extended
rrr3 = psnr(Mean,xf);
rrr4 = psnr(Weighted_average,xf);
rrr5 = psnr(Product,xf);
rrr6 = psnr(Min,xf);
%rrr7 = ssim(Median,xf);
rrr8 = psnr(Max,xf);
rrr9 = psnr(fuzzy1,xf);
rrr10 = psnr(fuzzy2,xf);
PSNR = [rrr1,rrr2,rrr3,rrr4,rrr5,rrr6,rrr8,rrr9,rrr10]';

Method = [{'CS'}, {'Coarse'}, {'Mean'}, {'Weighted_average'}, {'Product'}, {'Min'}, {'Max'},{'Fuzzy'}, {'Fuzzy2'}]';

T = table(Method, SSIM, PSNR, Accuracy_Error);
T(1:9,:)