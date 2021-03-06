
clear
clc
close all

%% parameter set
p = 185;
Maxtheta = 179;
    theta = 0:3.6:Maxtheta; % 2.4,3.6,7.2 (Coarse reconstruction failed at 7)
    nA = numel(theta); % Number of view-angles
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

%% noise 
rng(1);
aaa= 0;
bbb= 0; %choose accordingly 0.02, 0.05 0.1
ccc = aaa+randn(size(yf))*bbb;
yf = yf+ccc;

yc = yc+ccc;

%%
% % Construction of adjusted coarse A_collapsed
[A_collapsed] = system_shrink(N,n,Af);
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
    % noise
    yt1 = yt1+ccc;
    F = transpose(A_collapsed);
    G = (F*A_collapsed)\F; %(A^T*A)^-1*A^T
Xc_t1 = G*yt1;
Condition_Number = cond((F*A_collapsed));

%
figure(211)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(Xc_t1,n,n)), xlabel('Coarse'),set(gca,'FontSize',14)
colormap gray

figure(221)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(Xc_t1,100),title('Coarse'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
%% I = x;
p= N;
q =N;

% Run TVAL3
clear opts
opts.mu = [2^10];
opts.beta = [2^9];
opts.mu0 = 2^10;      % trigger continuation shceme
opts.beta0 = 2^0;    % trigger continuation scheme
opts.tol_inn = 1e-4;
opts.tol = 1e-5;
opts.maxit = 1000;
opts.TVnorm = 1;
opts.nonneg = true;
%opts.disp = true;
%opts.Ut = I;


t = cputime;
[U,out] = TVAL3(Af,yf,p,q,opts);
t = cputime - t;
RRMSE = (norm(xf-U(:))/norm(xf))*100;
TV_y = Af*U(:);

M_RRMSE = (norm(TV_y-yf)/norm(yf))*100;

figure(14)
subplot(1,2,2),imagesc(reshape(U,N,N)), xlabel('TV-regularization');
title(['RRMSE (%) = ', num2str(RRMSE)])
subplot(1,2,1),imagesc(reshape(xf,N,N)),xlabel('Origional');
colormap gray


figure(15)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(U(:),100),title('TV'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)

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
 
lu = sort(xf(:));
mu = sort(U(:));
su = sort(x(:));


figure(1126)
stairs(lu);
hold on
stairs(mu);
legend('Actual', 'TV')

figure(1127)
stairs(lu);
hold on
stairs(su);
legend('Actual', 'Enlarged Coarse')

figure(1128)
stairs(su);
hold on
stairs(mu);
legend('Enlarged Coarse', 'TV') 

figure(1130)
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

fis2 = readfis('SingleFIS9'); % DoD0.28 for DoD 0.28 only
input =  [x1(:)];
fuzzy2 = evalfis(input,fis2);

y_fuzzy2 = Af*fuzzy2;


figure(31)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Actual'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(fuzzy2,N,N)), xlabel('Fuzzy2'),set(gca,'FontSize',14)
colormap gray

figure(32)
subplot(1,2,1), histogram(xf,100), title('Actual'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(fuzzy2,100),title('Fuzzy2'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)

%
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
                                            % Measurement Error

Mean_y = Af*Mean;  
Weighted_average_y = Af*Weighted_average;
Product_y = Af*Product;
Min_y = Af*Min;
Median_y = Af*Median;
Max_y = Af*Max;

M1 = (norm(yf-TV_y)/norm(yf))*100;
M2 = (norm(yf-y_coarse)/norm(yf))*100;
M3 = (norm(yf-Mean_y)/norm(yf))*100;                                           
M4 = (norm(yf-Weighted_average_y)/norm(yf))*100; 
M5 = (norm(yf-Product_y)/norm(yf))*100;
M6 = (norm(yf-Min_y)/norm(yf))*100;
% M7 = (norm(yf-Median_y)/norm(yf))*100;
M8 = (norm(yf-Max_y)/norm(yf))*100;
M9 = (norm(yf-y_fuzzy1)/norm(yf))*100;
M10= (norm(yf-y_fuzzy2)/norm(yf))*100;

%

figure(21)
subplot(2,5,1), imshow((reshape(xf,[N,N]))), xlabel('Actual'),set(gca,'FontSize',14)
title(['DoD =' , num2str(DoD2)])
ylabel(['Angles =', num2str(nA)])
subplot(2,5,2), imshow((reshape(U(:),[N,N]))), xlabel('TV-regularization'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A1)])
ylabel(['Measurement =', num2str(M1)])
subplot(2,5,3), imshow((reshape(Xc_t1,[n,n]))), xlabel('Enlarged Coarse'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A2)])
ylabel(['Measurement = ', num2str(M2)])
subplot(2,5,4), imshow((reshape(Mean,[N,N]))), xlabel('Mean'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A3)])
ylabel(['Measurement = ', num2str(M3)])
subplot(2,5,5), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A4)])
ylabel(['Measurement = ', num2str(M4)])
subplot(2,5,6), imshow((reshape(Product,[N,N]))), xlabel('Product'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A5)])
ylabel(['Measurement = ', num2str(M5)])
subplot(2,5,7), imshow((reshape(Min,[N,N]))), xlabel('Min'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A6)])
ylabel(['Measurement = ', num2str(M6)])
subplot(2,5,8), imshow((reshape(Max,[N,N]))), xlabel('Max'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A8)])
ylabel(['Measurement = ', num2str(M8)])
subplot(2,5,9), imshow((reshape(fuzzy1,[N,N]))), xlabel('Fuzzy 1'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A9)])
ylabel(['Measurement = ', num2str(M9)])
subplot(2,5,10), imshow((reshape(fuzzy2,[N,N]))), xlabel('Fuzzy 2'),set(gca,'FontSize',14)
title(['Accuracy Error = ' , num2str(A10)])
ylabel(['Measurement = ', num2str(M10)])
colormap gray


figure(22)
subplot(2,5,1), imshow((reshape(xf,[N,N]))), xlabel('Actual'),set(gca,'FontSize',16)
%title(['DoD =' , num2str(DoD2)])
subplot(2,5,2), imshow((reshape(U(:),[N,N]))), xlabel('TV-regularization'),set(gca,'FontSize',16)
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

Method = [{'TV'}, {'Coarse'}, {'Mean'}, {'Weighted_average'}, {'Product'}, {'Min'}, {'Max'},{'Fuzzy'}, {'Fuzzy2'}]';

T = table(Method, SSIM, PSNR, Accuracy_Error);
T(1:9,:)

% Coarse image metric

yyy_accuracy = (norm(xc-Xc_t1)/norm(xc))*100
yyy_ssim = ssim(Xc_t1,xc)
yyy_psnr = psnr(Xc_t1,xc)


% %%
% 
% 
% suu = sort(fuzzy2(:));
% 
% 
% % figure(1226)
% % stairs(lu);
% % hold on
% % stairs(mu);
% % legend('Actual', 'TV')
% 
% figure(1227)
% stairs(suu);
% hold on
% stairs(lu);
% legend('Fuzzy2', 'Actual')
% 
% figure(1228)
% stairs(suu);
% hold on
% stairs(mu);
% legend('Fuzzy2', 'TV') 
%% Stop here

% 
% 
% 
% %%                                           FUSION
% x1 = U(:);
% % create a matrix for fusion 
% B = [x,x1];
%     
% % different fusion methods 
% Mean = mean(B,2);
% Weighted_average = 0.99*(B(:,1)) + 0.01*(B(:,2));
% Product = real(sqrt(dot(B(:,1),B(:,2),2)));
% Min = min(B,[],2);
% Median = median(B,2); 
% Max = max(B,[],2);
% 
% % fuzzy Fusion 
% fis = readfis('Triangle_FUSION');
% input =  [x(:),x1(:)];
% fuzzy = evalfis(input,fis);
% 
% 
% 
% 
%                                             % Accuracy Error 
% j1 = (norm(xf-Mean)/norm(xf))*100;
% j2 = (norm(xf-Weighted_average)/norm(xf))*100;
% j3 = (norm(xf-Product)/norm(xf))*100;
% j4 = (norm(xf-Min)/norm(xf))*100;
% j5 = (norm(xf-Median)/norm(xf))*100;
% j6 = (norm(xf-Max)/norm(xf))*100;
% j7 = (norm(xf-fuzzy)/norm(xf))*100;
% 
% 
%                                             % Measurement Error
% 
% Mean_y = Af*Mean;  
% Weighted_average_y = Af*Weighted_average;
% Product_y = Af*Product;
% Min_y = Af*Min;
% Median_y = Af*Median;
% Max_y = Af*Max;
% fuzzy_y = Af*fuzzy;
% 
% k1 = (norm(yf-Mean_y)/norm(yf))*100;                                           
% k2 = (norm(yf-Weighted_average_y)/norm(yf))*100; 
% k3 = (norm(yf-Product_y)/norm(yf))*100;
% k4 = (norm(yf-Min_y)/norm(yf))*100;
% k5 = (norm(yf-Median_y)/norm(yf))*100;
% k6 = (norm(yf-Max_y)/norm(yf))*100;
% k7 = (norm(yf-fuzzy_y)/norm(yf))*100;
% 
% 
% % figure(1)
% % subplot(1,2,1), imagesc((reshape(x,[n,n]))), xlabel('Original'), axis image
% % subplot(1,2,2), imagesc((reshape(x1,[n,n]))), xlabel('CS'), axis image
% % title(['Filtere Back-projection DoD =' , num2str(DoD) ' and RRMSE ', num2str(rrmse)])
% % colormap gray
% 
% 
% 
%                                                %  figures
%                                                
% figure(12)
% subplot(2,4,1), imagesc((reshape(xf,[N,N]))), xlabel('Original')
% title(['DoD =' , num2str(DoD2)])
% subplot(2,4,2), imagesc((reshape(U(:),[N,N]))), xlabel('TV')
% title(['Accuracy Error = ' , num2str(RRMSE)])
% ylabel(['Measurement Error =', num2str(M_RRMSE)])
% subplot(2,4,3), imagesc((reshape(Mean,[N,N]))), xlabel('Mean')
% title(['Accuracy Error = ' , num2str(j1)])
% ylabel(['Measurement Error = ', num2str(k1)])
% subplot(2,4,4), imagesc((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average')
% title(['Accuracy Error = ' , num2str(j2)])
% ylabel(['Measurement Error = ', num2str(k2)])
% subplot(2,4,5), imagesc((reshape(Product,[N,N]))), xlabel('Product')
% title(['Accuracy Error = ' , num2str(j3)])
% ylabel(['Measurement Error = ', num2str(k3)])
% subplot(2,4,6), imagesc((reshape(Min,[N,N]))), xlabel('Min')
% title(['Accuracy Error = ' , num2str(j4)])
% ylabel(['Measurement Error = ', num2str(k4)])
% 
% % Use only if merge more than two images 
% % subplot(2,4,6), imagesc((reshape(Median,[n,n]))), xlabel('Median') 
% % title([' RRMSE ', num2str(j5)]) 
% 
% subplot(2,4,7), imagesc((reshape(Max,[N,N]))), xlabel('Max')
% title(['Accuracy Error = ' , num2str(j6)])
% ylabel(['Measurement Error = ', num2str(k6)])
% subplot(2,4,8), imagesc((reshape(fuzzy,[N,N]))), xlabel('Fuzzy')
% title(['Accuracy Error = ' , num2str(j7)])
% ylabel(['Measurement Error = ', num2str(k7)])
% colormap gray
% 
% 
%                                 % Mean and Standard deviation 
% M0 = mean(xf); S0 = std(xf); % Original
% M01 = mean(x1); S01 = std(x1); % TV-regularization
% M1 = mean(Mean); S1 = std(Mean);
% M2 = mean(Weighted_average); S2 = std(Weighted_average);
% M3 = mean(Product); S3 = std(Product);
% M4 = mean(Min); S4 = std(Min);
% M5 = mean(Median); S5 = std(Median);
% M6 = mean(Max); S6 = std(Max);
% M7 = mean(fuzzy); S7 = std(fuzzy);                                   
%                                         
% 
%                                         % Histograms
% figure(13)
% subplot(2,4,1), histogram(xf), xlabel('Origional'), legend('Origional')
% xlim([-1 1])
% title([' Mean ', num2str(M0),' and SD ', num2str(S0)])
% subplot(2,4,2), histogram(xf), hold on, histogram(x1), xlabel('TV'), legend('Origional','TV')
% title([' Mean ', num2str(M01),' and SD ', num2str(S01)])
% subplot(2,4,3), histogram(xf), hold on, histogram(Mean), xlabel('Mean'), legend('Origional', 'Mean')
% title([' Mean ', num2str(M1),' and SD ', num2str(S1)])
% subplot(2,4,4), histogram(xf), hold on, histogram(Weighted_average), xlabel('Weighted Average'),legend('Origional', 'Weighted average')
% title([' Mean ', num2str(M2),' and SD ',num2str(S2)])
% subplot(2,4,5), histogram(xf), hold on, histogram(Product), xlabel('Product'),legend('Origional', 'Product')
% title([' Mean ', num2str(M3),' and SD ',num2str(S3)])
% subplot(2,4,6), histogram(xf), hold on, histogram(Min), xlabel('Min'),legend('Origional', 'Min')
% title([' Mean ', num2str(M4),' and SD ',num2str(S4)])
% 
% % Use only if merge more than two images
% % subplot(2,4,6), histogram(x), hold on, histogram(Median), xlabel('Median'),legend('Origional', 'Median')
% % title([' Mean ', num2str(M5),' and SD ',num2str(S5)])
% 
% subplot(2,4,7), histogram(xf), hold on, histogram(Max), xlabel('Max'),legend('Origional', 'Max')
% title([' Mean ', num2str(M6),' and SD ',num2str(S6)])
% subplot(2,4,8), histogram(xf), hold on, histogram(fuzzy), xlabel('fuzzy'),legend('Origional', 'fuzzy')
% title([' Mean ', num2str(M7),' and SD ',num2str(S7)])
% 
% 
% %% figure for thesis
% 
% figure(21)
% subplot(2,4,1), imshow((reshape(xf,[N,N]))), xlabel('Original'),set(gca,'FontSize',14)
% %title(['DoD =' , num2str(DoD2)])
% subplot(2,4,2), imshow((reshape(U(:),[N,N]))), xlabel('TV'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(RRMSE)])
% %ylabel(['Measurement Error =', num2str(M_RRMSE)])
% subplot(2,4,3), imshow((reshape(Mean,[N,N]))), xlabel('Mean'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(j1)])
% %ylabel(['Measurement Error = ', num2str(k1)])
% subplot(2,4,4), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(j2)])
% %ylabel(['Measurement Error = ', num2str(k2)])
% subplot(2,4,5), imshow((reshape(Product,[N,N]))), xlabel('Product'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(j3)])
% %ylabel(['Measurement Error = ', num2str(k3)])
% subplot(2,4,6), imshow((reshape(Min,[N,N]))), xlabel('Min'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(j4)])
% %ylabel(['Measurement Error = ', num2str(k4)])
% 
% % Use only if merge more than two images 
% % subplot(2,4,6), imagesc((reshape(Median,[n,n]))), xlabel('Median') 
% % title([' RRMSE ', num2str(j5)]) 
% 
% subplot(2,4,7), imshow((reshape(Max,[N,N]))), xlabel('Max'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(j6)])
% %ylabel(['Measurement Error = ', num2str(k6)])
% subplot(2,4,8), imshow((reshape(fuzzy,[N,N]))), xlabel('Fuzzy'),set(gca,'FontSize',14)
% %title(['Accuracy Error = ' , num2str(j7)])
% %ylabel(['Measurement Error = ', num2str(k7)])
% colormap gray
% %%
% figure(22)
% subplot(2,4,1), histogram(xf), title('Original'), legend('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M0),' and SD ', num2str(S0)])
% subplot(2,4,2), histogram(xf), hold on, histogram(x1), title('TV'), legend('Original','TV'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M01),' and SD ', num2str(S01)])
% subplot(2,4,3), histogram(xf), hold on, histogram(Mean), title('Mean'), legend('Original', 'Mean'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M1),' and SD ', num2str(S1)])
% subplot(2,4,4), histogram(xf), hold on, histogram(Weighted_average), title('Weighted Average'),legend('Original', 'Weighted average'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M2),' and SD ',num2str(S2)])
% subplot(2,4,5), histogram(xf), hold on, histogram(Product), title('Product'),legend('Original', 'Product'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M3),' and SD ',num2str(S3)])
% subplot(2,4,6), histogram(xf), hold on, histogram(Min), title('Min'),legend('Original', 'Min'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M4),' and SD ',num2str(S4)])
% 
% % Use only if merge more than two images
% % subplot(2,4,6), histogram(x), hold on, histogram(Median), xlabel('Median'),legend('Origional', 'Median')
% % title([' Mean ', num2str(M5),' and SD ',num2str(S5)])
% 
% subplot(2,4,7), histogram(xf), hold on, histogram(Max), title('Max'),legend('Original', 'Max'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
% %title([' Mean ', num2str(M6),' and SD ',num2str(S6)])
% subplot(2,4,8), histogram(xf), hold on, histogram(fuzzy), title('Fuzzy'),legend('Original', 'Fuzzy'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M7),' and SD ',num2str(S7)])