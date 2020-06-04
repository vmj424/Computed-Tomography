
clear
clc
close all
 
%% parameter set

p = 300;
Maxtheta = 179;
    theta = 0:5:Maxtheta;
    nA = numel(theta); % Number of angles
    nM = p*nA; % Number of measurments
    gamma = 0.01;
    % Output variables and related parameters
    n = 128; % Coarse
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


% noise 
rng(1);
aaa= 0;
bbb= 0; %choose accordingly 0.02, 0.05 0.1
ccc = aaa+randn(size(yf))*bbb;
yf = yf+ccc;

yc = yc+ccc;
%% load stored A_collapsed matrix 
%load('CT_128by128.mat')
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
    yt1 = yt1+ccc;
    F = transpose(A_collapsed);
    G = (F*A_collapsed)\F; %(A^T*A)^-1*A^T
Xc_t1 = G*yt1;

%
figure(001)
histogram(Xc_t1,100)
%%
figure(01)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(Xc_t1,n,n)), xlabel('Coarse'),set(gca,'FontSize',14)
colormap gray

figure(02)
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

figure(11)
subplot(1,2,2),imagesc(reshape(U,N,N)), xlabel('TV-regularization'),set(gca,'FontSize',16);
%title(['RRMSE (%) = ', num2str(RRMSE)])%, xlabel('10643 elements in A with 29 projections')
subplot(1,2,1),imagesc(reshape(xf,N,N)),xlabel('Actual'),set(gca,'FontSize',16);
colormap gray

 

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
 %%                                 DATA VISUALATION
 
mu = sort(U(:));
lu = sort(x(:));
figure(1126)
ylim([-0.2 1.2])
stairs(lu);
hold on
stairs(mu);
legend('Enlarged Coarse', 'TV')

figure(21)
subplot(1,2,1), histogram(U(:),100), xlabel('TV-regularization')
subplot(1,2,2), histogram(x,100), xlabel('Enlarged coarse')

%%                                           FUSION
x1 = U(:);
% create a matrix for fusion 
B = [x,x1];
    
% different fusion methods 
Mean = mean(B,2);
Weighted_average = 0.5*(B(:,1)) + 0.5*(B(:,2));
Product = real(sqrt(dot(B(:,1),B(:,2),2)));
Min = min(B,[],2);
Median = median(B,2); 
Max = max(B,[],2);

% fuzzy Fusion 
fis = readfis('Triangle_FUSION2');
input =  [x(:),x1(:)];
fuzzy = evalfis(input,fis);

%% fuzzy test 
fis2 = readfis('SingleFIS2');
input =  [x1(:)];
fuzzy2 = evalfis(input,fis2);

figure(31)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(fuzzy2,N,N)), xlabel('Fuzzy2'),set(gca,'FontSize',14)
colormap gray

figure(32)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(fuzzy2,100),title('Fuzzy2'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)

%%
                                            % Accuracy Error 
j1 = (norm(xf-Mean)/norm(xf))*100;
j2 = (norm(xf-Weighted_average)/norm(xf))*100;
j3 = (norm(xf-Product)/norm(xf))*100;
j4 = (norm(xf-Min)/norm(xf))*100;
j5 = (norm(xf-Median)/norm(xf))*100;
j6 = (norm(xf-Max)/norm(xf))*100;
j7 = (norm(xf-fuzzy)/norm(xf))*100;
j8 = (norm(xf-fuzzy2)/norm(xf))*100;

                                            % Measurement Error

Mean_y = Af*Mean;  
Weighted_average_y = Af*Weighted_average;
Product_y = Af*Product;
Min_y = Af*Min;
Median_y = Af*Median;
Max_y = Af*Max;
fuzzy_y = Af*fuzzy;
fuzzy2_y = Af*fuzzy2;

k1 = (norm(yf-Mean_y)/norm(yf))*100;                                           
k2 = (norm(yf-Weighted_average_y)/norm(yf))*100; 
k3 = (norm(yf-Product_y)/norm(yf))*100;
k4 = (norm(yf-Min_y)/norm(yf))*100;
k5 = (norm(yf-Median_y)/norm(yf))*100;
k6 = (norm(yf-Max_y)/norm(yf))*100;
k7 = (norm(yf-fuzzy_y)/norm(yf))*100;
k8 = (norm(yf-fuzzy2_y)/norm(yf))*100;

% figure(1)
% subplot(1,2,1), imagesc((reshape(x,[n,n]))), xlabel('Original'), axis image
% subplot(1,2,2), imagesc((reshape(x1,[n,n]))), xlabel('CS'), axis image
% title(['Filtere Back-projection DoD =' , num2str(DoD) ' and RRMSE ', num2str(rrmse)])
% colormap gray



                                               %  figures
                                               
figure(12)
subplot(2,4,1), imagesc((reshape(xf,[N,N]))), xlabel('Original')
title(['DoD =' , num2str(DoD2)])
subplot(2,4,2), imagesc((reshape(U(:),[N,N]))), xlabel('TV')
title(['Accuracy Error = ' , num2str(RRMSE)])
ylabel(['Measurement Error =', num2str(M_RRMSE)])
subplot(2,4,3), imagesc((reshape(Mean,[N,N]))), xlabel('Mean')
title(['Accuracy Error = ' , num2str(j1)])
ylabel(['Measurement Error = ', num2str(k1)])
subplot(2,4,4), imagesc((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average')
title(['Accuracy Error = ' , num2str(j2)])
ylabel(['Measurement Error = ', num2str(k2)])
subplot(2,4,5), imagesc((reshape(Product,[N,N]))), xlabel('Product')
title(['Accuracy Error = ' , num2str(j3)])
ylabel(['Measurement Error = ', num2str(k3)])
subplot(2,4,6), imagesc((reshape(Min,[N,N]))), xlabel('Min')
title(['Accuracy Error = ' , num2str(j4)])
ylabel(['Measurement Error = ', num2str(k4)])

% Use only if merge more than two images 
% subplot(2,4,6), imagesc((reshape(Median,[n,n]))), xlabel('Median') 
% title([' RRMSE ', num2str(j5)]) 

subplot(2,4,7), imagesc((reshape(Max,[N,N]))), xlabel('Max')
title(['Accuracy Error = ' , num2str(j6)])
ylabel(['Measurement Error = ', num2str(k6)])
subplot(2,4,8), imagesc((reshape(fuzzy,[N,N]))), xlabel('Fuzzy')
title(['Accuracy Error = ' , num2str(j7)])
ylabel(['Measurement Error = ', num2str(k7)])
colormap gray


                                % Mean and Standard deviation 
M0 = mean(xf); S0 = std(xf); % Original
M01 = mean(x1); S01 = std(x1); % TV-regularization
M1 = mean(Mean); S1 = std(Mean);
M2 = mean(Weighted_average); S2 = std(Weighted_average);
M3 = mean(Product); S3 = std(Product);
M4 = mean(Min); S4 = std(Min);
M5 = mean(Median); S5 = std(Median);
M6 = mean(Max); S6 = std(Max);
M7 = mean(fuzzy); S7 = std(fuzzy);                                   
                                        

                                        % Histograms
figure(13)
subplot(2,4,1), histogram(xf), xlabel('Origional'), legend('Origional')
xlim([-1 1])
title([' Mean ', num2str(M0),' and SD ', num2str(S0)])
subplot(2,4,2), histogram(xf), hold on, histogram(x1), xlabel('TV'), legend('Origional','TV')
title([' Mean ', num2str(M01),' and SD ', num2str(S01)])
subplot(2,4,3), histogram(xf), hold on, histogram(Mean), xlabel('Mean'), legend('Origional', 'Mean')
title([' Mean ', num2str(M1),' and SD ', num2str(S1)])
subplot(2,4,4), histogram(xf), hold on, histogram(Weighted_average), xlabel('Weighted Average'),legend('Origional', 'Weighted average')
title([' Mean ', num2str(M2),' and SD ',num2str(S2)])
subplot(2,4,5), histogram(xf), hold on, histogram(Product), xlabel('Product'),legend('Origional', 'Product')
title([' Mean ', num2str(M3),' and SD ',num2str(S3)])
subplot(2,4,6), histogram(xf), hold on, histogram(Min), xlabel('Min'),legend('Origional', 'Min')
title([' Mean ', num2str(M4),' and SD ',num2str(S4)])

% Use only if merge more than two images
% subplot(2,4,6), histogram(x), hold on, histogram(Median), xlabel('Median'),legend('Origional', 'Median')
% title([' Mean ', num2str(M5),' and SD ',num2str(S5)])

subplot(2,4,7), histogram(xf), hold on, histogram(Max), xlabel('Max'),legend('Origional', 'Max')
title([' Mean ', num2str(M6),' and SD ',num2str(S6)])
subplot(2,4,8), histogram(xf), hold on, histogram(fuzzy), xlabel('fuzzy'),legend('Origional', 'fuzzy')
title([' Mean ', num2str(M7),' and SD ',num2str(S7)])


%% figure for thesis

figure(14)
subplot(2,4,1), imshow((reshape(xf,[N,N]))), xlabel('Original'),set(gca,'FontSize',14)
%title(['DoD =' , num2str(DoD2)])
subplot(2,4,2), imshow((reshape(U(:),[N,N]))), xlabel('TV'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(RRMSE)])
%ylabel(['Measurement Error =', num2str(M_RRMSE)])
subplot(2,4,3), imshow((reshape(Mean,[N,N]))), xlabel('Mean'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j1)])
%ylabel(['Measurement Error = ', num2str(k1)])
subplot(2,4,4), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j2)])
%ylabel(['Measurement Error = ', num2str(k2)])
subplot(2,4,5), imshow((reshape(Product,[N,N]))), xlabel('Product'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j3)])
%ylabel(['Measurement Error = ', num2str(k3)])
subplot(2,4,6), imshow((reshape(Min,[N,N]))), xlabel('Min'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j4)])
%ylabel(['Measurement Error = ', num2str(k4)])

% Use only if merge more than two images 
% subplot(2,4,6), imagesc((reshape(Median,[n,n]))), xlabel('Median') 
% title([' RRMSE ', num2str(j5)]) 

subplot(2,4,7), imshow((reshape(Max,[N,N]))), xlabel('Max'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j6)])
%ylabel(['Measurement Error = ', num2str(k6)])
subplot(2,4,8), imshow((reshape(fuzzy,[N,N]))), xlabel('Fuzzy'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j7)])
%ylabel(['Measurement Error = ', num2str(k7)])
colormap gray
%%
figure(22)
subplot(2,4,1), histogram(xf), title('Original'), legend('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M0),' and SD ', num2str(S0)])
subplot(2,4,2), histogram(xf), hold on, histogram(x1), title('TV'), legend('Original','TV'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M01),' and SD ', num2str(S01)])
subplot(2,4,3), histogram(xf), hold on, histogram(Mean), title('Mean'), legend('Original', 'Mean'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M1),' and SD ', num2str(S1)])
subplot(2,4,4), histogram(xf), hold on, histogram(Weighted_average), title('Weighted Average'),legend('Original', 'Weighted average'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M2),' and SD ',num2str(S2)])
subplot(2,4,5), histogram(xf), hold on, histogram(Product), title('Product'),legend('Original', 'Product'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M3),' and SD ',num2str(S3)])
subplot(2,4,6), histogram(xf), hold on, histogram(Min), title('Min'),legend('Original', 'Min'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M4),' and SD ',num2str(S4)])

% Use only if merge more than two images
% subplot(2,4,6), histogram(x), hold on, histogram(Median), xlabel('Median'),legend('Origional', 'Median')
% title([' Mean ', num2str(M5),' and SD ',num2str(S5)])

subplot(2,4,7), histogram(xf), hold on, histogram(Max), title('Max'),legend('Original', 'Max'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M6),' and SD ',num2str(S6)])
subplot(2,4,8), histogram(xf), hold on, histogram(fuzzy), title('Fuzzy'),legend('Original', 'Fuzzy'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',12)
%title([' Mean ', num2str(M7),' and SD ',num2str(S7)])

%%%% figure for thesis

figure(21)
subplot(2,4,1), imshow((reshape(xf,[N,N]))), xlabel('Original'),set(gca,'FontSize',14)
%title(['DoD =' , num2str(DoD2)])
subplot(2,4,2), imshow((reshape(U(:),[N,N]))), xlabel('TV'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(RRMSE)])
%ylabel(['Measurement Error =', num2str(M_RRMSE)])
subplot(2,4,3), imshow((reshape(Mean,[N,N]))), xlabel('Mean'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j1)])
%ylabel(['Measurement Error = ', num2str(k1)])
subplot(2,4,4), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j2)])
%ylabel(['Measurement Error = ', num2str(k2)])
subplot(2,4,5), imshow((reshape(Product,[N,N]))), xlabel('Product'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j3)])
%ylabel(['Measurement Error = ', num2str(k3)])
subplot(2,4,6), imshow((reshape(Min,[N,N]))), xlabel('Min'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j4)])
%ylabel(['Measurement Error = ', num2str(k4)])

% Use only if merge more than two images 
% subplot(2,4,6), imagesc((reshape(Median,[n,n]))), xlabel('Median') 
% title([' RRMSE ', num2str(j5)]) 

subplot(2,4,7), imshow((reshape(Max,[N,N]))), xlabel('Max'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j6)])
%ylabel(['Measurement Error = ', num2str(k6)])
subplot(2,4,8), imshow((reshape(fuzzy,[N,N]))), xlabel('Fuzzy'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j7)])
%ylabel(['Measurement Error = ', num2str(k7)])
colormap gray
%%
figure(210)
subplot(2,5,1), imshow((reshape(xf,[N,N]))), xlabel('Original'),set(gca,'FontSize',14)
%title(['DoD =' , num2str(DoD2)])
subplot(2,5,2), imshow((reshape(U(:),[N,N]))), xlabel('TV'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(RRMSE)])
%ylabel(['Measurement Error =', num2str(M_RRMSE)])
subplot(2,5,3), imshow((reshape(Xc_t1,[n,n]))), xlabel('Coarse'),set(gca,'FontSize',14)
subplot(2,5,4), imshow((reshape(Mean,[N,N]))), xlabel('Mean'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j1)])
%ylabel(['Measurement Error = ', num2str(k1)])
subplot(2,5,5), imshow((reshape(Weighted_average,[N,N]))), xlabel('Weighted Average'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j2)])
%ylabel(['Measurement Error = ', num2str(k2)])
subplot(2,5,6), imshow((reshape(Product,[N,N]))), xlabel('Product'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j3)])
%ylabel(['Measurement Error = ', num2str(k3)])
subplot(2,5,7), imshow((reshape(Min,[N,N]))), xlabel('Min'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j4)])
%ylabel(['Measurement Error = ', num2str(k4)])

% Use only if merge more than two images 
% subplot(2,4,6), imagesc((reshape(Median,[n,n]))), xlabel('Median') 
% title([' RRMSE ', num2str(j5)]) 

subplot(2,5,8), imshow((reshape(Max,[N,N]))), xlabel('Max'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j6)])
%ylabel(['Measurement Error = ', num2str(k6)])
subplot(2,5,9), imshow((reshape(fuzzy,[N,N]))), xlabel('Fuzzy 1'),set(gca,'FontSize',14)
%title(['Accuracy Error = ' , num2str(j7)])
%ylabel(['Measurement Error = ', num2str(k7)])
subplot(2,5,10), imshow((reshape(fuzzy2,[N,N]))), xlabel('Fuzzy 2'),set(gca,'FontSize',14)
colormap gray

%%
 
sss = ssim(U(:),xf);
sss1 = ssim(fuzzy,xf);
sss2 = ssim(fuzzy2,xf);
SSIM = [sss, sss1, sss2]';

rr = immse(U(:),xf);
rr1 = immse(fuzzy,xf);
rr2 = immse(fuzzy2,xf);
IMMSE = [rr, rr1, rr2]';

ssrr = psnr(U(:),xf);
ssrr1 = psnr(fuzzy,xf);
ssrr2 = psnr(fuzzy2,xf);
PSNR = [ssrr, ssrr1, ssrr2]';

Method = [{'TV'}, {'Fuzzy'}, {'Fuzzy2'}]';

T = table(Method, SSIM, PSNR, IMMSE);
T(1:3,:)

%% 

aaa = sort(xf);
bbb = sort(U(:));
ggg = sort(fuzzy);
hhh = sort(fuzzy2);

figure(111)
plot(aaa);
hold on 
plot(bbb);
hold on 
plot(ggg)
hold on 
plot(hhh);
legend('actual','tv','fuzzy1','fuzzy2')
