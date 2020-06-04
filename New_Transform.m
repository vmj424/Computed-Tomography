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
    N = 128; % Fine
    k = 2;
    n = N/k; % Coarse

    d_c = round(sqrt(2)*n); % first to last beam distance in coarse 
    d_f = d_c*2;
    %d_f = round(sqrt(2)*N); % first to last beam distance in fine 
    
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
load('CT_128by128.mat') 


%[A_collapsed] = system_shrink(N,n,Af);


% Image adjustments

% Enlarge coarse to obtain adjusted fine 
x_enlarged = kron(reshape(xc,n,n),ones(k)); % Fine image from actual coarse phantom
x_enlarged = x_enlarged(:);

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
figure(01)
subplot(1,2,1),imshow(reshape(xf,N,N)), xlabel('Original'),set(gca,'FontSize',14)
subplot(1,2,2),imshow(reshape(Xc_t1,n,n)), xlabel('Coarse'),set(gca,'FontSize',14)
colormap gray

figure(02)
subplot(1,2,1), histogram(xf,100), title('Original'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(Xc_t1,100),title('Coarse'), xlabel('Grayscale value'), ylabel('Number of pixels'),set(gca,'FontSize',14)
% I = x;
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

M_RRMSE = (norm(yf-TV_y)/norm(yf))*100;

figure(11)
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
%                                 DATA VISUALATION
 
mu = sort(U(:));
lu = sort(x(:));
figure(1126)
ylim([-0.2 1.2])
stairs(lu);
hold on
stairs(mu);
legend('Enlarged Coarse', 'TV')

figure(212)
subplot(1,2,1), histogram(U(:),100), xlabel('TV-regularization'),set(gca,'FontSize',14)
subplot(1,2,2), histogram(x,100), xlabel('Enlarged coarse'),set(gca,'FontSize',14)


%% 
%%
T = fft(xf(:));

T_x = fft(x(:));
 
T_x1 = (fft(U(:)));

real_T = real(T);
imag_T = imag(T);

real_T_x = real(T_x);
imag_T_x = imag(T_x);

real_T_x1 = real(T_x1);
imag_T_x1 = imag(T_x1);

t = sort(real_T); tt = sort(real_T_x); tt1 = sort(real_T_x1);
g = sort(imag_T); gt = sort(imag_T_x); gt1 = sort(imag_T_x1);

figure(543)
subplot(1,3,1), stairs(tt), xlabel('Coarse')
subplot(1,3,2), stairs(tt1), xlabel('TV solution')
subplot(1,3,3), stairs(t), xlabel('Actual')
sgtitle('Real Part')

figure(5432)
subplot(1,3,1), stairs(gt), xlabel('Coarse')
subplot(1,3,2), stairs(gt1), xlabel('TV solution')
subplot(1,3,3), stairs(g), xlabel('Actual')
sgtitle('Imaginary Part')

figure(544)
subplot(1,3,1), hist(real_T_x,100), xlabel('Coarse')
subplot(1,3,2), hist(real_T_x1,100), xlabel('TV solution')
subplot(1,3,3), hist(real_T,100), xlabel('Actual')
sgtitle('Real Part')

figure(5442)
subplot(1,3,1), hist(imag_T_x,100), xlabel('Coarse')
subplot(1,3,2), hist(imag_T_x1,100), xlabel('TV solution')
subplot(1,3,3), hist(imag_T,100), xlabel('Actual')
sgtitle('Imaginary Part')

figure(1212)
subplot(1,2,2), stairs(g), hold on , stairs(gt), hold on, stairs(gt1), xlabel('Imagenry'), legend('Original','coarse','TV')
subplot(1,2,1), stairs(t), hold on , stairs(tt), hold on, stairs(tt1), xlabel('Real'), legend('Original','coarse','TV')

%%
figure(54412)
subplot(1,3,1), hist(real_T_x,100), xlabel('Coarse')
xlim([-150 150])
subplot(1,3,2), hist(real_T_x1,100), xlabel('TV solution')
xlim([-150 150])
subplot(1,3,3), hist(real_T,100), xlabel('Actual')
xlim([-150 150])
sgtitle('Real Part')
%
figure(54422)
subplot(1,3,1), hist(imag_T_x,100), xlabel('Coarse')
xlim([-30 30])
subplot(1,3,2), hist(imag_T_x1,100), xlabel('TV solution')
xlim([-30 30])
subplot(1,3,3), hist(imag_T,100), xlabel('Actual')
xlim([-30 30])
sgtitle('Imaginary Part')



%%                                                    FUSION

%%
T = [T_x,T_x1];
T_Max = max(T,[],2);
T_Min = min(T,[],2);
T_Mean = mean(T,2);

T_Product = real(sqrt(dot(T(:,1),T(:,2),2)));

Weighted_average = 0.01*(T(:,1)) + 0.99*(T(:,2));

F_T_Max = ifft(T_Max);
F_T_Min = ifft(T_Min);
F_T_Mean = ifft(T_Mean);
F_T_Product = ifft(T_Product);
F_T_WA = ifft(Weighted_average);


EE = (F_T_Min+F_T_Max)/1.86;
EE_error = (norm(xf(:)-EE)/norm(EE))*100;
figure(11)
imshow(reshape(EE(:),N,N)), xlabel(['Accuracy error = ', num2str(EE_error)]),set(gca,'Fontsize',14)

A1 = (norm(xf-U(:))/norm(xf))*100; % TV
A2 = (norm(xf-x)/norm(xf))*100;  % Coarse
A3 = (norm(xf-F_T_Max)/norm(xf))*100; % Max
A4 = (norm(xf-F_T_Min)/norm(xf))*100; % Min
A5 = (norm(xf-F_T_Mean)/norm(xf))*100; % Mean
A6 = (norm(xf-F_T_Product)/norm(xf))*100; % PRoduct
A7 = (norm(xf-F_T_WA)/norm(xf))*100; % WA

%%                                                   % FUZZY FUSION

rfis = readfis('Real_FUSION');
Real_inputs = [real_T_x, real_T_x1];
fused_real = evalfis(Real_inputs,rfis);


ifis  = readfis('Imag_FUSION');
Imag_inputs = [imag_T_x,imag_T_x1];
fused_imag = evalfis(Imag_inputs,ifis);

fused = complex(fused_real,fused_imag);
F_T_FUSED = real(ifft(fused));

A8 = (norm(xf-F_T_FUSED)/norm(xf))*100;



%%
                                                          %Refinement 
% Real 
fis = readfis('Real_refine.fis');

input_real =  [real_T_x1(:)];
real_refine_part = evalfis(input_real,fis);

% Imaginary
fi = readfis('Imag_refine.fis');

input_imag = [imag_T_x1(:)];
imag_refine_part = evalfis(input_imag,fi);
% 
% figure(1754)
% subplot(1,2,1),hist(real_refine_part,100), xlabel('Real part after refinement')
% xlim([-150 150])
% subplot(1,2,2), hist(imag_refine_part,100), xlabel('Imaginary part after refinement')
% xlim([-30 30])
% sgtitle('Histograms of Fuzzy Refinement solution')


resu = complex(real_refine_part, imag_refine_part);



F_T_Refinement = real(ifft(resu));
A9 = (norm(xf(:)-F_T_Refinement(:))/norm(xf(:)))*100;


%%                                                       FIGURES


figure(546)
subplot(2,5,1), imshow(reshape(xf(:),N,N)), xlabel('Actual')
subplot(2,5,2), imshow(reshape(U(:),N,N)), xlabel('TV-regularization'),title(['Accuracy Error = ' , num2str(A1)])
subplot(2,5,3), imshow(reshape(x,N,N)), xlabel('Enlarged Coarse'),title(['Accuracy Error = ' , num2str(A2)])
subplot(2,5,4), imshow(reshape(F_T_Mean,N,N)), xlabel('Mean'),title(['Accuracy Error = ' , num2str(A5)])
subplot(2,5,5), imshow(reshape(F_T_WA,N,N)), xlabel('Weighted Average'),title(['Accuracy Error = ' , num2str(A7)])
subplot(2,5,6), imshow(reshape(F_T_Product,N,N)), xlabel('Geometric Mean'),title(['Accuracy Error = ' , num2str(A6)])
subplot(2,5,7), imshow(reshape(F_T_Min,N,N)), xlabel('Minimum'),title(['Accuracy Error = ' , num2str(A4)])
subplot(2,5,8), imshow(reshape(F_T_Max,N,N)), xlabel('Maximum'),title(['Accuracy Error = ' , num2str(A3)])
subplot(2,5,9), imshow(reshape(F_T_FUSED,N,N)), xlabel('Fuzzy'),title(['Accuracy Error = ' , num2str(A8)])
subplot(2,5,10), imshow(reshape(F_T_Refinement,N,N)),xlabel('Image Refinement'), title(['Accuracy Error = ' , num2str(A9)])
colormap gray




figure(548)
subplot(2,5,1), imshow(reshape(xf(:),N,N)), xlabel('Actual'),set(gca,'FontSize',14)
subplot(2,5,2), imshow(reshape(U(:),N,N)), xlabel('TV-regularization'),set(gca,'FontSize',14)
subplot(2,5,3), imshow(reshape(x,N,N)), xlabel('Enlarged Coarse'),set(gca,'FontSize',14)
subplot(2,5,4), imshow(reshape(F_T_Mean,N,N)), xlabel('Mean'),set(gca,'FontSize',14)
subplot(2,5,5), imshow(reshape(F_T_WA,N,N)), xlabel('Weighted Average'),set(gca,'FontSize',14)
subplot(2,5,6), imshow(reshape(F_T_Product,N,N)), xlabel('Geometric Mean'),set(gca,'FontSize',14)
subplot(2,5,7), imshow(reshape(F_T_Min,N,N)), xlabel('Minimum'),set(gca,'FontSize',14)
subplot(2,5,8), imshow(reshape(F_T_Max,N,N)), xlabel('Maximum'),set(gca,'FontSize',14)
subplot(2,5,9), imshow(reshape(F_T_FUSED,N,N)), xlabel('Fuzzy'),set(gca,'FontSize',14)
subplot(2,5,10), imshow(reshape(F_T_Refinement,N,N)),xlabel('Image Refinement'),set(gca,'FontSize',14)
colormap gray


%%
sss1 = ssim(U(:),xf);% TV
sss2 = ssim(x(:),xf);% coarse extended
sss8 = ssim(F_T_Max,xf);
sss6 = ssim(F_T_Min,xf);
sss3 = ssim(F_T_Mean,xf);
sss5 = ssim(F_T_Product,xf);
%sss7 = ssim(Median,xf);
sss4 = ssim(F_T_WA,xf);
sss9 = ssim(F_T_FUSED,xf);
sss10 = ssim(F_T_Refinement,xf);
SSIM = [sss1,sss2,sss3,sss4,sss5,sss6,sss8,sss9,sss10]';


rrr1 = psnr(U(:),xf);% TV
rrr2 = psnr(x(:),xf);% coarse extended
rrr3 = psnr(F_T_Mean,xf);
rrr4 = psnr(F_T_WA,xf);
rrr5 = psnr(F_T_Product,xf);
rrr6 = psnr(F_T_Min,xf);
%rrr7 = ssim(Median,xf);
rrr8 = psnr(F_T_Max,xf);
rrr9 = psnr(F_T_FUSED,xf);
rrr10 = psnr(F_T_Refinement,xf);
PSNR = [rrr1,rrr2,rrr3,rrr4,rrr5,rrr6,rrr8,rrr9,rrr10]';

Accuracy_Error = [A1,A2,A5,A7,A6,A4,A3,A8,A9]';
Method = [{'TV'}, {'Coarse'}, {'Mean'},{'WA'},{'Product'}, {'Min'} ,{'Max'},{'Fuzzy'},{'Refinement'}]';

T = table(Method, Accuracy_Error,SSIM,PSNR);
T(1:9,:)