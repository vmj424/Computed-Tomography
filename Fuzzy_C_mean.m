% a = [1,200,204,5,198,4,4,203,2,195,3,200]';
% [cc, sp, ofn] = fcm(a,2);
% b=(sp>0.5);
% 
% 
% [a,sp',b']
% 
% cc


%%
% fuzzy c-mean clustering
clc
clear all
a = phantom(128);
a = a(:);
% a = a';
d = sort(a);
histogram(a,100)

% fuzzy C-mean clustering
options = 3; % default value is 2 (Exponent of the fuzzy partition)

[cc, sp, ofn] = fcm(a,6,options);