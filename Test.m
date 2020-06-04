clc
clear all
 
y = [1 2 3  4 5 6 7 8 9];
x = [3 2 1 6 5 4 9 7 8];



for ii = 1:1:length(x)
    if x(ii) >= y(ii)
        new(ii) = x(ii);
    else
        new(ii) = y(ii);
       ii = ii+1;
    end
    %ii = ii+1;
end
T_x = x;

T_x1 = y;


for kk = 1:1:length(xf)
    if T_x1(kk) >= T_x(kk)
    new(kk) = T_x1;
    else
    new(kk) = T_x(kk);
    kk = kk+1;
    end
end
