% N = 50;
% a = zeros(N);
% a(min()N/4:N/1)=1;
% imshow(a)
% ph =MyRectanglePhantom(50)

% function ph = MyRectanglePhantom(N)
% % Initialize the result
% ph = zeros(N);
% 
% % Add rectangles
% ph(max(2,round(N/10)):min(end-1,round(9*N/10)),round(N/8):round(N/3)) = 1;
% ph(max(2,round(N/3)):min(end-1,round(6*N/10)),round(6*N/10):round(8*N/10)) = 1;


function ph = MyRectanglePhantom(N)

% Initialize the result
ph = zeros(N);
%ph((5:15),(5:6))=1;
%ph((round(6*N/10)):min(end-1,round(5*N/10))= 1;
% Add rectangles
% ph(max(2,round(N/10)):min(end-1,round(9*N/10)),round(N/8):round(N/3)) = 1;
%ph(max(2,round(N/3)):min(end-1,round(6*N/10)),round(6*N/10):round(8*N/10)) = 1

ph(max(2,round(2*N/10)):min(end-1,round(5*N/10)),round(N/6):round(N/3)) = 0.65;
ph(max(2,round(6*N/10)):min(end-1,round(8*N/10)),round(N/5):round(N/2)) = 1;
ph(max(2,round(N/3)):min(end-1,round(6*N/10)),round(6*N/10):round(8*N/10)) = 0.25;
