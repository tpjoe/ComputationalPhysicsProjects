clear
clc
% Define varibles
N = 40;
w = ones(1,N);
cal_u = w*9999;
h = 1/(N-1);
x = 0:h:1;
f = h.^2*100*exp(-10*x);
A = triDiaMat(N,2,-1,-1);
exact_u = 1-(1-exp(-10)).*x-exp(-10*x);
Matlab_cal_u = (f/A);
%%
% Forward Subsitute
for i = 2:size(A,1)
    ratio = -(A(i,i-1)/A(i-1,i-1));
    A(i,:) = A(i,:)+A(i-1,:)*ratio;
    %w(i) = w(i)+w(i-1)*ratio;
    f(i) = f(i)+f(i-1)*ratio;
end
%Matlab_cal_u = (f/A);
% Backward Subsitue
w(:) = 1;
cal_u(N) = f(N)/w(N)/A(N,N);
for i = (size(A,1)-1):-1:1
    cal_u(i) = (f(i)-cal_u(i+1)*w(i+1)*A(i,i+1))/(A(i,i)*w(i));
end

plot(x,cal_u,x,exact_u,x,Matlab_cal_u)
legend('cal','exact','Matlab');