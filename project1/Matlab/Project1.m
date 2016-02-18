clear
clc

% Define test cases
testCase = [11];
%testCase = 1000:1000:5000;
error = testCase;
MatlabTime = testCase;
GaussianTime = testCase;
LU_time = testCase;

% Exact Solution
figure
hold on
title('Solution of Poisson Equation using different number of points');
x = 0:0.001:1;
x(1) = [];
exact_u = 1-(1-exp(-10)).*x-exp(-10*x);
plot(x,exact_u);
legendText = 'Exact Solution';

% Numerical Solution
for j = 1:length(testCase)
% Define varibles
N = testCase(j)
w = ones(1,N);
cal_u = w*9999;
h = 1/(N-1);
x = 0:h:1
x(1) = [];
f = h.^2*100*exp(-10*x)
f0 = f;
A = triDiaMat(N-1,2,-1,-1);
exact_u = 1-(1-exp(-10)).*x-exp(-10*x);
%% Matlab slash solver
tic
Matlab_cal_u = (f/A);
MatlabTime(j) = toc;
plot(x,Matlab_cal_u,'o');
%% Gaussian
tic
% Forward Subsitute
for i = 2:size(A,1)
    ratio = -(A(i,i-1)/A(i-1,i-1));
    A(i,:) = A(i,:)+A(i-1,:)*ratio;
    %w(i) = w(i)+w(i-1)*ratio;
    f(i) = f(i)+f(i-1)*ratio;
end

%Matlab_cal_u = (f/A);
% Backward Subsitue
%w(:) = 1;
cal_u(N) = f(N)/w(N)/A(N,N);
for i = (size(A,1)-1):-1:1
    cal_u(i) = (f(i)-cal_u(i+1)*w(i+1)*A(i,i+1))/(A(i,i)*w(i));
end
GaussianTime(j) = toc;

relativeError = (cal_u(1:end-2)-exact_u(2:end-1))./exact_u(2:end-1);
error(j) = max(log10(abs(relativeError)));
plot(x,cal_u,'.')
legendText = char(legendText,[num2str(N) ' points']);

%% LU

% Generate the matrix
fullMatrix = 2*diag(ones(N,1)) -diag(ones(N-1,1),1) -diag(ones(N-1,1),-1);
% Solve the problem
tic
[L,U] = lu(fullMatrix);
LU_results = f0/L/U;
LU_time(j) = toc;

end
legend(legendText);

figure
plot(log10(testCase),error,'o-');
title('Maximum relative error');
xlabel('Number of points(log)');

figure
plot(log10(testCase),log10(MatlabTime),log10(testCase),log10(GaussianTime),log10(testCase),log10(LU_time));
title('Time used for different method');
legend('Matlab slash slover','Gaussian o(n) method','Standard LU o(n^3) method')