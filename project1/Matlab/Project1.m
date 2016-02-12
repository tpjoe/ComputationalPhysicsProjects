clear
clc

% Define test cases
testCase = [10, 100, 1000, 10000];
testCase = testCase + 1;
error = zeros(size(testCase));
MatlabTime = testCase;
GaussianTime = testCase;
LU_time = testCase;

% Cslculate exact solution
figure
hold on
title('Solution of Poisson Equation using different number of points');
x = 0:0.001:1;
exact_u = 1-(1-exp(-10)).*x-exp(-10*x);
plot(x,exact_u);
legendText = 'Exact Solution';

% Numerical Solution
for j = 1:length(testCase)
    % Define varibles
    N = testCase(j)
    w = ones(1,N-2);
    cal_u = w*9999;
    h = 1/(N-1);
    x = h:h:1-h;
    f = h.^2*100*exp(-10*x);
    f0 = f;
    A = triDiaMat(N-2,2,-1,-1);
    exact_u = 1-(1-exp(-10)).*x-exp(-10*x);
    %% Matlab slash solver
    tic
    Matlab_cal_u = (f/A);
    MatlabTime(j) = toc;
    %% Gaussian simplified method
    tic
    % Forward Subsitute
    for i = 2:size(A,1)
        ratio = -(A(i,i-1)/A(i-1,i-1));
        A(i,:) = A(i,:)+A(i-1,:)*ratio;
        f(i) = f(i)+f(i-1)*ratio;
    end

    % Backward Subsitue
    %w(:) = 1;
    cal_u(N-2) = f(N-2)/w(N-2)/A(N-2,N-2);
    for i = (size(A,1)-1):-1:1
        cal_u(i) = (f(i)-cal_u(i+1)*w(i+1)*A(i,i+1))/(A(i,i)*w(i));
    end
    GaussianTime(j) = toc;

    relativeError = log10(abs((cal_u-exact_u)./exact_u));
    error(j) = max(relativeError);
    %plot(x,relativeError)
    plot(x,cal_u,'.')
    legendText = char(legendText,[num2str(N-1) ' points']);

    %% LU

    % Generate the matrix
    fullMatrix = 2*diag(ones(N-2,1)) -diag(ones(N-3,1),1) -diag(ones(N-3,1),-1);
    % Solve the problem
    tic
    [L,U] = lu(fullMatrix);
    LU_results = (f0/L)/U;
    LU_time(j) = toc;

end
legend(legendText);

figure
plot(log10(testCase),error,'o-');
title('Maximum relative error');
xlabel('Number of points(log)');

figure
plot(log10(testCase),(MatlabTime),log10(testCase),(GaussianTime),log10(testCase),(LU_time));
title('Time used for different method');
xlabel('Number of points(log)');
ylabel('Running time/s');
legend('Matlab slash slover','Gaussian o(n) method','Standard LU o(n^3) method','Location','NorthWest')