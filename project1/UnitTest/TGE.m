function [ cal_u ] =  TGE(N, d, e1, e2, f)
% Solve system of linear equations using TGE method
% function [ cal_u ] =  TGE(N, d, e1, e2, f)
    N = N+2;
    w = ones(1,N-2);
    cal_u = w*9999;
    %h = 1/(N-1);
    %x = h:h:1-h;
    %f = h.^2*100*exp(-10*x);
    A = triDiaMat(N-2,d,e1,e2);
    %exact_u = 1-(1-exp(-10)).*x-exp(-10*x);
    for i = 2:size(A,1)
        ratio = -(A(i,i-1)/A(i-1,i-1));
        A(i,:) = A(i,:)+A(i-1,:)*ratio;
        f(i) = f(i)+f(i-1)*ratio;
    end

    % Backward Subsitue
    cal_u(N-2) = f(N-2)/A(N-2,N-2);
    for i = (size(A,1)-1):-1:1
        cal_u(i) = (f(i)-cal_u(i+1)*A(i,i+1))/(A(i,i));
    end
    
end