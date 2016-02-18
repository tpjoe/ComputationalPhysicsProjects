function [ A ] = triDiaMat( n, d, e1, e2 )
%TRIDIAMAT generate tridiagonal matrix
%   [ A ] = triDiaMat( n, d, e1, e2 )

A = d*eye(n);
for i = 1:n-1
    A(i,i+1) = e1;
    A(i+1,i) = e2;
end

end

