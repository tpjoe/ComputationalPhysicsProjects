% test TGE



%% Test 1: unity matrix test
N = 3;
d = 1;
e1 = 0;
e2 = 0;
f = [ 1 2 3 ];
result = f;

u = TGE(N, d, e1, e2, f);
for i = 1:length(f)
    assert(u(i) == result(i))
end
 
%% Test 2: 2x2 Matrix test
N = 2;
d = 2;
e1 = 1;
e2 = 1;
f = [ 1 2 ];
result = [ 0 1 ];

u = TGE(N, d, e1, e2, f);
for i = 1:length(f)
    assert(u(i) == result(i))
end

%% Test 3: 10 points Poisson's Equation
N = 9;
d = 2;
e1 = -1;
e2 = -1;
h = 0.1;
x = h:h:1-h;
f = h.^2*100*exp(-10*x);
result = [0.4899 0.6119 0.5986 0.5356 0.4542 0.3660 0.2754 0.1839 0.0920];
err = 1E-4;

u = TGE(N, d, e1, e2, f);
for i = 1:length(f)
    assert((u(i)-result(i))<err)
end