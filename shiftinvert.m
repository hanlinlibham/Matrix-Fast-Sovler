function [lambda,its,Error] = shiftinvert(n,s,tol,itmax)
epsilon = 1e-2; c=1;
N=(n-2)^2;
f = zeros(N,1);
[u,K,M] = fpsq(epsilon,c,f);
its = 0; error = 1;
w = ones(N,1);w = w/norm(w); Error = [];
while error>tol & its<itmax
    its = its+1;
    w = (K-s*M)\(M*w);
    w = w/norm(w);
    lambda = (w'*K*w)/(w'*M*w);
    error = norm(K*w-lambda*M*w);
    Error = [Error error];
end
    