function T=tridiags(v,n)

e = ones(n,1);
T = spdiags(kron(v,e),-1:1, n,n);