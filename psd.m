function [x,k,resvec] = psd(A,b,x,tol,kmax)
epsilon = 1e-2; c = 1;
k = 0;
n = length(b);
r = b-A*x;
res = norm(r);resvec=res;
[u,K,M] = fpsq(epsilon,c,b);
while k<kmax & res>tol
    z=K\r;
    alpha=(r'*z)/(z'*A*z);
    x=x+alpha*z;
    r=b-A*x;
    res=norm(r);resvec=[resvec res];
    k=k+1;
end