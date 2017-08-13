clear all;load nrdiff;
tol = 1e-6; kmax = 1e5;
x0_1 = zeros(length(K1.RHS),1);
x0_2 = zeros(length(K2.RHS),1);
x0_3 = zeros(length(K3.RHS),1);
[x1,k1,resvec1] = psd(K1.A,K1.RHS,x0_1,tol,kmax);k1
[x2,k2,resvec2] = psd(K2.A,K2.RHS,x0_2,tol,kmax);k2
[x3,k3,resvec3] = psd(K3.A,K3.RHS,x0_3,tol,kmax);k3