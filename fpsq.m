function [u,K,M]=fpsq(epsilon,c,f)
%setup
n=sqrt(length(f)); E0=speye(n); h=1/(n+1);
% step 0:  generate eigenvector matrix V for T
    [i,j]=meshgrid(1:n,1:n);V=sqrt(2/(n+1))* sin(i.*j*pi/(n+1));
% generate eigenvalue matrix D for T;
    ME1 = c*h*h/36; ME2 = 4*c*h*h/36;
    MT1 = 4*c*h*h/36; MT2 = 16*c*h*h/36;    
    E1 = ME1 - epsilon/3; E2 = ME2 - epsilon/3;
    T1 = MT1 - epsilon/3; T2 = MT2 + 8*epsilon/3;
    E = tridiags([E1,E2,E1],n);T = tridiags([T1,T2,T1],n);
    ME = tridiags([ME1,ME2,ME1],n);
    MT = tridiags([MT1,MT2,MT1],n);
    Z = tridiags([0,0,1], n);
    D0 = V*T*V;
    D1 = V*E*V;
% generate M,K,Tcal 
    M = kron(E0,MT)+kron(Z,ME)+kron(Z,ME)';
    K = kron(E0,T)+kron(Z,E)+kron(Z,E)';
    Tcal=kron(E0,D0)+kron(Z,D1)+kron(Z,D1)';
% step 1:  pack f into a matrix G, compute G=VF, unpack G
    F=reshape(f,n,n);
    G=V*F;g=G(:);
% step 2: solve for z and pack it into Z 
    z=Tcal\g;Z=reshape(z,n,n);
% step 3: pack Z and solve for U, unpack u
    U=V*Z;u=U(:);