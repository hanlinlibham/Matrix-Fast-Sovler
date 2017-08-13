clear all;
lambda = 1+pi^2/10;
[lambda1,its1,Error1] = shiftinvert(64,lambda,1e-8,100);
[lambda2,its2,Error2] = shiftinvert(128,lambda,1e-8,100);
[lambda3,its3,Error3] = shiftinvert(256,lambda,1e-8,100);
ERROR = abs(lambda - [lambda1 lambda2 lambda3])
convergence = ERROR(2)/ERROR(1)