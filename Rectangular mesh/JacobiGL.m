function [x,w] = JacobiGL(alpha,beta,N)

% function [x] = JacobiGL(alpha,beta,N)
% Purpose: Compute the N'th order Gauss Lobatto quadrature 
%          points, x,and weights,w, associated with the Jacobi polynomial,
%          of type (alpha,beta) > -1 ( <> -0.5). 

x = zeros(N+1,1);
w = zeros(N+1,1);
if (N==1) x(1)=-1.0; x(2)=1.0; w(1)=1;w(2)=1; return; end

[xint,w1] = JacobiGQ(alpha+1,beta+1,N-2);
wint = (2*N+1)/(N*(N+1)) ./ (JacobiP(xint,0,0,N)).^2 ;

x = [-1, xint', 1]';
w = [2/(N)/(N+1),wint',2/(N)/(N+1)]';
return;
