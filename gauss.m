function [nodes,weights] = gauss(n,varargin)

% Compute Gauss-Legendre quadrature on compact interval [a,b]
%
% Usage: [NODES,WEIGHTS] = gauss(N [,A,B])
%
% Here, N is the length of the Gaussian quadrature rule, and the
% optional scalars A and B determine the compact interval. By default,
% there holds A = -1 and b = +1.
%
% The function returns an (N x 1) column vector WEIGHTS containing the
% quadrature weights and an (1 x N) row vector NODES containing the
% corresponding quadrature nodes.
%
% Example: Suppose we aim to approximate the integral 
%
%    int_a^b f dx
%
% Assume that F is a Matlab function, which takes a column vector of
% evaluation points X and returns a column vector Y of function values,
% where Y(j) = F(X(j)). Then, the numerical integration reads
%
%    [NODES,WEIGHTS] = gauss(n,A,B);
%    INT = WEIGHTS*F(NODES);
%
% Author: Dirk Praetorius - 11.03.2010


beta = (1:n-1)./sqrt((2*(1:n-1)).^2-1);
A = diag(beta,-1)+diag(beta,1);
[eigenvector,nodes] = eig(A);
[nodes,idx] = sort(diag(nodes));
weights = 2*eigenvector(1,idx).^2;

if nargin >= 3
    a = varargin{1};
    b = varargin{2};
    weights = 0.5*abs(b-a)*weights;
    nodes = 0.5*( a+b + nodes*(b-a) );
end
