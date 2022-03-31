function [x,w] = GaussLegendre( N )
%% function [x,w] = GaussLegendre( N )
%
% compute the Gauss Legendre quadrature points and weights using
% Trefethen's method from Spectral Methods in MATLAB
%

  beta     = 0.5 ./ sqrt( 1.0 - ( 2.0 * (1:N-1) ).^(-2) );
  T        = diag( beta, 1 ) + diag( beta, -1 );
  [ V, D ] = eig( T );
  
  x     = diag( D );
  [x,i] = sort( x );
  w     = [ 2 * V(1,i).^2 ]';

return;
