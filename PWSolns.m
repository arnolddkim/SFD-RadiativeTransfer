function [ eval, evec ] = PWSolns( xk, yk, L, albedo, mu, phi, M )
%% function [ eval, evec ] = PWSolns( xk, yk, L, a, s, mu, wt, phi, dphi, M )
% 
% This function computes the plane wave solutions that are used to compute
% the homogeneous solution for the radiative transfer equation at each
% spatial frequency pair, ( xk, yk ).
%
% Inputs:
%     xk   := spatial frequency in x
%     yk   := spatial frequency in y
%     L    := scattering operator (2M^2 by 2M^2 matrix)
%     a    := absorption coefficient
%     s    := scattering coefficient
%     mu   := Gauss-Legendre quadrature points
%     phi  := periodic trapezoid rule points
%     M    := order of the Gauss-Legendre quadrature rule
%
% Outputs:
%     eval := 1 by M^2 vector of eigenvalues with positive real part
%     evec := 2M^2 by M^2 matrix of corresponding eigenvectors
%
% Written by Arnold D. Kim, March 6, 2017

%% compute the composite angle arrays

MU  = repelem( mu, 2 * M );
PHI = repmat( phi, M, 1 );

%% set up the eigenvalue problem

A = - eye( 2 * M^2 ) + albedo * L ...
    - 1i * xk * diag( sqrt( 1.0 - MU.^2 ) .* cos( PHI ) ) ...
    - 1i * yk * diag( sqrt( 1.0 - MU.^2 ) .* sin( PHI ) );

A = diag( 1.0 ./ MU ) * A;

%% solve the eigenvalue problem

[ V, D ] = eig( A );

%% sort the eigenvalues by their real part

eval = diag(D);

[ ~, indx ] = sort( real(diag(D)) );
eval        = eval(indx(M^2+1:2*M^2));
evec        = V(:,indx(M^2+1:2*M^2));
  
return;
