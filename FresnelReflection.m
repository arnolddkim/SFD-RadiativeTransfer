function Rbc = FresnelReflection( nrel, mu, M )
%% function Rbc = FresnelReflection( nrel, mu, M )
%
% This function computes the Fresnel reflection coefficient for use in
% prescribing boundary conditions for the RTE.
%
% Inputs:
%     nrel := relative refractive index
%     mu   := M by 1 vector of Gauss-Legendre quadrature points
%     M    := number of Gauss-Legendre quadrature points
%
% Output:
%     Rbc  := M^2 by 1 vector of the Fresnel reflection coefficients
%
% Written by Arnold D. Kim, March 3, 2017

% compute the critical angle (if it exists)

if nrel >= 1.0

    mu_c = sqrt( 1.0 - 1.0 / nrel^2 );

else
    
    mu_c = 1.0;
    
end;

% compute the polar angle dependent Fresnel reflection coefficients

Rbc = zeros( M^2, 1 );

for im = 1 : M / 2
    
    mu2 = mu(im+M/2);
    
    if mu2 > mu_c
        
        mu1 = sqrt( 1.0 - nrel^2 * ( 1.0 - mu2^2 ) );
    
        rtmp = 0.5 * abs( ( nrel * mu1 - mu2 ) / ( nrel * mu1 + mu2 ) ).^2 ...
             + 0.5 * abs( ( nrel * mu2 - mu1 ) / ( nrel * mu2 + mu1 ) ).^2;
        
    else
        
        rtmp = 1.0;
        
    end;
    
    i = ( im - 1 ) * 2 * M + ( 1 : 2 * M );
    Rbc(i) = rtmp;
    
end;

return;