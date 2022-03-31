%% SourceIntegrationMethod.m
%
% This code implements the source integration method for interpolating the
% angular dependence of the radiance exiting a half space.
%
% Written by A. D. Kim on 9/20/2021 and modified on 9/24/2021

function I = SourceIntegrationMethod( MU, PHI, WT, MUgrid, PHIgrid, ...
    mu0, V, lambda, c, Psi, albedo, g, xk )

    % compute useful index arrays

    NN = length( MUgrid );
    
    [ M, N ] = size( V );

    [ indx, qndx ] = ndgrid( (1:NN), (1:M) );
    [ mndx, jndx ] = ndgrid( (1:N), (1:NN) );
    [ qqdx, jjdx ] = ndgrid( (1:M), (1:NN) );

    % ----------------------- %
    % COMPUTE THE INTERPOLANT %
    % ----------------------- %

    % compute the phase function

    CosTheta = MUgrid(indx) .* MU(qndx) ...
        + sqrt( ( 1 - MUgrid(indx).^2 ) .* ( 1 - MU(qndx).^2 ) ) ...
        .* cos( PHIgrid(indx) - PHI(qndx) );

    P = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * CosTheta ).^3;

    % contribution from homogeneous solution

    gamma0 = 1 + 1i * xk * sqrt( 1 - mu0^2 );
    gamma  = 1 + 1i * xk * sqrt( 1 - MUgrid(jndx).^2 ) .* cos( PHIgrid(jndx) );

    DOM_H = V * ( c(mndx) ./ ( gamma - lambda(mndx) .* MUgrid(jndx) ) );

    % contribution from the particular solution

    gamma = 1 + 1i * xk * sqrt( 1 - MUgrid(jjdx).^2 ) .* cos( PHIgrid(jjdx) );

    DOM_P = Psi(qqdx) * mu0 ./ ( gamma * mu0 - gamma0 * MUgrid(jjdx) );

    % contribution from the nonhomogeneous source

    CosTheta0 = MUgrid * mu0 ...
        + sqrt( 1 - MUgrid.^2 ) * sqrt( 1 - mu0^2 ) .* cos( PHIgrid );

    p0 = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * CosTheta0 ).^3;

    gamma  = 1 + 1i * xk * sqrt( 1 - MUgrid.^2 ) .* cos( PHIgrid );

    SRC = mu0 * p0 ./ ( gamma * mu0 - gamma0 * MUgrid );

    % compute the interpolant

    I = albedo * ( diag( P * diag(WT) * ( DOM_H + DOM_P ) ) + SRC );

return

%     % compute meshgrid of angles
% 
%     [ MUgrid, PHIgrid ] = meshgrid( mu, phi );
% 
%     % stretch the meshgrids to vectors
% 
%     MUgrid  = MUgrid(:);
%     PHIgrid = PHIgrid(:);
% 
%     % compute useful index arrays
% 
%     NN = length( MUgrid );
% 
%     [ indx, qndx ] = ndgrid( (1:NN), (1:length(V)) );
%     [ mndx, jndx ] = ndgrid( (1:length(c)), (1:NN) );
% 
%     % compute the interpolant
% 
%     CosTheta = MUgrid(indx) .* MU(qndx) ...
%         + sqrt( ( 1 - MUgrid(indx).^2 ) .* ( 1 - MU(qndx).^2 ) ) ...
%         .* cos( PHIgrid(indx) - PHI(qndx) );
% 
%     P = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * CosTheta ).^3;
% 
%     gamma0 = 1 + 1i * xk * sqrt( 1 - mu0^2 );
%     gamma  = 1 + 1i * xk * sqrt( 1 - MUgrid(jndx).^2 ) .* cos( PHIgrid(jndx) );
% 
%     C = c(mndx) ./ ( gamma - lambda(mndx) .* MUgrid(jndx) );
% 
%     CosTheta0 = MUgrid * mu0 ...
%         + sqrt( 1 - MUgrid.^2 ) * sqrt( 1 - mu0^2 ) .* cos( PHIgrid );
% 
%     p0 = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * CosTheta0 ).^3;
% 
%     gammagrid  = 1 + 1i * xk * sqrt( 1 - MUgrid.^2 ) .* cos( PHIgrid );
% 
%     I = albedo * diag( P * diag(WT) * ( V(ineg,:) * C ...
%         + Psi( qndx ) * mu0 ./ ( gamma * mu0 - gamma0 * MUgrid(indx) ) ) ...
%         + albedo * mu0 * p0 ./ ( gammagrid * mu0 - gamma0 * MUgrid );
% 
% 
% %     gamma0 = 1 + 1i * xk * sqrt( 1 - mu0^2 );
% %     gamma  = 1 + 1i * xk * sqrt( 1 - mu^2 ) * cos( phi );
% %     
% %     IDOM = V * diag( 1 ./ ( gamma - lambda * mu ) ) * c ...
% %         + Psi * mu0 ./ ( gamma * mu0 - gamma0 * mu );
% % 
% %     CosTheta = mu * MU ...
% %         + sqrt( 1 - mu^2 ) * sqrt( 1 - MU.^2 ) ...
% %         .* cos( PHI - phi );
% % 
% %     p = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * CosTheta ).^3;
% % 
% %     CosTheta0 = mu * mu0 ...
% %         + sqrt( 1 - mu^2 ) * sqrt( 1 - mu0^2 ) * cos( phi );
% % 
% %     p0 = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * CosTheta0 ).^3;
% % 
% %     I = albedo * p.' * diag(WT) * IDOM ...
% %         + albedo * mu0 * p0 / ( gamma * mu0 - gamma0 * mu );