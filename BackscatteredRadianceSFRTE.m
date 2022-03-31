%% BackscatteredRadianceSFRTE.m
%
% This Matlab code solves the 3D RTE in the spatial frequency domain for 
% a fixed spatial frequency.
%
% A spatially modulated plane wave incident obliquely on a half space 
% composed of a uniform absorbing and scattering medium. Scattering is
% governed by the Henyey-Greenstein scattering model. Therefore, the
% optical properties of the medium are the albedo, anisotropy factor, and
% relative refractive index. The user must also specify the angle of
% incidence (on the xz-plane) and the spatial frequency.
%
% This code implements the discrete ordinate method (DOM) using a product
% quadrature rule.
%
% To interpolate the radiance with respect to angles, this code implements
% a source integration method.
%
% NOTE: The source is prescribed just inside the medium, so it does not 
% take into account specular reflection, nor the transmission coefficient 
% from outside of the medium to inside of it.
%
% Last updated by A. D. Kim on 2/21/2022

clear;

%% set the figure parameters

set(0,'defaultaxesfontsize',18,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% set the optical properties

albedo = 0.990; % albedo
g      = 0.800; % anisotropy factor
nrel   = 1.000; % relative refractive index of the half space

%% set the spatial frequency and wavenumber

f  = 4.0;
xk = 2 * pi * f;

%% set the incident direction of the plane wave

theta0 = 20.0;               % in degrees
theta0 = theta0 / 180 * pi;  % convert to radians
mu0    = cos( theta0 );      % compute the cosine of the polar angle

%% set the order of the quadrature rule

M = 32;

%% compute the quadrature rules used for the discrete ordinate method

% Gauss-Legendre Quadrature Rule

[ mu, wt ] = GaussLegendre( M );

% Periodic Trapezoid Rule

dphi = pi / M;
phi  = ( -pi : dphi : pi - dphi )';

%% compute the Fresnel reflection coefficient

Rbc = FresnelReflection( nrel, mu, M );

%% compute the composite angle arrays

MU  = repelem( mu, 2 * M );
PHI = repmat( phi, M, 1 );
WT  = repelem( wt, 2 * M ) * dphi;

%% compute the scattering operator

% compute the spherical harmonics

[ Ynm, nvec ] = ComputeSphericalHarmonics( M, MU, PHI );

% compute the discrete Henyey-Greenstein scattering operator using the
% product quadrature rule

L = Ynm * diag( g.^nvec ) * Ynm' * diag( WT );

%% compute the source term

% compute the scattering angle with respect to mu' = mu0 and phi' = 0

XI0  = MU * mu0 + sqrt( 1 - MU.^2 ) * sqrt( 1 - mu0^2 ) .* cos( PHI );

% compute the angular distribution of the source

src0 = 0.25 / pi * ( 1 - g^2 ) ./ sqrt( 1 + g^2 - 2 * g * XI0 ).^3;

%% compute index arrays for the forward and backward hemispheres

ipos = ( 1 : 2 * M^2 );      % forward hemisphere
ineg = ( 2 * M^2 : -1 : 1 ); % backward hemisphere

% reorder subvectors of ineg so that indices increase with phi

for jm = 1 : M
    
    ineg((jm-1)*2*M+(1:2*M)) = ineg((jm-1)*2*M+(2*M:-1:1));
    
end

%% DISCRETE ORDINATE METHOD

% compute the particular solution: alpha * exp( -gamma * z )

gamma = ( 1 + 1i * xk * sqrt( 1 - mu0^2 ) ) / mu0;  

Ap    = - gamma * diag( MU ) ...
        + 1i * xk * diag( sqrt( 1 - MU.^2 ) .* cos(PHI) ) ...
        + eye( 2*M^2 ) - albedo * L;

Psi = Ap \ ( albedo * src0 );

% compute the homogeneous solution

[ lambda, V ] = PWSolns( xk, 0, L, albedo, mu, phi, M );

% compute the expansion coefficients by imposing boundary condition

Ah = V(ineg(M^2+1:2*M^2),:) - diag( Rbc ) * V(ipos(M^2+1:2*M^2),:);
b  = - Psi(ipos(M^2+1:2*M^2)) + diag( Rbc ) * Psi(ineg(M^2+1:2*M^2));
c  = Ah \ b;

%% INTERPOLATION OVER THE BACK HEMISPHERE JUST INSIDE THE MEDIUM

% compute angle grids

Ngrid       = 101;
muout_grid  = linspace( -1, 0, Ngrid );
phi_grid    = linspace( pi/2, 3*pi/2, Ngrid );

% compute meshgrid of angles out of the medium

[ MUout_grid, PHI_grid ] = meshgrid( muout_grid, phi_grid );

% stretch the meshgrids to vectors

MUout_grid  = MUout_grid(:);
PHI_grid    = PHI_grid(:);

% compute the grid of angles in the medium

MUin_grid = -sqrt( 1 - ( 1 - MUout_grid.^2 ) / nrel^2 );

% compute the interpolant

I = SourceIntegrationMethod( MU, PHI, WT, MUin_grid, PHI_grid, ...
    mu0, V(ineg,:), lambda, c, Psi, albedo, g, xk );

% compute Fresnel transmission coefficient

FresnelT = ( 2 * nrel * MUin_grid ./ ( MUin_grid + nrel * MUout_grid ) ).^2 ...
    + ( 2 * nrel * MUin_grid ./ ( nrel * MUin_grid + MUout_grid ) ).^2;

FresnelT = 0.5 * FresnelT .* MUout_grid / nrel^3 ./ MUin_grid;

% multiply radiance by Fresnel transmission coefficient and reshape

I = reshape( FresnelT .* I, Ngrid, Ngrid );

%% compute mu-star_out

t           = linspace( pi/2, 3*pi/2, Ngrid );
t           = t(2:end-1);
theta_star  = atan( tan(theta0) ./ cos(t) );

% critical angle

indx = find( -sin(theta_star) <= 1/nrel );

t_star      = t(indx);
mu_star_out = -cos( asin( nrel * sin( theta_star(indx) ) ) );

%% plot the result

figure(6)
pcolor( phi_grid, muout_grid, abs( I )' );
xticks([pi/2 3*pi/4 pi 5*pi/4 3*pi/2]);
xticklabels({'$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$','$3\pi/2$'});
set(gca,'TickLabelInterpreter', 'latex')
shading flat;
colorbar;
hold on;
plot( t_star, mu_star_out, 'r:' );
hold off;
xlabel( '$\varphi$', 'Interpreter', 'LaTeX' );
ylabel( '$\mu_{\mbox{out}}$', 'Interpreter', 'LaTeX' );
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% axis tight