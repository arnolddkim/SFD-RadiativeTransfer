%% ComputeSphericalHarmonics.m
%
% This code computes a matrix whose columns are evaluations of the 
% spherical harmonics evaluated on the grid points supplied as inputs up 
% to order N - 1 with N supplied as an input.
%
% This code calls Compute_Ynm.m.
%
% Written by A. D. Kim on 4/15/2018. Modified most recently on 4/11/2021

function [ Ynm, nvec ] = ComputeSphericalHarmonics( N, MU, PHI )

% allocate memory for the output

nvec = zeros( N^2, 1 );
Ynm  = zeros( length(MU), N^2 );

% initialize the counter

j = 0;

for n = 0 : N - 1
    
    for m = -n : n

        j        = j + 1;
        nvec(j)  = n;
        Ynm(:,j) = Compute_Ynm( n, m, MU, PHI );
        
    end
    
end

return
