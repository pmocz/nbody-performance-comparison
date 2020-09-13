function [ a ] = getAcc_vec( pos, mass, G, softening )
%GETACC Calculate the acceleration on each particle due to Newton's Law 
%   pos  is an N x 3 matrix of positions
%   mass is an N x 1 vector of masses
%   G is Newton's Gravitational constant
%   softening is the softening length
%   a is N x 3 matrix of accelerations

% positions r = [x,y,z] for all particles
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);

% matrix that stores all pairwise particle separations: r_j - r_i
dx = x' - x;
dy = y' - y;
dz = z' - z;

% matrix that stores 1/r^3 for all particle pairwise particle separations 
inv_r3 = (dx.^2 + dy.^2 + dz.^2 + softening.^2).^(-3/2);

ax = G * (dx .* inv_r3) * mass;
ay = G * (dy .* inv_r3) * mass;
az = G * (dz .* inv_r3) * mass;

% pack together the acceleration components
a = [ax ay az];


end

