function [ a ] = getAcc_for( pos, mass, G, softening )
%GETACC Calculate the acceleration on each particle due to Newton's Law
%   pos  is an N x 3 matrix of positions
%   mass is an N x 1 vector of masses
%   G is Newton's Gravitational constant
%   softening is the softening length
%   a is N x 3 matrix of accelerations

N = size(pos,1);
a = zeros(N,3);

for i = 1:N
    for j = 1:N
        dx = pos(j,1) - pos(i,1);
        dy = pos(j,2) - pos(i,2);
        dz = pos(j,3) - pos(i,3);
        inv_r3 = (dx^2 + dy^2 + dz^2 + softening^2)^(-1.5);
        a(i,1) = a(i,1) + G * (dx * inv_r3) * mass(j);
        a(i,2) = a(i,2) + G * (dy * inv_r3) * mass(j);
        a(i,3) = a(i,3) + G * (dz * inv_r3) * mass(j);
    end
end


end

