function [Z,P] = orthonorm(Y)
%ORTHONORM orthonormalization by Gram-Schmidt process
%   output:
%         Z = matrix of orthonormal column vectors
%         P = upper triangle matrix such that Z = YP
%   input:
%         Y = matrix of linearly independent column vector solutions [y1,
%         y2, y3, ...]

    % properly index
[rows,cols] = size(Y);
    % preallocate
Z = zeros( size(Y) );
P = zeros(cols,cols);
W = zeros(1,cols);

w1 = sqrt( Y(:,1).' * conj(Y(:,1)) );
W(1) = w1;
Z(:,1) = Y(:,1) / w1;
for i = 2:cols
    proj = 0;
    ynew = Y(:,i);
    for j = 1:i-1
        z_old = Z(:,j);
        proj = proj + (ynew.' * conj(z_old)) * z_old;
    end
    t = ynew - proj; w = sqrt( t.' * conj(t)); z = t / w;
    Z(:,i) = z; W(i) = w;
end

% for i = 1:cols
%     for j = 1:cols
%         if i < j
%             P(i,j) = - ( Z(:,i).' * conj(Y(:,j))) / ( W(i) * W(j) );
%         elseif i == j
%             P(i,j) = 1 / W(j);
%         else
%             P(i,j) = 0;
%         end
%     end
% end
P = Y \ Z;
end