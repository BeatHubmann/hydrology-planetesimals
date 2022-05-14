function [C] = isapprox(A, B, atol)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
assert(all(size(A) == size(B)));
C = false(size(A));
for j=1:size(A, 2)
    for i=1:size(A, 1)
        C(i, j) = abs(A(i,j)-B(i,j)) < max(1e4*eps(min(abs(A(i,j)),abs(B(i,j)))), atol);
    end
end
end