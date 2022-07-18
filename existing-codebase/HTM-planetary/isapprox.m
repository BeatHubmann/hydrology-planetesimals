function [C] = isapprox(A, B, atol)
% isapprox Compare two matrices A, B for approximate identity
% from https://ch.mathworks.com/matlabcentral/answers/97069-how-can-i-compare-numbers-for-equality-within-a-tolerance-in-matlab-8-0-r2012b
% See also the 2015 function http://www.mathworks.com/help/matlab/ref/ismembertol.html
assert(all(size(A) == size(B)));
C = false(size(A));
for j=1:size(A, 2)
    for i=1:size(A, 1)
        C(i, j) = abs(A(i,j)-B(i,j)) < max(1e4*eps(min(abs(A(i,j)),abs(B(i,j)))), atol);
    end
end
end