function [correct] = check(grid1,grid2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
assert(all(size(grid1) == size(grid2)));
atol = 1e-4;
correct = sum(sum(bsxfun(@(x,y) isapprox(x,y,atol),grid1,grid2))) == numel(grid1);
avg_L1error = sum(sum(abs(grid1-grid2)))/numel(grid1);
max_L1error = max(max(abs(grid1-grid2)));
rel_errors = abs((grid1-grid2)./grid1);
[max_relerror, I] = max(rel_errors(:));
[I_row, I_col] = ind2sub(size(rel_errors), I);
disp(strcat("avg L1 error: ", num2str(avg_L1error)));
disp(strcat("max L1 error: ", num2str(max_L1error)));
disp(strcat("max rel error [%]: ", num2str(max_relerror*100), " at: ", num2str([I_row, I_col])));
disp([grid1(I_row, I_col) grid2(I_row, I_col)]);
end