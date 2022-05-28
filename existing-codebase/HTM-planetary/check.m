function [correct] = check(grid1,grid2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
assert(all(size(grid1)==size(grid2)));
[ysize, xsize] = size(grid1);
maximum = max(max(max(grid1)), max(max(grid2)));
eps_at_max = eps(maximum);
assert(all(size(grid1) == size(grid2)));
atol = 1e-4;
correct = sum(sum(bsxfun(@(x,y) isapprox(x,y,atol),grid1,grid2))) == numel(grid1);
L1_errors = abs(grid1-grid2);
avg_L1_error = sum(sum(L1_errors))/numel(grid1);
rel_errors = abs((grid1-grid2)./grid1);
rel_errors(isinf(rel_errors)) = 0;
[max_L1_error, I_L1] = max(L1_errors(:));
[I_row_L1, I_col_L1] = ind2sub(size(L1_errors), I_L1);
[max_rel_error, I_rel] = max(rel_errors(:));
[I_row_rel, I_col_rel] = ind2sub(size(rel_errors), I_rel);
disp(strcat("eps(", num2str(maximum),") = ", num2str(eps_at_max)));
disp(strcat("avg L1 error: ", num2str(avg_L1_error)));
disp(strcat("max L1 error: ", num2str(max_L1_error), " at: ", num2str([I_row_L1, I_col_L1])));
disp([grid1(I_row_L1, I_col_L1) grid2(I_row_L1, I_col_L1)]);
disp(strcat("max rel error [%]: ", num2str(max_rel_error*100), " at: ", num2str([I_row_rel, I_col_rel])));
disp([grid1(I_row_rel, I_col_rel) grid2(I_row_rel, I_col_rel)]);
ymin = max(1, I_row_L1-32);
ymax = min(ysize, I_row_L1+32);
xmin = max(1, I_col_L1-32);
xmax = min(xsize, I_col_L1+32);
if (1<ysize && 1<xsize)
    heatmap(abs(grid1(ymin:ymax,xmin:xmax)-grid2(ymin:ymax,xmin:xmax)),...
        'XData', xmin:xmax, 'YData', ymin:ymax);
    
end
end