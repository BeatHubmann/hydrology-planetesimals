function [result_1, result_2] = compare_mat(func, data, range)
%compare_mat compare two data series wrt to func over all steps in a range
    path_1 = '.\';
    name_1 = 'planetary-NGD_';
    path_2 = 'C:\Users\ich\outTest\';
    name_2 = 'M_1078_';
    result_1 = zeros(length(aaa));
    result_2 = zeros(length(aaa));
    for step = range
        root_file_1 = [path_1 name_1 num2str(step) '.mat'];
        root_file_2 = [path_2 name_2 num2str(step) '.jld2'];
        data_1 = load(root_file_1, data);
        data_2 = h5read(root_file_2, ['/' data];
        result_1 = func(data_1);
        result_2 = func(data_2);
    end
end