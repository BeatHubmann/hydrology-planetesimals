function [data_1, data_2] = load_mat(data, step)
%load_mat Load two versions of data from .m and .jld2 at step
    path_1 = '.\';
    name_1 = 'planetary-NGD_';
    path_2 = 'C:\Users\ich\outTest\';
    name_2 = 'M_1078_';
    root_file_1 = [path_1 name_1 num2str(step) '.mat'];
    root_file_2 = [path_2 name_2 num2str(step) '.jld2'];
    data_1 = load(root_file_1, data);
    data_2 = h5read(root_file_2, ['/' data]);
end