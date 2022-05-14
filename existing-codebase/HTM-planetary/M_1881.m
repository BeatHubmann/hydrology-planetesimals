function [phim_jl] = M_1881(step)
phim_jl = h5read(['/Users/z7717/Desktop/test/M_1881_' num2str(step) '.jld2'], '/phim')';
end