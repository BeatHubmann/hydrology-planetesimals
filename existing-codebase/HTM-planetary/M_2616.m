function [xm_jl,ym_jl,tm_jl,tkm_jl,phim_jl,sxxm_jl,sxym_jl,...
    etavpm_jl,marknum_jl,dt_jl,dtm_jl,timesum_jl] = M_2616(step)
xm_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/xm')';
ym_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/ym')';
tm_jl = double(h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/tm'))';
tkm_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/tkm')';
phim_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/phim')';
sxxm_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/sxxm')';
sxym_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/sxym')';
etavpm_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/etavpm')';
marknum_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/marknum')';
dt_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/dt');
dtm_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/dtm');
timesum_jl = h5read(['/Users/z7717/Desktop/test/M_2616_' num2str(step) '.jld2'], '/timesum');
end