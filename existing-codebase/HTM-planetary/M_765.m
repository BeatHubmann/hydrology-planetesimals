function [SP_jl,FI_jl,gx_jl,gy_jl,dt_jl,dtelastic_jl,ETA_jl,ETA00_jl,...
    YNY_jl,YNY00_jl] = M_765(step)
%M_680 Summary of this function goes here
%   Detailed explanation goes here
SP_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/SP');
FI_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/FI');
gx_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/gx');
gy_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/gy');
dt_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/dt');
dtelastic_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/dtelastic');
ETA_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/ETA');
ETA00_jl = h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/ETA00');
YNY_jl = double(h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/YNY'));
YNY00_jl = double(h5read(['/Users/z7717/Desktop/test/M_765_' num2str(step) '.jld2'], '/YNY00'));
end