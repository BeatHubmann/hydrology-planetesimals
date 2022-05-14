function [xm_jl,ym_jl,sxxm_jl,sxym_jl] = M_2235(step)
xm_jl = h5read(['/Users/z7717/Desktop/test/M_2235_' num2str(step) '.jld2'], '/xm')';
ym_jl = h5read(['/Users/z7717/Desktop/test/M_2235_' num2str(step) '.jld2'], '/ym')';
sxxm_jl = h5read(['/Users/z7717/Desktop/test/M_2235_' num2str(step) '.jld2'], '/sxxm')';
sxym_jl = h5read(['/Users/z7717/Desktop/test/M_2235_' num2str(step) '.jld2'], '/sxym')';
end