function [HS_jl,HA_jl] = M_1622(step)
HS_jl = h5read(['/Users/z7717/Desktop/test/M_1622_' num2str(step) '.jld2'], '/HS');
HA_jl = h5read(['/Users/z7717/Desktop/test/M_1622_' num2str(step) '.jld2'], '/HA');
end