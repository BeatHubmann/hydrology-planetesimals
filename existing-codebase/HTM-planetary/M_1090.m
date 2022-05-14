function [APHI_jl,pr_jl,pr0_jl,pf_jl,pf0_jl,PHI_jl,ETAPHI_jl,...
    BETTAPHI_jl,dt_jl,aphimax_jl] = M_1090(step)
APHI_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/APHI');
pr_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/pr');
pr0_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/pr0');
pf_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/pf');
pf0_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/pf0');
PHI_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/PHI');
ETAPHI_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/ETAPHI');
BETTAPHI_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/BETTAPHI');
dt_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/dt');
aphimax_jl = h5read(['/Users/z7717/Desktop/test/M_1090_' num2str(step) '.jld2'], '/aphimax');
end