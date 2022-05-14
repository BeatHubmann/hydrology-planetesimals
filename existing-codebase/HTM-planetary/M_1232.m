function [SXX_jl,APHI_jl,PHI_jl,pr_jl,pf_jl,ps_jl] = M_1232(step)
SXX_jl = h5read(['/Users/z7717/Desktop/test/M_1232_' num2str(step) '.jld2'], '/SXX');
APHI_jl = h5read(['/Users/z7717/Desktop/test/M_1232_' num2str(step) '.jld2'], '/APHI');
PHI_jl = h5read(['/Users/z7717/Desktop/test/M_1232_' num2str(step) '.jld2'], '/PHI');
pr_jl = h5read(['/Users/z7717/Desktop/test/M_1232_' num2str(step) '.jld2'], '/pr');
pf_jl = h5read(['/Users/z7717/Desktop/test/M_1232_' num2str(step) '.jld2'], '/pf');
ps_jl = h5read(['/Users/z7717/Desktop/test/M_1232_' num2str(step) '.jld2'], '/ps');
end