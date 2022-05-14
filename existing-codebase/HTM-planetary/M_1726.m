function [tk0_jl,tk1_jl,tk2_jl,DT_jl,DT0_jl,RHOCP_jl,KX_jl,KY_jl,HR_jl,...
    HA_jl,HS_jl,LT_jl,RT_jl,ST_jl,dtm_jl,dtt_jl] = M_1726(step, titer)
tk0_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/tk0');
tk1_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/tk1');
tk2_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/tk2');
DT_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/DT');
DT0_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/DT0');
RHOCP_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/RHOCP');
KX_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/KX');
KY_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/KY');
HR_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/HR');
HA_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/HA');
HS_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/HS');
LT_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/LT_d');
RT_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/RT');
ST_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/ST');
dtm_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/dtm');
dtt_jl = h5read(['/Users/z7717/Desktop/test/M_1726_' num2str(step) '_' num2str(titer) '.jld2'], '/dtt');
end