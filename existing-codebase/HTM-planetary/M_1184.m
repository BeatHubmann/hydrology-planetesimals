function [EXY_jl,SXY_jl,DSXY_jl,EXX_jl,SXX_jl,DSXX_jl,EII_jl,...
    SII_jl,ETAP_jl,GGGP_jl,SXX0_jl,ETA_jl,GGG_jl,SXY0_jl,vx_jl,vy_jl] = M_1184(step)
EXY_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/EXY');
SXY_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/SXY');
DSXY_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/DSXY');
EXX_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/EXX');
SXX_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/SXX');
DSXX_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/DSXX');
EII_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/EII');
SII_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/SII');
ETAP_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/ETAP');
GGGP_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/GGGP');
SXX0_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/SXX0');
ETA_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/ETA');
GGG_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/GGG');
SXY0_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/SXY0');
vx_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/vx');
vy_jl = h5read(['/Users/z7717/Desktop/test/M_1184_' num2str(step) '.jld2'], '/vy');
end