function [ETA0_jl,ETA_jl,YNY_jl,GGG_jl,SXY0_jl,COH_jl,TEN_jl,FRI_jl,...
    RHOX_jl,RHOFX_jl,KX_jl,PHIX_jl,RX_jl,RHOY_jl,RHOFY_jl,KY_jl,PHIY_jl,...
    RY_jl,GGGP_jl,SXX0_jl,RHO_jl,RHOCP_jl,ALPHA_jl,ALPHAF_jl,HR_jl,...
    PHI_jl,BETTAPHI_jl,tk1_jl,xm_jl,ym_jl,tm_jl,phim_jl,sxxm_jl,sxym_jl,...
    etafluidcur_inv_kphim_jl,ktotalm_jl,KXSUM_jl,KYSUM_jl,TKSUM_jl,...
    WTXSUM_jl,WTYSUM_jl,WTPSUM_jl] = M_680(step)
%M_680 Summary of this function goes here
%   Detailed explanation goes here
    ETA0_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/ETA0');
    ETA_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/ETA');
    YNY_jl = double(h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/YNY'));
    GGG_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/GGG');
    SXY0_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/SXY0');
    COH_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/COH');
    TEN_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/TEN');
    FRI_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/FRI');
    RHOX_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RHOX');
    RHOFX_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RHOFX');
    KX_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/KX');
    PHIX_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/PHIX');
    RX_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RX');
    RHOY_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RHOY');
    RHOFY_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RHOFY');
    KY_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/KY');
    PHIY_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/PHIY');
    RY_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RY');
    GGGP_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/GGGP');
    SXX0_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/SXX0');
    RHO_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RHO');
    RHOCP_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/RHOCP');
    ALPHA_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/ALPHA');
    ALPHAF_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/ALPHAF');
    HR_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/HR');
    PHI_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/PHI');
    BETTAPHI_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/BETTAPHI');
    tk1_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/tk1');
    xm_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/xm')';
    ym_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/ym')';
    tm_jl = double(h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/tm')');
    phim_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/phim')';
    sxxm_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/sxxm')';
    sxym_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/sxym')';
    etafluidcur_inv_kphim_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/etafluidcur_inv_kphim')';
    ktotalm_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/ktotalm')';
    KXSUM_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/KXSUM');
    KYSUM_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/KYSUM');
    TKSUM_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/TKSUM');
    WTXSUM_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/WTXSUM');
    WTYSUM_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/WTYSUM');
    WTPSUM_jl = h5read(['/Users/z7717/Desktop/test/M_680_' num2str(step) '.jld2'], '/WTPSUM');
end