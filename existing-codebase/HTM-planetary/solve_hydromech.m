function [vx,vy,pr,qxD,qyD,pf,S] = solve_hydromech(L,R,pscale)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
S=L\R; % Obtaining algebraic vector of solutions S()
Nx1 = sqrt(size(L,2)/6);
Ny1 = sqrt(size(L,1)/6);
vx = zeros(Ny1, Nx1);
vy = zeros(Ny1, Nx1);
pr = zeros(Ny1, Nx1);
qxD = zeros(Ny1, Nx1);
qyD = zeros(Ny1, Nx1);
pf = zeros(Ny1, Nx1);
% Reload solutions S() to vx(), vy(), p()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; % Vx solid
        kvy=kvx+1; % Vy solid
        kpm=kvx+2; % Ptotal
        kqx=kvx+3; % qx Darcy
        kqy=kvx+4; % qy Darcy
        kpf=kvx+5; % P fluid
        % Reload solution
        vx(i,j)=S(kvx);
        vy(i,j)=S(kvy);
        pr(i,j)=S(kpm)*pscale;
        qxD(i,j)=S(kqx);
        qyD(i,j)=S(kqy);
        pf(i,j)=S(kpf)*pscale;
    end
end
end