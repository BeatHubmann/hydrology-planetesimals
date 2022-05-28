function [] = display_hydromech_index(L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Nx1 = sqrt(size(L,2)/6);
Ny1 = sqrt(size(L,1)/6);
for j=1:1:Nx1
    for i=1:1:Ny1
        kvx=((j-1)*Ny1+i-1)*6+1; % Vx solid
        kvy=kvx+1; % Vy solid
        kpm=kvx+2; % Ptotal
        kqx=kvx+3; % qx Darcy
        kqy=kvx+4; % qy Darcy
        kpf=kvx+5; % Pfluid    
        disp(strcat("i=",num2str(i)," j=",num2str(j), ...
                    " / vx:", num2str(kvx),...
                    " vy:", num2str(kvy),...
                    " pr:", num2str(kpm),...
                    " qx:", num2str(kqx),...
                    " qy:", num2str(kqy),...
                    " pf:", num2str(kpf)));
    end
end
end