function [] = lookup_hydromech_index(ii, jj, L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Nx1 = sqrt(size(L,2)/6);
Ny1 = sqrt(size(L,1)/6);
switch mod(ii, 6)
    case 0
        disp(strcat("row: pf ",num2str(ii)));
        kvx=ii-5;
    case 1
        disp(strcat("row: vx ",num2str(ii)));
        kvx=ii;
    case 2
        disp(strcat("row: vy ",num2str(ii)));
        kvx=ii-1;
    case 3
        disp(strcat("row: pr ",num2str(ii)));
        kvx=ii-2;
    case 4
        disp(strcat("row: qx ",num2str(ii)));
        kvx=ii-3;
    case 5
        disp(strcat("row: qy ",num2str(ii)));
        kvx=ii-4;
end
switch mod(jj, 6)
    case 0
        disp(strcat("col: pf ",num2str(jj)));
    case 1
        disp(strcat("col: vx ",num2str(jj)));
    case 2
        disp(strcat("col: vy ",num2str(jj)));
    case 3
        disp(strcat("col: pr ",num2str(jj)));
    case 4
        disp(strcat("col: qx ",num2str(jj)));
    case 5
        disp(strcat("col: qy ",num2str(jj)));
end
jjj=-1;
for j=1:1:Nx1
    for i=1:1:Ny1
        if (kvx==((j-1)*Ny1+i-1)*6+1)
            jjj=j;
            break
        end
    end
    if jjj>0
        break
    end
end
disp(strcat("i=",num2str(i)," j=",num2str(j)));
end