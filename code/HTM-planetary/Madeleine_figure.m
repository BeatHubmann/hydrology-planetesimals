% Jpeg files
nnamefile='madcph_'; % jpeg filename Fig1
for timefile=100:100:1111450
namefile=[nnamefile,num2str(timefile)];
load(namefile);
nname1 =[nnamefile,'a_'];
nname2 =[nnamefile,'b_'];
nname3 =[nnamefile,'c_'];

MARKPERCELL=zeros(Nx-1,Ny-1);
for m=1:1:marknum
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    MARKPERCELL(i,j)=MARKPERCELL(i,j)+1;  
end


avrpor=0;
numpor=0;
for m=1:1:marknum
    if(tm(m)<3)
        j=fix((xm(m)-x(1))/dx)+1;
        i=fix((ym(m)-y(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        MARKPERCELL(i,j)=MARKPERCELL(i,j)+1;  
        wtm=1/MARKPERCELL(i,j);
        avrpor=avrpor+phim(m)*wtm;
        numpor=numpor+wtm;
    end
end
avrpor=avrpor/numpor

DIVV=zeros(Ny1,Nx1);
gmag=zeros(Ny1,Nx1);
for i=2:1:Ny
  for j=2:1:Nx
    % EXX, SXX, DSXX
    DIVV(i,j)=(vx(i,j)-vx(i,j-1))/dx+(vy(i,j)-vy(i-1,j))/dy;
    gmag(i,j)=(((gx(i,j)+gx(i,j-1))/2)^2+((gy(i,j)+gy(i-1,j))/2)^2)^0.5;
  end
end



% 8) Visialize RHO() P(), vx(), vy() 
figure(1);clf;colormap('Jet')
subplot(3,3,1)
pcolor(xp,yp,RHO);
colorbar; 
shading interp;
axis ij image;
title(['RHO    step = ',num2str(timestep),'   time = ',num2str((timesum+dt)/(365.25*24*3600*1e+6)),' Myr'])

subplot(3,3,2)
pcolor(x,y,log10(ETA));
colorbar; % 
shading interp; 
axis ij image; 
title(['log10 viscosity, Pa*s  step = ',num2str(timestep),'   dtelastic = ',num2str(dt/(365.25*24*3600)),' yr','   dt = ',num2str(dtm/(365.25*24*3600)),' yr'])

subplot(3,3,3)
pcolor(xvx,yvx,vx);
colorbar; 
shading interp;
axis ij image;
title('vxs, m/s')

subplot(3,3,4)
pcolor(xvy,yvy,vy); 
colorbar;
shading interp;
axis ij image;
title('vys, m/s')

subplot(3,3,5)
pcolor(xp/1000,yp/1000,log10(PHI./(1-PHI)))
shading interp;
axis ij image;
colorbar
title('logPHI/(1-PHI)')


subplot(3,3,6)
pcolor(xp(2:Nx),yp(2:Ny),log10(EII(2:Ny,2:Nx)));
colorbar;
shading interp;
axis ij image;
title('log10 EII, 1/s')

subplot(3,3,7)
pcolor(xp(2:Nx),yp(2:Ny),SII(2:Ny,2:Nx));
colorbar;
shading interp;
axis ij image;
title('SII, Pa')

subplot(3,3,8)
pcolor(xp(2:Nx),yp(2:Ny),gmag(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('gmag, m/s^2')

subplot(3,3,9)
pcolor(xp/1000,yp/1000,tk2)
shading interp;
axis ij image;
colorbar
title('T, K')


    namejpg    =  [nname1,num2str(timestep)];
    print ('-djpeg', '-r300',namejpg);

    
% 8) Visialize RHO() P(), vx(), vy() 
figure(2);clf;colormap('Jet')
subplot(4,4,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10 viscosity, Pa*s  step = ',num2str(timestep),'   dtelastic = ',num2str(dt/(365.25*24*3600)),' yr','   dt = ',num2str(dtm/(365.25*24*3600)),' yr'])

hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(4,4,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title(['Ptotal, Pa    step = ',num2str(timestep),'   time = ',num2str((timesum+dt)/(365.25*24*3600*1e+6)),' Myr'])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,3)
pcolor(xp/1000,yp/1000,vxp)
shading interp;
axis ij image;
colorbar
title('vx, m/s')
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,4)
pcolor(xp/1000,yp/1000,vyp)
shading interp;
axis ij image;
colorbar
title('vy, m/s')
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,5)
pcolor(xp/1000,yp/1000,HS)
shading interp;
axis ij image;
colorbar
title('HS, W/m^3')

subplot(4,4,6)
pcolor(xp/1000,yp/1000,HA)
shading interp;
axis ij image;
colorbar
title('HA, W/m^3')

subplot(4,4,7)
pcolor(xp/1000,yp/1000,pr-pf)
shading interp;
axis ij image;
colorbar
title('Ptotal-Pfluid, kg/m^3')

subplot(4,4,8)
pcolor(xvx/1000,yvx/1000,log10(KX))
shading interp;
axis ij image;
colorbar
title('logK, W/m/K')

subplot(4,4,9)
pcolor(xp/1000,yp/1000,tk2)
shading interp;
axis ij image;
colorbar
title('T, K')

subplot(4,4,10)
pcolor(xp/1000,yp/1000,log10(PHI./(1-PHI)))
shading interp;
axis ij image;
colorbar
title('logPHI/(1-PHI)')

subplot(4,4,11)
pcolor(xp/1000,yp/1000,pf)
shading interp;
axis ij image;
colorbar
title('Pfluid, Pa')

subplot(4,4,12)
pcolor(xvx/1000,yvx/1000,qxD)
shading interp;
axis ij image;
colorbar
title('qxD, m/s')

subplot(4,4,13)
pcolor(xvy/1000,yvy/1000,qyD)
shading interp;
axis ij image;
colorbar
title('qyD, m/s')

subplot(4,4,14)
pcolor(xvx/1000,yvx/1000,log10(RX))
shading interp;
axis ij image;
colorbar
title('logRX')

subplot(4,4,15)
pcolor(xp/1000,yp/1000,log10(ETAPHI))
shading interp;
axis ij image;
colorbar
title('logETAphi, Pa*s')

subplot(4,4,16)
pcolor(xp,yp,RHO);
colorbar; 
shading interp;
axis ij image;
title(['RHO, kg/m^3'])


    namejpg    =  [nname2,num2str(timestep)];
    print ('-djpeg', '-r300',namejpg);
    
    
% Add markers to empty areas
dxm=dxm/2;
dym=dym/2;
Nxm=Nxm*2;
Nym=Nym*2;
marknumold=marknum
mdis=1e30*ones(Nym,Nxm);
mnum=zeros(Nym,Nxm);
mtyp=zeros(Nym,Nxm);
mpor=zeros(Nym,Nxm);
xxm=dxm/2:dxm:xsize-dxm/2;
yym=dym/2:dym:ysize-dym/2;
for m=1:1:marknum
    
    % Check markers with the nearest nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xxm(1))/dxm)+1;
    i=fix((ym(m)-yym(1))/dym)+1;
    if(j<1)
        j=1;
    elseif(j>Nxm-1)
        j=Nxm-1;
    end
    if(i<1)
        i=1;
    elseif(i>Nym-1)
        i=Nym-1;
    end
    Nxm0=(Nx-1)*Nxmc;
    Nym0=(Ny-1)*Nymc;
    i1(m)=fix((m-fix(m/Nym0)*Nym0)/Nymc);
    j1(m)=fix(m/(Nym0*Nxmc));
    tm(m)=tm(m)+(i1(m)-fix(i1(m)/2)*2)-(j1(m)-fix(j1(m)/2)*2);
    % Check nodes
    % i,j Node
    % Compute distance
    dxmj=xm(m)-xxm(j);
    dymi=ym(m)-yym(i);
    dismij=(dxmj^2+dymi^2)^0.5;
    if(dismij<mdis(i,j))
        mdis(i,j)=dismij;
        mnum(i,j)=m;
        mtyp(i,j)=tm(m);
        mpor(i,j)=phim(m);
    end
    % i+1,j Node
    % Compute distance
    dxmj=xm(m)-xxm(j);
    dymi=ym(m)-yym(i+1);
    dismi1j=(dxmj^2+dymi^2)^0.5;
    if(dismi1j<mdis(i+1,j))
        mdis(i+1,j)=dismi1j;
        mnum(i+1,j)=m;
        mtyp(i+1,j)=tm(m);
        mpor(i+1,j)=phim(m);
    end
    % i,j+1 Node
    % Compute distance
    dxmj=xm(m)-xxm(j+1);
    dymi=ym(m)-yym(i);
    dismij1=(dxmj^2+dymi^2)^0.5;
    if(dismij1<mdis(i,j+1))
        mdis(i,j+1)=dismij1;
        mnum(i,j+1)=m;
        mtyp(i,j+1)=tm(m);
        mpor(i,j+1)=phim(m);
    end
    % i+1,j+1 Node
    % Compute distance
    dxmj=xm(m)-xxm(j+1);
    dymi=ym(m)-yym(i+1);
    dismi1j1=(dxmj^2+dymi^2)^0.5;
    if(dismi1j1<mdis(i+1,j+1))
        mdis(i+1,j+1)=dismi1j1;
        mnum(i+1,j+1)=m;
        mtyp(i+1,j+1)=tm(m);
        mpor(i+1,j+1)=phim(m);
    end
end

% 8) Visialize markers 
figure(3);clf;colormap('Jet')
pcolor(xxm,yym,mtyp);% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['Marker type  step = ',num2str(timestep),'   time = ',num2str((timesum+dt)/(365.25*24*3600*1e+6)),' Myr'])


end
