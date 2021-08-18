function Ellipsoid_Initialize(Ep,LO,LI,dL,dI,IRes,nProtons,gaussParam,noiseParam,a,b,B0)
% Initialize proton images from Gaussian ellipsoid of magnetic field

%Constants
c = 2.998*10^8;                 % Speed of light (m/s)
e = 1.602*10^-19;               % Elementary charge (C)
e0 = 8.854*10^-12;              % Vacuum permittivity (F/m)
me = 0.5110*10^6;               % Rest energy of electron (eV/c^2)
Me = 9.109*10^-31;
mp = 1.67e-27;      %proton mass (kg)
Mp = 938.3*10^6;    %proton rest mass (eV)

%Probe proton parameters
% Ep3 = 3e6;
% Ep14 = 14.7e6;
% g3 = Ep3/Mp+1;
% g14 = Ep14/Mp+1;
% v3 = c*sqrt(1-g3^-2);
% v14 = c*sqrt(1-g14^-2);

%System properties
% R = 0.38e-3;        %radius of wire (m)
% r0 = 0.0e-3;        %lateral position of proton source (m)
% B0 = 13.5;            %maximum magnetic field at wire surface (T)
% LO = 1e-2;          %distance from proton source to system center (m)
% LI = 15e-2;         %distance from proton source to image plane (m)
LOI = LI-LO;        %distance from system center to image plane (m)
% dL = 5e-3;          %length of the interaction region (m)
% dI = 10e-2;         %width of image plane (m)

%Define image properties
% IRes = 400;
x2D = linspace(-dI/2,dI/2,IRes);%2*length(x2)-1);
dx = x2D(2)-x2D(1);
xedges = [x2D-dx/2 dI/2+dx/2];
[xmesh,zmesh] = meshgrid(x2D,x2D);
[xInd,yInd] = ndgrid(x2D,x2D);
rmesh = sqrt(xmesh.^2+zmesh.^2);

% gaussParam = 2;

%Angles from proton source to cover
alphamin = atan2(0,LO); 	%Not sure what to do with this now
alphamax = atan2(dI/2,LI);  %Goes to the edge of the image
alphares = IRes;           %Number of protons to initialize for the true image
alpha0 = linspace(-alphamax,alphamax,alphares); %Initializing initial angles

%Azimuthal angle
theta0 = linspace(-pi,pi,360*2);
[alpha0mesh,theta0mesh] = meshgrid(alpha0,theta0);

x0 = linspace(-LO*tan(alphamax),LO*tan(alphamax),IRes);
z0 = linspace(-LO*tan(alphamax),LO*tan(alphamax),IRes);
[x0mesh,z0mesh] = meshgrid(x0,z0);
[x0nd,z0nd] = ndgrid(x0,z0);

r0mesh = sqrt(x0mesh.^2+z0mesh.^2);
alpha0mesh = atan2(r0mesh,LO);
theta0mesh = atan2(z0mesh,x0mesh);

r0nd = sqrt(x0nd.^2+z0nd.^2);
alpha0nd = atan2(r0nd,LO);
theta0nd = atan2(z0nd,x0nd);

%% Integrated magnetic field (Gaussian ellipsoid)

%Ellipsoid properties
% a = 0.75*1e-3;      %ellipsoid radius (m)
% b = 2e-3;           %ellipsoid length (m)
% B0 = -70;           %ellipsoid maximum magnetic field (may not actually be maximum?)
xc0 = 0*1e-3;       %ellipsoid center (x) relative to the object plane
yc0 = 0;            %ellipsoid center (y)
zc0 = 0*1e-3;       %ellipsoid center (z)
xphi = 0*pi/2+pi/2; %ellipsoid orientation (pi/2 to probe parallel)

%Define functional form of the vector potential
AZfun = @(x,y,z,xc,yc,zc,LO,a,b,B,xphi) ...
    B.*a/2.*exp(-((x-xc)./a).^2-(((y-yc-LO)*cos(xphi)+(z-zc)*sin(xphi))./a).^2-...
    ((-(y-yc-LO)*sin(xphi)+(z-zc)*cos(xphi))./b).^2);

%Define B structure
Bf = @(x,y,z,xc,yc,zc,a,b) B0*exp(-(x-xc).^2./a.^2-(y-yc).^2./a^2-(z-zc).^2./b.^2);

%Initialize zero integrated vector potential meshes
IY_Div = 0*x0nd;
IZ_Div = IY_Div;
IY_Coll = 0*x0nd;
IZ_Coll = IY_Coll;

xc = xc0;
yc = yc0*cos(xphi)+zc0*sin(xphi);
zc = -yc0*sin(xphi)+zc0*cos(xphi);

%Diverging
IZfun_Div = integral(@(y) cos(xphi)*AZfun(x0mesh.*y./LO,y,z0mesh.*y./LO,xc,yc,zc,LO,a,b,B0,xphi),...
    0,LO+dL,'ArrayValued',true);
IYfun_Div = integral(@(y) sin(xphi)*AZfun(x0mesh.*y./LO,y,z0mesh.*y./LO,xc,yc,zc,LO,a,b,B0,xphi),...
    0,LO+dL,'ArrayValued',true);

IY_Div = IY_Div+IYfun_Div;
IZ_Div = IZ_Div+IZfun_Div;

%Collimated
IZfun_Coll = integral(@(y) cos(xphi)*AZfun(x0nd,y,z0nd,xc,yc,zc,LO,a,b,B0,xphi),...
    0,LO+dL,'ArrayValued',true);
IYfun_Coll = integral(@(y) sin(xphi)*AZfun(x0nd,y,z0nd,xc,yc,zc,LO,a,b,B0,xphi),...
    0,LO+dL,'ArrayValued',true);

IY_Coll = IY_Coll+IYfun_Coll;
IZ_Coll = IZ_Coll+IZfun_Coll;

% end

figure(5)
clf
imagesc(x0*1e3,z0*1e3,-abs(IY_Div/1e-6))
axis equal tight
set(gca,'ydir','normal','FontSize',40,'PlotBoxAspectRatio',[1 1 1])
% caxis([0,max(abs(IY_Div(:)/1e-6))])
xticks(-4:1:4)
yticks(-4:1:4)
xlabel('$r_0$ (mm)','Interpreter','LaTeX','FontSize',40)
ylabel('$r_I$ (cm)','Interpreter','LaTeX','FontSize',40)
colormap pink
cbar = colorbar;
ylabel(cbar,'$\int A_y$ d$y$ (T mm$^2$)','FontSize',40,'Interpreter','LaTeX')
title(['$B_0=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(gcf,['EllipsoidVectorPotential_2D_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['EllipsoidVectorPotential_2D_' num2str(round(B0)) 'T.png'],'-dpng','-painters')

%Calculate deflections and final position for the point source
[dIYx_Div,dIYz_Div] = gradient(IY_Div,x0,z0);
[dIZx_Div,dIZz_Div] = gradient(IZ_Div,x0,z0);

dIr = sqrt((dIYx_Div+dIZx_Div).^2+(dIYz_Div+dIZz_Div).^2);

figure(6)
clf
imagesc(x0*1e3,z0*1e3,abs(dIr)*1e3)
axis equal tight
set(gca,'ydir','normal','FontSize',40);%'PlotBoxAspectRatio',[1 1 1]*0.9)
caxis([0,1e3*max(abs(dIr(:)))])
xticks(-4:1:4)
yticks(-4:1:4)
xlabel('$r_0$ (mm)','Interpreter','LaTeX','FontSize',40)
ylabel('$r_I$ (cm)','Interpreter','LaTeX','FontSize',40)
colormap pink
cbar = colorbar;
ylabel(cbar,'$\int B_r$ d$y$ (T mm)','FontSize',40,'Interpreter','LaTeX')
title(['$B_0=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(gcf,['EllipsoidMagneticField_2D_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['EllipsoidMagneticField_2D_' num2str(round(B0)) 'T.png'],'-dpng','-painters')

%Turn the deflection map into a functional form
% fdyx_Div = griddedInterpolant(alpha0nd,theta0nd,dIYx_Div);
fdyx_Div = griddedInterpolant(x0nd,z0nd,dIYx_Div');
fdyz_Div = griddedInterpolant(x0nd,z0nd,dIYz_Div');
fdzx_Div = griddedInterpolant(x0nd,z0nd,dIZx_Div');
fdzz_Div = griddedInterpolant(x0nd,z0nd,dIZz_Div');

figure(10)
clf
imagesc(fdyx_Div(x0nd,z0nd)+fdyz_Div(x0nd,z0nd))

%% Initialize protons and get deflections from the field function

%FOR A RANDOMIZED ARRAY OF PROTONS IN 2D
% nProtons = 2e6;

%Uniform distribution
IIxnd = -dI/2+2*dI/2.*rand(nProtons,1);
IIznd = -dI/2+2*dI/2.*rand(nProtons,1);

%Create absolute r position
IIrnd = sqrt(IIxnd.^2+IIznd.^2);

% II3 = hist3([IIxnd,IIynd]/Z,{IIx-0.0,IIy});
II3 = hist3([IIxnd(abs(IIxnd)<=dI/2 & abs(IIznd)<=dI/2),...
    IIznd(abs(IIxnd)<=dI/2 & abs(IIznd)<=dI/2)],{x2D',x2D'});
figure(3)
clf
imagesc(II3)
axis equal tight

%Convert to angles
alpha0nd = atan2(IIrnd,LI);
theta0nd = atan2(IIznd,IIxnd);

%Proton positions at object center
px0 = LO*tan(alpha0nd).*cos(theta0nd);
pz0 = LO*tan(alpha0nd).*sin(theta0nd);
pr0 = sqrt(px0.^2+pz0.^2);

%Proton positions at image plane
pxI0 = LO*tan(alpha0nd).*cos(theta0nd)+LOI*tan(alpha0nd).*cos(theta0nd);
pzI0 = LO*tan(alpha0nd).*sin(theta0nd)+LOI*tan(alpha0nd).*sin(theta0nd);
prI0 = sqrt(pxI0.^2+pzI0.^2);

prIO = sqrt((LO*tan(alpha0nd).*cos(theta0nd)).^2+(LO*tan(alpha0nd).*sin(theta0nd)).^2);

% pxI0 = pxI0(abs(pxI0)<dI/2 & abs(pyI0)<dI/2);

%Make image via histogram
pI0 = hist3([pxI0(abs(pxI0)<=dI/2 & abs(pzI0)<=dI/2),...
    pzI0(abs(pxI0)<=dI/2 & abs(pzI0)<=dI/2)],{x2D',x2D'});

%Display image
figure(1)
clf
imagesc(pI0)
set(gca,'ydir','normal')
axis equal tight

%Proton probe energies
% Ep = [3,14.7]*1e6; %Proton energy array (eV)

%Relativistic gamma factor (used for determining velocity)
gamma = Ep/Mp+1;    %Gamma array

%Resulting proton velocities
v0 = c*sqrt(1-gamma.^-2);   %Velocity array [low, high]

%Calculate ballistic proton time of flight
t1 = (LO-dL/2)./v0;     %time to start of interaction region (s)
t2 = (LO+dL/2)./v0;     %time to far end of interaction region
t3 = (LI)./v0;          %time to image plane from initial location
t23 = t3-t2;            %time from end of interaction region to image plane
M = (t3)/(t2);  %magnification

% alphacutoff = 5*pi/120;

%Calculating proton positions
rO = v0'.*sin(alpha0nd').*t1';
r1 = v0'.*sin(alpha0nd').*t2';    %ballistic positions at L1
r12 = v0'.*sin(alpha0nd').*(t3'); %ballistic positions at LI
vr1 = v0'.*sin(alpha0nd');      %initial proton velocities along r (perpendicular)
vx1 = v0'.*sin(alpha0nd').*cos(theta0nd');
vz1 = v0'.*sin(alpha0nd').*sin(theta0nd');
vy1 = v0'.*cos(alpha0nd');      %initial proton velocities along y (parallel)

%
Dvx_Div = e/mp*(fdyx_Div(px0,pz0)'+1*(vz1./vy1).*fdzx_Div(px0,pz0)');
Dvy_Div = e/mp*(-(vx1./vy1).*fdyx_Div(px0,pz0)'-(vz1./vy1).*fdyz_Div(px0,pz0)');
Dvz_Div = e/mp*(fdyz_Div(px0,pz0)'-1*(vz1./vy1).*fdzz_Div(px0,pz0)');

vx2 = vx1+Dvx_Div;
vz2 = vz1+Dvz_Div;
vy2 = sqrt(v0'.^2.*(vx1./vx1)-vx2.^2-vz2.^2);
vr2 = sqrt(vx2.^2+vz2.^2);
alpha2 = atan2(vr2,vy2)-alpha0nd';

px2 = ones(2,1).*px0'+vx2*LOI./vy2;
pz2 = ones(2,1).*pz0'+vz2*LOI./vy2;
pr2 = sqrt(px2.^2+pz2.^2);

pxI1 = px0'+vx2(1,:)*LOI./vy2(1,:);
pzI1 = pz0'+vz2(1,:)*LOI./vy2(1,:);
prI1 = sqrt(pxI1.^2+pzI1.^2);

pxI1 = pxI1';
pzI1 = pzI1';
prI1 = prI1';

pxI2 = px0'+vx2(2,:)*LOI./vy2(2,:);
pzI2 = pz0'+vz2(2,:)*LOI./vy2(2,:);
prI2 = sqrt(pxI2.^2+pzI2.^2);

pxI2 = pxI2';
pzI2 = pzI2';
prI2 = prI2';

%We also want to find the radial deflection in a single line
r0line = v0'.*sin(alpha0).*t1';
r1line = v0'.*sin(alpha0).*t2';    %ballistic positions at L1
r12line = v0'.*sin(alpha0).*(t3'); %ballistic positions at LI
vr1line = v0'.*sin(alpha0);      %initial proton velocities along r (perpendicular)
vx1line = v0'.*sin(alpha0).*cos(0);
vz1line = v0'.*sin(alpha0).*sin(0);
vy1line = v0'.*cos(alpha0);      %initial proton velocities along y (parallel)

%Proton positions at object center
px0line = ones(2,1)*LO*tan(alpha0).*cos(0);
pz0line = ones(2,1)*LO*tan(alpha0).*sin(0);
pr0line = sqrt(px0line.^2+pz0line.^2);

%Proton positions at image plane
pxI0line = ones(2,1)*(LO*tan(alpha0).*cos(0)+LOI*tan(alpha0).*cos(0));
pzI0line = ones(2,1)*(LO*tan(alpha0).*sin(0)+LOI*tan(alpha0).*sin(0));
prI0line = sqrt(pxI0line.^2+pzI0line.^2);

prIOline = ones(2,1)*sqrt((LO*tan(alpha0).*cos(0)).^2+(LO*tan(alpha0).*sin(0)).^2);


Dvx_line = e/mp*(fdyx_Div(px0line,pz0line)+1*(vz1line./vy1line).*fdzx_Div(px0line,pz0line));
Dvy_line = e/mp*(-(vx1line./vy1line).*fdyx_Div(px0line,pz0line)-(vz1line./vy1line).*fdyz_Div(px0line,pz0line));
Dvz_line = e/mp*(fdyz_Div(px0line,pz0line)-1*(vz1line./vy1line).*fdzz_Div(px0line,pz0line));

vx2line = vx1line+Dvx_line;
vz2line = vz1line+Dvz_line;
vy2line = sqrt(v0'.^2.*(vx1line./vx1line)-vx2line.^2-vz2line.^2);
vr2line = sqrt(vx2line.^2+vz2line.^2);
alpha2line = atan2(vx2line,vy2line)-ones(2,1)*alpha0;

px2line = px0line+vx2line*LOI./vy2line;
pz2line = pz0line+vz2line*LOI./vy2line;
pr2line = sqrt(px2line.^2+pz2line.^2);

%% Make images

%Make image via histcounts2 (faster?)
pI1 = histcounts2(pxI1(abs(pxI1)<=dI/2 & abs(pzI1)<=dI/2),...
    pzI1(abs(pxI1)<=dI/2 & abs(pzI1)<=dI/2),xedges',xedges');
pI1 = imgaussfilt(pI1,gaussParam);
mpI1 = max(pI1(:));
% pI1 = medfilt2(pI1,[1,1]*3);
pI1 = pI1/mpI1;
if noiseParam~=0
pI1(pI1>=0.01) = imnoise(pI1(pI1>=0.01),'gaussian',0,noiseParam);
end

pI2 = histcounts2(pxI2(abs(pxI2)<=dI/2 & abs(pzI2)<=dI/2),...
    pzI2(abs(pxI2)<=dI/2 & abs(pzI2)<=dI/2),xedges',xedges');
pI2 = imgaussfilt(pI2,gaussParam);
% pI2 = medfilt2(pI2,[1,1]*3);
pI2 = pI2/mpI1;
if noiseParam~=0
pI2(pI2>=0.02) = imnoise(pI2(pI2>=0.02),'gaussian',0,noiseParam);
end

% pI1 = Radiograph_Weight(pxI1(:),pzI1(:),dI/2,dI/IRes,IRes);

%Plot proton positions as a function of initial positions

%Plots scaled for images?
figure(1010)
clf
% scatter(r1,r2/M)
plot(px0line(1,:)*1e3,pr2line(1,:)*1e2,'b','linewidth',2)
hold on
plot(px0line(2,:)*1e3,pr2line(2,:)*1e2,'r','linewidth',2)
% scatter(r1(1,:)*1e3,abs(r1(1,:))*M*1e2,'--','color',[1 1 1]*0.5,'linewidth',2)
axis([-0,4,0,ceil(max(pr2line(1,:)*1e2))])
set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
xticks(-4:1:4)
yticks(0:1:10)
% yticklabels({[],0.8,[],1.0,[],1.2,[]})
% ytickformat('%.1f')
xlabel('$r_0$ (mm)','Interpreter','LaTeX','FontSize',40)
ylabel('$r_I$ (cm)','Interpreter','LaTeX','FontSize',40)
title(['$B_0=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
LEG = legend({['$E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
    ['$E_p=$' num2str(Ep(2)/1e6) ' MeV'],'Ballistic trajectory'});
set(LEG,'FontSize',40,'Interpreter','LaTeX','location','southeast')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(gcf,['EllipsoidTrajectory_2D_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['EllipsoidTrajectory_2D_' num2str(round(B0)) 'T.png'],'-dpng','-painters')

figure(1011)
clf
% scatter(r1,r2/M)
% scatter(r1(1,1:1000:end)*1e3,alpha2(1,1:1000:end),'b.','linewidth',2)
plot(px0line(1,:)*1e3,alpha2line(1,:),'b','linewidth',2)
hold on
% scatter(r1(2,1:1000:end)*1e3,alpha2(2,1:1000:end),'r.','linewidth',2)
plot(px0line(2,:)*1e3,alpha2line(2,:),'r','linewidth',2)
% plot(r1(1,:)*1e3,r1(1,:)*M*1e2,'--','color',[1 1 1]*0.5,'linewidth',2)
axis([0,4,0,ceil(max(alpha2line(1,:)*10))/10])
set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
xticks(0:1:4)
yticks(0:0.05:0.5)
% yticklabels({[],0.8,[],1.0,[],1.2,[]})
% ytickformat('%.1f')
xlabel('$r_0$ (mm)','Interpreter','LaTeX','FontSize',40)
ylabel('Deflection angle $\alpha$','Interpreter','LaTeX','FontSize',40)
title(['$B_0=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
LEG = legend({['$E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
    ['$E_p=$' num2str(Ep(2)/1e6) ' MeV'],'Ballistic trajectory'});
set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
% print(gcf,['EllipsoidDeflection_2D_Truncated_Lineout_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['EllipsoidDeflection_2D_Lineout_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['EllipsoidDeflection_2D_Lineout_' num2str(round(B0)) 'T.png'],'-dpng','-painters')

%Make images
figure(1021)
clf
imagesc(x2D*100,x2D*100,(1-(pI1/max(pI1(:)))))
set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
xticks(-5:1:5)
yticks(-5:1:5)
colormap(gray(1000))
caxis([0,1])
cbar = colorbar;
% yticklabels({[],0.8,[],1.0,[],1.2,[]})
% ytickformat('%.1f')
xlabel('$x_I$ (cm)','Interpreter','LaTeX','FontSize',40)
ylabel('$z_I$ (cm)','Interpreter','LaTeX','FontSize',40)
ylabel(cbar,'Relative proton incidence','Interpreter','LaTeX','FontSize',40)
title(['$B_0=' num2str(B0) '$ T, $E_p=$' num2str(Ep(1)/1e6) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(1021,['EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(Ep(1)/1e6) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
print(1021,['EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(Ep(1)/1e6) '_MeV.png'],'-dpng','-painters')

%Make images
figure(1022)
clf
imagesc(x2D*100,x2D*100,(1-(pI2/max(pI1(:)))))
set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
xticks(-5:1:5)
yticks(-5:1:5)
colormap(gray(1000))
caxis([0,1])
cbar = colorbar;
% yticklabels({[],0.8,[],1.0,[],1.2,[]})
% ytickformat('%.1f')
xlabel('$x_I$ (cm)','Interpreter','LaTeX','FontSize',40)
ylabel('$z_I$ (cm)','Interpreter','LaTeX','FontSize',40)
ylabel(cbar,'Relative proton incidence','Interpreter','LaTeX','FontSize',40)
title(['$B_0=' num2str(B0) '$ T, $E_p=$' num2str(Ep(2)/1e6) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(1022,['EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
print(1022,['EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.png'],'-dpng','-painters')

clearvars Data1 Data2

%Save data file
Data1.x = xmesh*1e2;
Data1.z = zmesh*1e2;
Data1.I = pI1;
Data2.x = xmesh*1e2;
Data2.z = zmesh*1e2;
Data2.I = pI1;
save(['Data1_EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(Ep(1)/1e6) '_MeV'],'Data1')
save(['Data2_EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV'],'Data2')
% writematrix(Data1,['Data1_EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])
% writematrix(Data2,['Data2_EllipsoidImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])



%% Create a lineout from the images

xcenter = IRes/2;
ycenter = IRes/2;

rres = 200;
r = linspace(0,IRes/2,rres);
rI = r*dI/IRes;
thetares = 80;
theta = linspace((1/2-1/2)*pi,(1/2+1/2)*pi,thetares);
xr = zeros(thetares,rres);
yr = xr;
IIr = xr;
IIr1 = xr;
IIr2 = xr;
for iii = 1:thetares
    xr(iii,:) = xcenter+r*cos(theta(iii));
    yr(iii,:) = ycenter+r*sin(theta(iii));
    IIr1(iii,:) = interp2(abs(pI1),xr(iii,:),yr(iii,:),'linear');
    IIr2(iii,:) = interp2(abs(pI2),xr(iii,:),yr(iii,:),'linear');
end
% [XInd,YInd]  = ndgrid(XI,XI);
pI1f = griddedInterpolant(xInd,yInd,abs(pI1),'linear','nearest');

IIr1(isnan(IIr1))=0;
IIr1sum = median(abs(IIr1));    %median prevents drop in intensity at last pixel
% IIr1sum = IIr1(1,:);
% IIr1sum = IIr1sum/m1;

IIr2(isnan(IIr2))=0;
IIr2sum = median(abs(IIr2));
% IIr2sum = IIr2(10,:);   %different index than I1 so noise is not correlated
% IIr2sum = IIr2sum/m1;

m1 = max(max(IIr1sum),max(IIr2sum));

figure(165)
clf
plot(rI*100,IIr1sum/m1,'b-','linewidth',2)
hold on
plot(rI*100,IIr2sum/m1,'r-','linewidth',2)
axis([0,max(x2D)*100,0,1.2])
    set(gca,'FontSize',40)%,'PlotBoxAspectRatio',[1 1 1])
    xticks(0:1:5.6)
    yticks(0:0.2:1.2)
    xlabel('Image position (cm)','fontsize',40,'interpreter','latex')
    ylabel('Proton fluence','fontsize',40,'interpreter','latex')
    LEG = legend({['$E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
        ['$E_p=$' num2str(Ep(2)/1e6) ' MeV'],...
        ['Reconstructed, $E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
        ['Reconstructed, $E_p=$' num2str(Ep(2)/1e6) ' MeV']});
    set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northoutside')
title(['$B_0=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(165,['EllipsoidImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
print(165,['EllipsoidImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.png'],'-dpng')

clearvars Data1 Data2

%Save the images as csv to be read in later
Data1 = [rI*1e2; IIr1sum];
Data2 = [rI*1e2; IIr2sum];
writematrix(Data1,['Data1_EllipsoidImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(Ep(1)/1e6) '_MeV.csv'])
writematrix(Data2,['Data2_EllipsoidImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])

%Save the radial deflection function
Deflection1 = [px0line(1,:); alpha2line(1,:)];
Deflection2 = [px0line(2,:); alpha2line(2,:)];
writematrix(Deflection1,['Deflection1_EllipsoidImage_2D_Lineout_'...
    num2str(round(B0)) 'T_' num2str(Ep(1)/1e6) '_MeV.csv'])
writematrix(Deflection2,['Deflection2_EllipsoidImage_2D_Lineout_'...
    num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])

end
