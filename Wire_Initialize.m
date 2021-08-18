function Wire_Initialize(Ep,LO,LI,dL,dI,IRes,nProtons,gaussParam,R,B0)
% Initialize proton images from the magnetic field around a
% current-carrying wire

%Constants
q = 1.6e-19;
mp = 1.67e-27;      %proton mass (kg)
Mp = 938.3*10^6;    %proton rest mass (eV)
c = 2.998e8;        %speed of light (m/s)

%System properties
% R = 0.38e-3;        %radius of wire (m)
r0 = 0.0e-3;        %lateral position of proton source (m)
% B0 = 5;            %maximum magnetic field at wire surface (T)
% LO = 1e-2;          %distance from proton source to system center (m)
% LI = 15e-2;         %distance from proton source to image plane (m)
LOI = LI-LO;        %distance from system center to image plane (m)
% dL = 4e-3;          %length of the interaction region (m)
% dI = 10e-2;         %width of image plane (m)

x2D = linspace(-dI/2,dI/2,IRes);%2*length(x2)-1);
dx = x2D(2)-x2D(1);
xedges = [x2D-dx/2 dI/2+dx/2];
[xmesh,ymesh] = meshgrid(x2D,x2D);
[xInd,yInd] = ndgrid(x2D,x2D);
rmesh = sqrt(xmesh.^2+ymesh.^2);

%Angles from proton source to cover
alphamin = atan2(-r0+R,LO); %Starts from where the wire object is
alphamax = atan2(dI/2,LI);  %Goes to the edge of the image
alphares = 10000;           %Number of protons to initialize for the true image
alpha0 = linspace(-alphamax,alphamax,alphares); %Initializing initial angles

%Azimuthal angle
theta0 = linspace(-pi,pi,360*2);

[alpha0mesh,theta0mesh] = meshgrid(alpha0,theta0);

%Image filter parameter
%NOTE: The Gaussian filter command imgaussfilt require the MATLAB image
%processing toolbox. If unavailable, try the filter function
% gaussParam = 2;

%FOR A RANDOMIZED ARRAY OF PROTONS IN 2D
%Random only
IIxnd = -dI/2+2*dI/2.*rand(nProtons,1);
IIynd = -dI/2+2*dI/2.*rand(nProtons,1);

%Create absolute r position
IIrnd = sqrt(IIxnd.^2+IIynd.^2);

II3 = hist3([IIxnd(abs(IIxnd)<=dI/2 & abs(IIynd)<=dI/2),...
    IIynd(abs(IIxnd)<=dI/2 & abs(IIynd)<=dI/2)],{x2D',x2D'});
figure(3)
clf
imagesc(II3)
axis equal tight

%Convert to angles
alpha0nd = atan2(IIrnd,LI);
theta0nd = atan2(IIynd,IIxnd);

%Project to image plane
pxI0 = LO*tan(alpha0nd).*cos(theta0nd)+LOI*tan(alpha0nd).*cos(theta0nd);
pyI0 = LO*tan(alpha0nd).*sin(theta0nd)+LOI*tan(alpha0nd).*sin(theta0nd);
prI0 = sqrt(pxI0.^2+pyI0.^2);

prIO = sqrt((LO*tan(alpha0nd).*cos(theta0nd)).^2+(LO*tan(alpha0nd).*sin(theta0nd)).^2);

%Make image via histogram
pI0 = hist3([pxI0(abs(pxI0)<=dI/2 & abs(pyI0)<=dI/2 & prIO>R),...
    pyI0(abs(pxI0)<=dI/2 & abs(pyI0)<=dI/2 & prIO>R)],{x2D',x2D'});

%Display image
figure(1)
clf
imagesc(pI0)
set(gca,'ydir','normal')
axis equal tight

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

%Calculating proton positions
rO = v0'.*sin(alpha0nd').*t1';
r1 = v0'.*sin(alpha0nd').*t2';    %ballistic positions at L1
r12 = v0'.*sin(alpha0nd').*(t3'); %ballistic positions at LI
vr1 = v0'.*sin(alpha0nd');      %initial proton velocities along r (perpendicular)
vx1 = v0'.*sin(alpha0nd').*cos(theta0nd');
vy1 = v0'.*sin(alpha0nd').*sin(theta0nd');
vz1 = v0'.*cos(alpha0nd');      %initial proton velocities along z (parallel)

%% Integrated magnetic field (assuming 1/r)
BInt = -(B0*R.*cot((alpha0nd')).*(log(v0'.*sin((alpha0nd')).*t2'+r0)-log(v0'.*sin((alpha0nd')).*t1'+r0)))...
    .*(abs(alpha0nd')<=2*alphamax).*(abs(alpha0nd')>alphamin);
vr2 = q/mp*BInt+v0'.*sin(alpha0nd');    %New r velocity
vz2 = sqrt(v0'.^2-vr2.^2);              %Calculate new z velocity
alpha2 = atan2(vr2,vz2)-alpha0nd';      %Calculate new alpha angle (assuming theta unchanged)

%Calculate x and y velocities
vx2 = v0'.*sin(alpha2+alpha0nd').*cos(theta0nd');
vy2 = v0'.*sin(alpha2+alpha0nd').*sin(theta0nd');
% vz2 = v0'.*cos(alpha2+alpha0nd');      %initial proton velocities along z (parallel)

r2 = v0'.*sin(alpha0nd').*t2'+r0+vr2.*t3';  %Deflected radial proton positions

%Energy 1
pxI1 = LO*tan(alpha0nd).*cos(theta0nd)+...
    LOI*tan(alpha0nd+alpha2(1,:)').*cos(theta0nd);
pyI1 = LO*tan(alpha0nd).*sin(theta0nd)+...
    LOI*tan(alpha0nd+alpha2(1,:)').*sin(theta0nd);

%Energy 2
pxI2 = LO*tan(alpha0nd).*cos(theta0nd)+...
    LOI*tan(alpha0nd+alpha2(2,:)').*cos(theta0nd);
pyI2 = LO*tan(alpha0nd).*sin(theta0nd)+...
    LOI*tan(alpha0nd+alpha2(2,:)').*sin(theta0nd);

%Make image via histcounts2 (faster?)
pI1 = histcounts2(pxI1(abs(pxI1)<=(dI/2+dx/2) & abs(pyI1)<=(dI/2+dx/2) & prI0>LI*tan(alphamin)),...
    pyI1(abs(pxI1)<=(dI/2+dx/2) & abs(pyI1)<=(dI/2+dx/2) & prI0>LI*tan(alphamin)),xedges',xedges');
pI1 = imgaussfilt(pI1,gaussParam,'Padding','circular');

pI2 = histcounts2(pxI2(abs(pxI2)<=(dI/2+dx/2) & abs(pyI2)<=(dI/2+dx/2) & prI0>LI*tan(alphamin)),...
    pyI2(abs(pxI2)<=(dI/2+dx/2) & abs(pyI2)<=(dI/2+dx/2) & prI0>LI*tan(alphamin)),xedges',xedges');
pI2 = imgaussfilt(pI2,gaussParam);

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

BIntline = -(B0*R.*cot((alpha0)).*(log(v0'.*sin((alpha0)).*t2'+r0)-log(v0'.*sin((alpha0)).*t1'+r0)))...
    .*(abs(alpha0)<=2*alphamax).*(abs(alpha0)>alphamin);

vr2line = q/mp*BIntline+v0'.*sin(alpha0);    %New r velocity
vz2line = sqrt(v0'.^2-vr2line.^2);              %Calculate new z velocity
alpha2line = atan2(vr2line,vz2line)-ones(2,1)*alpha0;    %Calculate new alpha angle (assuming theta unchanged)

pr2line = v0'.*sin(alpha0).*t2'+r0+vr2line.*t3';  %Deflected radial proton positions

%% Plots

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
title(['$B=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
LEG = legend({['$E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
    ['$E_p=$' num2str(Ep(2)/1e6) ' MeV'],'Ballistic trajectory'});
set(LEG,'FontSize',40,'Interpreter','LaTeX','location','southeast')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(gcf,['WireTrajectory_2D_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['WireTrajectory_2D_' num2str(round(B0)) 'T.png'],'-dpng','-painters')

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
title(['$B=' num2str(B0) '$ T'],'FontSize',40,'Interpreter','LaTeX')
LEG = legend({['$E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
    ['$E_p=$' num2str(Ep(2)/1e6) ' MeV'],'Ballistic trajectory'});
set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
% print(gcf,['FilamentDeflection_2D_Truncated_Lineout_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['WireDeflection_2D_Lineout_' num2str(round(B0)) 'T.pdf'],'-dpdf','-painters','-bestfit')
print(gcf,['WireDeflection_2D_Lineout_' num2str(round(B0)) 'T.png'],'-dpng','-painters')

%Make images
figure(1021)
clf
% imagesc(x2D*100,x2D*100,1-pI1./max(pI1(:)))
imagesc(x2D*100,x2D*100,1-pI1./900) %Normalized to 15 T
set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
colormap(gray(1000))
caxis([0,1])
cbar = colorbar;
xticks(-5:1:5)
yticks(-5:1:5)
% yticklabels({[],0.8,[],1.0,[],1.2,[]})
% ytickformat('%.1f')
xlabel('Position (cm)','Interpreter','LaTeX','FontSize',40)
% ylabel('Image intensity','Interpreter','LaTeX','FontSize',40)
title(['$B=' num2str(B0) '$ T, $E_p=$' num2str(Ep(1)/1e6) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(1021,['WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(1)/1e6)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
print(1021,['WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(1)/1e6)) '_MeV.png'],'-dpng')

%Make images
figure(1022)
clf
% imagesc(x2D*100,x2D*100,1-pI2./max(pI1(:)))
imagesc(x2D*100,x2D*100,1-pI2./900) %Normalized to 15 T
set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
colormap(gray(1000))
caxis([0,1])
cbar = colorbar;
xticks(-5:1:5)
yticks(-5:1:5)
% yticklabels({[],0.8,[],1.0,[],1.2,[]})
% ytickformat('%.1f')
xlabel('Position (cm)','Interpreter','LaTeX','FontSize',40)
% ylabel('Image intensity','Interpreter','LaTeX','FontSize',40)
title(['$B=' num2str(B0) '$ T, $E_p=$' num2str(Ep(2)/1e6) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(1022,['WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
print(1022,['WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.png'],'-dpng')

clearvars Data1 Data2

%Save data file
Data1.x = xmesh*1e2;
Data1.y = ymesh*1e2;
Data1.I = pI1;
Data2.x = xmesh*1e2;
Data2.y = ymesh*1e2;
Data2.I = pI1;
save(['Data1_WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(1)/1e6)) '_MeV'],'Data1')
save(['Data2_WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV'],'Data2')
% writematrix(Data1,['Data1_WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])
% writematrix(Data2,['Data2_WireImage_2D_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])

%% Create a lineout from the images

xcenter = IRes/2;
ycenter = IRes/2;

rres = 200;
r = linspace(0,IRes/2,rres);
rI = r*dI/IRes;
thetares = 200;
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
% IIf = scatteredInterpolant(((xr(:)-xcenter)*pix*1e-4/16)*1.0./Ffr(rrr(:)),...
%     -((yr(:)-ycenter)*pix*1e-4/16)*1.0./Ffr(rrr(:)),IIr(:));

IIr1(isnan(IIr1))=0;
IIr1sum = median(abs(IIr1));
m1 = max(IIr1sum);
% m1 = 857; %For 15 T normalization
% IIr1sum = IIr1sum/m1;

IIr2(isnan(IIr2))=0;
IIr2sum = median(abs(IIr2));
% IIr2sum = IIr2sum/m1;

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
    LEG = legend({['Data, $E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
        ['Data, $E_p=$' num2str(Ep(2)/1e6) ' MeV'],...
        ['Reconstructed, $E_p=$' num2str(Ep(1)/1e6) ' MeV'],...
        ['Reconstructed, $E_p=$' num2str(Ep(2)/1e6) ' MeV']});
    set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northoutside')
title(['$B=' num2str(B0) '$ T, $E_p=$' num2str(Ep(2)/1e6) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
set(gcf, 'PaperSize', [10 10])
set(gcf, 'PaperPosition', [0 0 10 10])
print(165,['WireImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
print(165,['WireImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.png'],'-dpng')

clearvars Data1 Data2

Data1 = [rI*1e2; IIr1sum];
Data2 = [rI*1e2; IIr2sum];
writematrix(Data1,['Data1_WireImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(1)/1e6)) '_MeV.csv'])
writematrix(Data2,['Data2_WireImage_2D_Lineout_' num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])

%Save the radial deflection function
Deflection1 = [px0line(1,:); alpha2line(1,:)];
Deflection2 = [px0line(2,:); alpha2line(2,:)];
writematrix(Deflection1,['Deflection1_WireImage_2D_Lineout_'...
    num2str(round(B0)) 'T_' num2str(round(Ep(1)/1e6)) '_MeV.csv'])
writematrix(Deflection2,['Deflection2_WireImage_2D_Lineout_'...
    num2str(round(B0)) 'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'])

end