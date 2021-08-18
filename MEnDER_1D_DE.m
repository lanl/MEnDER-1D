function MEnDER_1D_DE(filename1,filename2,Ep,R,LO,LI,B0,gaussParam,alphafres,alphafres2,...
    ftype,niterations,weights,resultnum,nruns)
%% Load data files
LASTN = maxNumCompThreads(4);   %Limits number of threads in MATLAB

%Save?
saveon = 1;

dfilename1 = ['Deflection' filename1(5:end)];
dfilename2 = ['Deflection' filename2(5:end)];

%Load the files
Data1 = load(filename1);
Data2 = load(filename2);

if exist(dfilename1,'file') > 0
    Deflection1 = load(dfilename1);
    Deflection2 = load(dfilename2);
else
    Deflection1 = 0*Data1;
    Deflection2 = 0*Data2;
end

%Set filter parameters (if filter using instead of Gaussian imgaussfilt)
dofilt = 0; %Whether to use this filter or not (1 = yes)
a = 1;
windowSize = 5;
sig = 1;
wdist = -((windowSize-1)/2):1:((windowSize-1)/2);
b = exp(-(wdist.^2)/sig.^2);

%For the mean filter
% b = (1/windowSize)*ones(1,windowSize);
% b = [0.25,0.5,0.25];

%Redefine the positions and intensities
x1 = Data1(1,:);
I1 = Data1(2,:);
%If pre-processing the image
% I1 = filter(b,a,I1);  %Alternate filter
% I1 = imgaussfilt(I1,1);

x2 = Data2(1,:);
I2 = Data2(2,:);
% I2 = filter(b,a,I2);
% I2 = imgaussfilt(I2,1);

%Find maximum intensity of the input images
m1 = max(max(I1),max(I2));
I1 = I1/m1;
I2 = I2/m1;

% I1(I1>=0.02) = imnoise(I1(I1>=0.02),'gaussian',0,0.01);
% I2(I2>=0.02) = imnoise(I2(I2>=0.02),'gaussian',0,0.01);

%Functional (interpolant) form of the images
fI1 = griddedInterpolant(x1,I1,'pchip','nearest');
fI2 = griddedInterpolant(x2,I2,'pchip','nearest');

dx = x1(2)-x1(1);
xmax = max(x2)*5/5; %To set maximum position in reconstruction
xnew = x1(1):dx:xmax;

xedges = [xnew-dx/2 xmax+dx/2];

%Create extrapolated images;
I1 = fI1(xnew);
I2 = fI2(xnew);

%Define physical parameters
E1 = Ep(1)/1e6;             %proton energy 1 (MeV)
E2 = Ep(2)/1e6;          %proton energy 2 (MeV)
R = R*100;        %radius of wire (cm)
r0 = 0.0e-1;        %lateral position of proton source (cm)
LO = LO*100;          %distance from proton source to system center (cm)
LI = LI*100;         %distance from proton source to image plane (cm)
LOI = LI-LO;        %distance from system center to image plane (cm)
% dL = 4e-1;          %length of the interaction region (cm)
dI = 2*max(xnew);         %width of image plane (cm)
M = LI/LO;

%Angles from proton source to cover
alphamin = atan2(R,LO); %Starts from where the wire object is
alphamax = atan2(dI/2,LI);  %Goes to the edge of the image
alphares = 5000;           %Number of protons to initialize for the true image

R0 = linspace((R*M),(dI/2),alphares);
alpha0 = atan2(R0,LI);

%Project to image plane
pxI0 = LO*tan(alpha0)+LOI*tan(alpha0);

Original.image1 = I1;
Original.image2 = I2;
Original.E1 = E1;
Original.E2 = E2;
Original.R = R;
Original.LO = LO;
Original.LOI = LOI;
Original.LI = LI;
Original.alpha0 = alpha0;
Original.xI0 = pxI0;

%% First pass
%Define locations of function anchor points
alphaA = linspace(alphamin,alphamax*1,alphafres); %anchor locations

%x (r) location of anchor points in image plane
xAnd = LI*tan(alphaA);

%Define the number of candidates in the population
if alphafres <= 30
    candidates = 100;%10*alphafres;
else
    candidates = 4*alphafres;
    %     candidates = 100;
end

%Set the bounds of the deflection funciton
vmin = -0.01*0;
vmax = 0.5;

wi = weights;

%Set mean filter parameters
% windowSize = 3;
% % b = (1/windowSize)*ones(1,windowSize);
% b = [0.25,0.5,0.25];
% a = 1;

stepPlots = 1; %Whether to plot the intermediate steps or not

clearvars Reconstruction

%% Loop for number of runs

for run = 1:nruns
    
    %Set error weights
    wi = weights;
    
    %Initial weights
    w1 = wi(1);
    w2 = wi(2);
    w3 = wi(3);

    R1 = NewPopulation(I1,I2,E1,E2,LO,LOI,alphares,alpha0,pxI0,dI,xedges,...
        alphafres,alphaA,candidates,vmax,wi,ftype,gaussParam,a,b,dofilt);

    [Emin,loc] = min(R1.candidateError);
    minError = Emin;
    Echange = 0;
    tests = 0;
    
    figure(2)
    clf
    plot(xnew,I1,'b-','linewidth',2)
    hold on
    plot(xnew,I2,'r-','linewidth',2)
    plot(xnew,R1.pI1{loc},'b--','linewidth',2)
    plot(xnew,R1.pI2{loc},'r--','linewidth',2)
    
    figure(3)
    clf
    plot(alpha0,R1.wfdeflection{loc}(alpha0))
    hold on
    plot(alpha0,R1.wfdeflection{loc}(alpha0)*sqrt(E1/E2))
    plot(alpha0,LI*alphamin./(LI.*tan(alpha0))*0.15)
    
    figure(50)
    clf
    scatter((R1.Err1+R1.Err2),R1.Err3,'.')
    
    figure(51)
    clf
    scatter3((R1.Err1+R1.Err2),R1.Err3,R1.candidateError,'filled')
    
    %% Start evolution algorithm
    
    CR0 = 0.6;
    F0 = 0.5;
    step1 = 0;
    step2 = 0;
    
    %Start loop
    nochange = 0;
    downshift = 0;
    
    %Error limits (stopping criteria)
%     Elim = [0.1,0.075,0.05]*4; %MRE
    Elim = [0.5,0.25,0.1]/1; %MRE
    
    while downshift < 3 && tests < niterations && ...
            ((weights(1)~=0)*R1.Err1(loc)+(weights(2)~=0)*R1.Err2(loc)) > Elim(3)*1.0
        %     for iter = 1:niterations
        tests = tests+1;
        %         if tests == 1
        if downshift == 0
            CR = CR0;
            F = F0;
            %         elseif tests == 1e6
            %             CR = CR0/2;
            %             F = F0/2;
        end
        
        %Use this to reduce the mutation rate and thus speed up final
        %evolution steps
        if ((weights(1)~=0)*R1.Err1(loc)+(weights(2)~=0)*R1.Err2(loc)) <= Elim(1) && step1 == 0
            CR = CR0/2;
            downshift = 1;
            step1 = 1;
        elseif ((weights(1)~=0)*R1.Err1(loc)+(weights(2)~=0)*R1.Err2(loc)) <= Elim(2) && step2 == 0
            CR = CR0/2;
            downshift = 2;
            step2 = 1;
        end
        
        %Only let it go for so long before restarting if error is high
        if ((weights(1)~=0)*R1.Err1(loc)+(weights(2)~=0)*R1.Err2(loc)) > Elim(1) && tests >= niterations/2
%             downshift = 3;
            R1 = NewPopulation(I1,I2,E1,E2,LO,LOI,alphares,alpha0,pxI0,dI,xedges,...
                alphafres,alphaA,candidates,vmax,wi,ftype,gaussParam,a,b,dofilt);
            
            [Emin,loc] = min(R1.candidateError);
            minError = Emin;
            Echange = 0;
            tests = 0;
            
            w1 = wi(1);
            w2 = wi(2);
            w3 = wi(3);

            CR0 = 0.6;
            F0 = 0.5;
            step1 = 0;
            step2 = 0;
            nochange = 0;
            downshift = 0;
        end
        
        if ((weights(1)~=0)*R1.Err1(loc)+(weights(2)~=0)*R1.Err2(loc)) < Elim(2) && tests >= niterations/2
            downshift = 3;
        end
        
        %     F = 1-min(0.8*tests/1.5e6,0.8);
        if rem(tests,1e4)==0
            disp(['Iteration ' num2str(tests) ', minimum Error = ' num2str(Emin)])
            %         CR = CR0-0.0125*floor(tests/1e5);
            %         CR = max(CR,0.5);
            %         F = F0-0.05*floor(tests/1e5);
            %         F = max(F,0.2);
            
            %makePlots
            if stepPlots == 1
                figure(2)
                clf
                plot(xnew,I1,'b:','linewidth',2)
                hold on
                plot(xnew,I2,'r:','linewidth',2)
                plot(xnew,R1.pI1{loc},'b-','linewidth',2)
                plot(xnew,R1.pI2{loc},'r-','linewidth',2)
                axis([0,max(x1),0,1.2])
                set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
                xticks(0:1:5.6)
                yticks(0:0.1:1.2)
                xlabel('Image position (cm)','fontsize',40,'interpreter','latex')
                ylabel('Proton fluence','fontsize',40,'interpreter','latex')
                % LEG = legend({['Data, $E_p=$' num2str(E1) ' MeV'],...
                %     ['Data, $E_p=$' num2str(E2) ' MeV'],...
                %     ['Reconstructed, $E_p=$' num2str(E1) ' MeV'],...
                %     ['Reconstructed, $E_p=$' num2str(E2) ' MeV']});
                % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
                set(gcf, 'PaperSize', [10 10])
                set(gcf, 'PaperPosition', [0 0 10 10])
                
                figure(3)
                clf
%                 plot(LO*tan(alpha0)*10,LI*alphamin./(LI.*tan(alpha0))*0.106,'k-','linewidth',2)
                plot(Deflection1(1,:)*1e3,Deflection1(2,:),'b:','linewidth',2)
                hold on
                plot(Deflection2(1,:)*1e3,Deflection2(2,:),'r:','linewidth',2)
                plot(LO*tan(alpha0)*10,R1.wfdeflection{loc}(alpha0),'b-','linewidth',2)
                plot(LO*tan(alpha0)*10,R1.wfdeflection{loc}(alpha0)*sqrt(E1/E2),'r-','linewidth',2)
                scatter(LO*tan(alphaA)*10,R1.defmatrix{loc},'b+')
                scatter(LO*tan(alphaA)*10,R1.defmatrix{loc}*sqrt(E1/E2),'r+')
                scatter(LO*tan(alphaA)*10,mode(cell2mat(R1.defmatrix'),2),'b','filled')
                axis([0,max(x1),vmin,vmax])
                set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
                % xticks(0:1:4)
                % yticks(0:1:5)
                xlabel('Object position (mm)','fontsize',40,'interpreter','latex')
                ylabel('Deflection angle','fontsize',40,'interpreter','latex')
                % LEG = legend({['Reference $f \propto 1/r$'],...
                %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV'],...
                %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV']});
                % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
                set(gcf, 'PaperSize', [10 10])
                set(gcf, 'PaperPosition', [0 0 10 10])
                
                figure(4)
                clf
                plot(Echange,log(minError))
                
                figure(41)
                clf
                hist(R1.candidateError,100)
                
                figure(50)
                clf
                scatter((R1.Err1+R1.Err2),R1.Err3,'.')
                
            end
        end
        %Choose random candidate and test agents
        cx = randperm(candidates,4);
        
        %Choose random anchor index to change regardless
        Ri = randi(alphafres,1);
        
        %Get random number for each anchor index
        r = rand(alphafres,1);
        
        %DE/rand/1
        if R1.candidateError(cx(1)) < R1.candidateError(cx(2)) || cx(2)==loc
            %         if Err3(cx(1)) < Err3(cx(2))
            newdef = (r>=CR).*R1.defmatrix{cx(1)}+(r<CR).*R1.defmatrix{cx(1)}+...
                (r<CR).*F*(1+rand).*(R1.defmatrix{cx(2)}-R1.defmatrix{cx(3)});
            newdef(Ri) = R1.defmatrix{cx(1)}(Ri)+F*(1+rand).*(R1.defmatrix{cx(2)}(Ri)-...
                R1.defmatrix{cx(3)}(Ri))+0.0*randn*vmax;
        else
            newdef = (r>=CR).*R1.defmatrix{cx(2)}+(r<CR).*R1.defmatrix{cx(2)}+...
                (r<CR).*F*(1+rand).*(R1.defmatrix{cx(1)}-R1.defmatrix{cx(3)});
            newdef(Ri) = R1.defmatrix{cx(1)}(Ri)+F*(1+rand).*(R1.defmatrix{cx(2)}(Ri)-...
                R1.defmatrix{cx(3)}(Ri))+0.0*randn*vmax;
        end
        
        newdef(newdef<vmin) = vmin+abs(vmin-newdef(newdef<vmin));
        newdef(newdef>vmax) = rand*0.1*vmax;
        %To enfore 0 deflection at r = 0 (by symmetry?)
        if R == 0
            newdef(1) = 0;
        end
        
        if rem(tests,1e3)==0
%             newdef = R1.defmatrix{loc};
%             newdef = median(cell2mat(R1.defmatrix'),2);
            newdef = harmmean(cell2mat(R1.defmatrix'),2);
            
            %Some weird median filter
            newmesh = meshgrid(newdef,newdef);
            newmesh = medfilt2(newmesh,[1,1]*round(alphafres/10*(1+rand)))';
            newdef = newmesh(:,round(alphafres/2));
            
%             newdef = median(cell2mat(R1.defmatrix'),2);
        end
        
        if rem(tests,1e4)==0
            newdef = mode(cell2mat(R1.defmatrix'),2);
        end
        
        %Get error of new agent
        defAnchors = newdef;
        wfdeflectionTest = griddedInterpolant(alphaA,defAnchors,ftype); %functional weights
        alpha1Test = wfdeflectionTest(alpha0);
        alpha2Test = alpha1Test*sqrt(E1/E2);
        
        pxI1Test = LO*tan(alpha0)+LOI*tan(alpha1Test+alpha0);
        % pI1Test = hist(pxI1Test(pxI1Test<max(x1)),x1);
        %         pI1Test = histcounts(pxI1Test(pxI1Test<max(xnew)),xedges);
        pI1Test = radialweightedHist(pxI0,pxI1Test,xedges);
        % pI1Test = histcounts(pxI1Test(abs(pxI1Test)<=max(x2) pxI0>LI*tan(alphamin)),xedges');
        
        pxI2Test = LO*tan(alpha0)+LOI*tan(alpha2Test+alpha0);
        % pI2Test = hist(pxI2Test(pxI2Test<max(x2)),x2);
        %         pI2Test = histcounts(pxI2Test(pxI2Test<max(xnew)),xedges);
        pI2Test = radialweightedHist(pxI0,pxI2Test,xedges);
        
        if gaussParam > 0
            pI1Test = imgaussfilt(pI1Test,gaussParam);
            pI2Test = imgaussfilt(pI2Test,gaussParam);
        end
        
        %Apply mean filter
        if dofilt == 1
            pI1Test = filter(b,a,pI1Test);
            pI2Test = filter(b,a,pI2Test);
        end
        
        pI2Test = pI2Test*sum(I2)/sum(pI2Test);
        pI1Test = pI1Test*sum(I1)/sum(pI1Test);
        
        %Calculate error heuristics
        
        %Mean relative error of squares
        Err1Test = mean(abs(pI1Test.^2-I1.^2)./(I1+2e-2).^2);
        Err2Test = mean(abs(pI2Test.^2-I2.^2)./(I2+2e-2).^2);

        %add gradient smoothness?
        gga = gradient(alpha1Test,10*LO*tan(alpha0));
        ggda = sqrt((10*LO*tan(alpha0(2:end))-10*LO*tan(alpha0(1:end-1))).^2+(gga(2:end)-gga(1:end-1)).^2);
        Err3Test = (sum(ggda)-(10*LO*tan(alphamax)-10*LO*tan(alphamin)))/...
            (10*LO*tan(alphamax)-10*LO*tan(alphamin));
        
        %Total candidate error
        candidateErrorTest = (w1*Err1Test+w2*Err2Test)*(1+w3*Err3Test);

        if candidateErrorTest <= 1.25*R1.candidateError(cx(1))
            if candidateErrorTest <= 1.0*R1.candidateError(cx(1))
                R1.defmatrix{cx(1)} = newdef;
                R1.wfdeflection{cx(1)} = wfdeflectionTest;

                R1.pI1{cx(1)} = pI1Test;

                R1.pI2{cx(1)} = pI2Test;
                
                R1.Err1(cx(1)) = Err1Test;
                R1.Err2(cx(1)) = Err2Test;
                R1.Err3(cx(1)) = Err3Test;
                R1.candidateError(cx(1)) = candidateErrorTest;
                
            elseif rand<0.1 && cx(1)~=loc
                R1.defmatrix{cx(1)} = newdef;
                R1.wfdeflection{cx(1)} = wfdeflectionTest;

                R1.pI1{cx(1)} = pI1Test;

                R1.pI2{cx(1)} = pI2Test;
                
                R1.Err1(cx(1)) = Err1Test;
                R1.Err2(cx(1)) = Err2Test;
                R1.Err3(cx(1)) = Err3Test;
                R1.candidateError(cx(1)) = candidateErrorTest;
            end
            
        end
        
        [Emin,loc] = min(R1.candidateError);
        if Emin < minError(end)
            minError = [minError Emin];
            Echange = [Echange tests];
            if minError(end-1)/Emin-1 > 0.01/(downshift+1)
                nochange = 0;
            else
                nochange = nochange+1;
            end
        else
            nochange = nochange+1;
        end
        
        if nochange == 0.5e5
            nochange = 0;
            CR = CR/2;
            downshift = downshift+1;
        end
        
    end
    
    [Emin,loc] = min(R1.candidateError);
    % loc = cx(1);
    
    % makePlots
    
    figure(2)
    clf
    plot(xnew,I1,'b:','linewidth',2)
    hold on
    plot(xnew,I2,'r:','linewidth',2)
    plot(xnew,R1.pI1{loc},'b-','linewidth',2)
    plot(xnew,R1.pI2{loc},'r-','linewidth',2)
    axis([0,max(x1),0,1.2])
    set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
    xticks(0:1:5.6)
    yticks(0:0.1:1.2)
    xlabel('Image position (cm)','fontsize',40,'interpreter','latex')
    ylabel('Proton fluence','fontsize',40,'interpreter','latex')
    % LEG = legend({['Data, $E_p=$' num2str(E1) ' MeV'],...
    %     ['Data, $E_p=$' num2str(E2) ' MeV'],...
    %     ['Reconstructed, $E_p=$' num2str(E1) ' MeV'],...
    %     ['Reconstructed, $E_p=$' num2str(E2) ' MeV']});
    % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    
    figure(3)
    clf
%     plot(LO*tan(alpha0)*10,LI*alphamin./(LI.*tan(alpha0))*0.106,'k-','linewidth',2)
    plot(Deflection1(1,:)*1e3,Deflection1(2,:),'b:','linewidth',2)
                hold on
    plot(Deflection2(1,:)*1e3,Deflection2(2,:),'r:','linewidth',2)
    plot(LO*tan(alpha0)*10,R1.wfdeflection{loc}(alpha0),'b-','linewidth',2)
    plot(LO*tan(alpha0)*10,R1.wfdeflection{loc}(alpha0)*sqrt(E1/E2),'r-','linewidth',2)
    axis([0,max(x1),vmin,vmax])
    set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
    % xticks(0:1:4)
    % yticks(0:1:5)
    xlabel('Object position (mm)','fontsize',40,'interpreter','latex')
    ylabel('Deflection angle','fontsize',40,'interpreter','latex')
    % LEG = legend({['Reference $f \propto 1/r$'],...
    %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV'],...
    %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV']});
    % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    
    figure(4)
    clf
    plot(Echange,minError)
    
    figure(41)
    clf
    hist(R1.candidateError,20)
    
    figure(50)
    clf
    scatter((R1.Err1+R1.Err2),R1.Err3,'.')
    
    %% Second pass
    %Now we are going to take the results from the first pass and make a finer
    %mesh of nodes (by a factor of 2). We will initialize the candidates by interpolating values
    %from the existing candidates.
    
    %Define NEW locations of nodes by adding one between each set
    alphaA2 = linspace(alphamin,alphamax*1,alphafres2); %new node locations

    xA2nd = LI*tan(alphaA2);    %x location of node on image
    
    
    %Set error weights
    wi = weights;
    
    %Set mean filter parameters
%     windowSize = 3;
    % b = (1/windowSize)*ones(1,windowSize);
%     b = [0.25,0.5,0.25];
%     a = 1;
    
    stepPlots = 1; %Whether to plot the intermediate steps or not
    
    clearvars Reconstruction2
    
    %% Define the new candidates
    
    %Initial weights
    w1 = wi(1);
    w2 = wi(2);
    w3 = wi(3);
    
    %Initialize the anchor values for all new candidates by interpolation
    R2.defmatrix = cell(candidates,1);
    for can = 1:candidates
        R2.defmatrix{can} = interp1(alphaA',R1.defmatrix{can},alphaA2',ftype); %interpolation
        
        %To add some randonmess
        if can ~= loc
           R2.defmatrix{can} = R2.defmatrix{can}.*(1+0.02*randn(alphafres2,1)); 
        end
    end
    
    %Initialize error array
    R2.candidateError = zeros(candidates,1);%,vres);
    
    %Initialize candidate function array
    R2.wfdeflection = cell(candidates,1);
    
    %Initialize other arrays for each candidate
    R2.alpha1 = zeros(alphares,1);
    R2.alpha2 = zeros(alphares,1);
    R2.pxI1 = zeros(alphares,1);
    R2.pxI2 = zeros(alphares,1);
    
    R2.pI1 = cell(candidates,1);
    R2.pI2 = cell(candidates,1);
    R2.Err1 = zeros(candidates,1);
    R2.Err2 = zeros(candidates,1);
    R2.Err3 = zeros(candidates,1);
    
    %Calculate images for each candidate
    for iii = 1:candidates
        defAnchors = R2.defmatrix{iii};
        R2.wfdeflection{iii} = griddedInterpolant(alphaA2,defAnchors,ftype); %functional weights
        
        R2.alpha1 = R2.wfdeflection{iii}(alpha0);
        R2.alpha2 = R2.alpha1*sqrt(E1/E2);
        
        R2.pxI1 = LO*tan(alpha0)+LOI*tan(R2.alpha1+alpha0);     %x locations of protons
        R2.pI1{iii} = radialweightedHist(pxI0,R2.pxI1,xedges);  %intensity counts
        
        R2.pxI2 = LO*tan(alpha0)+LOI*tan(R2.alpha2+alpha0);
        R2.pI2{iii} = radialweightedHist(pxI0,R2.pxI2,xedges);
        
        if gaussParam > 0
            R2.pI1{iii} = imgaussfilt(R2.pI1{iii},gaussParam);
            R2.pI2{iii} = imgaussfilt(R2.pI2{iii},gaussParam);
        end
        
        R2.pI2{iii} = R2.pI2{iii}*sum(I2)/sum(R2.pI2{iii});
        R2.pI1{iii} = R2.pI1{iii}*sum(I1)/sum(R2.pI1{iii});
        
        %Mean relative error of squares
        R2.Err1(iii) = mean(abs(R2.pI1{iii}.^2-I1.^2)./(I1+1e-2).^2);
        R2.Err2(iii) = mean(abs(R2.pI2{iii}.^2-I2.^2)./(I2+1e-2).^2);

        gga = gradient(R2.alpha1,10*LO*tan(alpha0));
        ggda = sqrt((10*LO*tan(alpha0(2:end))-10*LO*tan(alpha0(1:end-1))).^2+(gga(2:end)-gga(1:end-1)).^2);
        R2.Err3(iii) = (sum(ggda)-(10*LO*tan(alphamax)-10*LO*tan(alphamin)))/...
            (10*LO*tan(alphamax)-10*LO*tan(alphamin));
        
        R2.candidateError(iii) = (w1*R2.Err1(iii)+w2*R2.Err2(iii))*(1+w3*R2.Err3(iii));
    end
    
    [Emin,loc] = min(R2.candidateError);
    minError = [minError,Emin];
    tests1 = tests+1;
    Echange = [Echange,tests1];
    tests = 0;
    
    figure(2)
    clf
    plot(xnew,I1,'b-','linewidth',2)
    hold on
    plot(xnew,I2,'r-','linewidth',2)
    plot(xnew,R2.pI1{loc},'b--','linewidth',2)
    plot(xnew,R2.pI2{loc},'r--','linewidth',2)
    
    figure(3)
    clf
    plot(alpha0,R2.wfdeflection{loc}(alpha0))
    hold on
    plot(alpha0,R2.wfdeflection{loc}(alpha0)*sqrt(E1/E2))
    plot(alpha0,LI*alphamin./(LI.*tan(alpha0))*0.15)
    
    figure(50)
    clf
    scatter((R2.Err1+R2.Err2),R2.Err3,'.')
    
    %% Start evolution algorithm
    
    niterations2 = 1.0*niterations;
    
    CR0 = 0.2;
    F0 = 0.5;
    step1 = 0;
    step2 = 0;
    
    %Start loop
    nochange = 0;
    downshift = 0;

    Elim2 = [0.04,0.03,0.015].*1; %MRE
    
    while downshift < 3 && tests < niterations2 && (R2.Err1(loc)+R2.Err2(loc)) > Elim2(3)
        tests = tests+1;
        if downshift == 0
            CR = CR0;
            F = F0;
            %         elseif tests == 1e6
            %             CR = CR0/2;
            %             F = F0/2;
        end
        
        %Use this to reduce the mutation rate and thus speed up final
        %evolution steps
        if (R2.Err1(loc)+R2.Err2(loc)) <= Elim2(1) && step1 == 0
            CR = CR0/2;
            downshift = 1;
            step1 = 1;
        elseif (R2.Err1(loc)+R2.Err2(loc)) <= Elim2(1) && step2 == 0
            downshift = 2;
            step2 = 1;
        end
        
%         %If the error has reached a suitable stage
%         if (R2.Err1(loc)+R2.Err2(loc)) < 0.05 && tests >= niterations2/2
%             downshift = 3;
%         end
        
        if rem(tests,1e4)==0
            disp(['Iteration ' num2str(tests) ', minimum Error = ' num2str(Emin)])
            
            %makePlots
            if stepPlots == 1
                figure(2)
                clf
                plot(xnew,I1,'b:','linewidth',2)
                hold on
                plot(xnew,I2,'r:','linewidth',2)
                plot(xnew,R2.pI1{loc},'b-','linewidth',2)
                plot(xnew,R2.pI2{loc},'r-','linewidth',2)
                axis([0,max(x1),0,1.2])
                set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
                xticks(0:1:5.6)
                yticks(0:0.1:1.2)
                xlabel('Image position (cm)','fontsize',40,'interpreter','latex')
                ylabel('Proton fluence','fontsize',40,'interpreter','latex')
                % LEG = legend({['Data, $E_p=$' num2str(E1) ' MeV'],...
                %     ['Data, $E_p=$' num2str(E2) ' MeV'],...
                %     ['Reconstructed, $E_p=$' num2str(E1) ' MeV'],...
                %     ['Reconstructed, $E_p=$' num2str(E2) ' MeV']});
                % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
                set(gcf, 'PaperSize', [10 10])
                set(gcf, 'PaperPosition', [0 0 10 10])
                
                figure(3)
                clf
%                 plot(LO*tan(alpha0)*10,LI*alphamin./(LI.*tan(alpha0))*0.106,'k-','linewidth',2)
                plot(Deflection1(1,:)*1e3,Deflection1(2,:),'b:','linewidth',2)
                hold on
                plot(Deflection2(1,:)*1e3,Deflection2(2,:),'r:','linewidth',2)
                plot(LO*tan(alpha0)*10,R2.wfdeflection{loc}(alpha0),'b-','linewidth',2)
                plot(LO*tan(alpha0)*10,R2.wfdeflection{loc}(alpha0)*sqrt(E1/E2),'r-','linewidth',2)
                scatter(LO*tan(alphaA2)*10,R2.defmatrix{loc},'b+')
                scatter(LO*tan(alphaA2)*10,R2.defmatrix{loc}*sqrt(E1/E2),'r+')
                scatter(LO*tan(alphaA2)*10,mean(cell2mat(R2.defmatrix'),2),'filled')
                axis([0,max(x1),vmin,vmax])
                set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
                % xticks(0:1:4)
                % yticks(0:1:5)
                xlabel('Object position (mm)','fontsize',40,'interpreter','latex')
                ylabel('Deflection angle','fontsize',40,'interpreter','latex')
                % LEG = legend({['Reference $f \propto 1/r$'],...
                %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV'],...
                %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV']});
                % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
                set(gcf, 'PaperSize', [10 10])
                set(gcf, 'PaperPosition', [0 0 10 10])
                
                figure(4)
                clf
                plot(Echange,log(minError))
                
                figure(41)
                clf
                hist(R2.candidateError,100)
                
                figure(50)
                clf
                scatter((R2.Err1+R2.Err2),R2.Err3,'.')
            end
        end
        %Choose random candidate and test agents
        cx = randperm(candidates,4);
        
        %Choose random anchor index to change regardless
        Ri = randi(alphafres2,1);
        
        %Get random number for each anchor index
        r = rand(alphafres2,1);
        
        %DE/rand/1
        if R2.candidateError(cx(1)) < R2.candidateError(cx(2)) || cx(2)==loc
            %         if Err3(cx(1)) < Err3(cx(2))
            newdef = (r>=CR).*R2.defmatrix{cx(1)}+(r<CR).*R2.defmatrix{cx(1)}+...
                (r<CR).*F*(1+rand).*(R2.defmatrix{cx(2)}-R2.defmatrix{cx(3)})+...
                (r<CR).*0.00.*randn(alphafres2,1)*vmax;
            newdef(Ri) = R2.defmatrix{cx(1)}(Ri)+F*(1+rand).*...
                (R2.defmatrix{cx(2)}(Ri)-R2.defmatrix{cx(3)}(Ri))...
                +0.01*randn*vmax;
        else
            newdef = (r>=CR).*R2.defmatrix{cx(2)}+(r<CR).*R2.defmatrix{cx(2)}+...
                (r<CR).*F*(1+rand).*(R2.defmatrix{cx(1)}-R2.defmatrix{cx(3)})+...
                (r<CR).*0.00.*randn(alphafres2,1)*vmax;
            newdef(Ri) = R2.defmatrix{cx(1)}(Ri)+F*(1+rand).*...
                (R2.defmatrix{cx(2)}(Ri)-R2.defmatrix{cx(3)}(Ri))...
                +0.01*randn*vmax;
        end

        newdef(newdef<vmin) = vmin+abs(vmin-newdef(newdef<vmin));
        newdef(newdef>vmax) = rand*vmax;
        %To enfore 0 deflection at r = 0 (by symmetry?)
        if R == 0
            newdef(1) = 0;
        end
        
        if rem(tests,1e4)==0
            
            %Try mean
            newdef = mean(cell2mat(R2.defmatrix'),2);

            % A weird median filter
            newmesh = meshgrid(newdef,newdef);
            newmesh = medfilt2(newmesh,[1,1]*round(alphafres2/5*(1+rand)))';
            newdef = newmesh(:,round(alphafres2/2));
        end
        
        %Get error of new agent
        defAnchors = newdef;
        wfdeflectionTest = griddedInterpolant(alphaA2,defAnchors,ftype); %functional weights
        alpha1Test = wfdeflectionTest(alpha0);
        alpha2Test = alpha1Test*sqrt(E1/E2);
        
        pxI1Test = LO*tan(alpha0)+LOI*tan(alpha1Test+alpha0);
        pI1Test = radialweightedHist(pxI0,pxI1Test,xedges);
        
        pxI2Test = LO*tan(alpha0)+LOI*tan(alpha2Test+alpha0);
        pI2Test = radialweightedHist(pxI0,pxI2Test,xedges);
        
        if gaussParam > 0
            pI1Test = imgaussfilt(pI1Test,gaussParam);
            pI2Test = imgaussfilt(pI2Test,gaussParam);
        end
        
        %Apply mean filter
        if dofilt == 1
            pI1Test = filter(b,a,pI1Test);
            pI2Test = filter(b,a,pI2Test);
        end
        
        pI2Test = pI2Test*sum(I2)/sum(pI2Test);
        pI1Test = pI1Test*sum(I1)/sum(pI1Test);
        
        %Mean relative error of squares
        Err1Test = mean(abs(pI1Test.^2-I1.^2)./(I1+2e-2).^2);
        Err2Test = mean(abs(pI2Test.^2-I2.^2)./(I2+2e-2).^2);

        %add gradient smoothness?
        gga = gradient(alpha1Test,10*LO*tan(alpha0));
        ggda = sqrt((10*LO*tan(alpha0(2:end))-10*LO*tan(alpha0(1:end-1))).^2+(gga(2:end)-gga(1:end-1)).^2);
        Err3Test = (sum(ggda)-(10*LO*tan(alphamax)-10*LO*tan(alphamin)))/...
            (10*LO*tan(alphamax)-10*LO*tan(alphamin));
        
        %Total candidate error
        candidateErrorTest = (w1*Err1Test+w2*Err2Test)*(1+w3*Err3Test);

        if candidateErrorTest <= 1.25*R2.candidateError(cx(1))
            if candidateErrorTest <= 1.0*R2.candidateError(cx(1))
                R2.defmatrix{cx(1)} = newdef;
                R2.wfdeflection{cx(1)} = wfdeflectionTest;

                R2.pI1{cx(1)} = pI1Test;

                R2.pI2{cx(1)} = pI2Test;
                
                R2.Err1(cx(1)) = Err1Test;
                R2.Err2(cx(1)) = Err2Test;
                R2.Err3(cx(1)) = Err3Test;
                R2.candidateError(cx(1)) = candidateErrorTest;

            elseif rand<0.1 && cx(1)~=loc
                R2.defmatrix{cx(1)} = newdef;
                R2.wfdeflection{cx(1)} = wfdeflectionTest;

                R2.pI1{cx(1)} = pI1Test;

                R2.pI2{cx(1)} = pI2Test;
                
                R2.Err1(cx(1)) = Err1Test;
                R2.Err2(cx(1)) = Err2Test;
                R2.Err3(cx(1)) = Err3Test;
                R2.candidateError(cx(1)) = candidateErrorTest;
            end
            
        end
        
        [Emin,loc] = min(R2.candidateError);
        if Emin < minError(end)
            minError = [minError Emin];
            Echange = [Echange tests1+tests];
            if minError(end-1)/Emin-1 > 0.01/(downshift+1)
                nochange = 0;
            else
                nochange = nochange+1;
            end
        else
            nochange = nochange+1;
        end
        
        if nochange == 0.5e5
            nochange = 0;
            CR = CR/2;
            downshift = downshift+1;
        end
        
    end
    
    [Emin,loc] = min(R2.candidateError);
    % loc = cx(1);
    
    % makePlots
    
    figure(2)
    clf
    plot(xnew,I1,'b:','linewidth',2)
    hold on
    plot(xnew,I2,'r:','linewidth',2)
    plot(xnew,R2.pI1{loc},'b-','linewidth',2)
    plot(xnew,R2.pI2{loc},'r-','linewidth',2)
    axis([0,max(x1),0,1.2])
    set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
    xticks(0:1:5.6)
    yticks(0:0.1:1.2)
    xlabel('Image position (cm)','fontsize',40,'interpreter','latex')
    ylabel('Proton fluence','fontsize',40,'interpreter','latex')
    % LEG = legend({['Data, $E_p=$' num2str(E1) ' MeV'],...
    %     ['Data, $E_p=$' num2str(E2) ' MeV'],...
    %     ['Reconstructed, $E_p=$' num2str(E1) ' MeV'],...
    %     ['Reconstructed, $E_p=$' num2str(E2) ' MeV']});
    % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    
    figure(3)
    clf
%     plot(LO*tan(alpha0)*10,LI*alphamin./(LI.*tan(alpha0))*0.106,'k-','linewidth',2)
    plot(Deflection1(1,:)*1e3,Deflection1(2,:),'b:','linewidth',2)
    hold on
    plot(Deflection2(1,:)*1e3,Deflection2(2,:),'r:','linewidth',2)
    plot(LO*tan(alpha0)*10,R2.wfdeflection{loc}(alpha0),'b-','linewidth',2)
    plot(LO*tan(alpha0)*10,R2.wfdeflection{loc}(alpha0)*sqrt(E1/E2),'r-','linewidth',2)
    scatter(LO*tan(alphaA2)*10,R2.defmatrix{loc},'b+')
    scatter(LO*tan(alphaA2)*10,R2.defmatrix{loc}*sqrt(E1/E2),'r+')
    axis([0,max(x1),vmin,vmax])
    set(gca,'FontSize',40,'PlotBoxAspectRatio',[1 1 1])
    % xticks(0:1:4)
    % yticks(0:1:5)
    xlabel('Object position (mm)','fontsize',40,'interpreter','latex')
    ylabel('Deflection angle','fontsize',40,'interpreter','latex')
    % LEG = legend({['Reference $f \propto 1/r$'],...
    %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV'],...
    %     ['Reconstructed $f$, $E_p=$' num2str(E2) ' MeV']});
    % set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northeast')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    
    figure(4)
    clf
    plot(Echange,log10(minError))
    
    figure(41)
    clf
    hist(R2.candidateError,20)
    
    figure(50)
    clf
    scatter((R2.Err1+R2.Err2),R2.Err3,'.')
    
    
    Reconstruction(run).anchors = R2.defmatrix{loc};
    Reconstruction(run).alphanodes = alphaA2;
    Reconstruction(run).xnodes = LO*tan(alphaA2);
    Reconstruction(run).ftype = ftype;
    Reconstruction(run).function = R2.wfdeflection{loc};
    Reconstruction(run).alpha0 = alpha0;
    Reconstruction(run).alpha1 = R2.wfdeflection{loc}(alpha0);
    Reconstruction(run).alpha2 = R2.wfdeflection{loc}(alpha0)*sqrt(E1/E2);
    Reconstruction(run).image1 = R2.pI1{loc};
    Reconstruction(run).data1 = I1;
    Reconstruction(run).image2 = R2.pI2{loc};
    Reconstruction(run).data2 = I2;
    Reconstruction(run).error1 = R2.Err1(loc);
    Reconstruction(run).error2 = R2.Err2(loc);
    Reconstruction(run).error3 = R2.Err3(loc);
    Reconstruction(run).error = R2.candidateError(loc);
    Reconstruction(run).echange = minError;
    Reconstruction(run).history = Echange;
    
end

%% Save results
% saveon = 1;
% resultnum = 3;
if saveon == 1
    % Reconstruction.anchors = R2.defmatrix(loc,:);
    % Reconstruction.ftype = ftype;
    % Reconstruction.function = wf{loc};
    % Reconstruction.alpha1 = alpha1{loc};
    % Reconstruction.alpha2 = alpha2{loc};
    % Reconstruction.image1 = pxI1{loc};
    % Reconstruction.data1 = pI1{loc};
    % Reconstruction.image2 = pxI2{loc};
    % Reconstruction.data2 = pI2{loc};
    % Reconstruction.error1 = Err1{loc};
    % Reconstruction.error2 = Err2{loc};
    % Reconstruction.error = candidateError(loc);
    
    x0 = Reconstruction(1).xnodes;
    x02 = linspace(x0(1),x0(end),length(Reconstruction(1).alpha1));
    
    save([filename1(1:end-4) '_Multi-Step_Reconstruction_' num2str(alphafres2) '_nodes_' ...
        num2str(resultnum)],'Reconstruction','Original')
    
    % writematrix(Reconstruction,[filename1(1:end-5) 'Reconstruction.csv'])
    
%     [Emin,loc] = min([Reconstruction.error]);
    [Emin,loc] = min([Reconstruction.error1]+[Reconstruction.error2]);
    
    figure(2)
    clf
    plot(xnew,I1/max(max(I1),max(I2)),'b:','linewidth',2)
    hold on
    plot(xnew,I2/max(max(I1),max(I2)),'r:','linewidth',2)
    plot(xnew,Reconstruction(loc).image1/max(max(I1),max(I2)),'b-','linewidth',2)
    plot(xnew,Reconstruction(loc).image2/max(max(I1),max(I2)),'r-','linewidth',2)
    axis([0,max(x1),0,1.2])
    set(gca,'FontSize',40)%,'PlotBoxAspectRatio',[1 1 1])
    xticks(0:1:5.6)
    yticks(0:0.1:1.2)
    xlabel('Image position (cm)','fontsize',40,'interpreter','latex')
    ylabel('Proton fluence','fontsize',40,'interpreter','latex')
    LEG = legend({['Data, $E_p=$' num2str(E1) ' MeV'],...
        ['Data, $E_p=$' num2str(E2) ' MeV'],...
        ['Reconstruction, $E_p=$' num2str(E1) ' MeV'],...
        ['Reconstruction, $E_p=$' num2str(E2) ' MeV']});
    set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northoutside')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    print(2,[filename1(1:end-4) '_Multi-Step_ReconstructedImage_' num2str(alphafres) '_nodes_' ...
        num2str(resultnum) '.pdf'],'-dpdf','-painters')
    print(2,[filename1(1:end-4) '_Multi-Step_ReconstructedImage_' num2str(alphafres) '_nodes_' ...
        num2str(resultnum) '.png'],'-dpng')
    
    figure(3)
    clf
%     plot(LO*tan(alpha0)*10,LI*alphamin./(LI.*tan(alpha0))*0.2324,'b--','linewidth',2)
%     hold on
%     plot(LO*tan(alpha0)*10,LI*alphamin./(LI.*tan(alpha0))*0.2324*sqrt(E1/E2),'r--','linewidth',2)
    plot(Deflection1(1,:)*1e3,Deflection1(2,:),'b:','linewidth',2)
    hold on
    plot(Deflection2(1,:)*1e3,Deflection2(2,:),'r:','linewidth',2)
    plot(x02*10,Reconstruction(loc).alpha1,'b-','linewidth',2)
    plot(x02*10,Reconstruction(loc).alpha2,'r-','linewidth',2)
    scatter(x0*10,Reconstruction(loc).anchors,'b+')
    scatter(x0*10,Reconstruction(loc).anchors*sqrt(E1/E2),'r+')
%     boxplot(LO*tan(alphaA2)*10,[Reconstruction(:).anchors])
%     axis([0,max(x0*10),-0.01,vmax])%ceil(max(Reconstruction(loc).anchors(:))*10)/10])
    axis([0,max(x0*10),-0.01,ceil(max(1.2*Reconstruction(loc).anchors(:))*10)/10])
    set(gca,'FontSize',40)%,'PlotBoxAspectRatio',[1 1 1])
    xticks(0:1:4)
    % yticks(0:1:5)
    xlabel('Object position (mm)','fontsize',40,'interpreter','latex')
    ylabel('Deflection angle','fontsize',40,'interpreter','latex')
    LEG = legend({['Data, $E_p=$' num2str(E1) ' MeV'],...
        ['Data, $E_p=$' num2str(E2) ' MeV'],...
        ['Reconstruction, $E_p=$' num2str(E1) ' MeV'],...
        ['Reconstruction, $E_p=$' num2str(E2) ' MeV']});
    set(LEG,'FontSize',40,'Interpreter','LaTeX','location','northoutside')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    print(3,[filename1(1:end-4) '_Multi-Step_ReconstructedFunction_' num2str(alphafres) '_nodes_' ...
        num2str(resultnum) '.pdf'],'-dpdf','-painters')
    print(3,[filename1(1:end-4) '_Multi-Step_ReconstructedFunction_' num2str(alphafres) '_nodes_' ...
        num2str(resultnum) '.png'],'-dpng')
    
    %Plot 2D representation
    alpha0 = Reconstruction(loc).alpha0;
    alphamax = max(alpha0);
    rI = linspace(0,max(LI*tan(alpha0)*1),200);
    
    %Define grid positions to plot
    xlin = linspace(-max(LI*tan(alpha0)*1),max(LI*tan(alpha0)*1),400);
    ylin = linspace(-max(LI*tan(alpha0)*1),max(LI*tan(alpha0)*1),400);
    [xmesh,ymesh] = meshgrid(xlin,ylin);
    
    %Convert to input for functions
    thetamesh = atan2(ymesh,xmesh);
    rmesh = sqrt(xmesh.^2+ymesh.^2);
    
    %Create interpolants
    II1F = griddedInterpolant(rI,Reconstruction(loc).image1,'linear','nearest');
    II2F = griddedInterpolant(rI,Reconstruction(loc).image2,'linear','nearest');
    
    %Calculate images
    IIr1 = II1F(rmesh);
    IIr2 = II2F(rmesh);
    
    %Make images
    figure(1021)
    clf
    % imagesc(x2D*100,x2D*100,1-pI1./max(pI1(:)))
    imagesc(xlin,ylin,1-IIr1) %Normalized to 15 T
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
    title(['$B=' num2str(B0) '$ T, $E_p=$' num2str(E1) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    print(1021,['ReconstructedImage_2D_' num2str(round(B0)) 'T_' num2str(round(E1)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
    print(1021,['ReconstructedImage_2D_' num2str(round(B0)) 'T_' num2str(round(E1)) '_MeV.png'],'-dpng')

    %Make images
    figure(1022)
    clf
    % imagesc(x2D*100,x2D*100,1-pI2./max(pI1(:)))
    imagesc(xlin,ylin,1-IIr2) %Normalized to 15 T
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
    title(['$B=' num2str(B0) '$ T, $E_p=$' num2str(E2) ' MeV'],'FontSize',40,'Interpreter','LaTeX')
    set(gcf, 'PaperSize', [10 10])
    set(gcf, 'PaperPosition', [0 0 10 10])
    print(1022,['ReconstructedImage_2D_' num2str(round(B0)) 'T_' num2str(round(E2)) '_MeV.pdf'],'-dpdf','-painters','-bestfit')
    print(1022,['ReconstructedImage_2D_' num2str(round(B0)) 'T_' num2str(round(E2)) '_MeV.png'],'-dpng')
    
end
end