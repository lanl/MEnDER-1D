function R = NewPopulation(I1,I2,E1,E2,LO,LOI,alphares,alpha0,pxI0,dI,xedges,...
    alphafres,alphaA,candidates,vmax,wi,ftype,gaussParam,a,b,dofilt)

%Initial weights
w1 = wi(1);
w2 = wi(2);
w3 = wi(3);

%Randomly initialize the anchor values for all candidates
R.defmatrix = cell(candidates,1);
for can = 1:candidates
    R.defmatrix{can} = 0.1*vmax*rand(alphafres,1).*double(alphaA~=0)';
end

%Initialize error array
R.candidateError = zeros(candidates,1);%,vres);

%Initialize candidate function array
R.wfdeflection = cell(candidates,1);

%Initialize other arrays for each candidate
R.alpha1 = zeros(alphares,1);
R.alpha2 = zeros(alphares,1);
R.pxI1 = zeros(alphares,1);
R.pxI2 = zeros(alphares,1);

R.pI1 = cell(candidates,1);
R.pI2 = cell(candidates,1);
R.Err1 = zeros(candidates,1);
R.Err2 = zeros(candidates,1);
R.Err3 = zeros(candidates,1);
% R.Err4 = zeros(candidates,1);
% R.Err5 = zeros(candidates,1);

%Calculate images for each candidate
for iii = 1:candidates
    defAnchors = R.defmatrix{iii};
    R.wfdeflection{iii} = griddedInterpolant(alphaA,defAnchors,ftype); %functional weights
    
    R.alpha1 = R.wfdeflection{iii}(alpha0);
    R.alpha2 = R.alpha1*sqrt(E1/E2);
    
    R.pxI1 = LO*tan(alpha0)+LOI*tan(R.alpha1+alpha0);
    %         R.pI1{iii} = histcounts(R.pxI1(R.pxI1<max(x1)),xedges);
    R.pI1{iii} = radialweightedHist(pxI0,R.pxI1,xedges);
    
    R.pxI2 = LO*tan(alpha0)+LOI*tan(R.alpha2+alpha0);
    %         R.pI2{iii} = histcounts(R.pxI2(R.pxI2<max(x2)),xedges);
    R.pI2{iii} = radialweightedHist(pxI0,R.pxI2,xedges);
    
    if gaussParam > 0
        R.pI1{iii} = imgaussfilt(R.pI1{iii},gaussParam);
        R.pI2{iii} = imgaussfilt(R.pI2{iii},gaussParam);
    end
    
    if dofilt == 1
        R.pI1{iii} = filter(b,a,R.pI1{iii});
        R.pI2{iii} = filter(b,a,R.pI2{iii});
    end
    
    R.pI2{iii} = R.pI2{iii}*sum(I2)/sum(R.pI2{iii});
    R.pI1{iii} = R.pI1{iii}*sum(I1)/sum(R.pI1{iii});
    
    %Mean relative error of squares
    R.Err1(iii) = mean(abs(R.pI1{iii}.^2-I1.^2)./(I1+2e-2).^2);
    R.Err2(iii) = mean(abs(R.pI2{iii}.^2-I2.^2)./(I2+2e-2).^2);
    
    %Arc length of deflection gradient
    gga = gradient(R.alpha1,10*LO*tan(alpha0));
    ggda = sqrt((10*LO*tan(alpha0(2:end))-10*LO*tan(alpha0(1:end-1))).^2+(gga(2:end)-gga(1:end-1)).^2);
    R.Err3(iii) = (sum(ggda)-(10*LO*tan(alpha0(end))-10*LO*tan(alpha0(1))))/...
        (10*LO*tan(alpha0(end))-10*LO*tan(alpha0(1)));
    
    %Total candidate error
    R.candidateError(iii) = (w1*R.Err1(iii)+w2*R.Err2(iii))*(1+w3*R.Err3(iii));

end