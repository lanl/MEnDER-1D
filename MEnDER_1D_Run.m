% MEnDER-1D run script
% This script will guide the user through using the 1D field reconstruction
% algorithm for the test cases presented in: J. M. Levesque and L. J.
% Beesley, "A new method for Reconstructing magnetic deflections from sets
% of proton images," (2021).

% Need to call the type of test (Gaussian, wire, wire+)
% type = 'Gaussian';
% type = 'wire';
type = 'twire';

% Probe proton parameters
Ep = [3,14.7]*1e6; %Proton energy array (eV)

% System properties
LO = 1e-2;          %distance from proton source to system center (m)
LI = 15e-2;         %distance from proton source to image plane (m)
LOI = LI-LO;        %distance from system center to image plane (m)
dL = 4e-3;          %length of the interaction region (m)

% Image parameters
dI = 10e-2;         %width of image plane (m)
IRes = 400;         %image size
nProtons = 2e7;     %number of protons per 2D image
gaussParam = 2;     %Gaussian blur
noiseParam = 0.0;   %Gaussian noise variance

% Basic (exact duplicate of paper) or configurable (user gives parameters)

% Give parameters for the field
switch type
    case 'Gaussian'
        % Gaussian ellipsoid
        R = 0.0e-3;         %radius of wire (m) (also where to start reconstruction)
        a = 0.75*1e-3;      %ellipsoid radius (m)
        b = 2e-3;           %ellipsoid length (m)
        B0 = -20;             %maximum magnetic field at wire surface (T)
    case 'wire'
        % Wire field
        R = 0.38e-3;        %radius of wire (m) (if used)
        B0 = -5;             %maximum magnetic field at wire surface (T)
    case 'twire'
        % Truncated wire field
        R = 0.38e-3;        %radius of wire (m) (if used)
        B0 = -5;             %maximum magnetic field at wire surface (T)
end

%% Run field initialization
init = 1;   %Whether to create new file or just get filenames

switch type
    case 'Gaussian'
        % Gaussian ellipsoid
        if init == 1
            Ellipsoid_Initialize(Ep,LO,LI,dL,dI,IRes,nProtons,gaussParam,noiseParam,a,b,B0)
        end
        filename1 = ['Data1_EllipsoidImage_2D_Lineout_' num2str(round(B0))...
            'T_' num2str(round(Ep(1)/1e6)) '_MeV.csv'];
        filename2 = ['Data2_EllipsoidImage_2D_Lineout_' num2str(round(B0))...
            'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'];
    case 'wire'
        % Wire field
        if init == 1
            Wire_Initialize(Ep,LO,LI,dL,dI,IRes,nProtons,gaussParam,R,B0)
        end
        filename1 = ['Data1_WireImage_2D_Lineout_' num2str(round(B0))...
            'T_' num2str(round(Ep(1)/1e6)) '_MeV.csv'];
        filename2 = ['Data2_WireImage_2D_Lineout_' num2str(round(B0))...
            'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'];
    case 'twire'
        % Truncated wire field
        if init == 1
            Truncated_Wire_Initialize(Ep,LO,LI,dL,dI,IRes,nProtons,gaussParam,R,B0)
        end
        filename1 = ['Data1_TruncatedWireImage_2D_Lineout_' num2str(round(B0))...
            'T_' num2str(round(Ep(1)/1e6)) '_MeV.csv'];
        filename2 = ['Data2_TruncatedWireImage_2D_Lineout_' num2str(round(B0))...
            'T_' num2str(round(Ep(2)/1e6)) '_MeV.csv'];
        
end

%% Run reconstruction algorithm
% Parameters for the reconstruction method
alphafres = 21;     %number of nodes of reconstructed deflection function
alphafres2 = 21;%alphafres+(alphafres-1);    %number of nodes aftre reinitialization
ftype = 'pchip';   %interpolant function type, pchip is preferred
iterations = 150e3;  %Maximum number of evolution iterations
weights = [1,2,0.5];  %Initial heuristic weights; if weight 1 or 2 = 0, will only use the one image
resultnum = 1;      %A number identifier for saving output
nruns = 1;         %Number of separate runs to initialize

MEnDER_1D_DE(filename1,filename2,Ep,R,LO,LI,B0,gaussParam,alphafres,alphafres2,...
    ftype,iterations,weights,resultnum,nruns)

