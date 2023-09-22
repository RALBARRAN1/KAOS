clc
clear all
close all
format long
tic

% Note: Run this code to load phase-space grid for kinetic code in FORTRAN.

% -------------------------------------------------------

% SET PHASE-SPACE GRID PARAMETERS:

electronflag= 1;
RE= 638e4; % Earth radius [m]
GG= 667e-13; % Universal gravitational constant [N m^2/kg^2]
ME= 598e22; % Earth mass [kg]
kB= 138e-25; % Boltzmann constant [m^2 kg s^-2 K^-1]
mO= (16e0)*(1.66054e-27); % O+ mass [kg]
Stot= 1e0; % Total number of particle species
IONEXPORTflag= 1 % Set Ion grid Export flag
ENAEXPORTflag= 0 % Set ENA grid Export flag
IONVPERPVECflag= 1 % Set to 1 for 2D ion Vperp grid
HEATINGflag= 1;
SymVparGridflag= 1;

NqGpF= 28e0 ; %33e0; %73e0; % Number of Q Grid Cells (including lower and upper boundary ghost cells i.e. +3)
NVparGpF= 30e0 ; %42e0; % Number of Vpar Grid Cells (even (div by 2 odd) for +/- log10 Vpar grid) (+ 3)
if IONVPERPVECflag == 1
   NVperp1GpF= 30e0 ; %42e0; % Number of Vperp1 Grid Cells (even (div by 2 odd) for +/- log10 Vpar grid) (+ 3)
   NVperp2GpF= 30e0 ; %42e0; % Number of Vperp2 Grid Cells (even (div by 2 odd) for +/- log10 Vpar grid) (+ 3)
   
   if HEATINGflag == 0
       Vperp12sigmaFac= 4; % 8 for thermal 45 for heating Sigma factor with linear grid to resolve thermal core of MB distrib.
       VparsigmaFac= 4; % 8 for thermal 45 for heating
       
       Vperp12NlinRange= (NVperp1GpF)/2e0; % 15 for thermal, 3 for heating
       if (SymVparGridflag == 1)
           VparNlinRange= (NVparGpF)/2e0; % 20 for thermal, 3 for heating
       end
       if (SymVparGridflag == 0)
           VparNlinRange= NVparGpF; % 20 for thermal, 3 for heating
       end
   end
   if HEATINGflag == 1
       Vperp12sigmaFac= 25; % 25 for W0F, W0F1, W0F2, 14 for W0, 15 for VV0Ft, 7 for VV0F, 11 for VV0, 8 for LA Sigma factor with linear grid to resolve thermal core of MB distrib.
       VparsigmaFac= 25; % 25 for W0F, W0F1, W0F2, 14 for W0, 15 for VV0Ft, 7 for VV0F, 11 for VV0, 8 for LA
       
       Vperp12NlinRange= (NVperp1GpF)/2e0; % 15 for thermal, 3 for heating
       if (SymVparGridflag == 1)
           VparNlinRange= (NVparGpF)/2e0; % 20 for thermal, 3 for heating
       end
       if (SymVparGridflag == 0)
           VparNlinRange= NVparGpF; % 20 for thermal, 3 for heating
       end
   end
else
   NVperpGpF= NVparGpF/2e0; % Number of Vperp Grid Cells
end
NVpGpF= NVperp1GpF; % Number of Vp Grid Cells (even (div by 2 odd) for +/- log10 Vp grid) (+ 3)
NVqGpF= NVparGpF; % Number of Vq Grid Cells (even (div by 2 odd) for +/- log10 Vq grid) (+ 3)
NVphiGpF= NVperp2GpF; % Number of Vphi Grid Cells (even (div by 2 odd) for +/- log10 Vphi grid) (+ 3)

if HEATINGflag == 0
    VpphisigmaFac= Vperp12sigmaFac; % 8 for thermal 75 for heating Sigma factor with linear grid to resolve thermal core of MB distrib.
    VqsigmaFac= VparsigmaFac; % 8 for thermal 75 for heating 
    
    VpphiNlinRange= (NVpGpF)/2e0; % 5 for thermal, 3 for heating
    if (SymVparGridflag == 1)
        VqNlinRange= (NVqGpF)/2e0; % 5 for thermal, 3 for heating
    end
    if (SymVparGridflag == 0)
        VqNlinRange= NVqGpF; % 5 for thermal, 3 for heating
    end
end
if HEATINGflag == 1
    VpphisigmaFac= Vperp12sigmaFac; % 8 for thermal 75 for heating Sigma factor with linear grid to resolve thermal core of MB distrib.
    VqsigmaFac= VparsigmaFac; % 8 for thermal 75 for heating 
    
    VpphiNlinRange= (NVpGpF)/2e0; % 5 for thermal, 3 for heating
    if (SymVparGridflag == 1)
        VqNlinRange= (NVqGpF)/2e0; % 5 for thermal, 3 for heating
    end
    if (SymVparGridflag == 0)
        VqNlinRange= NVqGpF; % 5 for thermal, 3 for heating
    end
end

% Define directory for data export:
datadir= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOGridW0F2/'; % CS
% For VISIONS-1 CS and Parametric Study (Ti= 500K, 1000K, 5000K), (Te= 500K, 5000K, 10000K):
dT= [0e0 200e0 1e3];

Ti= 2.5e3+ dT(1); % Ion initialization temperature [K] W0, W0F, W0F1, W0F2 sims
Te= 2.5e3; % Electron initialization temperature [K]

% Initialize as isotropic (in pitch-angle) distributions
TsPerp= Ti; % Ti= (1/3)*TsPar+ (2/3)*TsPerp
TsPar= Ti;
Vperp12sigma= sqrt(kB*TsPerp/mO); % MB sigma for Vperp12
Vparsigma= sqrt(kB*TsPar/mO); % MB sigma for Vpar
Vperp12mean= sqrt(2e0*kB*TsPerp/(pi*mO)); % Mean 1D MB speed
Vparmean= 0e0; % sqrt(2e0*kB*TsPar/(pi*mO)); % Mean 1D MB speed

mNeut= mO; % O neutral mass [kg]
TNeut= 848d0; % O neutral temperature [K] from NRLMSISE-00

Vpphisigma= sqrt(kB*TNeut/mNeut); % MB sigma for Vp, Vphi
Vqsigma= sqrt(kB*TNeut/mNeut); % MB sigma for Vq
Vpphimean= sqrt(2e0*kB*TNeut/(pi*mNeut)); % Mean 1D MB speed
Vqmean= sqrt(2e0*kB*TNeut/(pi*mNeut)); % Mean 1D MB speed

for s= 1:1:Stot
    
    if s == 1 % (O+)
        Specie(s).Nf= 1e0; % Number of dipole flux-tubes
        Specie(s).ms= mO; % Mass [kg]
        Specie(s).qs= 1.602e-19; % Charge [C]
        Specie(s).rads= 152e-12; % Atomic radius [m]
    end
    if s == 2 % (H+)
        Specie(s).Nf= 1e0;
        Specie(s).ms= mO;
        Specie(s).qs= 1.602e-19;
        Specie(s).rads= 53e-12;
    end
    
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot;
    
    if numel(find(isnan(Specie(s).Nf)== 1))~= 0 | numel(find(isinf(Specie(s).Nf)== 1))~= 0 | ...
            numel(find(isempty(Specie(s).Nf)== 1))~= 0
        disp(horzcat('ERROR: Nf has NaN, Inf, or empty array element for particle specie= ', ...
            num2str(s)))
    end
    
    if numel(find(isnan(Specie(s).ms)== 1))~= 0 | numel(find(isinf(Specie(s).ms)== 1))~= 0 | ...
            numel(find(isempty(Specie(s).ms)== 1))~= 0
        disp(horzcat('ERROR: ms has NaN, Inf, or empty array element for particle specie= ', ...
            num2str(s)))
    end
    
    if numel(find(isnan(Specie(s).qs)== 1))~= 0 | numel(find(isinf(Specie(s).qs)== 1))~= 0 | ...
            numel(find(isempty(Specie(s).qs)== 1))~= 0
        disp(horzcat('ERROR: qs has NaN, Inf, or empty array element for particle specie= ', ...
            num2str(s)))
    end
    
end

% -------------------------------------------------------

% SET PRELIMINARY NUMBER OF CONFIGURATION-SPACE GRID CELLS PER PARTICLE SPECIES AND
% FLUX TUBE:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        
        % Total number of q grid values (Make even number to avoid equatorial
        % grid boundary):
        if s == 1 & f < Specie(s).Nf/2
            Specie(s).FluxTube(f).NqGp= NqGpF;
        else if s == 1 & f >= Specie(s).Nf/2
                Specie(s).FluxTube(f).NqGp= NqGpF;
            end
        end
        if s == 2 & f < Specie(s).Nf/2
            Specie(s).FluxTube(f).NqGp= NqGpF;
        else if s == 2 & f >= Specie(s).Nf/2
                Specie(s).FluxTube(f).NqGp= NqGpF;
            end
        end
        
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        
        if numel(find(isnan(Specie(s).FluxTube(f).NqGp)== 1))~= 0 | ...
                numel(find(isinf(Specie(s).FluxTube(f).NqGp)== 1))~= 0 | ...
                numel(find(isempty(Specie(s).FluxTube(f).NqGp)== 1))~= 0
            disp(horzcat('ERROR: NqGp has NaN, Inf, or empty array element for particle specie= ', ...
                num2str(s), ' and flux tube= ', num2str(f)))
        end
        
    end
end

% -------------------------------------------------------

% SET MAGNETIC L-SHELL VALUES AND PRELIMINARY FIELD-ALIGNED GRID CELLS FOR EACH PARTICLE SPECIES
% AND FLUX-TUBE:

% Geomagnetic coords of VISONS 1 latICR, lonICR:
% 58.84N (S)= 90- 58.84= 31.16S	89.54W

% Geomagnetic coords of VISONS 1 from BBELF data latICR, lonICR:
% 68.915 S= 90+ 68.915= 158.92 S 89.54W

% Geomagnetic coords of (Wu '99) latICR, lonICR:
% 73.9 N = 90- 73.9= 16.09S 179.2E
% thetaICRM= 16.09;

thetaICRM= 1.587845981205766e+02; % theta (co-latitude) of max BBELF wave power from VISIONS 1 flight

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        Specie(s).FluxTube(f).pG= 15d0; %7.88d0; % Set grid L-shell
        Specie(s).FluxTube(f).qGA= -0.135d0; % Set lower boundary q value

        if Specie(s).FluxTube(f).qGA <= 0
            SMagHemFlag= 1;
            Specie(s).FluxTube(f).qGB= -0.046e0; % Set upper boundary q value
        end
        if Specie(s).FluxTube(f).qGA > 0
            SMagHemFlag= 0;
            Specie(s).FluxTube(f).qGB= -0.0605e0;
        end
        Specie(s).FluxTube(f).dq= 1e0; % dq length of config cells
        Specie(s).FluxTube(f).NqG= Specie(s).FluxTube(f).NqGp/Specie(s).FluxTube(f).dq- 1e0; % Number of config cells
        
    end
end

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        
        Specie(s).FluxTube(f).qG(:)= linspace(Specie(s).FluxTube(f).qGA, Specie(s).FluxTube(f).qGB, ...
            Specie(s).FluxTube(f).NqGp); % Set preliminary range of q-coord
        
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        
        if numel(find(isnan(Specie(s).FluxTube(f).pG)== 1))~= 0 | ...
                numel(find(isinf(Specie(s).FluxTube(f).pG)== 1))~= 0 | ...
                numel(find(isempty(Specie(s).FluxTube(f).pG)== 1))~= 0
            disp(horzcat('ERROR: pG has NaN, Inf, or empty array element for particle specie= ', ...
                num2str(s), ' and flux tube= ', num2str(f)))
        end
        
        if numel(find(isnan(Specie(s).FluxTube(f).dq)== 1))~= 0 | ...
                numel(find(isinf(Specie(s).FluxTube(f).dq)== 1))~= 0 | ...
                numel(find(isempty(Specie(s).FluxTube(f).dq)== 1))~= 0
            disp(horzcat('ERROR: dq has NaN, Inf, or empty array element for particle specie= ', ...
                num2str(s), ' and flux tube= ', num2str(f)))
        end
        
        if numel(find(isnan(Specie(s).FluxTube(f).NqG)== 1))~= 0 | ...
                numel(find(isinf(Specie(s).FluxTube(f).NqG)== 1))~= 0 | ...
                numel(find(isempty(Specie(s).FluxTube(f).NqG)== 1))~= 0
            disp(horzcat('ERROR: NqG has NaN, Inf, or empty array element for particle specie= ', ...
                num2str(s), ' and flux tube= ', num2str(f)))
        end
        
        if numel(find(isnan(Specie(s).FluxTube(f).qGA)== 1))~= 0 | ...
                numel(find(isinf(Specie(s).FluxTube(f).qGA)== 1))~= 0 | ...
                numel(find(isempty(Specie(s).FluxTube(f).qGA)== 1))~= 0
            disp(horzcat('ERROR: qGA has NaN, Inf, or empty array element for particle specie= ', ...
                num2str(s), ' and flux tube= ', num2str(f)))
        end
        
        if numel(find(isnan(Specie(s).FluxTube(f).qGB)== 1))~= 0 | ...
                numel(find(isinf(Specie(s).FluxTube(f).qGB)== 1))~= 0 | ...
                numel(find(isempty(Specie(s).FluxTube(f).qGB)== 1))~= 0
            disp(horzcat('ERROR: qGB has NaN, Inf, or empty array element for particle specie= ', ...
                num2str(s), ' and flux tube= ', num2str(f)))
        end
        
    end
end

% -------------------------------------------------------

% SET FINAL RANGE OF FIELD-ALINGED GRID AND INTIAL TEMPERATURE PER PARTICLE SPECIES AND FLUX TUBE:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        for Qind= 1:1:Specie(s).FluxTube(f).NqGp;
            
            % Set final range of q-coord
            Specie(s).FluxTube(f).QCell(Qind).qG= Specie(s).FluxTube(f).qG(Qind);
            
        end
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqGp;
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).qG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).qG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).qG)== 1))~= 0
                disp(horzcat('ERROR: qG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
        end
    end
end

% -------------------------------------------------------

% RUN QUARTIC DIPOLE POLYNOMIAL ROOT-FINDER FOR TOTAL FA CARTESIAN VALUES OF GRID LIMITS:

for s= 1:1:Stot;
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqGp;
            
            Specie(s).FluxTube(f).QCell(Qind).Ap= 0e0; % From quartic form
            Specie(s).FluxTube(f).QCell(Qind).Bp= 0e0;
            Specie(s).FluxTube(f).QCell(Qind).Cp= 1e0/(Specie(s).FluxTube(f).pG* ...
                Specie(s).FluxTube(f).QCell(Qind).qG^2e0);
            Specie(s).FluxTube(f).QCell(Qind).Dp= -1e0/(Specie(s).FluxTube(f).QCell(Qind).qG^2e0);
            
            Specie(s).FluxTube(f).QCell(Qind).Ab= -Specie(s).FluxTube(f).QCell(Qind).Bp; % From resolvent cubic form
            Specie(s).FluxTube(f).QCell(Qind).Bb= Specie(s).FluxTube(f).QCell(Qind).Ap* ...
                Specie(s).FluxTube(f).QCell(Qind).Cp- 4e0*Specie(s).FluxTube(f).QCell(Qind).Dp;
            Specie(s).FluxTube(f).QCell(Qind).Cb= 4e0*Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Dp- Specie(s).FluxTube(f).QCell(Qind).Cp^2e0 ...
                - Specie(s).FluxTube(f).QCell(Qind).Ap^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Dp;
            
            Specie(s).FluxTube(f).QCell(Qind).D3= Specie(s).FluxTube(f).QCell(Qind).Bb^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ab^2e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Cb* ...
                Specie(s).FluxTube(f).QCell(Qind).Ab^3e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Bb^3e0 ...
                + 18e0*Specie(s).FluxTube(f).QCell(Qind).Ab* ...
                Specie(s).FluxTube(f).QCell(Qind).Bb* ...
                Specie(s).FluxTube(f).QCell(Qind).Cb ...
                - 27e0*Specie(s).FluxTube(f).QCell(Qind).Cb^2e0; % Resolvent cubic discriminant
            
            Specie(s).FluxTube(f).QCell(Qind).D4= (Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Cp^3e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^3e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp^3e0 ...
                + 18e0*Specie(s).FluxTube(f).QCell(Qind).Cp^3e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap ...
                - 27e0*Specie(s).FluxTube(f).QCell(Qind).Cp^4e0+ ...
                256e0*Specie(s).FluxTube(f).QCell(Qind).Dp^3e0) ...
                + Specie(s).FluxTube(f).QCell(Qind).Dp* ...
                (-4e0*Specie(s).FluxTube(f).QCell(Qind).Bp^3e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                + 18e0*Specie(s).FluxTube(f).QCell(Qind).Cp* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^3e0 ...
                + 16e0*Specie(s).FluxTube(f).QCell(Qind).Bp^4e0- ...
                80e0*Specie(s).FluxTube(f).QCell(Qind).Cp* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp^2e0 ...
                *Specie(s).FluxTube(f).QCell(Qind).Ap- ...
                6e0*Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                + 144e0*Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp)+ ...
                Specie(s).FluxTube(f).QCell(Qind).Dp^2e0 ...
                *(-27e0*Specie(s).FluxTube(f).QCell(Qind).Ap^4e0+ ...
                144e0*Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                - 128e0*Specie(s).FluxTube(f).QCell(Qind).Bp^2e0- ...
                192e0*Specie(s).FluxTube(f).QCell(Qind).Cp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap); % Original quartic disrcriminant
            
            Specie(s).FluxTube(f).QCell(Qind).Abb= Specie(s).FluxTube(f).QCell(Qind).Bb- ...
                Specie(s).FluxTube(f).QCell(Qind).Ab^2e0/3e0;
            Specie(s).FluxTube(f).QCell(Qind).Bbb= 2e0*Specie(s).FluxTube(f).QCell(Qind).Ab^3e0/27e0- ...
                Specie(s).FluxTube(f).QCell(Qind).Ab ...
                *Specie(s).FluxTube(f).QCell(Qind).Bb/3e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).Cb;
            Specie(s).FluxTube(f).QCell(Qind).Delta= -Specie(s).FluxTube(f).QCell(Qind).Bbb/2e0;
            Specie(s).FluxTube(f).QCell(Qind).Epsilon= Specie(s).FluxTube(f).QCell(Qind).Abb/3e0;
            
            Specie(s).FluxTube(f).QCell(Qind).Thetap= acos(Specie(s).FluxTube(f).QCell(Qind).Delta ...
                /sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon^3e0));
            Specie(s).FluxTube(f).QCell(Qind).sigma21= ...
                2e0*sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon) ...
                *cos(Specie(s).FluxTube(f).QCell(Qind).Thetap/3e0) ...
                - Specie(s).FluxTube(f).QCell(Qind).Ab/3e0; % Resolvent cubic roots for D3< 0
            Specie(s).FluxTube(f).QCell(Qind).sigma22= ...
                2e0*sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon) ...
                *cos((Specie(s).FluxTube(f).QCell(Qind).Thetap+ 2e0*pi)/3e0) ...
                - Specie(s).FluxTube(f).QCell(Qind).Ab/3e0;
            Specie(s).FluxTube(f).QCell(Qind).sigma23= ...
                2e0*sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon)...
                *cos((Specie(s).FluxTube(f).QCell(Qind).Thetap+ 4e0*pi)/3e0) ...
                - Specie(s).FluxTube(f).QCell(Qind).Ab/3e0;
            Specie(s).FluxTube(f).QCell(Qind).sigma2= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).sigma23)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).sigma23)^2e0); % Choose real cubic root
            Specie(s).FluxTube(f).QCell(Qind).sigma= ...
                Specie(s).FluxTube(f).QCell(Qind).sigma2;
            
            Specie(s).FluxTube(f).QCell(Qind).mu= ...
                sqrt(Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                Specie(s).FluxTube(f).QCell(Qind).Bp+ ...
                Specie(s).FluxTube(f).QCell(Qind).sigma);
            
            if Specie(s).FluxTube(f).QCell(Qind).mu ~= 0;
                Specie(s).FluxTube(f).QCell(Qind).nu= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    Specie(s).FluxTube(f).QCell(Qind).mu^2e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp+ ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).Ap* ...
                    Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    8e0*Specie(s).FluxTube(f).QCell(Qind).Cp- ...
                    Specie(s).FluxTube(f).QCell(Qind).Ap^3e0)/ ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).mu));
                Specie(s).FluxTube(f).QCell(Qind).pii= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    Specie(s).FluxTube(f).QCell(Qind).mu^2e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).Ap* ...
                    Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    8e0*Specie(s).FluxTube(f).QCell(Qind).Cp- ...
                    Specie(s).FluxTube(f).QCell(Qind).Ap^3e0)/ ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).mu));
            else
                Specie(s).FluxTube(f).QCell(Qind).nu= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp+ ...
                    2e0*sqrt(Specie(s).FluxTube(f).QCell(Qind).sigma^2e0- ...
                    4e0*Specie(s).FluxTube(f).QCell(Qind).Dp));
                Specie(s).FluxTube(f).QCell(Qind).pii= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    2e0*sqrt(Specie(s).FluxTube(f).QCell(Qind).sigma^2e0- ...
                    4e0*Specie(s).FluxTube(f).QCell(Qind).Dp));
            end
            
            Specie(s).FluxTube(f).QCell(Qind).gamma1p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).nu/2e0; % Original quartic roots
            Specie(s).FluxTube(f).QCell(Qind).gamma2p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0- ...
                Specie(s).FluxTube(f).QCell(Qind).nu/2e0;
            Specie(s).FluxTube(f).QCell(Qind).gamma3p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0- ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).pii/2e0;
            Specie(s).FluxTube(f).QCell(Qind).gamma4p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0- ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0- ...
                Specie(s).FluxTube(f).QCell(Qind).pii/2e0;
            
            Specie(s).FluxTube(f).QCell(Qind).gamma1= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma1p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma1p)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).gamma2= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma2p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma2p)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).gamma3= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma3p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma3p)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).gamma4= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma4p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma4p)^2e0);
            
            Specie(s).FluxTube(f).QCell(Qind).r1= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma1*RE; % Revert quartic roots with (q, p) into (r, theta)
            Specie(s).FluxTube(f).QCell(Qind).r2= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma2*RE;
            Specie(s).FluxTube(f).QCell(Qind).r3= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma3*RE;
            Specie(s).FluxTube(f).QCell(Qind).r4= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma4*RE;
            Specie(s).FluxTube(f).QCell(Qind).theta1= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r1/ ...
                (RE*Specie(s).FluxTube(f).pG)));
            Specie(s).FluxTube(f).QCell(Qind).theta2= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r2/...
                (RE*Specie(s).FluxTube(f).pG)));
            Specie(s).FluxTube(f).QCell(Qind).theta3= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r3/ ...
                (RE*Specie(s).FluxTube(f).pG)));
            Specie(s).FluxTube(f).QCell(Qind).theta4= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r4/ ...
                (RE*Specie(s).FluxTube(f).pG)));
            
            Specie(s).FluxTube(f).QCell(Qind).qtest1= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta1)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r1^2e0); % Revert back to (q, p)
            Specie(s).FluxTube(f).QCell(Qind).qtest2= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta2)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r2^2e0);
            Specie(s).FluxTube(f).QCell(Qind).qtest3= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta3)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r3^2e0);
            Specie(s).FluxTube(f).QCell(Qind).qtest4= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta4)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r4^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest1= ...
                Specie(s).FluxTube(f).QCell(Qind).r1/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta1))^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest2= ...
                Specie(s).FluxTube(f).QCell(Qind).r2/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta2))^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest3= ...
                Specie(s).FluxTube(f).QCell(Qind).r3/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta3))^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest4= ...
                Specie(s).FluxTube(f).QCell(Qind).r4/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta4))^2e0);
            
            % Note gamma 3 is real positive root. For q< 0 (q> 0), phase-shift theta by
            % pi (0) to get correct sign of q. (above and below dipole equator)
            
            Specie(s).FluxTube(f).QCell(Qind).rfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).r3; % Get final (r, theta) values
            
            if Specie(s).FluxTube(f).QCell(Qind).qG< 0 | ...
                    Specie(s).FluxTube(f).QCell(Qind).qG== 0e0; % Phase shift root 3 solution
                Specie(s).FluxTube(f).QCell(Qind).thetafinalG= ...
                    pi- Specie(s).FluxTube(f).QCell(Qind).theta3;
            else if Specie(s).FluxTube(f).QCell(Qind).qG> 0;
                    Specie(s).FluxTube(f).QCell(Qind).thetafinalG= ...
                        Specie(s).FluxTube(f).QCell(Qind).theta3;
                end
            end
            
            % Get final (x, y, z) values
            
            Specie(s).FluxTube(f).QCell(Qind).phifinalG= pi/2e0;
            Specie(s).FluxTube(f).QCell(Qind).xfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)* ...
                cos(Specie(s).FluxTube(f).QCell(Qind).phifinalG);
            Specie(s).FluxTube(f).QCell(Qind).yfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).phifinalG);
            Specie(s).FluxTube(f).QCell(Qind).zfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG* ...
                cos(Specie(s).FluxTube(f).QCell(Qind).thetafinalG);
            
            %             Specie(s).FluxTube(f).QCell(Qind).phifinalG= ...
            %                 atan2(Specie(s).FluxTube(f).QCell(Qind).yfinalG, ...
            %                 Specie(s).FluxTube(f).QCell(Qind).xfinalG);
            
            % Get final (q, p) value to compare with initial input:
            
            Specie(s).FluxTube(f).QCell(Qind).qfinalG= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).rfinalG^2e0);
            Specie(s).FluxTube(f).QCell(Qind).pfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).thetafinalG))^2e0);
            
            Specie(s).FluxTube(f).QCell(Qind).qfinalp= ...
                isnan(real(Specie(s).FluxTube(f).QCell(Qind).qfinalG)); % Set NaN values to zero
            Specie(s).FluxTube(f).QCell(Qind).pfinalp= ...
                isnan(real(Specie(s).FluxTube(f).QCell(Qind).pfinalG)); % Set NaN values to zero
            if Specie(s).FluxTube(f).QCell(Qind).qfinalp == 1;
                Specie(s).FluxTube(f).QCell(Qind).qfinalG= 0e0;
            end
            if Specie(s).FluxTube(f).QCell(Qind).pfinalp == 1;
                Specie(s).FluxTube(f).QCell(Qind).pfinalG= 0e0;
            end
            
            if abs(Specie(s).FluxTube(f).QCell(Qind).qG- ...
                    Specie(s).FluxTube(f).QCell(Qind).qfinalG) > 1e-14
                disp('BAD Q')
            end
            
            if abs(Specie(s).FluxTube(f).pG- Specie(s).FluxTube(f).QCell(Qind).pfinalG) > 1e-13
                disp('BAD P')
            end
            
        end
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqGp;
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).rfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).rfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).rfinalG)== 1))~= 0
                disp(horzcat('ERROR: rfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)== 1))~= 0
                disp(horzcat('ERROR: thetafinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).phifinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).phifinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).phifinalG)== 1))~= 0
                disp(horzcat('ERROR: phifinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).xfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).xfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).xfinalG)== 1))~= 0
                disp(horzcat('ERROR: xfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).yfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).yfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).yfinalG)== 1))~= 0
                disp(horzcat('ERROR: yfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).zfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).zfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).zfinalG)== 1))~= 0
                disp(horzcat('ERROR: zfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).qfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).qfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).qfinalG)== 1))~= 0
                disp(horzcat('ERROR: qfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).pfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).pfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).pfinalG)== 1))~= 0
                disp(horzcat('ERROR: pfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
        end
    end
end

% ---------------------------------------------

% SET CONFIGURATION SPACE GRID LIMITS:

for s= 1:1:Stot;
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqG;
            
            Specie(s).FluxTube(f).QCell(Qind).qGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ 1).qfinalG; % lower q limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).qGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).qfinalG; % upper q limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).qGC= ...
                (Specie(s).FluxTube(f).QCell(Qind).qGL+ ...
                Specie(s).FluxTube(f).QCell(Qind).qGH)/2e0; % center q values of FA cells
                        
            Specie(s).FluxTube(f).qGL(Qind)= Specie(s).FluxTube(f).QCell(Qind).qGL;
            Specie(s).FluxTube(f).qGH(Qind)= Specie(s).FluxTube(f).QCell(Qind).qGH;
            Specie(s).FluxTube(f).qGC(Qind)= Specie(s).FluxTube(f).QCell(Qind).qGC;
            
            Specie(s).FluxTube(f).QCell(Qind).pGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ 1).pfinalG; % lower p limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).pGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).pfinalG; % upper p limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).pGC= ...
                (Specie(s).FluxTube(f).QCell(Qind).pGL+ ...
                Specie(s).FluxTube(f).QCell(Qind).pGH)/2e0; % center p values of FA cells
            
            Specie(s).FluxTube(f).pGC(Qind)= Specie(s).FluxTube(f).QCell(Qind).pGC;
            
            Specie(s).FluxTube(f).QCell(Qind).phiGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ 1).phifinalG; % lower phid limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).phiGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).phifinalG; % upper phid limits of FA cells
            
            Specie(s).FluxTube(f).QCell(Qind).rGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ 1).rfinalG; % lower r limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).rGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).rfinalG; % upper r limits of FA cells
            
            Specie(s).FluxTube(f).rGL(Qind)= Specie(s).FluxTube(f).QCell(Qind).rGL;
            Specie(s).FluxTube(f).rGH(Qind)= Specie(s).FluxTube(f).QCell(Qind).rGH;
            
            Specie(s).FluxTube(f).QCell(Qind).thetaGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ 1).thetafinalG; % lower theta limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).thetaGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).thetafinalG; % upper theta limits of FA cells
            
            Specie(s).FluxTube(f).QCell(Qind).yGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ 1).yfinalG; % lower y limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).yGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).yfinalG; % upper y limits of FA cells
            
            Specie(s).FluxTube(f).QCell(Qind).zGL= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)* ...
                Specie(s).FluxTube(f).dq+ 1).zfinalG; % lower z limits of FA cells
            Specie(s).FluxTube(f).QCell(Qind).zGH= ...
                Specie(s).FluxTube(f).QCell((Qind- 1)*Specie(s).FluxTube(f).dq+ ...
                1+ Specie(s).FluxTube(f).dq).zfinalG; % upper z limits of FA cells
            
            Specie(s).FluxTube(f).QCell(Qind).xGL= 0e0; % Let phi= pi/2 s.t. all FA cells have x= 0
            Specie(s).FluxTube(f).QCell(Qind).xGH= 0e0;
            
        end
    end
end

% -------------------------------------------------------

% RUN QUARTIC DIPOLE POLYNOMIAL ROOT-FINDER FOR TOTAL FA CARTESIAN VALUES OF GRID CENTER VALUES:

for s= 1:1:Stot;
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqG;
            
            Specie(s).FluxTube(f).QCell(Qind).Ap= 0e0; % From quartic form
            Specie(s).FluxTube(f).QCell(Qind).Bp= 0e0;
            Specie(s).FluxTube(f).QCell(Qind).Cp= 1e0/(Specie(s).FluxTube(f).QCell(Qind).pGC* ...
                Specie(s).FluxTube(f).QCell(Qind).qGC^2e0);
            Specie(s).FluxTube(f).QCell(Qind).Dp= -1e0/(Specie(s).FluxTube(f).QCell(Qind).qGC^2e0);
            
            Specie(s).FluxTube(f).QCell(Qind).Ab= -Specie(s).FluxTube(f).QCell(Qind).Bp; % From resolvent cubic form
            Specie(s).FluxTube(f).QCell(Qind).Bb= Specie(s).FluxTube(f).QCell(Qind).Ap* ...
                Specie(s).FluxTube(f).QCell(Qind).Cp- 4e0*Specie(s).FluxTube(f).QCell(Qind).Dp;
            Specie(s).FluxTube(f).QCell(Qind).Cb= 4e0*Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Dp- Specie(s).FluxTube(f).QCell(Qind).Cp^2e0 ...
                - Specie(s).FluxTube(f).QCell(Qind).Ap^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Dp;
            
            Specie(s).FluxTube(f).QCell(Qind).D3= Specie(s).FluxTube(f).QCell(Qind).Bb^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ab^2e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Cb* ...
                Specie(s).FluxTube(f).QCell(Qind).Ab^3e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Bb^3e0 ...
                + 18e0*Specie(s).FluxTube(f).QCell(Qind).Ab* ...
                Specie(s).FluxTube(f).QCell(Qind).Bb* ...
                Specie(s).FluxTube(f).QCell(Qind).Cb ...
                - 27e0*Specie(s).FluxTube(f).QCell(Qind).Cb^2e0; % Resolvent cubic discriminant
            
            Specie(s).FluxTube(f).QCell(Qind).D4= (Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Cp^3e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^3e0 ...
                - 4e0*Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp^3e0 ...
                + 18e0*Specie(s).FluxTube(f).QCell(Qind).Cp^3e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap ...
                - 27e0*Specie(s).FluxTube(f).QCell(Qind).Cp^4e0+ ...
                256e0*Specie(s).FluxTube(f).QCell(Qind).Dp^3e0) ...
                + Specie(s).FluxTube(f).QCell(Qind).Dp* ...
                (-4e0*Specie(s).FluxTube(f).QCell(Qind).Bp^3e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                + 18e0*Specie(s).FluxTube(f).QCell(Qind).Cp* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^3e0 ...
                + 16e0*Specie(s).FluxTube(f).QCell(Qind).Bp^4e0- ...
                80e0*Specie(s).FluxTube(f).QCell(Qind).Cp* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp^2e0 ...
                *Specie(s).FluxTube(f).QCell(Qind).Ap- ...
                6e0*Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                + 144e0*Specie(s).FluxTube(f).QCell(Qind).Cp^2e0* ...
                Specie(s).FluxTube(f).QCell(Qind).Bp)+ ...
                Specie(s).FluxTube(f).QCell(Qind).Dp^2e0 ...
                *(-27e0*Specie(s).FluxTube(f).QCell(Qind).Ap^4e0+ ...
                144e0*Specie(s).FluxTube(f).QCell(Qind).Bp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap^2e0 ...
                - 128e0*Specie(s).FluxTube(f).QCell(Qind).Bp^2e0- ...
                192e0*Specie(s).FluxTube(f).QCell(Qind).Cp* ...
                Specie(s).FluxTube(f).QCell(Qind).Ap); % Original quartic disrcriminant
            
            Specie(s).FluxTube(f).QCell(Qind).Abb= Specie(s).FluxTube(f).QCell(Qind).Bb- ...
                Specie(s).FluxTube(f).QCell(Qind).Ab^2e0/3e0;
            Specie(s).FluxTube(f).QCell(Qind).Bbb= 2e0*Specie(s).FluxTube(f).QCell(Qind).Ab^3e0/27e0- ...
                Specie(s).FluxTube(f).QCell(Qind).Ab ...
                *Specie(s).FluxTube(f).QCell(Qind).Bb/3e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).Cb;
            Specie(s).FluxTube(f).QCell(Qind).Delta= -Specie(s).FluxTube(f).QCell(Qind).Bbb/2e0;
            Specie(s).FluxTube(f).QCell(Qind).Epsilon= Specie(s).FluxTube(f).QCell(Qind).Abb/3e0;
            
            Specie(s).FluxTube(f).QCell(Qind).Thetap= acos(Specie(s).FluxTube(f).QCell(Qind).Delta ...
                /sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon^3e0));
            Specie(s).FluxTube(f).QCell(Qind).sigma21= ...
                2e0*sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon) ...
                *cos(Specie(s).FluxTube(f).QCell(Qind).Thetap/3e0) ...
                - Specie(s).FluxTube(f).QCell(Qind).Ab/3e0; % Resolvent cubic roots for D3< 0
            Specie(s).FluxTube(f).QCell(Qind).sigma22= ...
                2e0*sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon) ...
                *cos((Specie(s).FluxTube(f).QCell(Qind).Thetap+ 2e0*pi)/3e0) ...
                - Specie(s).FluxTube(f).QCell(Qind).Ab/3e0;
            Specie(s).FluxTube(f).QCell(Qind).sigma23= ...
                2e0*sqrt(-Specie(s).FluxTube(f).QCell(Qind).Epsilon)...
                *cos((Specie(s).FluxTube(f).QCell(Qind).Thetap+ 4e0*pi)/3e0) ...
                - Specie(s).FluxTube(f).QCell(Qind).Ab/3e0;
            Specie(s).FluxTube(f).QCell(Qind).sigma2= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).sigma23)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).sigma23)^2e0); % Choose real cubic root
            Specie(s).FluxTube(f).QCell(Qind).sigma= ...
                Specie(s).FluxTube(f).QCell(Qind).sigma2;
            
            Specie(s).FluxTube(f).QCell(Qind).mu= ...
                sqrt(Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                Specie(s).FluxTube(f).QCell(Qind).Bp+ ...
                Specie(s).FluxTube(f).QCell(Qind).sigma);
            
            if Specie(s).FluxTube(f).QCell(Qind).mu ~= 0;
                Specie(s).FluxTube(f).QCell(Qind).nu= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    Specie(s).FluxTube(f).QCell(Qind).mu^2e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp+ ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).Ap* ...
                    Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    8e0*Specie(s).FluxTube(f).QCell(Qind).Cp- ...
                    Specie(s).FluxTube(f).QCell(Qind).Ap^3e0)/ ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).mu));
                Specie(s).FluxTube(f).QCell(Qind).pii= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    Specie(s).FluxTube(f).QCell(Qind).mu^2e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).Ap* ...
                    Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    8e0*Specie(s).FluxTube(f).QCell(Qind).Cp- ...
                    Specie(s).FluxTube(f).QCell(Qind).Ap^3e0)/ ...
                    (4e0*Specie(s).FluxTube(f).QCell(Qind).mu));
            else
                Specie(s).FluxTube(f).QCell(Qind).nu= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp+ ...
                    2e0*sqrt(Specie(s).FluxTube(f).QCell(Qind).sigma^2e0- ...
                    4e0*Specie(s).FluxTube(f).QCell(Qind).Dp));
                Specie(s).FluxTube(f).QCell(Qind).pii= ...
                    sqrt(3e0*Specie(s).FluxTube(f).QCell(Qind).Ap^2e0/4e0- ...
                    2e0*Specie(s).FluxTube(f).QCell(Qind).Bp- ...
                    2e0*sqrt(Specie(s).FluxTube(f).QCell(Qind).sigma^2e0- ...
                    4e0*Specie(s).FluxTube(f).QCell(Qind).Dp));
            end
            
            Specie(s).FluxTube(f).QCell(Qind).gamma1p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).nu/2e0; % Original quartic roots
            Specie(s).FluxTube(f).QCell(Qind).gamma2p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0- ...
                Specie(s).FluxTube(f).QCell(Qind).nu/2e0;
            Specie(s).FluxTube(f).QCell(Qind).gamma3p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0- ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0+ ...
                Specie(s).FluxTube(f).QCell(Qind).pii/2e0;
            Specie(s).FluxTube(f).QCell(Qind).gamma4p= ...
                -Specie(s).FluxTube(f).QCell(Qind).Ap/4e0- ...
                Specie(s).FluxTube(f).QCell(Qind).mu/2e0- ...
                Specie(s).FluxTube(f).QCell(Qind).pii/2e0;
            
            Specie(s).FluxTube(f).QCell(Qind).gamma1= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma1p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma1p)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).gamma2= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma2p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma2p)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).gamma3= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma3p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma3p)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).gamma4= ...
                sqrt(real(Specie(s).FluxTube(f).QCell(Qind).gamma4p)^2e0+ ...
                imag(Specie(s).FluxTube(f).QCell(Qind).gamma4p)^2e0);
            
            Specie(s).FluxTube(f).QCell(Qind).r1= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma1*RE; % Revert quartic roots with (q, p) into (r, theta)
            Specie(s).FluxTube(f).QCell(Qind).r2= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma2*RE;
            Specie(s).FluxTube(f).QCell(Qind).r3= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma3*RE;
            Specie(s).FluxTube(f).QCell(Qind).r4= ...
                Specie(s).FluxTube(f).QCell(Qind).gamma4*RE;
            Specie(s).FluxTube(f).QCell(Qind).theta1= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r1/ ...
                (RE*Specie(s).FluxTube(f).pG)));
            Specie(s).FluxTube(f).QCell(Qind).theta2= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r2/...
                (RE*Specie(s).FluxTube(f).pG)));
            Specie(s).FluxTube(f).QCell(Qind).theta3= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r3/ ...
                (RE*Specie(s).FluxTube(f).pG)));
            Specie(s).FluxTube(f).QCell(Qind).theta4= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r4/ ...
                (RE*Specie(s).FluxTube(f).pG)));
            
            Specie(s).FluxTube(f).QCell(Qind).qtest1= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta1)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r1^2e0); % Revert back to (q, p)
            Specie(s).FluxTube(f).QCell(Qind).qtest2= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta2)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r2^2e0);
            Specie(s).FluxTube(f).QCell(Qind).qtest3= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta3)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r3^2e0);
            Specie(s).FluxTube(f).QCell(Qind).qtest4= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).theta4)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).r4^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest1= ...
                Specie(s).FluxTube(f).QCell(Qind).r1/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta1))^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest2= ...
                Specie(s).FluxTube(f).QCell(Qind).r2/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta2))^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest3= ...
                Specie(s).FluxTube(f).QCell(Qind).r3/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta3))^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ptest4= ...
                Specie(s).FluxTube(f).QCell(Qind).r4/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).theta4))^2e0);
            
            % Note gamma 3 is real positive root. For q< 0 (q> 0), phase-shift theta by
            % pi (0) to get correct sign of q. (above and below dipole equator)
            
            Specie(s).FluxTube(f).QCell(Qind).rfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).r3; % Get final (r, theta) values
            
            if Specie(s).FluxTube(f).QCell(Qind).qG< 0 | ...
                    Specie(s).FluxTube(f).QCell(Qind).qG== 0e0; % Phase shift root 3 solution
                Specie(s).FluxTube(f).QCell(Qind).thetafinalG= ...
                    pi- Specie(s).FluxTube(f).QCell(Qind).theta3;
            else if Specie(s).FluxTube(f).QCell(Qind).qG> 0;
                    Specie(s).FluxTube(f).QCell(Qind).thetafinalG= ...
                        Specie(s).FluxTube(f).QCell(Qind).theta3;
                end
            end
            
            % Get final (x, y, z) values
            
            Specie(s).FluxTube(f).QCell(Qind).phifinalG= pi/2e0;
            Specie(s).FluxTube(f).QCell(Qind).xfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)* ...
                cos(Specie(s).FluxTube(f).QCell(Qind).phifinalG);
            Specie(s).FluxTube(f).QCell(Qind).yfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).phifinalG);
            Specie(s).FluxTube(f).QCell(Qind).zfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG* ...
                cos(Specie(s).FluxTube(f).QCell(Qind).thetafinalG);
            
            %             Specie(s).FluxTube(f).QCell(Qind).phifinalG= ...
            %                 atan2(Specie(s).FluxTube(f).QCell(Qind).yfinalG, ...
            %                 Specie(s).FluxTube(f).QCell(Qind).xfinalG);
            
            % Get final (q, p) value to compare with initial input:
            
            Specie(s).FluxTube(f).QCell(Qind).qfinalG= ...
                (RE^2e0)*cos(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)/ ...
                (Specie(s).FluxTube(f).QCell(Qind).rfinalG^2e0);
            Specie(s).FluxTube(f).QCell(Qind).pfinalG= ...
                Specie(s).FluxTube(f).QCell(Qind).rfinalG/ ...
                (RE*(sin(Specie(s).FluxTube(f).QCell(Qind).thetafinalG))^2e0);
            
            Specie(s).FluxTube(f).QCell(Qind).qfinalp= ...
                isnan(real(Specie(s).FluxTube(f).QCell(Qind).qfinalG)); % Set NaN values to zero
            Specie(s).FluxTube(f).QCell(Qind).pfinalp= ...
                isnan(real(Specie(s).FluxTube(f).QCell(Qind).pfinalG)); % Set NaN values to zero
            if Specie(s).FluxTube(f).QCell(Qind).qfinalp == 1;
                Specie(s).FluxTube(f).QCell(Qind).qfinalG= 0e0;
            end
            if Specie(s).FluxTube(f).QCell(Qind).pfinalp == 1;
                Specie(s).FluxTube(f).QCell(Qind).pfinalG= 0e0;
            end
            
            if abs(Specie(s).FluxTube(f).QCell(Qind).qGC- ...
                    Specie(s).FluxTube(f).QCell(Qind).qfinalG) > 1e-14
                disp('BAD Q')
            end
            
            if abs(Specie(s).FluxTube(f).QCell(Qind).pGC- Specie(s).FluxTube(f).QCell(Qind).pfinalG) > 1e-13
                disp('BAD P')
            end
            
            % SET CONFIGURATION SPACE GRID CENTER VALUES, METRIC FACTORS, AND VOLUMES:
            
            Specie(s).FluxTube(f).QCell(Qind).rGC= Specie(s).FluxTube(f).QCell(Qind).rfinalG;
            Specie(s).FluxTube(f).QCell(Qind).thetaGC= Specie(s).FluxTube(f).QCell(Qind).thetafinalG;
            Specie(s).FluxTube(f).QCell(Qind).phiGC= Specie(s).FluxTube(f).QCell(Qind).phifinalG;
            Specie(s).FluxTube(f).QCell(Qind).xGC= Specie(s).FluxTube(f).QCell(Qind).xfinalG;
            Specie(s).FluxTube(f).QCell(Qind).yGC= Specie(s).FluxTube(f).QCell(Qind).yfinalG;
            Specie(s).FluxTube(f).QCell(Qind).zGC= Specie(s).FluxTube(f).QCell(Qind).zfinalG;
            
            Specie(s).FluxTube(f).phiGC(Qind)= Specie(s).FluxTube(f).QCell(Qind).phiGC;
            
            Specie(s).FluxTube(f).rGCp(Qind)= Specie(s).FluxTube(f).QCell(Qind).rGC;
            
            Specie(s).FluxTube(f).QCell(Qind).ellGC= ...
                (1e0+ 3e0*cos(Specie(s).FluxTube(f).QCell(Qind).thetaGC)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ellGL= ...
                (1e0+ 3e0*cos(Specie(s).FluxTube(f).QCell(Qind).thetaGL)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).ellGH= ...
                (1e0+ 3e0*cos(Specie(s).FluxTube(f).QCell(Qind).thetaGH)^2e0);
                        
            Specie(s).FluxTube(f).QCell(Qind).hqC= ...
                abs(Specie(s).FluxTube(f).QCell(Qind).rGC^3e0/ ...
                (RE^2e0*sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC))); % q metric factor centered in FA cells
            Specie(s).FluxTube(f).QCell(Qind).hpC= ...
                abs(RE*sin(Specie(s).FluxTube(f).QCell(Qind).thetaGC)^3e0/ ...
                (sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC))); % p metric factor centered in FA cells
            Specie(s).FluxTube(f).QCell(Qind).hphiC= ...
                abs(Specie(s).FluxTube(f).QCell(Qind).rGC* ...
                sin(Specie(s).FluxTube(f).QCell(Qind).thetaGC)); % phi metric factor centered in FA cells
            
            % Note: Make dpC equal to the Lshell drift limit over entire simulation and dphidC equal to the
            % phid drift limit in kinetic solver f90 subroutine.
            
            Specie(s).FluxTube(f).QCell(Qind).dqC= ...
                abs(Specie(s).FluxTube(f).QCell(Qind).qGH- ...
                Specie(s).FluxTube(f).QCell(Qind).qGL); % dq across FA cells
            Specie(s).FluxTube(f).QCell(Qind).dpC= 1e-3;  % dp across FA cells [RE]
            Specie(s).FluxTube(f).QCell(Qind).dphiC= pi/180e0; % dphid across FA cells [1 deg= pi/180 rads]
            
            Specie(s).FluxTube(f).QCell(Qind).d3xC= ...
                Specie(s).FluxTube(f).QCell(Qind).hqC* ...
                Specie(s).FluxTube(f).QCell(Qind).hpC* ...
                Specie(s).FluxTube(f).QCell(Qind).hphiC* ...
                Specie(s).FluxTube(f).QCell(Qind).dqC* ...
                Specie(s).FluxTube(f).QCell(Qind).dpC* ...
                Specie(s).FluxTube(f).QCell(Qind).dphiC;
            
            if Specie(s).FluxTube(f).QCell(Qind).d3xC == 0e0
                Specie(s).FluxTube(f).QCell(Qind).d3xC= ...
                    1e-9;
            end
            
            if Specie(s).FluxTube(f).QCell(Qind).d3xC < 0e0
                disp(horzcat('ERROR: Negative d3xC for particle species= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ...
                    ', and Qind= ', num2str(Qind)))
            end
            Specie(s).FluxTube(f).d3xC(Qind)= Specie(s).FluxTube(f).QCell(Qind).d3xC;
        end
        
        % Compute lower (non-ghost cell) boundary flux-tube cross section [m^2]
        Specie(s).FluxTube(f).hpLB= ...
            abs(RE*sin(Specie(s).FluxTube(f).QCell(2).thetaGL)^3e0/ ...
                (sqrt(Specie(s).FluxTube(f).QCell(2).ellGL)));
        Specie(s).FluxTube(f).hphiLB= ...
            abs(Specie(s).FluxTube(f).QCell(2).rGL* ...
            sin(Specie(s).FluxTube(f).QCell(2).thetaGL));
        Specie(s).FluxTube(f).dpLB= 1e-3;
        Specie(s).FluxTube(f).dphiLB= pi/180e0;
                    
        % Compute upper boundary flux-tube cross section [m^2]
        Specie(s).FluxTube(f).hpUB= ...
            abs(RE*sin(Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).thetaGH)^3e0/ ...
                (sqrt(Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).ellGH)));
        Specie(s).FluxTube(f).hphiUB= ...
            abs(Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).rGH* ...
            sin(Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).thetaGH));
        Specie(s).FluxTube(f).dpUB= 1e-3;
        Specie(s).FluxTube(f).dphiUB= pi/180e0;
            
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqGp;
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).rfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).rfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).rfinalG)== 1))~= 0
                disp(horzcat('ERROR: rfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).thetafinalG)== 1))~= 0
                disp(horzcat('ERROR: thetafinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).phifinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).phifinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).phifinalG)== 1))~= 0
                disp(horzcat('ERROR: phifinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).xfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).xfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).xfinalG)== 1))~= 0
                disp(horzcat('ERROR: xfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).yfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).yfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).yfinalG)== 1))~= 0
                disp(horzcat('ERROR: yfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).zfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).zfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).zfinalG)== 1))~= 0
                disp(horzcat('ERROR: zfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).qfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).qfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).qfinalG)== 1))~= 0
                disp(horzcat('ERROR: qfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).pfinalG)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).pfinalG)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).pfinalG)== 1))~= 0
                disp(horzcat('ERROR: pfinalG has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
        end
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqG;
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).qGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).qGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).qGL)== 1))~= 0
                disp(horzcat('ERROR: qGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).qGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).qGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).qGH)== 1))~= 0
                disp(horzcat('ERROR: qGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).qGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).qGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).qGC)== 1))~= 0
                disp(horzcat('ERROR: qGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).pGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).pGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).pGL)== 1))~= 0
                disp(horzcat('ERROR: pGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).pGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).pGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).pGH)== 1))~= 0
                disp(horzcat('ERROR: pGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).pGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).pGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).pGC)== 1))~= 0
                disp(horzcat('ERROR: pGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).phiGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).phiGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).phiGL)== 1))~= 0
                disp(horzcat('ERROR: phiGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).phiGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).phiGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).phiGH)== 1))~= 0
                disp(horzcat('ERROR: phiGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).phiGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).phiGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).phiGC)== 1))~= 0
                disp(horzcat('ERROR: phiGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).rGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).rGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).rGL)== 1))~= 0
                disp(horzcat('ERROR: rGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).rGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).rGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).rGH)== 1))~= 0
                disp(horzcat('ERROR: rGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).rGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).rGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).rGC)== 1))~= 0
                disp(horzcat('ERROR: rGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).thetaGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).thetaGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).thetaGL)== 1))~= 0
                disp(horzcat('ERROR: thetaGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).thetaGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).thetaGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).thetaGH)== 1))~= 0
                disp(horzcat('ERROR: thetaGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).thetaGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).thetaGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).thetaGC)== 1))~= 0
                disp(horzcat('ERROR: thetaGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).ellGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).ellGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).ellGC)== 1))~= 0
                disp(horzcat('ERROR: ellGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).yGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).yGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).yGL)== 1))~= 0
                disp(horzcat('ERROR: yGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).yGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).yGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).yGH)== 1))~= 0
                disp(horzcat('ERROR: yGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).yGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).yGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).yGC)== 1))~= 0
                disp(horzcat('ERROR: yGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).zGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).zGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).zGL)== 1))~= 0
                disp(horzcat('ERROR: zGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).zGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).zGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).zGH)== 1))~= 0
                disp(horzcat('ERROR: zGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).zGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).zGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).zGC)== 1))~= 0
                disp(horzcat('ERROR: zGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).xGL)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).xGL)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).xGL)== 1))~= 0
                disp(horzcat('ERROR: xGL has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).xGH)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).xGH)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).xGH)== 1))~= 0
                disp(horzcat('ERROR: xGH has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).xGC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).xGC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).xGC)== 1))~= 0
                disp(horzcat('ERROR: xGC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).hqC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).hqC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).hqC)== 1))~= 0
                disp(horzcat('ERROR: hqC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).hpC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).hpC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).hpC)== 1))~= 0
                disp(horzcat('ERROR: hpC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).hphiC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).hphiC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).hphiC)== 1))~= 0
                disp(horzcat('ERROR: hphiC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).dqC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).dqC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).dqC)== 1))~= 0
                disp(horzcat('ERROR: dqC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).dpC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).dpC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).dpC)== 1))~= 0
                disp(horzcat('ERROR: dpC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).dphiC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).dphiC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).dphiC)== 1))~= 0
                disp(horzcat('ERROR: dphiC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).d3xC)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).d3xC)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).d3xC)== 1))~= 0
                disp(horzcat('ERROR: d3xC has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Q grid cell= ', num2str(Qind)))
            end
            
        end
    end
end

% -------------------------------------------------------

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        Specie(s).FluxTube(f).d3xCLB= Specie(s).FluxTube(f).QCell(1).d3xC; % LB Ghost cell volume
        Specie(s).FluxTube(f).d3xCUB= Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).d3xC; % UB Ghost cell volume
    end
end

% -------------------------------------------------------

% Compute low altitude reference density and altitude from PFISR data:
% load('/Users/robertalbarran/Dropbox/datasets/data_mat/07Feb2013_1.1min_rawdata.mat')
% launcht= 0; %launch 2/7/13 8:21 UTC
% ICRt= 591; %9.85min so 8:21+ 9.85= 30.85= 31
% landt= 901.8; %15.03min so 8:21+ 15.03= 36.03= 36
%
% beamind= 13;
%
% launchtind= 4*410+ 74;
% ICRtind= 4*410+ 83;
% landtind= 4*410+ 88;
% tofA= linspace(launcht, landt, 88-74);
% altA= squeeze(beamalt(1:17, 13))*1e-3;
% tofAICRind= 9;
%
% isTeA= squeeze(iste(:, beamind, 74:88));
% isTiA= squeeze(isti(:, beamind, 74:88));
% isneA= squeeze(isne(:, beamind, 74:88));
% isviA= squeeze(isvi(:, beamind, 74:88));
%
% Laltinds= find(altA < (Specie(s).FluxTube(f).QCell(1).rGL*1e-3- RE*1e-3));
% LaltAref= altA(Laltinds(end));
% LisneAref= isneA((tofAICRind), Laltinds(end));

% -------------------------------------------------------
horzcat('Number of Particle Species= ', num2str(Stot))
horzcat('Number of Flux Tubes= ', num2str(Specie(1).Nf))
horzcat('NqG= ', num2str(Specie(s).FluxTube(f).NqG- 2))
horzcat('L shell [RE]= ', num2str(Specie(1).FluxTube(1).pG))
horzcat('L shell drift limit [RE]= ', num2str(Specie(s).FluxTube(f).QCell(Qind).dpC))
horzcat('qGL= ', num2str(Specie(s).FluxTube(f).QCell(1).qGL))
horzcat('qGH= ', num2str(Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).qGH))
horzcat('rGL [km]= ', num2str((Specie(s).FluxTube(f).QCell(1).rGL- RE)*1e-3))
horzcat('rGH [km]= ', num2str((Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).rGH- RE)*1e-3))
horzcat('rGH [$R_E$]= ', num2str((Specie(s).FluxTube(f).QCell(Specie(s).FluxTube(f).NqG).rGH- RE)./RE))
% horzcat('Lower grid reference altitude [km]= ', num2str(LaltAref))
% horzcat('Lower grid reference density [m^-3]= ', num2str(LisneAref))

% zns0= RE+ 370.78e3; % For VISIONS-1 Case Study
% ns0= 6.678e10;
% Qindns0= 3e0;

zns0= (Specie(s).FluxTube(f).QCell(1).rGC- RE); % For VISIONS-1 Case Study
ns0= 6.678e10;
Qindns0= 1e0;

% zns0= RE+ 7656e3; % For (Barghouthi '97 and Crew '90) VS1 VS2
% ns0= 5e9;

% zns0= RE+ 10846e3; % For (Barghouthi '94) VS3
% ns0= 1e8;

nsnormfac= 1e17;

% Compute altitude dependent gravitational acceleration, integration over ds
for s= 1:1:Stot
    Specie(s).Qindns0= Qindns0;
    for f= 1:1:Specie(s).Nf;
        for Qind= 1:1:Specie(s).FluxTube(f).NqG;
            Specie(s).FluxTube(f).QCell(Qind).TsPerp= TsPerp; 
            Specie(s).FluxTube(f).QCell(Qind).TsPar= TsPar; 
            Specie(s).FluxTube(f).QCell(Qind).Te= Te; 
            Specie(s).FluxTube(f).TsPerp(Qind)= Specie(s).FluxTube(f).QCell(Qind).TsPerp;
            Specie(s).FluxTube(f).TsPar(Qind)= Specie(s).FluxTube(f).QCell(Qind).TsPar;
            Specie(s).FluxTube(f).Te(Qind)= Specie(s).FluxTube(f).QCell(Qind).Te;
                                    
            Specie(s).FluxTube(f).rGC(Qind)= Specie(s).FluxTube(f).QCell(Qind).rGC;
            Specie(s).FluxTube(f).thetaGC(Qind)= Specie(s).FluxTube(f).QCell(Qind).thetaGC;
            Specie(s).FluxTube(f).ellGC(Qind)= Specie(s).FluxTube(f).QCell(Qind).ellGC;
            Specie(s).FluxTube(f).hqC(Qind)= Specie(s).FluxTube(f).QCell(Qind).hqC;
            Specie(s).FluxTube(f).dqC(Qind)= Specie(s).FluxTube(f).QCell(Qind).dqC;
            Specie(s).FluxTube(f).hpC(Qind)= Specie(s).FluxTube(f).QCell(Qind).hpC;
            Specie(s).FluxTube(f).dpC(Qind)= Specie(s).FluxTube(f).QCell(Qind).dpC;
            Specie(s).FluxTube(f).dspC(Qind)= Specie(s).FluxTube(f).hpC(Qind)* ...
                Specie(s).FluxTube(f).dpC(Qind);
            Specie(s).FluxTube(f).hphiC(Qind)= Specie(s).FluxTube(f).QCell(Qind).hphiC;
            Specie(s).FluxTube(f).dphiC(Qind)= Specie(s).FluxTube(f).QCell(Qind).dphiC;
            Specie(s).FluxTube(f).dsphiC(Qind)= Specie(s).FluxTube(f).hphiC(Qind)* ...
                Specie(s).FluxTube(f).dphiC(Qind);
        end
        for Qind= 1:1:Specie(s).FluxTube(f).NqG;
            Ibbbp(Qind)= 2e0*GG*ME*((Specie(s).FluxTube(f).hqC(Qind)* ...
                Specie(s).FluxTube(f).dqC(Qind)* ...
                cos(Specie(s).FluxTube(f).thetaGC(Qind)))/ ...
                ((Specie(s).FluxTube(f).rGC(Qind)^2e0)* ...
                (sqrt(Specie(s).FluxTube(f).ellGC(Qind)))));
        end
        I0bbb= sum(Ibbbp(1:Qindns0)); % Assume summation does not include ghost cell
        for Qind= 1:1:Specie(s).FluxTube(f).NqG;
            % Altitude dependent gravitational acceleration, integration over ds
            Ibbb(Qind)= sum(Ibbbp(1:1:Qind)); % Assume summation does not include ghost cell
            gCbbb(Qind)= Ibbb(Qind)- I0bbb;

            if (electronflag == 1)
                argCbbb(Qind)= (Specie(s).ms*gCbbb(Qind))/((kB/1e0)*(Specie(s).FluxTube(f).QCell(Qind).TsPar+ Specie(s).FluxTube(f).QCell(Qind).Te));
            end
            if (electronflag == 0)
                argCbbb(Qind)= (Specie(s).ms*gCbbb(Qind))/((kB/1e0)*(Specie(s).FluxTube(f).QCell(Qind).TsPar));
            end
            
            Specie(s).FluxTube(f).HiCp(Qind)= abs(1e0/argCbbb(Qind)); % Compute ion scale height [m]
            Specie(s).FluxTube(f).dsC(Qind)= ...
                Specie(s).FluxTube(f).QCell(Qind).dqC*Specie(s).FluxTube(f).QCell(Qind).hqC;
            Specie(s).FluxTube(f).HiCratio(Qind)= ...
                (Specie(s).FluxTube(f).QCell(Qind).dqC*Specie(s).FluxTube(f).QCell(Qind).hqC)/ ...
                Specie(s).FluxTube(f).HiCp(Qind);
            
            nsCbbb(Qind)= ns0*exp(argCbbb(Qind)); % [m^-3]
            
            nsnormCbbb(Qind)= (nsCbbb(Qind)/nsnormfac)*(Specie(s).FluxTube(f).QCell(Qind).d3xC/360e0);

            Specie(s).FluxTube(f).d3xCp(Qind)= Specie(s).FluxTube(f).QCell(Qind).d3xC;
            Specie(s).FluxTube(f).rGCp(Qind)= Specie(s).FluxTube(f).QCell(Qind).rGC;
        end
        Specie(s).FluxTube(f).LBNominalDensity= nsCbbb(1);
        Specie(s).FluxTube(f).UBNominalDensity= nsCbbb(Specie(s).FluxTube(f).NqG);
    end
end

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        Specie(s).FluxTube(f).HiCratioMean= mean(Specie(s).FluxTube(f).HiCratio);
    end
end

horzcat('Mean ratio of FA grid length to ion scale height= ', num2str(Specie(s).FluxTube(f).HiCratioMean))

% fignum= 1; % Assign figure number
% fig(fignum)= figure(fignum);
% set(fig(fignum), 'Position', [10 1000 1200 500])
% 
% subplot(1, 2, 1)
% semilogx(Specie(s).FluxTube(f).HiCp.*1e-3, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'kx-', 'MarkerSize', 20, 'LineWidth', 2)
% hold on
% semilogx(Specie(s).FluxTube(f).dsC.*1e-3, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'mx-', 'MarkerSize', 20, 'LineWidth', 2)
% xlabel('[km]', 'interpreter', 'latex', 'FontSize', 25)
% ylabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25)
% LH(1)= plot(nan, nan, char('kx-'), 'MarkerSize', 20, 'LineWidth', 1);
% L{1}= horzcat('Plasma Scale Height');
% LH(2)= plot(nan, nan, char('mx-'), 'MarkerSize', 20, 'LineWidth', 1);
% L{2}= horzcat('Field-Aligned Grid Length');
% legend(LH, L, 'interpreter', 'latex', 'FontSize', 20, 'location', 'North')
% grid on
% set(gca, 'FontSize', 25)
% set(gcf, 'color', 'white');
% hold off
% 
% subplot(1, 2, 2)
% plot(Specie(s).FluxTube(f).HiCp.*1e-3, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'kx-', 'MarkerSize', 20, 'LineWidth', 2)
% hold on
% plot(Specie(s).FluxTube(f).dsC.*1e-3, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'mx-', 'MarkerSize', 20, 'LineWidth', 2)
% xlabel('[km]', 'interpreter', 'latex', 'FontSize', 25)
% ylabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25)
% LH(1)= plot(nan, nan, char('kx-'), 'MarkerSize', 20, 'LineWidth', 1);
% L{1}= horzcat('Plasma Scale Height');
% LH(2)= plot(nan, nan, char('mx-'), 'MarkerSize', 20, 'LineWidth', 1);
% L{2}= horzcat('Field-Aligned Grid Length');
% legend(LH, L, 'interpreter', 'latex', 'FontSize', 20, 'location', 'North')
% grid on
% set(gca, 'FontSize', 25)
% set(gcf, 'color', 'white');
% hold off
% 
% fignum= fignum+ 1e0; % Assign figure number
% fig(fignum)= figure(fignum);
% set(fig(fignum), 'Position', [10 1000 1200 500])
% 
% subplot(1, 3, 1)
% semilogx(nsCbbb, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'k', 'MarkerSize', 20, 'LineWidth', 2)
% xlabel('$n_C$ [m$^{-3}$]', 'interpreter', 'latex', 'FontSize', 25)
% ylabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25)
% grid on
% set(gca, 'FontSize', 25)
% set(gcf, 'color', 'white');
% hold off
% 
% subplot(1, 3, 2)
% plot(Specie(s).FluxTube(f).d3xCp/nsnormfac, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'k', 'MarkerSize', 20, 'LineWidth', 2)
% xlabel('$d^3x/\mu$ [m$^3$]', 'interpreter', 'latex', 'FontSize', 25)
% ylabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25)
% grid on
% set(gca, 'FontSize', 25)
% set(gcf, 'color', 'white');
% hold off
% 
% subplot(1, 3, 3)
% semilogx(nsnormCbbb, (Specie(s).FluxTube(f).rGCp- RE)./RE, 'k', 'MarkerSize', 20, 'LineWidth', 2)
% xlabel('$|n_C|$', 'interpreter', 'latex', 'FontSize', 25)
% ylabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25)
% grid on
% set(gca, 'FontSize', 25)
% set(gcf, 'color', 'white');
% hold off

disp('ION AND ENA CONFIG-SPACE GRID GENERATION COMPLETE')

%% -------------------------------------------------------

% SET PRELIMINARY NUMBER OF VELOCITY-SPACE GRID CELLS PER ION PARTICLE SPECIES AND
% FLUX TUBE:

clc
close all
tic

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                        
            if IONVPERPVECflag == 1
                
                % Total number of Vperp1 grid values:
                if s == 1 & f < Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp= NVperp1GpF;
                else if s == 1 & f >= Specie(s).Nf/2
                        Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp= NVperp1GpF;
                    end
                end
                if s == 2 & f < Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp= NVperp1GpF;
                else if s == 2 & f >= Specie(s).Nf/2
                        Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp= NVperp1GpF;
                    end
                end
                
                % Total number of Vperp2 grid values:
                if s == 1 & f < Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp= NVperp2GpF;
                else if s == 1 & f >= Specie(s).Nf/2
                        Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp= NVperp2GpF;
                    end
                end
                if s == 2 & f < Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp= NVperp2GpF;
                else if s == 2 & f >= Specie(s).Nf/2
                        Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp= NVperp2GpF;
                    end
                end
                
            else
                
               % Total number of Vperp grid values:
                if s == 1 & f < Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVperpGp= NVperpGpF;
                else if s == 1 & f >= Specie(s).Nf/2
                        Specie(s).FluxTube(f).QCell(Qind).NVperpGp= NVperpGpF;
                    end
                end
                if s == 2 & f < Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVperpGp= NVperpGpF;
                else if s == 2 & f >= Specie(s).Nf/2
                        Specie(s).FluxTube(f).QCell(Qind).NVperpGp= NVperpGpF;
                    end
                end
                
            end

            % Total number of Vpar grid values:
            if s == 1 & f < Specie(s).Nf/2
                Specie(s).FluxTube(f).QCell(Qind).NVparGp= NVparGpF;
            else if s == 1 & f >= Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVparGp= NVparGpF;
                end
            end
            if s == 2 & f < Specie(s).Nf/2
                Specie(s).FluxTube(f).QCell(Qind).NVparGp= NVparGpF;
            else if s == 2 & f >= Specie(s).Nf/2
                    Specie(s).FluxTube(f).QCell(Qind).NVparGp= NVparGpF;
                end
            end
            
        end
    end
end

% -------------------------------------------------------

% SET DIAGNOSTIC FLAGS FOR NaN, INFINITE, AND EMPTY ARRAY ELEMENTS:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf;
        for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
            
            if IONVPERPVECflag == 1
                
                if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp)== 1))~= 0 | ...
                        numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp)== 1))~= 0 | ...
                        numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp)== 1))~= 0
                    disp(horzcat('ERROR: NVperp1Gp has NaN, Inf, or empty array element for particle specie= ', ...
                        num2str(s), ', flux tube= ', num2str(f), ', and Qind= ', num2str(Qind)))
                end
                
                if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp)== 1))~= 0 | ...
                        numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp)== 1))~= 0 | ...
                        numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp)== 1))~= 0
                    disp(horzcat('ERROR: NVperp2Gp has NaN, Inf, or empty array element for particle specie= ', ...
                        num2str(s), ', flux tube= ', num2str(f), ', and Qind= ', num2str(Qind)))
                end
            
            else
                
                if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).NVperpGp)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).NVperpGp)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).NVperpGp)== 1))~= 0
                disp(horzcat('ERROR: NVperpGp has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Qind= ', num2str(Qind)))
                end
                
            end
            
            if numel(find(isnan(Specie(s).FluxTube(f).QCell(Qind).NVparGp)== 1))~= 0 | ...
                    numel(find(isinf(Specie(s).FluxTube(f).QCell(Qind).NVparGp)== 1))~= 0 | ...
                    numel(find(isempty(Specie(s).FluxTube(f).QCell(Qind).NVparGp)== 1))~= 0
                disp(horzcat('ERROR: NVparGp has NaN, Inf, or empty array element for particle specie= ', ...
                    num2str(s), ', flux tube= ', num2str(f), ', and Qind= ', num2str(Qind)))
            end
            
        end
    end
end

% -------------------------------------------------------

% CREATE PRELIMINARY FIELD-ALIGNED ION VELOCITY SPACE GRID:

for s= 1:1:Stot;
    for f= 1:1:Specie(s).Nf;
        for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
            
            % Set range of Eulerian Vperp1 coords.
            Specie(s).FluxTube(f).QCell(Qind).dVperp1= 1e0; 
            Specie(s).FluxTube(f).QCell(Qind).NVperp1G= (Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp/ ...
                Specie(s).FluxTube(f).QCell(Qind).dVperp1)- 1e0; 
            
            % Set range of Eulerian Vperp2 coords.
            Specie(s).FluxTube(f).QCell(Qind).dVperp2= 1e0; 
            Specie(s).FluxTube(f).QCell(Qind).NVperp2G= (Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp/ ...
                Specie(s).FluxTube(f).QCell(Qind).dVperp2)- 1e0; 
                        
            % Positive Vperp1 values
            Specie(s).FluxTube(f).QCell(Qind).Vperp1GA= 0e0; 
            Specie(s).FluxTube(f).QCell(Qind).Vperp1GB= Vperp12sigmaFac*Vperp12sigma;
            Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp1= ...
                linspace(Specie(s).FluxTube(f).QCell(Qind).Vperp1GA, ...
                Specie(s).FluxTube(f).QCell(Qind).Vperp1GB, Vperp12NlinRange);
                        
            % Negative Vperp1 values
            Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp2= ...
                flip(-Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp1);
            
            % Put negative and positive Vperp1 values together
            Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp= ...
                [Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp2(1:end) ...
                Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp1(2:end)];
            
            Specie(s).FluxTube(f).QCell(Qind).Vperp2Gp= ...
                Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp;
            
            Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp= Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp- 1;
            Specie(s).FluxTube(f).QCell(Qind).NVperp1G= Specie(s).FluxTube(f).QCell(Qind).NVperp1G- 1;
            
            Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp= Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp- 1;
            Specie(s).FluxTube(f).QCell(Qind).NVperp2G= Specie(s).FluxTube(f).QCell(Qind).NVperp2G- 1;
                        
            for Vperp1ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp1G;
                Specie(s).FluxTube(f).QCell(Qind).dVperp1G(Vperp1ind)= ...
                    Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp(Vperp1ind+ 1)- ...
                    Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp(Vperp1ind);
            end
            for Vperp2ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp2G;
                Specie(s).FluxTube(f).QCell(Qind).dVperp2G(Vperp2ind)= ...
                    Specie(s).FluxTube(f).QCell(Qind).Vperp2Gp(Vperp2ind+ 1)- ...
                    Specie(s).FluxTube(f).QCell(Qind).Vperp2Gp(Vperp2ind);
            end
            
            % Set range of Eulerian Vpar coords.
            Specie(s).FluxTube(f).QCell(Qind).dVpar= 1e0; 
            Specie(s).FluxTube(f).QCell(Qind).NVparG= (Specie(s).FluxTube(f).QCell(Qind).NVparGp/ ...
                Specie(s).FluxTube(f).QCell(Qind).dVpar)- 1e0; 
            
            if (SymVparGridflag == 1)
                % Positive Vpar values
                Specie(s).FluxTube(f).QCell(Qind).VparGA= 0e0; 
                Specie(s).FluxTube(f).QCell(Qind).VparGB= VparsigmaFac*Vparsigma;
                Specie(s).FluxTube(f).QCell(Qind).VparGp1= ...
                    linspace(Specie(s).FluxTube(f).QCell(Qind).VparGA, ...
                    Specie(s).FluxTube(f).QCell(Qind).VparGB, VparNlinRange);
                           
                % Negative Vpar values
                Specie(s).FluxTube(f).QCell(Qind).VparGp2= ...
                    flip(-Specie(s).FluxTube(f).QCell(Qind).VparGp1);
            
                % Put negative and positive Vpar values together
                Specie(s).FluxTube(f).QCell(Qind).VparGp= ...
                    [Specie(s).FluxTube(f).QCell(Qind).VparGp2(1:end) ...
                    Specie(s).FluxTube(f).QCell(Qind).VparGp1(2:end)];
            
                Specie(s).FluxTube(f).QCell(Qind).NVparGp= Specie(s).FluxTube(f).QCell(Qind).NVparGp- 1;
                Specie(s).FluxTube(f).QCell(Qind).NVparG= Specie(s).FluxTube(f).QCell(Qind).NVparG- 1;
            
            end
            if (SymVparGridflag == 0)
                % Positive Vpar values
                Specie(s).FluxTube(f).QCell(Qind).VparGA= -4e0*Vparsigma; 
                Specie(s).FluxTube(f).QCell(Qind).VparGB= VparsigmaFac*Vparsigma;
                Specie(s).FluxTube(f).QCell(Qind).VparGp1= ...
                    linspace(Specie(s).FluxTube(f).QCell(Qind).VparGA, ...
                    Specie(s).FluxTube(f).QCell(Qind).VparGB, VparNlinRange);

                Specie(s).FluxTube(f).QCell(Qind).NVparGp= Specie(s).FluxTube(f).QCell(Qind).NVparGp- 1;
                Specie(s).FluxTube(f).QCell(Qind).NVparG= Specie(s).FluxTube(f).QCell(Qind).NVparG- 1;
                if (SMagHemFlag == 0)  
                    Specie(s).FluxTube(f).QCell(Qind).VparGp= ...
                        Specie(s).FluxTube(f).QCell(Qind).VparGp1;
                end
                if (SMagHemFlag == 1)  
                    Specie(s).FluxTube(f).QCell(Qind).VparGp= ...
                        flip(-Specie(s).FluxTube(f).QCell(Qind).VparGp1);
                end
            end
                        
            for Vparind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVparG;
                Specie(s).FluxTube(f).QCell(Qind).dVparG(Vparind)= ...
                    Specie(s).FluxTube(f).QCell(Qind).VparGp(Vparind+ 1)- ...
                    Specie(s).FluxTube(f).QCell(Qind).VparGp(Vparind);
            end
                                        
            for Vperp1ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp;
                for Vperp2ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp;
                    for Vparind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVparGp;

                        Specie(s).FluxTube(f).QCell(Qind).Vperp1G(Vperp1ind, Vperp2ind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).Vperp1Gp(Vperp1ind);

                        Specie(s).FluxTube(f).QCell(Qind).Vperp2G(Vperp1ind, Vperp2ind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).Vperp2Gp(Vperp2ind);

                        Specie(s).FluxTube(f).QCell(Qind).VparG(Vperp1ind, Vperp2ind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VparGp(Vparind);

                    end
                end
            end
                            
        end
    end
end

% -------------------------------------------------------

% CREATE FINAL FIELD-ALIGNED VELOCITY SPACE GRID:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
            
            if IONVPERPVECflag == 1
                
                for Vperp1ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp;
                    for Vperp2ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp;
                        for Vparind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVparGp;

                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1G(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).Vperp1G(Vperp1ind, Vperp2ind, Vparind);
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2G(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).Vperp2G(Vperp1ind, Vperp2ind, Vparind);
                            Specie(s).FluxTube(f).QCell(Qind).V2Perp2Cell(Vperp1ind, Vperp2ind, Vparind).VparG(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VparG(Vperp1ind, Vperp2ind, Vparind);
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparG(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VparG(Vperp1ind, Vperp2ind, Vparind);

                        end
                    end
                end
            
            else
                
                for Vperpind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperpGp;
                    for Vparind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVparGp;

                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpG(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VperpG(Vperpind, Vparind);
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparG(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VparG(Vperpind, Vparind);

                    end
                end
                
            end
            
        end
    end
end

% ---------------------------------------------

% GET VELOCITY SPACE GRID BOUNDARIES, CENTER VALUES, METRIC FACTORS, AND VOLUMES:

for s= 1:1:Stot
    for f= 1:1:Specie(s).Nf
        for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
            
            if IONVPERPVECflag == 1
                
                for Vperp1ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp1G;
                    for Vperp2ind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperp2G;
                        for Vparind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVparG;

                            % Lower Vperp1 limits of Vperp1 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GL(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell((Vperp1ind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVperp1+ 1, Vperp2ind, 1).Vperp1G(1);
                            % Lower Vperp2 limits of Vperp1 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GL(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, (Vperp2ind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVperp2+ 1, 1).Vperp2G(1);
                            % Upper Vperp1 limits of Vperp1 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GH(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell((Vperp1ind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVperp1+ 1+ ...
                                Specie(s).FluxTube(f).QCell(Qind).dVperp1, Vperp2ind, 1).Vperp1G(1);
                            % Upper Vperp2 limits of Vperp2 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GH(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, (Vperp2ind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVperp2+ 1+ ...
                                Specie(s).FluxTube(f).QCell(Qind).dVperp2, 1).Vperp2G(1);
                            % Center Vperp1 values of Vperp1 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GC(1)= ...
                                (Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GL(1)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GH(1))/2e0;
                            Specie(s).FluxTube(f).QCell(Qind).Vperp1GCp(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GC(1);
                            % Center Vperp2 values of Vperp2 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GC(1)= ...
                                (Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GL(1)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GH(1))/2e0;
                            Specie(s).FluxTube(f).QCell(Qind).Vperp2GCp(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GC(1);
           
                            Specie(s).FluxTube(f).QCell(Qind).Vperp1GL(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GL(1);
                            Specie(s).FluxTube(f).QCell(Qind).Vperp2GL(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GL(1);
                            Specie(s).FluxTube(f).QCell(Qind).Vperp1GH(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GH(1);
                            Specie(s).FluxTube(f).QCell(Qind).Vperp2GH(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GH(1);
                            Specie(s).FluxTube(f).QCell(Qind).Vperp1GC(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GC(1);
                            Specie(s).FluxTube(f).QCell(Qind).Vperp2GC(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GC(1);
                            
                            % Lower Vpar limits of Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGL(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(1, 1, (Vparind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVpar+ 1).VparG(1);
                            % Upper Vpar limits of Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGH(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(1, 1, (Vparind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVpar+ ...
                                1+ Specie(s).FluxTube(f).QCell(Qind).dVpar).VparG(1);
                            % Center Vpar values of Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGC(1)= ...
                                (Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGL(1)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGH(1))/2e0;
                            Specie(s).FluxTube(f).QCell(Qind).VparGCp(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGC(1);
    
                            Specie(s).FluxTube(f).QCell(Qind).VparGL(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGL(1);
                            Specie(s).FluxTube(f).QCell(Qind).VparGH(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGH(1);
                            Specie(s).FluxTube(f).QCell(Qind).VparGC(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGC(1);
                            
                            % Vel-space metric factors (local FA Cartesian coords., where Vperp1= Vy= Vp, Vperp2= Vx= Vphid)
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVperp1C(1)= 1e0;
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVperp2C(1)= 1e0;
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVparC(1)= 1e0;

                            Specie(s).FluxTube(f).QCell(Qind).hVperp1Cp(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVperp1C(1);
                            Specie(s).FluxTube(f).QCell(Qind).hVperp2Cp(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVperp2C(1);
                            Specie(s).FluxTube(f).QCell(Qind).hVparCp(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVparC(1);

                            % dVperp1 across Vperp1 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVperp1C(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GH(1)- ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp1GL(1));
                            Specie(s).FluxTube(f).QCell(Qind).dVperp1C(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVperp1C(1);
                            % dVperp2 across Vperp2 cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVperp2C(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GH(1)- ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).Vperp2GL(1));
                            Specie(s).FluxTube(f).QCell(Qind).dVperp2C(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVperp2C(1);
                            % dVpar across Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVparC(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGH(1)- ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).VparGL(1));
                            Specie(s).FluxTube(f).QCell(Qind).dVparC(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVparC(1);

                            Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).d3vC(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVperp1C(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVperp2C(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).hVparC(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVperp1C(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVperp2C(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).dVparC(1);
                            
                            if Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).d3vC(1) == 0e0
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).d3vC(1)= ...
                                    1e-9;
                            end

                            Specie(s).FluxTube(f).QCell(Qind).d3vC(Vperp1ind, Vperp2ind, Vparind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).d3vC(1);

                            if Specie(s).FluxTube(f).QCell(Qind).V2PerpCell(Vperp1ind, Vperp2ind, Vparind).d3vC(1) < 0e0
                                disp(horzcat('ERROR: Negative d3vC for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vperp1ind= ', num2str(Vperp1ind), ...
                                    ', Vperp2ind= ', num2str(Vperp2ind), ', and Vparind= ', num2str(Vparind)))
                            end
                        end
                    end
                end
            
            else
                
                for Vperpind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVperpG;
                    for Vparind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVparG;

                        % Lower Vperp limits of Vperp cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGL(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell((Vperpind- 1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).dVperp+ 1, 1).VperpG(1);
                        % Upper Vperp limits of Vperp cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGH(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell((Vperpind- 1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).dVperp+ 1+ ...
                            Specie(s).FluxTube(f).QCell(Qind).dVperp, 1).VperpG(1);
                        % Center Vperp values of Vperp cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGC(1)= ...
                            (Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGL(1)+ ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGH(1))/2e0;
                        Specie(s).FluxTube(f).QCell(Qind).VperpGCp(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGC(1);
                        
                        Specie(s).FluxTube(f).QCell(Qind).VperpGL(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGL(1);
                        Specie(s).FluxTube(f).QCell(Qind).VperpGH(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGH(1);
                        Specie(s).FluxTube(f).QCell(Qind).VperpGC(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGC(1);
    
                        % Lower Vpar limits of Vpar cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGL(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(1, (Vparind- 1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).dVpar+ 1).VparG(1);
                        % Upper Vpar limits of Vpar cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGH(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(1, (Vparind- 1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).dVpar+ ...
                            1+ Specie(s).FluxTube(f).QCell(Qind).dVpar).VparG(1);
                        % Center Vpar values of Vpar cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGC(1)= ...
                            (Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGL(1)+ ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGH(1))/2e0;
                        Specie(s).FluxTube(f).QCell(Qind).VparGCp(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGC(1);
    
                        Specie(s).FluxTube(f).QCell(Qind).VparGL(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGL(1);
                        Specie(s).FluxTube(f).QCell(Qind).VparGH(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGH(1);
                        Specie(s).FluxTube(f).QCell(Qind).VparGC(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGC(1);
                        
                        % Vel-space metric factors (local FA cylindrical coords., where Vperp1= Vy= Vp, Vperp2= Vx= Vphid)
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVperpC(1)= 1e0;
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVparC(1)= 1e0;
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVthetaC(1)= ...
                            abs(Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGC(1));

                        Specie(s).FluxTube(f).QCell(Qind).hVperpCp(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVperpC(1);

                        Specie(s).FluxTube(f).QCell(Qind).hVparCp(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVparC(1);

                        Specie(s).FluxTube(f).QCell(Qind).hVthetaCp(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVthetaC(1);

                        % dVperp across Vperp cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVperpC(1)= ...
                            abs(Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGH(1)- ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGL(1));
                        Specie(s).FluxTube(f).QCell(Qind).dVperpC(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVperpC(1);
                        % dVpar across Vpar cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVparC(1)= ...
                            abs(Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGH(1)- ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGL(1));
                        Specie(s).FluxTube(f).QCell(Qind).dVparC(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVparC(1);
                        % dVtheta across Vtheta cells
                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVthetaC(1)= 2e0*pi;

                        Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC(1)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVperpC(1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVparC(1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).hVthetaC(1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVperpC(1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVparC(1)* ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).dVthetaC(1);

                        if Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC(1) == 0e0
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC(1)= ...
                                1e-9;
                        end

                        Specie(s).FluxTube(f).QCell(Qind).d3vC(Vperpind, Vparind)= ...
                            Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC(1);

                        if Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC(1) < 0e0
                            disp(horzcat('ERROR: Negative d3vC for particle species= ', ...
                                num2str(s), ', flux tube= ', num2str(f), ...
                                ', Qind= ', num2str(Qind), ', Vperpind= ', num2str(Vperpind), ...
                                ', and Vparind= ', num2str(Vparind)))
                        end
                    end
                end
            
            end
            
        end
    end
end

if IONVPERPVECflag == 1
    
    horzcat('NVperp1G= ', num2str(Specie(1).FluxTube(1).QCell(1).NVperp1G))
    horzcat('NVperp2G= ', num2str(Specie(1).FluxTube(1).QCell(1).NVperp2G))
    horzcat('NVparG= ', num2str(Specie(1).FluxTube(1).QCell(1).NVparG))
    horzcat('Vperp1GL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V2PerpCell(1, 1, 1).Vperp1GL(1)*1e-3))
    horzcat('Vperp1GH [km/s]= ', ...
        num2str(Specie(1).FluxTube(1).QCell(1).V2PerpCell(Specie(1).FluxTube(1).QCell(1).NVperp1G, 1, 1).Vperp1GH(1)*1e-3))
    horzcat('Vperp2GL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V2PerpCell(1, 1, 1).Vperp2GL(1)*1e-3))
    horzcat('Vperp2GH [km/s]= ', ...
        num2str(Specie(1).FluxTube(1).QCell(1).V2PerpCell(1, Specie(1).FluxTube(1).QCell(1).NVperp2G, 1).Vperp2GH(1)*1e-3))
    horzcat('VparGL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V2PerpCell(1, 1, 1).VparGL(1)*1e-3))
    horzcat('VparGH [km/s]= ', ...
        num2str(Specie(1).FluxTube(1).QCell(1).V2PerpCell(1, 1, Specie(1).FluxTube(1).QCell(1).NVparG).VparGH(1)*1e-3))
    
    Vperp12Gridlinspace= (Specie(s).FluxTube(f).QCell(Qind).Vperp1G(end, 1, 1)- ...
        Specie(s).FluxTube(f).QCell(Qind).Vperp1G(end- 1, 1, 1));
    VparGridlinspace= (Specie(s).FluxTube(f).QCell(Qind).VparG(1, 1, end)- ...
        Specie(s).FluxTube(f).QCell(Qind).VparG(1, 1, end- 1));
    
    horzcat('Vperp12 linear grid to resolve MB thermal core with resolution [km/s]= ', ...
        num2str(Vperp12Gridlinspace*1e-3))
    horzcat('Vpar linear grid to resolve MB thermal core with resolution [km/s]= ', ...
        num2str(VparGridlinspace*1e-3))
    horzcat('Vperp12 sigma [km/s]= ', num2str(Vperp12sigma*1e-3), ...
        ', Vperp12 sigma factor= ', num2str(Vperp12sigmaFac))
    horzcat('Vpar sigma [km/s]= ', num2str(Vparsigma*1e-3), ...
        ', Vpar sigma factor= ', num2str(VparsigmaFac))
    
else
    
    horzcat('NVperpG= ', num2str(Specie(1).FluxTube(1).QCell(1).NVperpG))
    horzcat('NVparG= ', num2str(Specie(1).FluxTube(1).QCell(1).NVparG))
    horzcat('VperpGL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).VCell(1, 1).VperpGL(1)*1e-3))
    horzcat('VperpGH [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).VCell(Specie(1).FluxTube(1).QCell(1).NVperpG, ...
        Specie(1).FluxTube(1).QCell(1).NVparG).VperpGH(1)*1e-3))
    horzcat('VparGL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).VCell(1, 1).VparGL(1)*1e-3))
    horzcat('VparGH [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).VCell(Specie(1).FluxTube(1).QCell(1).NVperpG, ...
        Specie(1).FluxTube(1).QCell(1).NVparG).VparGH(1)*1e-3))

end

toc

disp('ION PHASE-SPACE GRID GENERATION COMPLETE')

%% -------------------------------------------------------

% SET PRELIMINARY NUMBER OF VELOCITY-SPACE GRID CELLS PER ENA PARTICLE SPECIES AND
% FLUX TUBE:

clc
close all
tic

if ENAEXPORTflag == 1
    
    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf
            for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                
                % Total number of Vp, Vq, Vphi grid values (Make odd number for 2D distrib fnc
                % moment integration):
                
                Specie(s).FluxTube(f).QCell(Qind).NVpGp= NVpGpF;
                Specie(s).FluxTube(f).QCell(Qind).NVqGp= NVqGpF;
                Specie(s).FluxTube(f).QCell(Qind).NVphiGp= NVphiGpF;
                
            end
        end
    end
    
    % -------------------------------------------------------
    
    % CREATE PRELIMINARY DIPOLE ENA VELOCITY SPACE GRID:
    
    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf;
            for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                
                % Set range of Eulerian Vp coords.
                Specie(s).FluxTube(f).QCell(Qind).dVp= 1e0; 
                Specie(s).FluxTube(f).QCell(Qind).NVpG= (Specie(s).FluxTube(f).QCell(Qind).NVpGp/ ...
                    Specie(s).FluxTube(f).QCell(Qind).dVp)- 1e0; 
                
                % Set range of Eulerian Vphi coords.
                Specie(s).FluxTube(f).QCell(Qind).dVphi= 1e0; 
                Specie(s).FluxTube(f).QCell(Qind).NVphiG= (Specie(s).FluxTube(f).QCell(Qind).NVphiGp/ ...
                    Specie(s).FluxTube(f).QCell(Qind).dVphi)- 1e0; 
    
                % Positive Vp values
                Specie(s).FluxTube(f).QCell(Qind).VpGA= 0e0; 
                Specie(s).FluxTube(f).QCell(Qind).VpGB= VpphisigmaFac*Vpphisigma;
                Specie(s).FluxTube(f).QCell(Qind).VpGp1= ...
                    linspace(Specie(s).FluxTube(f).QCell(Qind).VpGA, ...
                    Specie(s).FluxTube(f).QCell(Qind).VpGB, VpphiNlinRange);

                % Negative Vp values
                Specie(s).FluxTube(f).QCell(Qind).VpGp2= ...
                    flip(-Specie(s).FluxTube(f).QCell(Qind).VpGp1);

                % Put negative and positive Vp values together
                Specie(s).FluxTube(f).QCell(Qind).VpGp= ...
                    [Specie(s).FluxTube(f).QCell(Qind).VpGp2(1:end- 1) ...
                    Specie(s).FluxTube(f).QCell(Qind).VpGp1(2:end)];
    
                Specie(s).FluxTube(f).QCell(Qind).VphiGp= ...
                    Specie(s).FluxTube(f).QCell(Qind).VpGp;
                
                Specie(s).FluxTube(f).QCell(Qind).NVpGp= Specie(s).FluxTube(f).QCell(Qind).NVpGp- 2;
                Specie(s).FluxTube(f).QCell(Qind).NVpG= Specie(s).FluxTube(f).QCell(Qind).NVpG- 2;
                
                Specie(s).FluxTube(f).QCell(Qind).NVphiGp= Specie(s).FluxTube(f).QCell(Qind).NVphiGp- 2;
                Specie(s).FluxTube(f).QCell(Qind).NVphiG= Specie(s).FluxTube(f).QCell(Qind).NVphiG- 2;

                for Vpind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVpG;
                    Specie(s).FluxTube(f).QCell(Qind).dVpG(Vpind)= ...
                        Specie(s).FluxTube(f).QCell(Qind).VpGp(Vpind+ 1)- ...
                        Specie(s).FluxTube(f).QCell(Qind).VpGp(Vpind);
                end
                
                for Vphiind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVphiG;
                    Specie(s).FluxTube(f).QCell(Qind).dVphiG(Vphiind)= ...
                        Specie(s).FluxTube(f).QCell(Qind).VphiGp(Vphiind+ 1)- ...
                        Specie(s).FluxTube(f).QCell(Qind).VphiGp(Vphiind);
                end
                
                % Set range of Eulerian Vq coords.
                Specie(s).FluxTube(f).QCell(Qind).dVq= 1e0; 
                Specie(s).FluxTube(f).QCell(Qind).NVqG= (Specie(s).FluxTube(f).QCell(Qind).NVqGp/ ...
                    Specie(s).FluxTube(f).QCell(Qind).dVq)- 1e0; 
    
                
                
                if (SymVparGridflag == 1)
                    % Positive Vq values
                    Specie(s).FluxTube(f).QCell(Qind).VqGA= 0e0; 
                    Specie(s).FluxTube(f).QCell(Qind).VqGB= VqsigmaFac*Vqsigma;
                    Specie(s).FluxTube(f).QCell(Qind).VqGp1= ...
                        linspace(Specie(s).FluxTube(f).QCell(Qind).VqGA, ...
                        Specie(s).FluxTube(f).QCell(Qind).VqGB, VqNlinRange);

                    % Negative Vq values
                    Specie(s).FluxTube(f).QCell(Qind).VqGp2= ...
                        flip(-Specie(s).FluxTube(f).QCell(Qind).VqGp1);

                    % Put negative and positive Vq values together
                    Specie(s).FluxTube(f).QCell(Qind).VqGp= ...
                        [Specie(s).FluxTube(f).QCell(Qind).VqGp2(1:end- 1) ...
                        Specie(s).FluxTube(f).QCell(Qind).VqGp1(2:end)];

                    Specie(s).FluxTube(f).QCell(Qind).NVqGp= Specie(s).FluxTube(f).QCell(Qind).NVqGp- 2;
                    Specie(s).FluxTube(f).QCell(Qind).NVqG= Specie(s).FluxTube(f).QCell(Qind).NVqG- 2;
            
                end
                if (SymVparGridflag == 0)
                    % Positive Vq values
                    Specie(s).FluxTube(f).QCell(Qind).VqGA= -4e0*Vqsigma; 
                    Specie(s).FluxTube(f).QCell(Qind).VqGB= VqsigmaFac*Vqsigma;
                    Specie(s).FluxTube(f).QCell(Qind).VqGp1= ...
                        linspace(Specie(s).FluxTube(f).QCell(Qind).VqGA, ...
                        Specie(s).FluxTube(f).QCell(Qind).VqGB, VqNlinRange);

                    Specie(s).FluxTube(f).QCell(Qind).NVqGp= Specie(s).FluxTube(f).QCell(Qind).NVqGp- 1;
                    Specie(s).FluxTube(f).QCell(Qind).NVqG= Specie(s).FluxTube(f).QCell(Qind).NVqG- 1;
                    if (SMagHemFlag == 0)  
                        Specie(s).FluxTube(f).QCell(Qind).VqGp= ...
                            Specie(s).FluxTube(f).QCell(Qind).VqGp1;
                    end
                    if (SMagHemFlag == 1)  
                        Specie(s).FluxTube(f).QCell(Qind).VqGp= ...
                            flip(-Specie(s).FluxTube(f).QCell(Qind).VqGp1);
                    end
                end
            
                for Vqind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVqG;
                    Specie(s).FluxTube(f).QCell(Qind).dVqG(Vqind)= ...
                        Specie(s).FluxTube(f).QCell(Qind).VqGp(Vqind+ 1)- ...
                        Specie(s).FluxTube(f).QCell(Qind).VqGp(Vqind);
                end
                            
                for Vpind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVpGp;
                    for Vqind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVqGp;
                        for Vphiind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVphiGp;
                            
                            Specie(s).FluxTube(f).QCell(Qind).VpG(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VpGp(Vpind);
                            Specie(s).FluxTube(f).QCell(Qind).VqG(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VqGp(Vqind);
                            Specie(s).FluxTube(f).QCell(Qind).VphiG(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VphiGp(Vphiind);
                            
                        end
                    end
                end
            end
        end
    end
    
    % -------------------------------------------------------
    
    % CREATE FINAL DIPOLE ENA VELOCITY SPACE GRID:
    
    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf
            for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                for Vpind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVpGp;
                    for Vqind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVqGp;
                        for Vphiind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVphiGp;
                            
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpG(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VpG(Vpind, Vqind, Vphiind);
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqG(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VqG(Vpind, Vqind, Vphiind);
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiG(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).VphiG(Vpind, Vqind, Vphiind);
                            
                        end
                    end
                end
            end
        end
    end
    
    % ---------------------------------------------
    
    % GET ENA VELOCITY SPACE GRID BOUNDARIES, CENTER VALUES, METRIC FACTORS, AND VOLUMES:
    
    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf
            for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                for Vpind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVpG;
                    for Vqind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVqG;
                        for Vphiind= 1:1:Specie(s).FluxTube(f).QCell(Qind).NVphiG;
                            
                            % Lower Vp limits of Vp cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGL(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell((Vpind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVp+ 1, 1, 1).VpG(1);
                            % Upper Vp limits of Vp cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGH(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell((Vpind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVp+ 1+ ...
                                Specie(s).FluxTube(f).QCell(Qind).dVp, 1, 1).VpG(1);
                            % Center Vp values of Vp cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC(1)= ...
                                (Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGL(1)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGH(1))/2e0;
                            Specie(s).FluxTube(f).QCell(Qind).VpGCp(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC(1);
                            
                            Specie(s).FluxTube(f).QCell(Qind).VpGL(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGL(1);
                            Specie(s).FluxTube(f).QCell(Qind).VpGH(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGH(1);
                            Specie(s).FluxTube(f).QCell(Qind).VpGC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC(1);
                            
                            % Lower Vq limits of Vq cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGL(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(1, (Vqind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVq+ 1, 1).VqG(1);
                            % Upper Vq limits of Vq cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGH(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(1, (Vqind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVq+ ...
                                1+ Specie(s).FluxTube(f).QCell(Qind).dVq, 1).VqG(1);
                            % Center Vq values of Vq cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC(1)= ...
                                (Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGL(1)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGH(1))/2e0;
                            Specie(s).FluxTube(f).QCell(Qind).VqGCp(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC(1);
                            
                            Specie(s).FluxTube(f).QCell(Qind).VqGL(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGL(1);
                            Specie(s).FluxTube(f).QCell(Qind).VqGH(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGH(1);
                            Specie(s).FluxTube(f).QCell(Qind).VqGC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC(1);
                            
                            % Lower Vphi limits of Vphi cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGL(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(1, 1, (Vphiind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVphi+ 1).VphiG(1);
                            % Upper Vphi limits of Vphi cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGH(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(1, 1, (Vphiind- 1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).dVphi+ ...
                                1+ Specie(s).FluxTube(f).QCell(Qind).dVphi).VphiG(1);
                            % Center Vphi values of Vphi cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGC(1)= ...
                                (Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGL(1)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGH(1))/2e0;
                            Specie(s).FluxTube(f).QCell(Qind).VphiGCp(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGC(1);
                            
                            Specie(s).FluxTube(f).QCell(Qind).VphiGL(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGL(1);
                            Specie(s).FluxTube(f).QCell(Qind).VphiGH(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGH(1);
                            Specie(s).FluxTube(f).QCell(Qind).VphiGC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGC(1);
                            
                            % Spherical velocity components
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VrGC(1)= ...
                                2e0*Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC* ...
                                cos(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC* ...
                                sin(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC);
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VthetaGC(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC* ...
                                sin(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC)- ...
                                2e0*Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC* ...
                                cos(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC);
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VellGC(1)= ...
                                (1e0+ 3e0*((cos(Specie(s).FluxTube(f).QCell(Qind).thetaGC))^2e0));
                            
                            Specie(s).FluxTube(f).QCell(Qind).VrGCp(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VrGC(1);
                            
                            Specie(s).FluxTube(f).QCell(Qind).VthetaGCp(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VthetaGC(1);
                            
                            Specie(s).FluxTube(f).QCell(Qind).VellGCp(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VellGC(1);
                            
                            % Test dipole velocity components
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGCtest(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VrGC* ...
                                sin(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/ ...
                                sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC)- ...
                                2e0*Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VthetaGC* ...
                                cos(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/ ...
                                sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC);
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGCtest(1)= ...
                                2e0*Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VrGC* ...
                                cos(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/ ...
                                sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC)+ ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VthetaGC* ...
                                sin(Specie(s).FluxTube(f).QCell(Qind).thetaGC)/ ...
                                sqrt(Specie(s).FluxTube(f).QCell(Qind).ellGC);
                            
                            if abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC(1)- ...
                                    Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGCtest(1)) > 1e-6
                                disp('BAD VP')
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC(1)- ...
                                    Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGCtest(1))
                            end
                            
                            if abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC(1)- ...
                                    Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGCtest(1)) > 1e-6
                                disp('BAD VQ')
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC(1)- ...
                                    Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGCtest(1))
                            end
                            
                            % ENA Vel-space metric factors
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVpC(1)= ...
                                abs(RE*sin(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VthetaGC)^3e0/ ...
                                (sqrt(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VellGC)));
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVqC(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VrGC^3e0/ ...
                                (RE^2e0*sqrt(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VellGC)));
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVphiC(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VrGC* ...
                                sin(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VthetaGC));
                            
                            Specie(s).FluxTube(f).QCell(Qind).hVpC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVpC(1);
                            Specie(s).FluxTube(f).QCell(Qind).hVqC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVqC(1);
                            Specie(s).FluxTube(f).QCell(Qind).hVphiC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVphiC(1);
                            
                            % dVp across Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVpC(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGH(1)- ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGL(1));
                            % dVq across Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVqC(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGH(1)- ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGL(1));
                            % dVphi across Vpar cells
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVphiC(1)= ...
                                abs(Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGH(1)- ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGL(1));
                            
                            Specie(s).FluxTube(f).QCell(Qind).dVpC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVpC(1);
                            Specie(s).FluxTube(f).QCell(Qind).dVqC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVqC(1);
                            Specie(s).FluxTube(f).QCell(Qind).dVphiC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVphiC(1);
                            
                            Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).d33vC(1)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVpC(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVqC(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).hVphiC(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVpC(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVqC(1)* ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).dVphiC(1);
                            
                            if Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).d33vC(1) == 0e0
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).d33vC(1)= ...
                                    1e-9;
                            end
                            
                            Specie(s).FluxTube(f).QCell(Qind).d33vC(Vpind, Vqind, Vphiind)= ...
                                Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).d33vC(1);
                            
                            if Specie(s).FluxTube(f).QCell(Qind).hVpC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative hVpCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                            if Specie(s).FluxTube(f).QCell(Qind).hVqC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative hVqCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                            if Specie(s).FluxTube(f).QCell(Qind).hVphiC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative hVphiCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                            if Specie(s).FluxTube(f).QCell(Qind).dVpC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative dVpCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                            if Specie(s).FluxTube(f).QCell(Qind).dVqC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative dVqCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                            if Specie(s).FluxTube(f).QCell(Qind).dVphiC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative dVphiCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                            if Specie(s).FluxTube(f).QCell(Qind).d33vC(Vpind, Vqind, Vphiind) < 0e0
                                disp(horzcat('ERROR: Negative d33vCp for particle species= ', ...
                                    num2str(s), ', flux tube= ', num2str(f), ...
                                    ', Qind= ', num2str(Qind), ', Vpind= ', num2str(Vpind), ', Vqind= ', num2str(Vqind), ...
                                    ', Vphiind= ', num2str(Vphiind)))
                            end
                            
                        end
                    end
                end
            end
        end
    end
        
    horzcat('NVpG= ', num2str(Specie(1).FluxTube(1).QCell(1).NVpG))
    horzcat('NVqG= ', num2str(Specie(1).FluxTube(1).QCell(1).NVqG))
    horzcat('NVphiG= ', num2str(Specie(1).FluxTube(1).QCell(1).NVphiG))
    horzcat('VpGL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V3Cell(1, 1, 1).VpGL(1)*1e-3))
    horzcat('VpGH [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V3Cell(Specie(1).FluxTube(1).QCell(1).NVpG, ...
        Specie(1).FluxTube(1).QCell(1).NVqG, Specie(1).FluxTube(1).QCell(1).NVphiG).VpGH(1)*1e-3))
    horzcat('VqGL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V3Cell(1, 1, 1).VqGL(1)*1e-3))
    horzcat('VqGH [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V3Cell(Specie(1).FluxTube(1).QCell(1).NVpG, ...
        Specie(1).FluxTube(1).QCell(1).NVqG, Specie(1).FluxTube(1).QCell(1).NVphiG).VqGH(1)*1e-3))
    horzcat('VphiGL [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V3Cell(1, 1, 1).VphiGL(1)*1e-3))
    horzcat('VphiGH [km/s]= ', num2str(Specie(1).FluxTube(1).QCell(1).V3Cell(Specie(1).FluxTube(1).QCell(1).NVpG, ...
        Specie(1).FluxTube(1).QCell(1).NVqG, Specie(1).FluxTube(1).QCell(1).NVphiG).VphiGH(1)*1e-3))
    
    VpphiGridlinspace= (Specie(s).FluxTube(f).QCell(Qind).VpG(end, 1, 1)- ...
        Specie(s).FluxTube(f).QCell(Qind).VpG(end- 1, 1, 1));
    
    VqGridlinspace= (Specie(s).FluxTube(f).QCell(Qind).VqG(1, end, 1)- ...
        Specie(s).FluxTube(f).QCell(Qind).VqG(1, end- 1, 1));
    
    horzcat('Vp, Vphi linear grid to resolve MB thermal core with resolution [km/s]= ', ...
        num2str(VpphiGridlinspace*1e-3))
    horzcat('Vq linear grid to resolve MB thermal core with resolution [km/s]= ', ...
        num2str(VqGridlinspace*1e-3))
    horzcat('Vp, Vphi sigma [km/s]= ', num2str(Vpphisigma*1e-3), ...
        ', Vp, Vphi sigma factor= ', num2str(VpphisigmaFac))
    horzcat('Vq sigma [km/s]= ', num2str(Vqsigma*1e-3), ...
        ', Vq sigma factor= ', num2str(VqsigmaFac))
        
end

toc
    
disp('ENA PHASE-SPACE GRID GENERATION COMPLETE')

%% -------------------------------------------------------

clc
close all
tic

% GLOBAL EARTH MAGNETIC DIPOLE FIELD MESH GRID:

BXdim= 2;
BYdim= 2;
BZdim= 2;
Bdim= 8;

[Y2, Z2]= meshgrid(-BYdim*RE:(BYdim*RE/Bdim):BYdim*RE, ...
    -BZdim*RE:(BZdim*RE/Bdim):BZdim*RE);

[X3, Y3, Z3]= meshgrid(-BXdim*RE:(BXdim*RE/Bdim):BXdim*RE, -BYdim*RE:(BYdim*RE/Bdim):BYdim*RE, ...
    -BZdim*RE:(BZdim*RE/Bdim):BZdim*RE);

r= RE;
th= 0:pi/50:2*pi;
xunit= r*cos(th);
yunit= r*sin(th);

% -------------------------------------------------------

% EARTH MAGNETIC DIPOLE FIELD VECTOR:

m= 8.06e15; % Earth mag dipole moment [T m^3]
R3= sqrt(X3.^2+ Y3.^2+ Z3.^2);
R2= sqrt(Y2.^2+ Z2.^2);
THETA3= acos(Z3./R3);
THETA2= acos(Z2./R2);
PHI3= atan2(Y3, X3);
PHI2= atan2(Y2, 0e0);
Lsh3= R3./(RE.*(sin(THETA3)).^2);
Lsh2= R2./(RE.*(sin(THETA2)).^2);
ell3= 1+ 3*cos(THETA3).^2;
ell2= 1+ 3*cos(THETA2).^2;
Bmag3= m.*sqrt(ell3)./(R3).^3;
Bmag2= m.*sqrt(ell2)./(R2).^3;
indzero3= find(R3<= RE);
indzero2= find(R2<= RE);
Bmag3(indzero3)= 0;
Bmag2(indzero2)= 0;
Lsh3(indzero3)= 0;
Lsh2(indzero2)= 0;

% fg= (Specie(s).qs.*Bmag2./Specie(s).ms)./(2e0*pi);

Bxvec3= 3.*Bmag3.*cos(THETA3).*sin(THETA3).*cos(PHI3)./sqrt(ell3); % Dipole Bfield in Cart coords.
Byvec3= 3.*Bmag3.*cos(THETA3).*sin(THETA3).*sin(PHI3)./sqrt(ell3);
Bzvec3= Bmag3.*(3.*cos(THETA3).^2- 1)./sqrt(ell3);

Byvec2= 3.*Bmag2.*cos(THETA2).*sin(THETA2).*sin(PHI2)./sqrt(ell2);
Bzvec2= Bmag2.*(3.*cos(THETA2).^2- 1)./sqrt(ell2);

rN= RE;
thetaN= pi;
phiN= 0;
xN= rN*cos(phiN)*sin(thetaN); % Earth N pole coords
yN= rN*sin(phiN)*sin(thetaN);
zN= rN*cos(thetaN);

rS= RE;
thetaS= 0;
phiS= 0;
xS= rS*cos(phiS)*sin(thetaS); % Earth S pole coords
yS= rS*sin(phiS)*sin(thetaS);
zS= rS*cos(thetaS);

% -------------------------------------------------------

% PLOT EARTH MAGNETIC DIPOLE FIELD AND CONFIGURATION-SPACE GRID:

% fignum= 1; % Assign figure number
% fig(fignum)= figure(fignum);
% set(fig(fignum), 'Position', [10 10 1200 1000])
%
% qvr3grid= quiver3(X3/RE, Y3/RE, Z3/RE, Bxvec3, Byvec3, Bzvec3, 3);
% hold on
% set(qvr3grid, 'Color', 'm', 'LineStyle', '-', 'LineWidth', 1)
%
% for s= 1:1:Stot;
%     for f= 1:1:Specie(s).Nf;
%         for Qind= 1:1:Specie(s).FluxTube(f).NqG;
%             plot3(Specie(s).FluxTube(f).QCell(Qind).xGL/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).yGL/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).zGL/RE, 'xk', 'MarkerSize', 50', 'LineWidth', 2.5)
%             plot3(Specie(s).FluxTube(f).QCell(Qind).xGH/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).yGH/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).zGH/RE, 'xk', 'MarkerSize', 50', 'LineWidth', 2.5)
%             plot3(Specie(s).FluxTube(f).QCell(Qind).xGC/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).yGC/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).zGC/RE, '.k', 'MarkerSize', 50', 'LineWidth', 2.5)
%         end
%     end
% end
%
% hold off
% grid on
% xlim([-BXdim BXdim])
% ylim([-BYdim BYdim])
% zlim([-BZdim BZdim])
% xlabel('$x$ [R$_E$]', 'interpreter', 'latex','FontSize', 25)
% ylabel('$y$ [R$_E$]', 'interpreter', 'latex', 'FontSize', 25)
% zlabel('$z$ [R$_E$]', 'interpreter', 'latex', 'FontSize', 25)
% legend('Dipole Field', 'Location', 'NE')
%
% % for s= 1:1:Stot;
% %     for f= 1:1:Specie(s).Nf;
% %         for Qind= 1:1:Specie(s).FluxTube(f).NqG;
% %             text(Specie(s).FluxTube(f).QCell(Qind).xGC/RE, ...
% %                 Specie(s).FluxTube(f).QCell(Qind).yGC/RE, ...
% %                 Specie(s).FluxTube(f).QCell(Qind).zGC/RE, ...
% %                 horzcat('$q_{GC}= ', num2str(Qind), '$'), ...
% %                 'interpreter', 'latex', 'FontSize', 20, 'color', 'red')
% %         end
% %     end
% % end
%
% set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
% set(gca, 'FontSize', 25)
% set(gcf, 'color','white')

% fignum= 2; % Assign figure number
% fig(fignum)= figure(fignum);
% set(fig(fignum), 'Position', [10 10 1000 1000])
%
% qvr2grid= quiver(Y2/RE, Z2/RE, Byvec2, Bzvec2, 2);
% hold on
% set(qvr2grid, 'Color', 'm', 'LineStyle', '-', 'LineWidth', 1)
% [starty startz]= meshgrid(-2:0.3:2, -2:0.3:2);
% sline= streamline(Y2/RE, Z2/RE, Byvec2, Bzvec2, starty, startz);
% set(sline, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1)
% plot(Byvec2, Bzvec2);
% fill(xunit/RE, yunit/RE, 'w');
%
% for s= 1:1:Stot;
%     for f= 1:1:Specie(s).Nf;
%         for Qind= 1:1:Specie(s).FluxTube(f).NqG;
%             plot(Specie(s).FluxTube(f).QCell(Qind).yGL/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).zGL/RE, 'xk', 'MarkerSize', 50', 'LineWidth', 2)
%             plot(Specie(s).FluxTube(f).QCell(Qind).yGH/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).zGH/RE, 'xk', 'MarkerSize', 50', 'LineWidth', 2)
%             plot(Specie(s).FluxTube(f).QCell(Qind).yGC/RE, ...
%                 Specie(s).FluxTube(f).QCell(Qind).zGC/RE, '.k', 'MarkerSize', 50', 'LineWidth', 2)
%         end
%     end
% end
%
% hold off
% grid on
% xlim([-BYdim BYdim])
% ylim([-BZdim BZdim])
% xlabel('$y$ [R$_E$]', 'interpreter', 'latex','FontSize', 25)
% ylabel('$z$ [R$_E$]', 'interpreter', 'latex', 'FontSize', 25)
% % legend('Dipole Field', 'Location', 'NE')
%
% set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
% set(gca, 'FontSize', 25)
% set(gcf, 'color','white')

% -------------------------------------------------------

% PLOT VELOCITY-SPACE GRID:

fignum= 2; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1000 500])

if IONVPERPVECflag == 1
    
    subplot(1, 3, 1)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vperp2ind= Specie(s).FluxTube(f).QCell(Qind).NVperp2G;
                    for Vparind= Specie(s).FluxTube(f).QCell(Qind).NVparG;
                        plot(Specie(s).FluxTube(f).QCell(Qind).Vperp1GC(:, Vperp2ind, Vparind)*1e-3, ...
                            'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                    end
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_{\perp 1}$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{{\perp 1}_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVperp1G]);
    % xlim([1 20])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')
    
    subplot(1, 3, 2)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vperp1ind= Specie(s).FluxTube(f).QCell(Qind).NVperp1G;
                    for Vparind= Specie(s).FluxTube(f).QCell(Qind).NVparG;
                        plot(Specie(s).FluxTube(f).QCell(Qind).Vperp2GC(Vperp1ind, :, Vparind)*1e-3, ...
                            'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                    end
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_{\perp 2}$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{{\perp 2}_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVperp2G]);
    % xlim([1 20])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')

    subplot(1, 3, 3)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vperp1ind= Specie(s).FluxTube(f).QCell(Qind).NVperp1G;
                    for Vperp2ind= Specie(s).FluxTube(f).QCell(Qind).NVperp2G;
                        plot(squeeze(Specie(s).FluxTube(f).QCell(Qind).VparGC(Vperp1ind, Vperp2ind, :))*1e-3, ...
                            'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                    end
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_{\parallel}$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{{\parallel}_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVparG])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')
    
    % expdir= '/Users/robertalbarran/Desktop/FULLDESKTOP/DISSFULL/DISSDOC/Figures/';
    % saveas(figure(fignum), [horzcat(expdir, 'logvelgrid.png')]);
    
else
    
    subplot(1, 2, 1)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vparind= Specie(s).FluxTube(f).QCell(Qind).NVparG;
                    plot(Specie(s).FluxTube(f).QCell(Qind).VperpGC(:, Vparind)*1e-3, ...
                        'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_{\perp}$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{{\perp}_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVperpG]);
    % xlim([1 20])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')

    subplot(1, 2, 2)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vperpind= Specie(s).FluxTube(f).QCell(Qind).NVperpG;
                    plot(Specie(s).FluxTube(f).QCell(Qind).VparGC(Vperpind, :)*1e-3, ...
                        'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_{\parallel}$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{{\parallel}_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVparG])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')
    
    % expdir= '/Users/robertalbarran/Desktop/FULLDESKTOP/DISSFULL/DISSDOC/Figures/';
    % saveas(figure(fignum), [horzcat(expdir, 'logvelgrid.png')]);
    
end

if ENAEXPORTflag == 1
    fignum= 3; % Assign figure number
    fig(fignum)= figure(fignum);
    set(fig(fignum), 'Position', [10 10 1000 500])

    subplot(1, 3, 1)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vqind= Specie(s).FluxTube(f).QCell(Qind).NVqG;
                    for Vphiind= Specie(s).FluxTube(f).QCell(Qind).NVphiG;
                        plot(Specie(s).FluxTube(f).QCell(Qind).VpGC(:, Vqind, Vphiind)*1e-3, ...
                            'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                    end
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_p$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{p_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVpG]);
    % xlim([1 20])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')

    subplot(1, 3, 2)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vpind= Specie(s).FluxTube(f).QCell(Qind).NVpG;
                    for Vphiind= Specie(s).FluxTube(f).QCell(Qind).NVphiG;
                        plot(Specie(s).FluxTube(f).QCell(Qind).VqGC(Vpind, :, Vphiind)*1e-3, ...
                            'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                    end
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_q$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{q_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVpG]);
    % xlim([1 20])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')

    subplot(1, 3, 3)
    for s= 1:1:Stot;
        for f= 1:1:Specie(s).Nf;
            for Qind= 1:1:Specie(s).FluxTube(f).NqG;
                for Vpind= Specie(s).FluxTube(f).QCell(Qind).NVpG;
                    for Vqind= Specie(s).FluxTube(f).QCell(Qind).NVqG;
                        plot(squeeze(Specie(s).FluxTube(f).QCell(Qind).VphiGC(Vpind, Vqind, :))*1e-3, ...
                            'x-k', 'MarkerSize', 20', 'LineWidth', 2)
                    end
                end
            end
        end
    end

    grid on
    xlabel('$\tilde v_\phi$', 'interpreter', 'latex','FontSize', 25)
    ylabel('$v_{\phi_{C}}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
    xlim([1 Specie(1).FluxTube(1).QCell(1).NVpG]);
    % xlim([1 20])
    set(gca, 'zminortick', 'on', 'yminortick', 'on','xminortick', 'on')
    set(gca, 'FontSize', 25)
    set(gcf, 'color','white')

end

toc

disp('PHASE-SPACE GRID VISUALIZATION COMPLETE')

%% -------------------------------------------------------

% EXPORT ALL PARAMETERS PER PARTICLE SPECIES, FLUX TUBE, AND PHASE-SPACE GRID CELL:

clc
close all
tic

if ENAEXPORTflag == 1
        
    VpphiNlinRangeID= fopen(horzcat(datadir, 'VpphiNlinRangemat.bin'), 'w');
    fwrite(VpphiNlinRangeID, VpphiNlinRange, 'real*8');
    fclose(VpphiNlinRangeID);
    
    VqNlinRangeID= fopen(horzcat(datadir, 'VqNlinRangemat.bin'), 'w');
    fwrite(VqNlinRangeID, VqNlinRange, 'real*8');
    fclose(VqNlinRangeID);
    
    VpphiGridlinspaceID= fopen(horzcat(datadir, 'VpphiGridlinspacemat.bin'), 'w');
    fwrite(VpphiGridlinspaceID, VpphiGridlinspace, 'real*8');
    fclose(VpphiGridlinspaceID);
    
    VpphisigmaID= fopen(horzcat(datadir, 'Vpphisigmamat.bin'), 'w');
    fwrite(VpphisigmaID, Vpphisigma, 'real*8');
    fclose(VpphisigmaID);
    
    VpphisigmaFacID= fopen(horzcat(datadir, 'VpphisigmaFacmat.bin'), 'w');
    fwrite(VpphisigmaFacID, VpphisigmaFac, 'real*8');
    fclose(VpphisigmaFacID);
    
    VqGridlinspaceID= fopen(horzcat(datadir, 'VqGridlinspacemat.bin'), 'w');
    fwrite(VqGridlinspaceID, VqGridlinspace, 'real*8');
    fclose(VqGridlinspaceID);
    
    VqsigmaID= fopen(horzcat(datadir, 'Vqsigmamat.bin'), 'w');
    fwrite(VqsigmaID, Vqsigma, 'real*8');
    fclose(VqsigmaID);
    
    VqsigmaFacID= fopen(horzcat(datadir, 'VqsigmaFacmat.bin'), 'w');
    fwrite(VqsigmaFacID, VqsigmaFac, 'real*8');
    fclose(VqsigmaFacID);    
end
    
if IONEXPORTflag == 1
    
    Vperp12NlinRangeID= fopen(horzcat(datadir, 'Vperp12NlinRangemat.bin'), 'w');
    fwrite(Vperp12NlinRangeID, Vperp12NlinRange, 'real*8');
    fclose(Vperp12NlinRangeID);
    
    VparNlinRangeID= fopen(horzcat(datadir, 'VparNlinRangemat.bin'), 'w');
    fwrite(VparNlinRangeID, VparNlinRange, 'real*8');
    fclose(VparNlinRangeID);
    
    Vperp12GridlinspaceID= fopen(horzcat(datadir, 'Vperp12Gridlinspacemat.bin'), 'w');
    fwrite(Vperp12GridlinspaceID, Vperp12Gridlinspace, 'real*8');
    fclose(Vperp12GridlinspaceID);
    
    Vperp12sigmaID= fopen(horzcat(datadir, 'Vperp12sigmamat.bin'), 'w');
    fwrite(Vperp12sigmaID, Vperp12sigma, 'real*8');
    fclose(Vperp12sigmaID);
    
    Vperp12sigmaFacID= fopen(horzcat(datadir, 'Vperp12sigmaFacmat.bin'), 'w');
    fwrite(Vperp12sigmaFacID, Vperp12sigmaFac, 'real*8');
    fclose(Vperp12sigmaFacID);
    
    VparGridlinspaceID= fopen(horzcat(datadir, 'VparGridlinspacemat.bin'), 'w');
    fwrite(VparGridlinspaceID, VparGridlinspace, 'real*8');
    fclose(VparGridlinspaceID);
    
    VparsigmaID= fopen(horzcat(datadir, 'Vparsigmamat.bin'), 'w');
    fwrite(VparsigmaID, Vparsigma, 'real*8');
    fclose(VparsigmaID);
    
    VparsigmaFacID= fopen(horzcat(datadir, 'VparsigmaFacmat.bin'), 'w');
    fwrite(VparsigmaFacID, VparsigmaFac, 'real*8');
    fclose(VparsigmaFacID);
    
    for s= 1:1:Stot
        NfID= fopen(horzcat(datadir, ...
            num2str(1), '_', 'Nfmat.bin'), 'w');
        fwrite(NfID, Specie(s).Nf, 'integer*8');
        fclose(NfID);
        
        for f= 1:1:Specie(s).Nf
                        
            NqGID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'NqGmat.bin'), 'w');
            fwrite(NqGID, Specie(s).FluxTube(f).NqG, 'integer*8');
            fclose(NqGID);
            HiCratioMeanID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'HiCratioMeanmat.bin'), 'w');
            fwrite(HiCratioMeanID, Specie(s).FluxTube(f).HiCratioMean, 'real*8');
            fclose(HiCratioMeanID);
            TeID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'Temat.bin'), 'w');
            fwrite(TeID, Specie(s).FluxTube(f).Te, 'real*8');
            fclose(TeID);            
                        
            for Qind= 1:1:1 %Specie(s).FluxTube(f).NqG;
                
                if IONVPERPVECflag == 1
                    
                    NVperp1GpID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVperp1Gpmat.bin'), 'w');
                    fwrite(NVperp1GpID, Specie(s).FluxTube(f).QCell(Qind).NVperp1Gp, 'integer*8');
                    fclose(NVperp1GpID);
                    NVperp1GID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVperp1Gmat.bin'), 'w');
                    fwrite(NVperp1GID, Specie(s).FluxTube(f).QCell(Qind).NVperp1G, 'integer*8');
                    fclose(NVperp1GID);
                    
                    NVperp2GpID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVperp2Gpmat.bin'), 'w');
                    fwrite(NVperp2GpID, Specie(s).FluxTube(f).QCell(Qind).NVperp2Gp, 'integer*8');
                    fclose(NVperp2GpID);
                    NVperp2GID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVperp2Gmat.bin'), 'w');
                    fwrite(NVperp2GID, Specie(s).FluxTube(f).QCell(Qind).NVperp2G, 'integer*8');
                    fclose(NVperp2GID);
                    
                else
                    
                    NVperpGpID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVperpGpmat.bin'), 'w');
                    fwrite(NVperpGpID, Specie(s).FluxTube(f).QCell(Qind).NVperpGp, 'integer*8');
                    fclose(NVperpGpID);
                    NVperpGID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVperpGmat.bin'), 'w');
                    fwrite(NVperpGID, Specie(s).FluxTube(f).QCell(Qind).NVperpG, 'integer*8');
                    fclose(NVperpGID);
                
                end
                
                NVparGpID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVparGpmat.bin'), 'w');
                fwrite(NVparGpID, Specie(s).FluxTube(f).QCell(Qind).NVparGp, 'integer*8');
                fclose(NVparGpID);
                NVparGID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVparGmat.bin'), 'w');
                fwrite(NVparGID, Specie(s).FluxTube(f).QCell(Qind).NVparG, 'integer*8');
                fclose(NVparGID);
                
                if ENAEXPORTflag == 1
                    NVpGID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVpGmat.bin'), 'w');
                    fwrite(NVpGID, Specie(s).FluxTube(f).QCell(Qind).NVpG, 'integer*8');
                    fclose(NVpGID);
                    NVqGID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVqGmat.bin'), 'w');
                    fwrite(NVqGID, Specie(s).FluxTube(f).QCell(Qind).NVqG, 'integer*8');
                    fclose(NVqGID);
                    NVphiGID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), '_', 'NVphiGmat.bin'), 'w');
                    fwrite(NVphiGID, Specie(s).FluxTube(f).QCell(Qind).NVphiG, 'integer*8');
                    fclose(NVphiGID);
                end
            end
        end
        Qindns0ID= fopen(horzcat(datadir, ...
            num2str(s), '_', 'Qindns0mat.bin'), 'w');
        fwrite(Qindns0ID, Specie(s).Qindns0, 'integer*8');
        fclose(Qindns0ID);
        msID= fopen(horzcat(datadir, ...
            num2str(s), '_', 'msmat.bin'), 'w');
        fwrite(msID, Specie(s).ms, 'real*8');
        fclose(msID);
        qsID= fopen(horzcat(datadir, ...
            num2str(s), '_', 'qsmat.bin'), 'w');
        fwrite(qsID, Specie(s).qs, 'real*8');
        fclose(qsID);
        
    end
            
    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf
            qGCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'qGCmat.bin'), 'w');
            fwrite(qGCID, Specie(s).FluxTube(f).qGC(:), 'real*8');
            fclose(qGCID);
            hqCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'hqCmat.bin'), 'w');
            fwrite(hqCID, Specie(s).FluxTube(f).hqC(:), 'real*8');
            fclose(hqCID);
            dpCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'dpCmat.bin'), 'w');
            fwrite(dpCID, Specie(s).FluxTube(f).dpC(:), 'real*8');
            fclose(dpCID);
            dqCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'dqCmat.bin'), 'w');
            fwrite(dqCID, Specie(s).FluxTube(f).dqC(:), 'real*8');
            fclose(dqCID);
            dphiCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'dphiCmat.bin'), 'w');
            fwrite(dphiCID, Specie(s).FluxTube(f).dphiC(:), 'real*8');
            fclose(dphiCID);
            TsPerpID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'TsPerpmat.bin'), 'w');
            fwrite(TsPerpID, Specie(s).FluxTube(f).TsPerp(:), 'real*8');
            fclose(TsPerpID);
            TsParID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'TsParmat.bin'), 'w');
            fwrite(TsParID, Specie(s).FluxTube(f).TsPar(:), 'real*8');
            fclose(TsParID);
            rGCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'rGCmat.bin'), 'w');
            fwrite(rGCID, Specie(s).FluxTube(f).rGC(:), 'real*8');
            fclose(rGCID);
            phiGCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'phiGCmat.bin'), 'w');
            fwrite(phiGCID, Specie(s).FluxTube(f).phiGC(:), 'real*8');
            fclose(phiGCID);                
            thetaGCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'thetaGCmat.bin'), 'w');
            fwrite(thetaGCID, Specie(s).FluxTube(f).thetaGC(:), 'real*8');
            fclose(thetaGCID);
            ellGCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'ellGCmat.bin'), 'w');
            fwrite(ellGCID, Specie(s).FluxTube(f).ellGC(:), 'real*8');
            fclose(ellGCID);
            qGLID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'qGLmat.bin'), 'w');
            fwrite(qGLID, Specie(s).FluxTube(f).qGL(:), 'real*8');
            fclose(qGLID);
            qGHID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'qGHmat.bin'), 'w');
            fwrite(qGHID, Specie(s).FluxTube(f).qGH(:), 'real*8');
            fclose(qGHID);
            pGCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'pGCmat.bin'), 'w');
            fwrite(pGCID, Specie(s).FluxTube(f).pGC(:), 'real*8');
            fclose(pGCID);
            d3xCID= fopen(horzcat(datadir, ...
                num2str(s), '_', num2str(f), '_', 'd3xCmat.bin'), 'w');
            fwrite(d3xCID, Specie(s).FluxTube(f).d3xC(:), 'real*8');
            fclose(d3xCID);
        end
    end
    
    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf
            for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                
                if IONVPERPVECflag == 1
                    
                    dVperp1CID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'dVperp1Cmat.bin'), 'w');
                    fwrite(dVperp1CID, Specie(s).FluxTube(f).QCell(Qind).dVperp1C(:, :, :), ...
                        'real*8');
                    Vperp1GLID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'Vperp1GLmat.bin'), 'w');
                    fwrite(Vperp1GLID, Specie(s).FluxTube(f).QCell(Qind).Vperp1GL(:, :, :), 'real*8');
                    fclose(Vperp1GLID);
                    Vperp1GHID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'Vperp1GHmat.bin'), 'w');
                    fwrite(Vperp1GHID, Specie(s).FluxTube(f).QCell(Qind).Vperp1GH(:, :, :), 'real*8');
                    fclose(Vperp1GHID);
                    Vperp1GCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'Vperp1GCmat.bin'), 'w');
                    fwrite(Vperp1GCID, Specie(s).FluxTube(f).QCell(Qind).Vperp1GC(:, :, :), 'real*8');
                    fclose(Vperp1GCID);
                    dVperp2CID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'dVperp2Cmat.bin'), 'w');
                    fwrite(dVperp2CID, Specie(s).FluxTube(f).QCell(Qind).dVperp2C(:, :, :), ...
                        'real*8');
                    fclose(dVperp2CID);
                    Vperp2GLID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'Vperp2GLmat.bin'), 'w');
                    fwrite(Vperp2GLID, Specie(s).FluxTube(f).QCell(Qind).Vperp2GL(:, :, :), 'real*8');
                    fclose(Vperp2GLID);
                    Vperp2GHID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'Vperp2GHmat.bin'), 'w');
                    fwrite(Vperp2GHID, Specie(s).FluxTube(f).QCell(Qind).Vperp2GH(:, :, :), 'real*8');
                    fclose(Vperp2GHID);
                    Vperp2GCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'Vperp2GCmat.bin'), 'w');
                    fwrite(Vperp2GCID, Specie(s).FluxTube(f).QCell(Qind).Vperp2GC(:, :, :), 'real*8');
                    fclose(Vperp2GCID);
                    dVparCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'dVparCmat.bin'), 'w');
                    fwrite(dVparCID, Specie(s).FluxTube(f).QCell(Qind).dVparC(:, :, :), ...
                        'real*8');
                    fclose(dVparCID);
                    VparGLID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VparGLmat.bin'), 'w');
                    fwrite(VparGLID, Specie(s).FluxTube(f).QCell(Qind).VparGL(:, :, :), 'real*8');
                    fclose(VparGLID);
                    VparGHID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VparGHmat.bin'), 'w');
                    fwrite(VparGHID, Specie(s).FluxTube(f).QCell(Qind).VparGH(:, :, :), 'real*8');
                    fclose(VparGHID);
                    VparGCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VparGCmat.bin'), 'w');
                    fwrite(VparGCID, Specie(s).FluxTube(f).QCell(Qind).VparGC(:, :, :), 'real*8');
                    fclose(VparGCID);
                    d3vCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'd3vCmat.bin'), 'w');
                    fwrite(d3vCID, Specie(s).FluxTube(f).QCell(Qind).d3vC(:, :, :), 'real*8');
                    fclose(d3vCID);
                    
                else
                   
                    dVperpCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'dVperpCmat.bin'), 'w');
                    fwrite(dVperpCID, Specie(s).FluxTube(f).QCell(Qind).dVperpC(:, :), ...
                        'real*8');
                    fclose(dVperpCID);
                    dVparCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'dVparCmat.bin'), 'w');
                    fwrite(dVparCID, Specie(s).FluxTube(f).QCell(Qind).dVparC(:, :), ...
                        'real*8');
                    fclose(dVparCID);
                    VperpGLID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        'VperpGLmat.bin'), 'w');
                    fwrite(VperpGLID, Specie(s).FluxTube(f).QCell(Qind).VperpGL(:, :), 'real*8');
                    fclose(VperpGLID);
                    VperpGHID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VperpGHmat.bin'), 'w');
                    fwrite(VperpGHID, Specie(s).FluxTube(f).QCell(Qind).VperpGH(:, :), 'real*8');
                    fclose(VperpGHID);
                    VperpGCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VperpGCmat.bin'), 'w');
                    fwrite(VperpGCID, Specie(s).FluxTube(f).QCell(Qind).VperpGC(:, :), 'real*8');
                    fclose(VperpGCID);
                    VparGLID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VparGLmat.bin'), 'w');
                    fwrite(VparGLID, Specie(s).FluxTube(f).QCell(Qind).VparGL(:, :), 'real*8');
                    fclose(VparGLID);
                    VparGHID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VparGHmat.bin'), 'w');
                    fwrite(VparGHID, Specie(s).FluxTube(f).QCell(Qind).VparGH(:, :), 'real*8');
                    fclose(VparGHID);
                    VparGCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'VparGCmat.bin'), 'w');
                    fwrite(VparGCID, Specie(s).FluxTube(f).QCell(Qind).VparGC(:, :), 'real*8');
                    fclose(VparGCID);
                    d3vCID= fopen(horzcat(datadir, ...
                        num2str(s), '_', num2str(f), '_', num2str(Qind), ...
                        '_', 'd3vCmat.bin'), 'w');
                    fwrite(d3vCID, Specie(s).FluxTube(f).QCell(Qind).d3vC(:, :), 'real*8');
                    fclose(d3vCID);
                    
                end
                
            end
        end
    end
end

if ENAEXPORTflag == 1
    
    mNeutID= fopen(horzcat(datadir, 'mNeutmat.bin'), 'w');
    fwrite(mNeutID, mNeut, 'real*8');
    fclose(mNeutID);

    TNeutID= fopen(horzcat(datadir, 'TNeutmat.bin'), 'w');
    fwrite(TNeutID, TNeut, 'real*8');
    fclose(TNeutID);

    for s= 1:1:Stot
        for f= 1:1:Specie(s).Nf
            for Qind= 1%:1:Specie(s).FluxTube(f).NqG;
                VpGLID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VpGLmat.bin'), 'w');
                fwrite(VpGLID, Specie(s).FluxTube(f).QCell(Qind).VpGL(:, :, :), 'real*8');
                fclose(VpGLID);
                VpGHID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VpGHmat.bin'), 'w');
                fwrite(VpGHID, Specie(s).FluxTube(f).QCell(Qind).VpGH(:, :, :), 'real*8');
                fclose(VpGHID);
                VpGCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VpGCmat.bin'), 'w');
                fwrite(VpGCID, Specie(s).FluxTube(f).QCell(Qind).VpGC(:, :, :), 'real*8');
                fclose(VpGCID);
                VqGLID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VqGLmat.bin'), 'w');
                fwrite(VqGLID, Specie(s).FluxTube(f).QCell(Qind).VqGL(:, :, :), 'real*8');
                fclose(VqGLID);
                VqGHID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VqGHmat.bin'), 'w');
                fwrite(VqGHID, Specie(s).FluxTube(f).QCell(Qind).VqGH(:, :, :), 'real*8');
                fclose(VqGHID);
                VqGCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VqGCmat.bin'), 'w');
                fwrite(VqGCID, Specie(s).FluxTube(f).QCell(Qind).VqGC(:, :, :), 'real*8');
                fclose(VqGCID);
                VphiGLID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VphiGLmat.bin'), 'w');
                fwrite(VphiGLID, Specie(s).FluxTube(f).QCell(Qind).VphiGL(:, :, :), 'real*8');
                fclose(VphiGLID);
                VphiGHID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VphiGHmat.bin'), 'w');
                fwrite(VphiGHID, Specie(s).FluxTube(f).QCell(Qind).VphiGH(:, :, :), 'real*8');
                fclose(VphiGHID);
                VphiGCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'VphiGCmat.bin'), 'w');
                fwrite(VphiGCID, Specie(s).FluxTube(f).QCell(Qind).VphiGC(:, :, :), 'real*8');
                fclose(VphiGCID);
                hVpCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'hVpCmat.bin'), 'w');
                fwrite(hVpCID, Specie(s).FluxTube(f).QCell(Qind).hVpC(:, :, :), 'real*8');
                fclose(hVpCID);
                hVqCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'hVqCmat.bin'), 'w');
                fwrite(hVqCID, Specie(s).FluxTube(f).QCell(Qind).hVqC(:, :, :), 'real*8');
                fclose(hVqCID);
                hVphiCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'hVphiCmat.bin'), 'w');
                fwrite(hVphiCID, Specie(s).FluxTube(f).QCell(Qind).hVphiC(:, :, :), 'real*8');
                fclose(hVphiCID);
                dVpCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'dVpCmat.bin'), 'w');
                fwrite(dVpCID, Specie(s).FluxTube(f).QCell(Qind).dVpC(:, :, :), 'real*8');
                fclose(dVpCID);
                dVqCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'dVqCmat.bin'), 'w');
                fwrite(dVqCID, Specie(s).FluxTube(f).QCell(Qind).dVqC(:, :, :), 'real*8');
                fclose(dVqCID);
                dVphiCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'dVphiCmat.bin'), 'w');
                fwrite(dVphiCID, Specie(s).FluxTube(f).QCell(Qind).dVphiC(:, :, :), 'real*8');
                fclose(dVphiCID);
                d33vCID= fopen(horzcat(datadir, ...
                    num2str(s), '_', num2str(f), '_', num2str(Qind), '_', ...
                    'd33vCmat.bin'), 'w');
                fwrite(d33vCID, Specie(s).FluxTube(f).QCell(Qind).d33vC(:, :, :), 'real*8');
                fclose(d33vCID);
            end
        end
    end
end

toc

disp('PHASE-SPACE GRID EXPORT COMPLETE')
