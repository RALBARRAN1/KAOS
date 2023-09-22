clc
clear all
close all

tic

% IMPORT MOMENTS

RE= 638e4; % Earth radius [m]
kB= 138e-25; % Boltzmann constant [m^2 kg s^-2 K^-1]
mion= (16e0)*(1.67e-27); % O+ ion mass [kg]
qion= 1.602e-19; % Charge [C]
NqICA= 1;
NqICB= 25; % 25 for W0 and W0F, 30 for LA and PC00, 19 for VV0 and VV0F

NNtT= 153;
ranksize= 2;
root= 1; r= 1; Stot= 1; s= 1; Rank(1).Specie(1).Nf= 1; f= 1;

DensityOutputflag= 0;
ICRflag= 0;
MomentFilterFlag= 1;
QExchangeFlag= 0;
MAfilterPt= 11e0;

NtT= 3e4;
Q0datfac= NtT/(NNtT);
ZeroNaNflag= 1;
NqLB= NqICA; NqUB= NqICB;

IONVPERPVECflag= 1;
if IONVPERPVECflag == 1
    if ICRflag == 0
        NVperp1G= 28; 
        NVperp2G= 28;
        NVparG= 28; 
    end
    if ICRflag == 1
        NVperp1G= 28; 
        NVperp2G= 28;
        NVparG= 28;
    end
else
    NVperpG= 9;
    NVparG= 19;
end

if QExchangeFlag == 1;
    NVpG= 30; 
    NVqG= 20;
    NVphiG= 30;
end

% simnum= 'PC00icr11scfINPUTa2'; % Assign simulation number
% datadirgrid= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOGridPC00/';

% simnum= 'VV0FicrV04epar2scINPUTa0'; % Assign simulation number
% datadirgrid= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOGridVV0F/';

simnum= 'W0icr00epar4scINPUTa1'; % Assign simulation number
datadirgrid= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOGridW0/';

% simnum= 'LA1ncscOUTPUTa0'; % Assign simulation number
% datadirgrid= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOGridLA1/';

datadir= horzcat('/Users/robertalbarran/Desktop/FULLDESKTOP/FFSIMS/', simnum, '/KMIOData', simnum, '/');
dataexpdir= horzcat('/Users/robertalbarran/Desktop/FULLDESKTOP/FFSIMS/', simnum, '/');

% datadir= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOData/';
% dataexpdir= '/Users/robertalbarran/Desktop/FULLDESKTOP/KMIOSims/';
 
qGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'qGCmat.bin'));
Rank(r).Specie(s).FluxTube(f).qGC= fread(qGCID, NqUB, 'real*8');
fclose(qGCID);
pGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'pGCmat.bin'));
Rank(r).Specie(s).FluxTube(f).pGC= fread(pGCID, NqUB, 'real*8');
fclose(pGCID);
hqCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'hqCmat.bin'));
Rank(r).Specie(s).FluxTube(f).hqC= fread(hqCID, NqUB, 'real*8');
fclose(hqCID);
dqCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'dqCmat.bin'));
Rank(r).Specie(s).FluxTube(f).dqC= fread(dqCID, NqUB, 'real*8');
fclose(dqCID);
qGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'qGLmat.bin'));
Rank(r).Specie(s).FluxTube(f).qGL= fread(qGLID, NqUB, 'real*8');
fclose(qGLID);
qGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'qGHmat.bin'));
Rank(r).Specie(s).FluxTube(f).qGH= fread(qGHID, NqUB, 'real*8');
fclose(qGHID);
d3xCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'd3xCmat.bin'));
Rank(r).Specie(s).FluxTube(f).d3xC= fread(d3xCID, NqUB, 'real*8');
fclose(d3xCID);
rGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'rGCmat.bin'));
Rank(r).Specie(s).FluxTube(f).rGC= fread(rGCID, NqUB, 'real*8');
fclose(rGCID);
thetaGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'thetaGCmat.bin'));
Rank(r).Specie(s).FluxTube(f).thetaGC= fread(thetaGCID, NqUB, 'real*8');
fclose(thetaGCID);
ellGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'ellGCmat.bin'));
Rank(r).Specie(s).FluxTube(f).ellGC= fread(ellGCID, NqUB, 'real*8');
fclose(ellGCID);
Rank(r).Specie(s).FluxTube(f).rGCp(:)= ...
    (Rank(r).Specie(s).FluxTube(f).rGC(:)- RE).*1e-3;
TsPerpID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'TsPerpmat.bin'));
Rank(r).Specie(s).FluxTube(f).TsPerp= fread(TsPerpID, NqUB, 'real*8');
fclose(TsPerpID);
Rank(r).Specie(s).FluxTube(f).TsPerpp(:)= Rank(r).Specie(s).FluxTube(f).TsPerp(:);
TsParID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', ...
    'TsParmat.bin'));
Rank(r).Specie(s).FluxTube(f).TsPar= fread(TsParID, NqUB, 'real*8');
fclose(TsParID);
Rank(r).Specie(s).FluxTube(f).TsParp(:)= Rank(r).Specie(s).FluxTube(f).TsPar(:);
Rank(r).Specie(s).FluxTube(f).Ts(:)= ...
    (1e0/3e0).*Rank(r).Specie(s).FluxTube(f).TsParp(:)+ ...
    (2e0/3e0).*Rank(r).Specie(s).FluxTube(f).TsPerpp(:);
Rank(r).Specie(s).FluxTube(f).Tsp(:)= Rank(r).Specie(s).FluxTube(f).Ts(:);
nsnormfacTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(s), '_', ...
    num2str(f), '_', 'nsnormfacTfort.bin'));
Rank(r).Specie(s).FluxTube(f).nsnormfacT= fread(nsnormfacTID, 1, 'real*8');
fclose(nsnormfacTID);

for nn= 1:1:NNtT+ 1
    TimeTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'TimeTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).TimeT(nn)= fread(TimeTID, 1, 'real*8');
    fclose(TimeTID);
    TeNTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'TeNTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).TeNT(nn)= fread(TimeTID, 1, 'real*8');
    fclose(TeNTID);
end

for Qind= 1
    if IONVPERPVECflag == 1
        
        dVperp1CID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'dVperp1Cmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp1C= ...
            fread(dVperp1CID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(dVperp1CID);
        dVperp2CID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'dVperp2Cmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp2C= ...
            fread(dVperp2CID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(dVperp2CID);
        dVparCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
           '_', 'dVparCmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVparC= ...
            fread(dVparCID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(dVparCID);
    
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp1C= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp1C, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp2C= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp2C, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVparC= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVparC, NVperp1G, NVperp2G, NVparG);
        
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp1Cp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp1C(:, :, :);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp2Cp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVperp2C(:, :, :);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVparCp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).dVparC(:, :, :);

        Vperp1GLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'Vperp1GLmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GL= ...
            fread(Vperp1GLID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(Vperp1GLID);
    
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GL= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GL, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GLp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GL(:, :, :);

        Vperp1GHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'Vperp1GHmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GH= ...
            fread(Vperp1GHID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(Vperp1GHID);
    
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GH= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GH, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GHp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GH(:, :, :);

        Vperp1GCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'Vperp1GCmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GC= ...
            fread(Vperp1GCID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(Vperp1GCID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GC= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GC, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GCp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp1GC(:, :, :);

        Vperp2GLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'Vperp2GLmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GL= ...
            fread(Vperp2GLID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(Vperp2GLID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GL= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GL, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GLp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GL(:, :, :);

        Vperp2GHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'Vperp2GHmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GH= ...
            fread(Vperp2GHID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(Vperp2GHID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GH= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GH, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GHp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GH(:, :, :);

        Vperp2GCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'Vperp2GCmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GC= ...
            fread(Vperp2GCID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(Vperp2GCID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GC= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GC, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GCp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).Vperp2GC(:, :, :);

        VparGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'VparGLmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGL= ...
            fread(VparGLID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(VparGLID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGL= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGL, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGLp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGL(:, :, :);

        VparGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'VparGHmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGH= ...
            fread(VparGHID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(VparGHID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGH= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGH, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGHp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGH(:, :, :);

        VparGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'VparGCmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGC= ...
            fread(VparGCID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(VparGCID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGC= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGC, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGCp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGC(:, :, :);

        d3vCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
            '_', 'd3vCmat.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).d3vC= ...
            fread(d3vCID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(d3vCID);

        Rank(r).Specie(s).FluxTube(f).QCell(Qind).d3vC= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).d3vC, NVperp1G, NVperp2G, NVparG);
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).d3vCp(:, :, :)= ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).d3vC(:, :, :);

    else
        
        for Vperpind= 1:1:NVperpG;
            for Vparind= 1:1:NVparG;
                
                VperpGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'VperpGLmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGL= ...
                    fread(VperpGLID, 1, 'real*8');
                fclose(VperpGLID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VperpGLp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGL;
                
                VperpGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'VperpGHmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGH= ...
                    fread(VperpGHID, 1, 'real*8');
                fclose(VperpGHID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VperpGHp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGH;
                
                VparGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'VparGLmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGL= ...
                    fread(VparGLID, 1, 'real*8');
                fclose(VparGLID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGLp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGL;
                
                VparGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'VparGHmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGH= ...
                    fread(VparGHID, 1, 'real*8');
                fclose(VparGHID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGHp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGH;
                
                VperpGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'VperpGCmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGC= ...
                    fread(VperpGCID, 1, 'real*8');
                fclose(VperpGCID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VperpGCp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VperpGC;
                
                VparGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'VparGCmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGC= ...
                    fread(VparGCID, 1, 'real*8');
                fclose(VparGCID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VparGCp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).VparGC;
                
                d3vCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                    '_', num2str(Vperpind), '_', num2str(Vparind), '_', 'd3vCmat.bin'));
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC= ...
                    fread(d3vCID, 1, 'real*8');
                fclose(d3vCID);
                
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).d3vCp(Vperpind, Vparind)= ...
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VCell(Vperpind, Vparind).d3vC;
            end
        end
    end
    
    if QExchangeFlag == 1
        for Vpind= 1:1:NVpG;
            for Vqind= 1:1:NVqG;
                for Vphiind= 1:1:NVphiG;
                    
                    VpGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VpGLmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGL= ...
                        fread(VpGLID, 1, 'real*8');
                    fclose(VpGLID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VpGLp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGL;
                    
                    VpGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VpGHmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGH= ...
                        fread(VpGHID, 1, 'real*8');
                    fclose(VpGHID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VpGHp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGH;
                    
                    VpGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VpGCmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC= ...
                        fread(VpGCID, 1, 'real*8');
                    fclose(VpGCID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VpGCp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VpGC;
                    
                    VqGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VqGLmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGL= ...
                        fread(VqGLID, 1, 'real*8');
                    fclose(VqGLID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VqGLp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGL;
                    
                    VqGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VqGHmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGH= ...
                        fread(VqGHID, 1, 'real*8');
                    fclose(VqGHID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VqGHp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGH;
                    
                    VqGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VqGCmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC= ...
                        fread(VqGCID, 1, 'real*8');
                    fclose(VqGCID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VqGCp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VqGC;
                    
                    VphiGLID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VphiGLmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGL= ...
                        fread(VphiGLID, 1, 'real*8');
                    fclose(VphiGLID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VphiGLp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGL;
                    
                    VphiGHID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VphiGHmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGH= ...
                        fread(VphiGHID, 1, 'real*8');
                    fclose(VphiGHID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VphiGHp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGH;
                    
                    VphiGCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'VphiGCmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGC= ...
                        fread(VphiGCID, 1, 'real*8');
                    fclose(VphiGCID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).VphiGCp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).VphiGC;
                    
                    d33vCID= fopen(horzcat(datadirgrid, num2str(s), '_', num2str(f), '_', num2str(1), ...
                        '_', num2str(Vpind), '_', num2str(Vqind), '_', num2str(Vphiind), '_', 'd33vCmat.bin'));
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).d33vC= ...
                        fread(d33vCID, 1, 'real*8');
                    fclose(d33vCID);
                    
                    Rank(r).Specie(s).FluxTube(f).QCell(Qind).d33vCp(Vpind, Vqind, Vphiind)= ...
                        Rank(r).Specie(s).FluxTube(f).QCell(Qind).V3Cell(Vpind, Vqind, Vphiind).d33vC;
                    
                end
            end
        end
    end
end

EGMagRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(s), '_', ...
    num2str(f), '_', 'EGMagRTfort.bin'));
Rank(r).Specie(s).FluxTube(f).EGMagRTr(:)= fread(EGMagRTID, ...
    (((NqUB- NqLB)+ 1)), 'real*8');
fclose(EGMagRTID);

for nn= 1:1:NNtT+ 1
    if QExchangeFlag == 1
        sigmaIonNeutRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'sigmaIonNeutRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).sigmaIonNeutRT(nn, :)= fread(sigmaIonNeutRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(sigmaIonNeutRTID);
        
        nuIonNeutRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'nuIonNeutRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).nuIonNeutRT(nn, :)= fread(nuIonNeutRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(nuIonNeutRTID);
    end
    M0phRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'M0phRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, :)= fread(M0phRTID, ...
        ((NqUB- NqLB)+ 1), 'real*8');
    fclose(M0phRTID);
        
    M1ParphRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'M1ParphRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, :)= fread(M1ParphRTID, ...
        ((NqUB- NqLB)+ 1), 'real*8');
    fclose(M1ParphRTID);
    
    if MomentFilterFlag == 1
        M0FiltAvrgRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M0FiltAvrgRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M0FiltAvrgRTr(nn, :)= fread(M0FiltAvrgRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M0FiltAvrgRTID);
        
        M1ParFiltAvrgRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M1ParFiltAvrgRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M1ParFiltAvrgRTr(nn, :)= fread(M1ParFiltAvrgRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M1ParFiltAvrgRTID);
    end
    
    if IONVPERPVECflag == 1
        
        M1Perp1phRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M1Perp1phRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, :)= fread(M1Perp1phRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M1Perp1phRTID);
        
        M2Perp1phRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M2Perp1phRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(nn, :)= fread(M2Perp1phRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M2Perp1phRTID);
        
        M1Perp2phRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M1Perp2phRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, :)= fread(M1Perp2phRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M1Perp2phRTID);
        
        M2Perp2phRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M2Perp2phRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(nn, :)= fread(M2Perp2phRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M2Perp2phRTID);
        
    else
        
        M1PerpphRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M1PerpphRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, :)= fread(M1PerpphRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M1PerpphRTID);
        
        M2PerpphRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', 'M2PerpphRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(nn, :)= fread(M2PerpphRTID, ...
            ((NqUB- NqLB)+ 1), 'real*8');
        fclose(M2PerpphRTID);
        
    end
    
    M2ParphRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'M2ParphRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).M2ParphRTr(nn, :)= fread(M2ParphRTID, ...
        ((NqUB- NqLB)+ 1), 'real*8');
    fclose(M2ParphRTID);
    
    M2phRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'M2phRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).M2phRTr(nn, :)= fread(M2phRTID, ...
        ((NqUB- NqLB)+ 1), 'real*8');
    fclose(M2phRTID);
    
    EAInertialRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'EAInertialRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).EAInertialRTr(nn, :)= fread(EAInertialRTID, ...
        (((NqUB- NqLB)+ 1)), 'real*8');
    fclose(EAInertialRTID);
    
    EAPressureRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'EAPressureRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).EAPressureRTr(nn, :)= fread(EAPressureRTID, ...
        (((NqUB- NqLB)+ 1)), 'real*8');
    fclose(EAPressureRTID);
    
    EAMagRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
        num2str(f), '_', 'EAMagRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).EAMagRTr(nn, :)= fread(EAMagRTID, ...
        (((NqUB- NqLB)+ 1)), 'real*8');
    fclose(EAMagRTID);
                       
    if (Rank(r).Specie(s).FluxTube(f).EAMagRTr(nn, :) ~= ...
        Rank(r).Specie(s).FluxTube(f).EAPressureRTr(nn, :)+ ...
        Rank(r).Specie(s).FluxTube(f).EAInertialRTr(nn, :))
        disp('Inconsistent Ambipolar magnitude summation')  
    end
    
    for Qind= NqLB:1:NqUB
        if Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind) == 0e0
            Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)= NaN;
            Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)= NaN;
            Rank(r).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind)= NaN;
            Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)= NaN;
            Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)= NaN;
            Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(nn, Qind)= NaN;
            Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(nn, Qind)= NaN;
        end 
    end
end
for nn= 1:1:NNtT+ 1
    for Qind= NqLB:1:NqUB;
        if Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind) == 0e0
            
            if ZeroNaNflag == 1
                Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).M2phRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)= 0e0;
                
                if IONVPERPVECflag == 1
                    
                    Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).Perp1Temp(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).Perp2Temp(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).Perp1Flux(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).Perp2Flux(nn, Qind)= 0e0;
                    PitchAngleM(nn, Qind)= 0e0;
                    GyroAngleM(nn, Qind)= 0e0;
                    
                else
                    
                    Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).PerpTemp(nn, Qind)= 0e0;
                    Rank(r).Specie(s).FluxTube(f).PerpFlux(nn, Qind)= 0e0;
                    
                end
                
                Rank(r).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).EAMagRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).EAPressureRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).EAInertialRTr(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).ParTemp(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).Temp(nn, Qind)= 0e0;
                Rank(r).Specie(s).FluxTube(f).ParFlux(nn, Qind)= 0e0;
                
            else
                
                Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).M2phRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)= NaN;
                
                if IONVPERPVECflag == 1
                    
                    Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).Perp1Temp(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).Perp2Temp(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).Perp1Flux(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).Perp2Flux(nn, Qind)= NaN;
                    PitchAngleM(nn, Qind)= NaN;
                    GyroAngleM(nn, Qind)= NaN;
                    
                else
                    
                    Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).PerpTemp(nn, Qind)= NaN;
                    Rank(r).Specie(s).FluxTube(f).PerpFlux(nn, Qind)= NaN;
                    
                end
                
                Rank(r).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).EAMagRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).EAPressureRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).EAInertialRTr(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).ParTemp(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).Temp(nn, Qind)= NaN;
                Rank(r).Specie(s).FluxTube(f).ParFlux(nn, Qind)= NaN;
                
            end
        end
        
        Rank(r).Specie(s).FluxTube(f).IonPressure(nn, Qind)= ...
            2e0*Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)* ...
            Rank(r).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind);
        
        Rank(r).Specie(s).FluxTube(f).ElecPressure(nn, Qind)= ...
            Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)*kB* ...
            Rank(r).Specie(s).FluxTube(f).TeNT(nn);
        
        Rank(r).Specie(s).FluxTube(f).PlasmaPressure(nn, Qind)= ...
            Rank(r).Specie(s).FluxTube(f).IonPressure(nn, Qind)+ ...
            Rank(r).Specie(s).FluxTube(f).ElecPressure(nn, Qind);
                
    end
end

if DensityOutputflag == 1
    datadirInput= horzcat('/Users/robertalbarran/Desktop/FULLDESKTOP/MSSims/Density', simnum, '/KMIODensity', simnum, '/');
    
    DensityOutputRTID= fopen(horzcat(datadirInput, num2str(s), '_', ...
        num2str(f), '_', 'DensityOutputRTfort.bin'));
    Specie(s).FluxTube(f).DensityOutputRTr(:)= fread(DensityOutputRTID, ...
        ((NqUB- NqLB)+ 1), 'real*8');
    fclose(DensityOutputRTID);
    
    nsnormCLBOutputID= fopen(horzcat(datadirInput, num2str(s), '_', ...
        num2str(f), '_', 'nsnormCLBOutputfort.bin'));
    Specie(s).FluxTube(f).nsnormCLBOutput= fread(nsnormCLBOutputID, ...
        1, 'integer*8');
    fclose(nsnormCLBOutputID);
    
    nsnormCUBOutputID= fopen(horzcat(datadirInput, num2str(s), '_', ...
        num2str(f), '_', 'nsnormCUBOutputfort.bin'));
    Specie(s).FluxTube(f).nsnormCUBOutput= fread(nsnormCUBOutputID, ...
        1, 'integer*8');
    fclose(nsnormCUBOutputID);  
    
    if Specie(s).FluxTube(f).DensityOutputRTr(:) ~= Rank(r).Specie(s).FluxTube(f).M0phRTr(end, :)
        disp('DensityOutputRTr not equal to ending M0phRtr value')
    end
end

Rank(r).Specie(s).FluxTube(f).Perp1Temp(:, :)= ...
    2e0.*Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(:, :)./kB; % [K]
Rank(r).Specie(s).FluxTube(f).Perp2Temp(:, :)= ...
    2e0.*Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(:, :)./kB; % [K]
Rank(r).Specie(s).FluxTube(f).PerpTemp(:, :)= ...
    (1e0/2e0).*Rank(r).Specie(s).FluxTube(f).Perp1Temp(:, :)+ ...
    (1e0/2e0).*Rank(r).Specie(s).FluxTube(f).Perp2Temp(:, :);
Rank(r).Specie(s).FluxTube(f).Perp1Flux(:, :)= ...
    Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
    Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(:, :); % [m^2/s]
Rank(r).Specie(s).FluxTube(f).Perp2Flux(:, :)= ...
    Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
    Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(:, :); % [m^2/s]
Rank(r).Specie(s).FluxTube(f).PerpFlux(:, :)= ...
    Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
    sqrt(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(:, :).^2e0+ ...
    Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(:, :).^2e0); % [m^2/s]

Rank(r).Specie(s).FluxTube(f).ParTemp(:, :)= ...
    2e0.*Rank(r).Specie(s).FluxTube(f).M2ParphRTr(:, :)./kB; % [K]
Rank(r).Specie(s).FluxTube(f).Temp(:, :)= ...
    (1e0/3e0)*Rank(r).Specie(s).FluxTube(f).ParTemp(:, :)+ ...
    (2e0/3e0)*Rank(r).Specie(s).FluxTube(f).PerpTemp(:, :);
Rank(r).Specie(s).FluxTube(f).ParFlux(:, :)= ...
    Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
    Rank(r).Specie(s).FluxTube(f).M1ParphRTr(:, :); % [m^2/s]

% Energy, pitch angle, gyroangle from moments:
VPerpM(:, :)= sqrt(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(:, :).^2e0+ ...
    Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(:, :).^2e0);
VParM(:, :)= -Rank(r).Specie(s).FluxTube(f).M1ParphRTr(:, :);
VelM(:, :)= sqrt(VPerpM(:, :).^2e0+ VParM(:, :).^2e0);
EnergyM(:, :)= (mion/2e0).*(VelM(:, :).^2e0);
PitchAngleM(:, :)= abs(90e0- atand(VPerpM(:, :)./VParM(:, :)));
GyroAngleM(:, :)= acosd(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(:, :)./VPerpM(:, :));

% for nn= 1:1:NNtT+ 1
%     for Qind= NqLB:1:NqUB
%         if Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind) == 0e0
%             if ZeroNaNflag == 1
%                 VPerpM(nn, Qind)= 0e0;
%                 VParM(nn, Qind)= 0e0;
%                 VelM(nn, Qind)= 0e0;
%                 EnergyM(nn, Qind)= 0e0;
%                 PitchAngleM(nn, Qind)= 0e0;
%                 GyroAngleM(nn, Qind)= 0e0;
%             else
%                 VPerpM(nn, Qind)= NaN;
%                 VParM(nn, Qind)= NaN;
%                 VelM(nn, Qind)= NaN;
%                 EnergyM(nn, Qind)= NaN;
%                 PitchAngleM(nn, Qind)= NaN;
%                 GyroAngleM(nn, Qind)= NaN;
%             end
%         end
%     end
% end

for Qind= NqLB:1:NqUB
    if (((2e0/3e0).*Rank(root).Specie(s).FluxTube(f).TsPerpp(Qind)+ ...
            (1e0/3e0).*Rank(root).Specie(s).FluxTube(f).TsParp(Qind)) ~= Rank(root).Specie(s).FluxTube(f).Tsp(Qind))
    end
      
    M0Tav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M0phRTr(:, Qind))/(NNtT+ 1e0);
    EAMagTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).EAMagRTr(:, Qind))/(NNtT+ 1e0);
    EAPressureTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).EAPressureRTr(:, Qind))/(NNtT+ 1e0);
    EAInertialTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).EAInertialRTr(:, Qind))/(NNtT+ 1e0);
    M1ParTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M1ParphRTr(:, Qind))/(NNtT+ 1e0);
    
    M1Perp1Tav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M1Perp1phRTr(:, Qind))/(NNtT+ 1e0);
    M1Perp2Tav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M1Perp2phRTr(:, Qind))/(NNtT+ 1e0);
    M2Perp1Tav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M2Perp1phRTr(:, Qind))/(NNtT+ 1e0);
    M2Perp2Tav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M2Perp2phRTr(:, Qind))/(NNtT+ 1e0);
    Perp1TempTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).Perp1Temp(:, Qind))/(NNtT+ 1e0);
    Perp2TempTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).Perp2Temp(:, Qind))/(NNtT+ 1e0);
    Perp1FluxTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).Perp1Flux(:, Qind))/(NNtT+ 1e0);
    Perp2FluxTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).Perp2Flux(:, Qind))/(NNtT+ 1e0);
    PitchAngleMTav(Qind)= sum(PitchAngleM(:, Qind))/(NNtT+ 1e0);
    GyroAngleMTav(Qind)= sum(GyroAngleM(:, Qind))/(NNtT+ 1e0);
    
    UPerpTav(Qind)= ...
        sum(VPerpM(:, Qind))/(NNtT+ 1e0);
    PerpTempTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).PerpTemp(:, Qind))/(NNtT+ 1e0);
    PerpFluxTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).PerpFlux(:, Qind))/(NNtT+ 1e0);
    
    M2ParTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M2ParphRTr(:, Qind))/(NNtT+ 1e0);
    
    if (MomentFilterFlag == 1)
        Rank(1).Specie(1).FluxTube(1).M0FiltAvrgTav(Qind)= ...
            sum(Rank(1).Specie(1).FluxTube(1).M0FiltAvrgRTr(:, Qind))/(NNtT+ 1e0);
        Rank(1).Specie(1).FluxTube(1).UParFiltAvrgTav(Qind)= ...
            sum(Rank(1).Specie(1).FluxTube(1).M1ParFiltAvrgRTr(:, Qind))/(NNtT+ 1e0);
    end
    
    M2Tav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).M2phRTr(:, Qind))/(NNtT+ 1e0);
    TempTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).Temp(:, Qind))/(NNtT+ 1e0);
    ParTempTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).ParTemp(:, Qind))/(NNtT+ 1e0);
    ParFluxTav(Qind)= ...
        sum(Rank(1).Specie(1).FluxTube(1).ParFlux(:, Qind))/(NNtT+ 1e0);
end

toc

disp('MOMENTS IMPORT COMPLETE')

%%

clc
close all
tic

% RUN QUARTIC DIPOLE POLYNOMIAL ROOT-FINDER BOUNDARY ALTITUDES:

for s= 1:1:Stot;
    for f= 1:1:1;
        pdip(:)= Rank(r).Specie(s).FluxTube(f).pGC(:); 
        qdip(:)= Rank(r).Specie(s).FluxTube(f).qGL(:);
        
        for Qind= NqLB:1:NqUB;
            
            Specie(s).FluxTube(f).QCell(Qind).Ap= 0e0; % From quartic form
            Specie(s).FluxTube(f).QCell(Qind).Bp= 0e0;
            Specie(s).FluxTube(f).QCell(Qind).Cp= 1e0/(pdip(Qind)* ...
                qdip(Qind)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).Dp= -1e0/(qdip(Qind)^2e0);
            
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
                (RE*pdip(Qind))));
            Specie(s).FluxTube(f).QCell(Qind).theta2= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r2/...
                (RE*pdip(Qind))));
            Specie(s).FluxTube(f).QCell(Qind).theta3= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r3/ ...
                (RE*pdip(Qind))));
            Specie(s).FluxTube(f).QCell(Qind).theta4= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r4/ ...
                (RE*pdip(Qind))));
            
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
            
            if qdip(Qind)< 0 | ...
                    qdip(Qind)== 0e0; % Phase shift root 3 solution
                Specie(s).FluxTube(f).QCell(Qind).thetafinalG= ...
                    pi- Specie(s).FluxTube(f).QCell(Qind).theta3;
            else if qdip(Qind)> 0;
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
            
            if abs(qdip(Qind)- ...
                    Specie(s).FluxTube(f).QCell(Qind).qfinalG) > 1e-14
                disp('BAD Q')
            end
            
            if abs(pdip(Qind)- Specie(s).FluxTube(f).QCell(Qind).pfinalG) > 1e-13
                disp('BAD P')
            end
            
            Specie(s).FluxTube(f).rfinalG(Qind)= Specie(s).FluxTube(f).QCell(Qind).rfinalG;
            
        end
        
        Rank(r).Specie(s).FluxTube(f).rGL(:)= Specie(s).FluxTube(f).rfinalG(:);
        
        pdip(:)= Rank(r).Specie(s).FluxTube(f).pGC(:); 
        qdip(:)= Rank(r).Specie(s).FluxTube(f).qGH(:);
        
        for Qind= NqLB:1:NqUB;
            
            Specie(s).FluxTube(f).QCell(Qind).Ap= 0e0; % From quartic form
            Specie(s).FluxTube(f).QCell(Qind).Bp= 0e0;
            Specie(s).FluxTube(f).QCell(Qind).Cp= 1e0/(pdip(Qind)* ...
                qdip(Qind)^2e0);
            Specie(s).FluxTube(f).QCell(Qind).Dp= -1e0/(qdip(Qind)^2e0);
            
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
                (RE*pdip(Qind))));
            Specie(s).FluxTube(f).QCell(Qind).theta2= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r2/...
                (RE*pdip(Qind))));
            Specie(s).FluxTube(f).QCell(Qind).theta3= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r3/ ...
                (RE*pdip(Qind))));
            Specie(s).FluxTube(f).QCell(Qind).theta4= ...
                asin(sqrt(Specie(s).FluxTube(f).QCell(Qind).r4/ ...
                (RE*pdip(Qind))));
            
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
            
            if qdip(Qind)< 0 | ...
                    qdip(Qind)== 0e0; % Phase shift root 3 solution
                Specie(s).FluxTube(f).QCell(Qind).thetafinalG= ...
                    pi- Specie(s).FluxTube(f).QCell(Qind).theta3;
            else if qdip(Qind)> 0;
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
            
            if abs(qdip(Qind)- ...
                    Specie(s).FluxTube(f).QCell(Qind).qfinalG) > 1e-14
                disp('BAD Q')
            end
            
            if abs(pdip(Qind)- Specie(s).FluxTube(f).QCell(Qind).pfinalG) > 1e-13
                disp('BAD P')
            end
            
            Specie(s).FluxTube(f).rfinalG(Qind)= Specie(s).FluxTube(f).QCell(Qind).rfinalG;
            
        end
        
        Rank(r).Specie(s).FluxTube(f).rGH(:)= Specie(s).FluxTube(f).rfinalG(:);
        
    end
end

toc

disp('ALTITUDE BOUNDARIES COMPLETE')

%%
clc
close all

tic

% FILTER MOMENTS
r= 1; % Root Rank
if MomentFilterFlag == 1
    for nn= 1:1:NNtT+ 1
        
        % FILTER TEST DENSITY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M0phRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M0FiltAvrgRTr(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M0FiltAvrgRTr(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M0FiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M0FiltAvrgRTr(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PARALLEL VELOCITY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M1ParphRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M1ParFiltAvrgRTr(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M1ParFiltAvrgRTr(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M1ParFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M1ParFiltAvrgRTr(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP1 VELOCITY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M1Perp1phRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M1Perp1FiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M1Perp1FiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M1Perp1FiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M1Perp1FiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP2 VELOCITY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M1Perp2phRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M1Perp2FiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M1Perp2FiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M1Perp2FiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M1Perp2FiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP1 ENERGY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M2Perp1phRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M2Perp1FiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M2Perp1FiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M2Perp1FiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M2Perp1FiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP2 ENERGY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M2Perp2phRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M2Perp2FiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M2Perp2FiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M2Perp2FiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M2Perp2FiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PARALLEL ENERGY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M2ParphRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M2ParFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M2ParFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M2ParFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M2ParFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER TOTAL ENERGY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).M2phRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).M2FiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).M2FiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).M2phRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            M2FiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).M2FiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER INERTIAL EAMB:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).EAInertialRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).EAInertialFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).EAInertialFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).EAInertialRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            EAInertialFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).EAInertialFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PRESSURE EAMB:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).EAPressureRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).EAPressureFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).EAPressureFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).EAPressureRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            EAPressureFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).EAPressureFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER MAGNITUDE EAMB:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).EAMagRTr(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).EAMAgFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).EAMagFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).EAMagRTr(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            EAMagFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).EAMagFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP1 TEMP MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).Perp1Temp(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).Perp1TempFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).Perp1TempFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).Perp1Temp(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            Perp1TempFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).Perp1TempFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP2 TEMP MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).Perp2Temp(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).Perp2TempFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).Perp2TempFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).Perp2Temp(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            Perp2TempFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).Perp2TempFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP TEMP MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).PerpTemp(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).PerpTempFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).PerpTempFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).PerpTemp(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            PerpTempFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).PerpTempFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP1 FLUX MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).Perp1Flux(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).Perp1FluxFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).Perp1FluxFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).Perp1Flux(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            Perp1FluxFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).Perp1FluxFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP2 FLUX MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).Perp2Flux(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).Perp2FluxFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).Perp2FluxFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).Perp2Flux(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            Perp2FluxFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).Perp2FluxFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PERP FLUX MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).PerpFlux(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).PerpFluxFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).PerpFluxFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).PerpFlux(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            PerpFluxFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).PerpFluxFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PAR TEMP MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).ParTemp(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).ParTempFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).ParTempFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).ParTemp(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            ParTempFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).ParTempFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER TOTAL TEMP MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).Temp(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).TempFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).TempFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).Temp(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            TempFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).TempFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PARALLEL FLUX MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= Rank(1).Specie(1).FluxTube(1).ParFlux(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                Rank(r).Specie(s).FluxTube(f).ParFluxFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                Rank(r).Specie(s).FluxTube(f).ParFluxFiltAvrg(nn, Qind)= ...
                    Rank(r).Specie(s).FluxTube(f).ParFlux(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            ParFluxFiltAvrgTav(Qind)= ...
                sum(Rank(r).Specie(s).FluxTube(f).ParFluxFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER TOTAL ENERGY MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= EnergyM(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                EnergyMFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                EnergyMFiltAvrg(nn, Qind)= ...
                    EnergyM(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            EnergyMFiltAvrgTav(Qind)= ...
                sum(EnergyMFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER PITCH ANGLE MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= PitchAngleM(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                PitchAngleMFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                PitchAngleMFiltAvrg(nn, Qind)= ...
                    PitchAngleM(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            PitchAngleMFiltAvrgTav(Qind)= ...
                sum(PitchAngleMFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
        % FILTER GYRO ANGLE MOMENT:
        for Qind= NqLB:1:NqUB
            MomentFiltIn(nn, Qind)= GyroAngleM(nn, Qind);
        end
        for Qind= NqLB:1:NqUB
            if ((Qind >= NqLB) & (Qind <= NqUB))
                
                % Filter low altitude moments
                if (Qind <= round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments below given cell
                    if (Qind == NqLB)
                        FiltSum2= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg1;
                    else if (Qind == NqLB+ 1e0)
                            FiltSum2= MomentFiltIn(nn, NqLB);
                            FiltAvrg2= FiltSum2;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum2= MomentFiltIn(nn, Qind- 1);
                            for MAfilterQind2= 2:1:round(real(Qind- NqLB))
                                FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                            end
                            FiltAvrg2= FiltSum2/round(real(Qind- NqLB));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
                
                % Filter mid altitude moments
                if ((Qind > round(NqLB- 1e0+ (MAfilterPt- 1e0)/2e0)) & ...
                        (Qind < (NqUB- round((MAfilterPt- 1e0)/2e0))))
                    % Average moments above (and including) given cell
                    FiltSum1= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind+ 1);
                    for MAfilterQind1= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                    end
                    % Average moments below given cell
                    FiltSum2= MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg1= FiltSum1/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    FiltAvrg2= FiltSum2/round((MAfilterPt- 1e0)/2e0);
                    MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                end
                
                % Filter high altitude moments
                if (Qind >= (NqUB- round((MAfilterPt- 1e0)/2e0)))
                    % Average moments below (and including) given cell
                    FiltSum2= MomentFiltIn(nn, Qind)+ MomentFiltIn(nn, Qind- 1);
                    for MAfilterQind2= 2:1:round((MAfilterPt- 1e0)/2e0)
                        FiltSum2= FiltSum2+ MomentFiltIn(nn, (Qind- MAfilterQind2));
                    end
                    FiltAvrg2= FiltSum2/(round((MAfilterPt- 1e0)/2e0)+ 1);
                    % Average moments above given cell
                    if (Qind == NqUB)
                        FiltSum1= 0e0;
                        MomentFiltOut(nn, Qind)= FiltAvrg2;
                    else if (Qind == NqUB- 1e0)
                            FiltSum1= MomentFiltIn(nn, NqUB);
                            FiltAvrg1= FiltSum1;
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        else
                            FiltSum1= MomentFiltIn(nn, Qind+ 1);
                            for MAfilterQind1= 2:1:round(real(NqUB- Qind))
                                FiltSum1= FiltSum1+ MomentFiltIn(nn, (Qind+ MAfilterQind1));
                            end
                            FiltAvrg1= FiltSum1/round(real(NqUB- Qind));
                            MomentFiltOut(nn, Qind)= (FiltAvrg1+ FiltAvrg2)/2e0;
                        end
                    end
                end
            else
                MomentFiltOut(nn, Qind)= 0e0;
            end
        end
        for Qind= NqLB:1:NqUB
            if (Qind > NqLB) & (Qind < NqUB)
                GyroAngleMFiltAvrg(nn, Qind)= MomentFiltOut(nn, Qind);
            else
                GyroAngleMFiltAvrg(nn, Qind)= ...
                    GyroAngleM(nn, Qind);
            end
        end
        for Qind= NqLB:1:NqUB;
            GyroAngleMFiltAvrgTav(Qind)= ...
                sum(GyroAngleMFiltAvrg(:, Qind))/(NNtT+ 1e0);
        end
        
    end
end

toc

disp('FILTER MOMENTS COMPLETE')

%%
tic
clc
close all

% PLOT ION MOMENTS:

clear L LH Q0datfacVec
Qindlen= ((NqUB- NqLB)+ 1);
for Qind= 1:1:Qindlen;
    rr(Qind)= Rank(r).Specie(s).FluxTube(f).rGCp(Qind);
    qq(Qind)= Rank(r).Specie(s).FluxTube(f).qGC(Qind);
end
nd(:)= Rank(r).Specie(s).FluxTube(f).TimeT(:);
[rrAB ndAB]= meshgrid(rr, nd);
[qqAB ndAB]= meshgrid(qq, nd);

iiNt= size(nd);

for ii= 1:1:iiNt(1)
    for Qind= NqLB:1:NqUB;
        if ii == 1
            Q0datfacVec(ii, Qind)= 0;
        else
            Q0datfacVec(ii, Qind)= Q0datfacVec(ii- 1, Qind)+ Q0datfac;
        end
    end
end

cbarloc= 'EastOutside'; logAEflag= 0; logTflag= 0; logWflag= 0; fignum= 1; figlabelflag= 0; FIGSAVEflag= 1;
pcolorflag= 1; cmapset= 'jet'; xfigsize= 750; yfigsize= 350; fsize= 20; linewidth= 2;
dxfigsize= 250; dyfigsize= 250; legloc= 'NorthEast';

logPitchAngleflag= 0; logGyroAngleflag= 0;

if IONVPERPVECflag == 1
    
    FPerp1plotflag= 0;
    FPerp2plotflag= 0;
    WPerp1plotflag= 0;
    WPerp2plotflag= 0;
    TPerp1plotflag= 0;
    TPerp2plotflag= 0;
    
else
    
    FPerpplotflag= 0;
    WPerpplotflag= 0;
    TPerpplotflag= 0;
    
end

LAflag= 0;
LATflag= 0;
PC00flag= 0;
VV0flag= 0;
W0flag= 1;
    
Ambplotflag= 1; Betaplotflag= 0;
Nplotflag= 1;
UParplotflag= 1; UPerp12plotflag= 1; UPerpplotflag= 0;
WParplotflag= 1; WPerp12plotflag= 1;
TParplotflag= 1; TPerp12plotflag= 1;
FParplotflag= 0; FPerp12plotflag= 0;
Wplotflag= 1; Tplotflag= 1;
PitchAngleplotflag= 0; GyroAngleplotflag= 0;

UPerplim1= 0e0; UPerplim2= 20e0; UPerplimflag= 0;
FluxPerplim1= 11.5e0; FluxPerplim2= 15e0; FluxPerplimflag= 0;
WPerplim1= 0e0; WPerplim2= 0.3e0; WPerplimflag= 1;
TPerplim1= 0e0; TPerplim2= 5.5e3; TPerplimflag= 1;

Aamblim1= 0e0; Aamblim2= 3e-6; Aamblimflag= 0;
if LAflag == 1
    Nlim1= 8; Nlim2= 11; Nlimflag= 1; % LA sims
end
if LATflag == 1
    Nlim1= 8; Nlim2= 11; Nlimflag= 1; % LAT sims
end
if PC00flag == 1
    Nlim1= 6; Nlim2= 8; Nlimflag= 1; % PC00 sims
end
if VV0flag == 1
    Nlim1= 6; Nlim2= 8; Nlimflag= 0; % VV0 sims
end
if W0flag == 1
    Nlim1= 6; Nlim2= 8; Nlimflag= 0; % W0 sims
end
Betalim1= -12e0; Betalim2= 0e0; Betalimflag= 0;
PitchAnglelim1= 0e0; PitchAnglelim2= 180e0; PitchAnglelimflag= 0;
GyroAnglelim1= 0e0; GyroAnglelim2= 100e0; GyroAnglelimflag= 0;

if LAflag == 1
    UPerp1lim1= -1e0; UPerp1lim2= 1e0; UPerp1limflag= 1;
    UPerp2lim1= -1e0; UPerp2lim2= 1e0; UPerp2limflag= 1;
    FluxPerp1lim1= 11.5e0; FluxPerp1lim2= 15e0; FluxPerp1limflag= 1;
    FluxPerp2lim1= 11.5e0; FluxPerp2lim2= 15e0; FluxPerp2limflag= 1;
    WPerp1lim1= 0e0; WPerp1lim2= 0.1e0; WPerp1limflag= 1;
    WPerp2lim1= 0e0; WPerp2lim2= 0.1e0; WPerp2limflag= 1;
    TPerp1lim1= 0e0; TPerp1lim2= 2e3; TPerp1limflag= 0;
    TPerp2lim1= 0e0; TPerp2lim2= 2e3; TPerp2limflag= 0;
    
    UParlim1= -1e0; UParlim2= 1; UParlimflag= 1;
    FluxParlim1= -15e0; FluxParlim2= 15e0; FluxParlimflag= 1;
    Wlim1= 0e0; Wlim2= 0.3e0; Wlimflag= 0;
    WParlim1= 0e0; WParlim2= 0.1e0; WParlimflag= 0;
    Tlim1= 0e0; Tlim2= 2e3; Tlimflag= 0;
    TParlim1= 0e0; TParlim2= 2e3; TParlimflag= 0;
end
if LATflag == 1
    UPerp1lim1= -1e0; UPerp1lim2= 1e0; UPerp1limflag= 1;
    UPerp2lim1= -1e0; UPerp2lim2= 1e0; UPerp2limflag= 1;
    FluxPerp1lim1= 11.5e0; FluxPerp1lim2= 15e0; FluxPerp1limflag= 1;
    FluxPerp2lim1= 11.5e0; FluxPerp2lim2= 15e0; FluxPerp2limflag= 1;
    WPerp1lim1= 0e0; WPerp1lim2= 0.1e0; WPerp1limflag= 1;
    WPerp2lim1= 0e0; WPerp2lim2= 0.1e0; WPerp2limflag= 1;
    TPerp1lim1= 0e0; TPerp1lim2= 2e3; TPerp1limflag= 1;
    TPerp2lim1= 0e0; TPerp2lim2= 2e3; TPerp2limflag= 1;
    
    UParlim1= -1e0; UParlim2= 1; UParlimflag= 1;
    FluxParlim1= -15e0; FluxParlim2= 15e0; FluxParlimflag= 1;
    Wlim1= 0e0; Wlim2= 0.3e0; Wlimflag= 1;
    WParlim1= 0e0; WParlim2= 0.1e0; WParlimflag= 1;
    Tlim1= 0e0; Tlim2= 2e3; Tlimflag= 1;
    TParlim1= 0e0; TParlim2= 2e3; TParlimflag= 1;
end
if PC00flag == 1
    UPerp1lim1= -1e0; UPerp1lim2= 1e0; UPerp1limflag= 1;
    UPerp2lim1= -1e0; UPerp2lim2= 1e0; UPerp2limflag= 1;
    FluxPerp1lim1= 8e0; FluxPerp1lim2= 16e0; FluxPerp1limflag= 0;
    FluxPerp2lim1= 8e0; FluxPerp2lim2= 16e0; FluxPerp2limflag= 0;
    WPerp1lim1= 0e0; WPerp1lim2= 6e0; WPerp1limflag= 0;
    WPerp2lim1= 0e0; WPerp2lim2= 6e0; WPerp2limflag= 0;
    TPerp1lim1= 0e0; TPerp1lim2= 1.3e5; TPerp1limflag= 0;
    TPerp2lim1= 0e0; TPerp2lim2= 1.3e5; TPerp2limflag= 0;
    
    UParlim1= -10e0; UParlim2= 10e0; UParlimflag= 1;
    FluxParlim1= -15e0; FluxParlim2= 25e0; FluxParlimflag= 0;
    Wlim1= 0e0; Wlim2= 12e0; Wlimflag= 0;
    WParlim1= 0e0; WParlim2= 1e0; WParlimflag= 0;
    Tlim1= 0e0; Tlim2= 1.2e5; Tlimflag= 0;
    TParlim1= 0e0; TParlim2= 1.5e4; TParlimflag= 0;
end
if VV0flag == 1
    UPerp1lim1= -1e0; UPerp1lim2= 1e0; UPerp1limflag= 1;
    UPerp2lim1= -1e0; UPerp2lim2= 1e0; UPerp2limflag= 1;
    FluxPerp1lim1= 11.5e0; FluxPerp1lim2= 15e0; FluxPerp1limflag= 1;
    FluxPerp2lim1= 11.5e0; FluxPerp2lim2= 15e0; FluxPerp2limflag= 1;
    WPerp1lim1= 0e0; WPerp1lim2= 1e0; WPerp1limflag= 1;
    WPerp2lim1= 0e0; WPerp2lim2= 1e0; WPerp2limflag= 1;
    TPerp1lim1= 0e0; TPerp1lim2= 2e4; TPerp1limflag= 1;
    TPerp2lim1= 0e0; TPerp2lim2= 2e4; TPerp2limflag= 1;
    
    UParlim1= -2e0; UParlim2= 2; UParlimflag= 1;
    FluxParlim1= -15e0; FluxParlim2= 15e0; FluxParlimflag= 1;
    Wlim1= 0e0; Wlim2= 2e0; Wlimflag= 1;
    WParlim1= 0e0; WParlim2= 0.3e0; WParlimflag= 1;
    Tlim1= 0e0; Tlim2= 2e4; Tlimflag= 1;
    TParlim1= 0e0; TParlim2= 0.7e4; TParlimflag= 1;
end
if W0flag == 1
    UPerp1lim1= -1e0; UPerp1lim2= 1e0; UPerp1limflag= 1;
    UPerp2lim1= -1e0; UPerp2lim2= 1e0; UPerp2limflag= 1;
    FluxPerp1lim1= 11.5e0; FluxPerp1lim2= 15e0; FluxPerp1limflag= 1;
    FluxPerp2lim1= 11.5e0; FluxPerp2lim2= 15e0; FluxPerp2limflag= 1;
    WPerp1lim1= 0e0; WPerp1lim2= 1e0; WPerp1limflag= 1;
    WPerp2lim1= 0e0; WPerp2lim2= 1e0; WPerp2limflag= 1;
    TPerp1lim1= 0e0; TPerp1lim2= 2e4; TPerp1limflag= 1;
    TPerp2lim1= 0e0; TPerp2lim2= 2e4; TPerp2limflag= 1;
    
    UParlim1= -5e0; UParlim2= 5; UParlimflag= 1;
    FluxParlim1= -15e0; FluxParlim2= 15e0; FluxParlimflag= 1;
    Wlim1= 0e0; Wlim2= 2e0; Wlimflag= 1;
    WParlim1= 0e0; WParlim2= 0.3e0; WParlimflag= 1;
    Tlim1= 0e0; Tlim2= 2e4; Tlimflag= 1;
    TParlim1= 0e0; TParlim2= 0.7e4; TParlimflag= 1;
    
%     UPerp1lim1= -1e0; UPerp1lim2= 1e0; UPerp1limflag= 1;
%     UPerp2lim1= -1e0; UPerp2lim2= 1e0; UPerp2limflag= 1;
%     FluxPerp1lim1= 11.5e0; FluxPerp1lim2= 15e0; FluxPerp1limflag= 1;
%     FluxPerp2lim1= 11.5e0; FluxPerp2lim2= 15e0; FluxPerp2limflag= 1;
%     WPerp1lim1= 0e0; WPerp1lim2= 17e0; WPerp1limflag= 1;
%     WPerp2lim1= 0e0; WPerp2lim2= 17e0; WPerp2limflag= 1;
%     TPerp1lim1= 0e0; TPerp1lim2= 4e5; TPerp1limflag= 1;
%     TPerp2lim1= 0e0; TPerp2lim2= 4e5; TPerp2limflag= 1;
%     
%     UParlim1= -15e0; UParlim2= 15; UParlimflag= 1;
%     FluxParlim1= -15e0; FluxParlim2= 15e0; FluxParlimflag= 1;
%     Wlim1= 0e0; Wlim2= 35e0; Wlimflag= 1;
%     WParlim1= 0e0; WParlim2= 5e0; WParlimflag= 1;
%     Tlim1= 0e0; Tlim2= 3e5; Tlimflag= 1;
%     TParlim1= 0e0; TParlim2= 5e4; TParlimflag= 1;
end

if LAflag == 1 | LATflag == 1
    ylim1= 350e0; % LA 350km
    ylim2= 1e3; % LA 1e3km
end
if PC00flag == 1
    ylim1= 12300e0; %rr(NqLB); % PC00
    ylim2= 14500e0; %rr(NqUB); % PC00
end
if VV0flag == 1
    ylim1= rr(NqLB); % VV0
    ylim2= rr(NqUB); % VV0
end
if W0flag == 1
    ylim1= rr(NqLB); % W0
    ylim2= rr(NqUB); % W0
end

fignum= 0;
% -------------------------------------------------------

if (Ambplotflag == 1) % Filtered Ambipolar Mag Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        if (logAEflag == 1)
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                log10(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr))));
        end
        if (logAEflag == 0)
            contourf(ndAB(:, 1)/60, rrAB(1, :)*1e3/RE, ...
                real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr)));
        end
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logAEflag == 1)
            xlabel(cbar, horzcat('log$_{10}(E_A)$ [V/m]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logAEflag == 0)
            xlabel(cbar, horzcat('$E_A$ [V/m]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Aamblimflag == 1)
            caxis([Aamblim1 Aamblim2]);
        end
        hold off
    else
        if (logAEflag == 1)
            pp= pcolor(ndAB/60, rrAB, ...
                log10(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr))));
        end
        if (logAEflag == 0)
            pp= pcolor(ndAB/60, rrAB, ...
                real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr)));
        end
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logAEflag == 1)
            xlabel(cbar, horzcat('log$_{10}(E_A)$ [V/m]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logAEflag == 0)
            xlabel(cbar, horzcat('$E_A$ [V/m]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Aamblimflag == 1)
            caxis([Aamblim1 Aamblim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    if (logAEflag == 1)
        plot(log10(real(EAMagTav)), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(log10(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(1, :)))), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot(log10(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(end, :)))), rrAB(1, :), 'r', 'LineWidth', linewidth)
    end
    if (logAEflag == 0)
        plot(real(abs(EAMagTav)), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(1, :))), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(end, :))), rrAB(1, :), 'r', 'LineWidth', linewidth)
    end
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    if (logAEflag == 1)
        L{1}= horzcat('log$_{10}(\bar E_A)$');
        L{2}= horzcat('log$_{10}(E_{Ai})$');
        L{3}= horzcat('log$_{10}(E_{Af})$');
    end
    if (logAEflag == 0)
        L{1}= horzcat('$\bar E_A$');
        L{2}= horzcat('$E_{Ai}$');
        L{3}= horzcat('$E_{Af}$');
    end
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    if (logAEflag == 1)
        xlabel(horzcat('log$_{10}(E_A)$ [V/m]'), ...
            'interpreter', 'latex', 'FontSize', 25)
    end
    if (logAEflag == 0)
        xlabel(horzcat('$E_A$ [V/m]'), ...
            'interpreter', 'latex', 'FontSize', 25)
    end
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (Aamblimflag == 1)
        caxis([Aamblim1 Aamblim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        if (logAEflag == 1)
            pp= pcolor(ndAB/60, rrAB, ...
                abs(log10(real(Rank(r).Specie(s).FluxTube(f).EAMagFiltAvrg))));
        end
        if (logAEflag == 0)
            pp= pcolor(ndAB/60, rrAB, ...
                abs(real(Rank(r).Specie(s).FluxTube(f).EAMagFiltAvrg)));
        end
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logAEflag == 1)
            xlabel(cbar, horzcat('log$_{10}(E_{Af})$ [V/m]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logAEflag == 0)
            xlabel(cbar, horzcat('$E_{Af}$ [V/m]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Aamblimflag == 1)
            caxis([Aamblim1 Aamblim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        if (logAEflag == 1)
            plot(log10(real(EAMagFiltAvrgTav)), rrAB, 'b', 'LineWidth', linewidth)
            hold on
            plot(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(1, :))), rrAB(1, :), 'g', 'LineWidth', linewidth)
        end
        if (logAEflag == 0)
            plot(abs(real(EAMagFiltAvrgTav)), rrAB, 'b', 'LineWidth', linewidth)
            hold on
            plot(real(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(1, :))), rrAB(1, :), 'g', 'LineWidth', linewidth)
        end
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        if (logAEflag == 1)
            L{1}= horzcat('log$_{10}(\bar E_{Af})$');
            L{2}= horzcat('log$_{10}(E_{Ai})$');
        end
        if (logAEflag == 0)
            L{1}= horzcat('$\bar E_{Af}$');
            L{2}= horzcat('$E_{Ai}$');
        end
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar E_{Af}$ [V/m]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Aamblimflag == 1)
            caxis([Aamblim1 Aamblim2]);
        end
        hold off
        
    end
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' E AMBIPOLAR MAG CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (PitchAngleplotflag == 1) % Pitch Angle Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    
    subplot(1, 2, 2)
    if pcolorflag == 0
        if (logPitchAngleflag == 1)
            contourf(ndAB(:, 1)/60, rrAB(1, :), real(log10(PitchAngleM')));
        end
        if (logPitchAngleflag == 0)
            contourf(ndAB(:, 1)/60, rrAB(1, :), real(PitchAngleM'));
        end
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logPitchAngleflag == 1)
            xlabel(cbar, horzcat('log$_{10}(\alpha)$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logPitchAngleflag == 0)
            xlabel(cbar, horzcat('$\alpha$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (PitchAnglelimflag == 1)
            caxis([PitchAnglelim1 PitchAnglelim2]);
        end
        hold off
    else
        if (logPitchAngleflag == 1)
            pp= pcolor(ndAB/60, rrAB, ...
                real(log10(PitchAngleM)));
        end
        if (logPitchAngleflag == 0)
            pp= pcolor(ndAB/60, rrAB, ...
                real(PitchAngleM));
        end
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar('EastOutside');
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if (logPitchAngleflag == 1)
            xlabel(cbar, horzcat('log$_{10}(\alpha)$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logPitchAngleflag == 0)
            xlabel(cbar, horzcat('$\alpha$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (PitchAnglelimflag == 1)
            caxis([PitchAnglelim1 PitchAnglelim2]);
        end
        hold off
    end
    
    subplot(1, 2, 1)
    if (logPitchAngleflag == 1)
        plot(real(log10(PitchAngleMTav)), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(real(log10(PitchAngleM(1, :))), rrAB, 'g', 'LineWidth', linewidth)
        plot(real(log10(PitchAngleM(end, :))), rrAB, 'r', 'LineWidth', linewidth)
    end
    if (logPitchAngleflag == 0)
        plot(real(PitchAngleMTav), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(real(PitchAngleM(1, :)), rrAB, 'g', 'LineWidth', linewidth)
        plot(real(PitchAngleM(end, :)), rrAB, 'r', 'LineWidth', linewidth)
    end
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    if (logPitchAngleflag == 1)
        xlabel(horzcat('log$_{10}(\alpha)$ [degs]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        L{1}= horzcat('log$_{10}(\bar \alpha)$');
        L{2}= horzcat('log$_{10}(\alpha_i)$');
        L{3}= horzcat('log$_{10}(\alpha_f)$');
    end
    if (logPitchAngleflag == 0)
        xlabel(horzcat('$\alpha$ [degs]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        L{1}= horzcat('$\bar \alpha$');
        L{2}= horzcat('$\alpha_i$');
        L{3}= horzcat('$\alpha_f$');
    end
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (PitchAnglelimflag == 1)
        xlim([PitchAnglelim1 PitchAnglelim2]);
    end
    hold off
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PITCH ANGLE MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (GyroAngleplotflag == 1) % Gyro Angle Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    
    subplot(1, 2, 2)
    if pcolorflag == 0
        if (logGyroAngleflag == 1)
            contourf(ndAB(:, 1)/60, rrAB(1, :), real(log10(GyroAngleM')));
        end
        if (logGyroAngleflag == 0)
            contourf(ndAB(:, 1)/60, rrAB(1, :), real(GyroAngleM'));
        end
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logPitchAngleflag == 1)
            xlabel(cbar, horzcat('log$_{10}(\varphi)$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logPitchAngleflag == 0)
            xlabel(cbar, horzcat('$\varphi$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (GyroAnglelimflag == 1)
            caxis([GyroAnglelim1 GyroAnglelim2]);
        end
        hold off
    else
        if (logGyroAngleflag == 1)
            pp= pcolor(ndAB/60, rrAB(1, :), real(log10(GyroAngleM)));
        end
        if (logGyroAngleflag == 0)
            pp= pcolor(ndAB/60, rrAB(1, :), real(GyroAngleM));
        end
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar('EastOutside');
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if (logGyroAngleflag == 1)
            xlabel(cbar, horzcat('log$_{10}(\varphi)$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        if (logGyroAngleflag == 0)
            xlabel(cbar, horzcat('$\varphi$ [degs]'), ...
                'interpreter', 'latex', 'FontSize', 25)
        end
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (GyroAnglelimflag == 1)
            caxis([GyroAnglelim1 GyroAnglelim2]);
        end
        hold off
    end
    
    subplot(1, 2, 1)
    if (logGyroAngleflag == 1)
        plot(real(log10(GyroAngleMTav)), rrAB(1, :), 'k', 'LineWidth', linewidth)
        hold on
        plot(real(log10(GyroAngleM(1, :))), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot(real(log10(GyroAngleM(end, :))), rrAB(1, :), 'r', 'LineWidth', linewidth)
    end
    if (logGyroAngleflag == 0)
        plot(real(GyroAngleMTav), rrAB(1, :), 'k', 'LineWidth', linewidth)
        hold on
        plot(real(GyroAngleM(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot(real(GyroAngleM(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
    end
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    if (logGyroAngleflag == 1)
        xlabel(horzcat('log$_{10}(\varphi)$ [degs]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        L{1}= horzcat('log$_{10}(\bar \varphi)$');
        L{2}= horzcat('log$_{10}(\varphi_i)$');
        L{3}= horzcat('log$_{10}(\varphi_f)$');
    end
    if (logGyroAngleflag == 0)
        xlabel(horzcat('$\varphi$ [degs]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        L{1}= horzcat('$\bar \varphi$');
        L{2}= horzcat('$\varphi_i$');
        L{3}= horzcat('$\varphi_f$');
    end
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (GyroAnglelimflag == 1)
        xlim([GyroAnglelim1 GyroAnglelim2]);
    end
    hold off
        
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION GYRO ANGLE MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (Nplotflag == 1) % Density Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), real((Rank(r).Specie(s).FluxTube(f).M0phRTr')));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('log$_{10}(n)$ [m$^{-3}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Nlimflag == 1)
            caxis([Nlim1 Nlim2]);
        end
        hold off
    else
%         pp= pcolor(ndAB/60, rrAB, ...
%            real(log10(Rank(r).Specie(s).FluxTube(f).NqRTp)));
        pp= pcolor(ndAB/60, rrAB, ...
           real(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr)));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar('EastOutside');
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbar, horzcat('log$_{10}(n)$ [m$^{-3}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Nlimflag == 1)
            caxis([Nlim1 Nlim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(log10(M0Tav), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('log$_{10}(\bar n)$');
    L{2}= horzcat('log$_{10}(n_i)$');
    L{3}= horzcat('log$_{10}(n_f)$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
    end
    xlabel(horzcat('log$_{10}(n)$ [m$^{-3}$]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (Nlimflag == 1)
        xlim([Nlim1 Nlim2]);
    end
    hold off
    
    clear LH
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real(log10(Rank(r).Specie(s).FluxTube(f).M0FiltAvrgRTr)));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar('EastOutside');
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbar, horzcat('log$_{10}(n_f)$ [m$^{-3}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Nlimflag == 1)
            caxis([Nlim1 Nlim2]);
        end
        hold off
        
        subplot(2, 2, 3)
        plot(log10(M0FiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('log$_{10}(\bar n_f)$');
        L{2}= horzcat('log$_{10}(n_i)$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('log$_{10}(\bar n_f)$ [m$^{-3}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Nlimflag == 1)
            xlim([Nlim1 Nlim2]);
        end
        hold off
        
    end
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION DENSITY MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (Wplotflag == 1) % Energy Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2phRTr'));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$w$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Wlimflag == 1)
            caxis([Wlim1 Wlim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2phRTr));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$w$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Wlimflag == 1)
            caxis([Wlim1 Wlim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real((6.242e18).*M2Tav), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2phRTr(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar w$');
    L{2}= horzcat('$w_{i}$');
    L{3}= horzcat('$w_{f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$w$ [eV]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (Wlimflag == 1)
        xlim([Wlim1 Wlim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2FiltAvrg));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$w_{f}$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Wlimflag == 1)
            caxis([Wlim1 Wlim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real((6.242e18).*M2FiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar w_{f}$');
        L{2}= horzcat('$w_{i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar w_{f}$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Wlimflag == 1)
            xlim([Wlim1 Wlim2]);
        end
        hold off
        
    end
            
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION ENERGY MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (Tplotflag == 1) % Temperature Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(Rank(r).Specie(s).FluxTube(f).Temp'));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$T$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Tlimflag == 1)
            caxis([Tlim1 Tlim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).Temp));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$T$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Tlimflag == 1)
            caxis([Tlim1 Tlim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real(TempTav), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real(Rank(r).Specie(s).FluxTube(f).Temp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real(Rank(r).Specie(s).FluxTube(f).Temp(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar T$');
    L{2}= horzcat('$T_{i}$');
    L{3}= horzcat('$T_{f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$T$ [K]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (Tlimflag == 1)
        xlim([Tlim1 Tlim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).TempFiltAvrg));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$T_{f}$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Tlimflag == 1)
            caxis([Tlim1 Tlim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real(TempFiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real(Rank(r).Specie(s).FluxTube(f).Temp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar T_{f}$');
        L{2}= horzcat('$T_{i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar T_{f}$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Tlimflag == 1)
            xlim([Tlim1 Tlim2]);
        end
        hold off
        
    end
        
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION TEMP MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (Betaplotflag == 1) % Beta Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    
    subplot(1, 2, 2)
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(log10(Rank(r).Specie(s).FluxTube(f).IonBeta')));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('log$_{10}(\beta)$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Betalimflag == 1)
            caxis([Betalim1 Betalim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(log10(Rank(r).Specie(s).FluxTube(f).IonBeta)));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(d)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbar, horzcat('log$_{10}(\beta)$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Betalimflag == 1)
            caxis([Betalim1 Betalim2]);
        end
        hold off
    end
    
    subplot(1, 2, 1)
    plot(log10(Rank(r).Specie(s).FluxTube(f).IonBetaTav), rrAB, 'k', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(d)', 'interpreter', 'latex', 'FontSize', 25)
    end
    xlabel(horzcat('log$_{10}(\bar \beta)$'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (Betalimflag == 1)
        xlim([Betalim1 Betalim2]);
    end
    hold off
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION BETA MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (UPerpplotflag == 1) % Perpendicular Velocity Plot
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    
    subplot(1, 2, 2)
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(VPerpM'*1e-3));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_\perp$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerplimflag == 1)
            caxis([UPerplim1 UPerplim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(VPerpM*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_\perp$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerplimflag == 1)
            caxis([UPerplim1 UPerplim2]);
        end
        hold off
    end
    
    subplot(1, 2, 1)
    plot(real(UPerpTav*1e-3), rrAB, 'k', 'LineWidth', linewidth)
    hold on
    plot(real(VPerpM(1, :)*1e-3), rrAB, 'g', 'LineWidth', linewidth)
    plot(real(VPerpM(end, :)*1e-3), rrAB, 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar u_\perp$');
    L{2}= horzcat('$u_{\perp i}$');
    L{3}= horzcat('$u_{\perp f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$u_\perp$ [km/s]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (UPerplimflag == 1)
        xlim([UPerplim1 UPerplim2]);
    end
    hold off
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP VEL MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (UPerp12plotflag == 1)  % Perpendicular12 Velocity Plots
        
    % Perpendicular1 Velocity Plots
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr'*1e-3));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\perp 1}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp1limflag == 1)
            caxis([UPerp1lim1 UPerp1lim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\perp 1}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp1limflag == 1)
            caxis([UPerp1lim1 UPerp1lim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real(M1Perp1Tav*1e-3), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(1, :)*1e-3), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(end, :)*1e-3), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar u_{\perp 1}$');
    L{2}= horzcat('$u_{\perp 1 i}$');
    L{3}= horzcat('$u_{\perp 1 f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$u_{\perp 1}$ [km/s]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (UPerp1limflag == 1)
        xlim([UPerp1lim1 UPerp1lim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).M1Perp1FiltAvrg*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\perp 1 f}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp1limflag == 1)
            caxis([UPerp1lim1 UPerp1lim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real(M1Perp1FiltAvrgTav*1e-3), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(1, :)*1e-3), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar u_{\perp 1 f}$');
        L{2}= horzcat('$u_{\perp 1 i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar u_{\perp 1 f}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp1limflag == 1)
            xlim([UPerp1lim1 UPerp1lim2]);
        end
        hold off
        
    end
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP1 VEL MOMENT CONTOUR.png')]);
    end
    
    % Perpendicular2 Velocity Plots
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr'*1e-3));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\perp 2}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp2limflag == 1)
            caxis([UPerp2lim1 UPerp2lim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\perp 2}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp2limflag == 1)
            caxis([UPerp2lim1 UPerp2lim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real(M1Perp2Tav*1e-3), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(1, :)*1e-3), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(end, :)*1e-3), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar u_{\perp 1}$');
    L{2}= horzcat('$u_{\perp 2 i}$');
    L{3}= horzcat('$u_{\perp 2 f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$u_{\perp 2}$ [km/s]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (UPerp2limflag == 1)
        xlim([UPerp2lim1 UPerp2lim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).M1Perp2FiltAvrg*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\perp 2 f}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp2limflag == 1)
            caxis([UPerp2lim1 UPerp2lim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real(M1Perp2FiltAvrgTav*1e-3), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(1, :)*1e-3), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar u_{\perp 2 f}$');
        L{2}= horzcat('$u_{\perp 2 i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar u_{\perp 1 f}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp2limflag == 1)
            xlim([UPerp2lim1 UPerp2lim2]);
        end
        hold off
        
    end
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP2 VEL MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (UParplotflag == 1) % Parallel Velocity Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr'*1e-3));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_\parallel$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            caxis([UParlim1 UParlim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_\parallel$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            caxis([UParlim1 UParlim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real(-M1ParTav*1e-3), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(1, :)*1e-3), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(end, :)*1e-3), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar u_\parallel$');
    L{2}= horzcat('$u_{\parallel i}$');
    L{3}= horzcat('$u_{\parallel f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$u_\parallel$ [km/s]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (UParlimflag == 1)
        xlim([UParlim1 UParlim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real(-Rank(r).Specie(s).FluxTube(f).M1ParFiltAvrgRTr*1e-3));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$u_{\parallel f}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            caxis([UParlim1 UParlim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real(-M1ParFiltAvrgTav*1e-3), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(1, :)*1e-3), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar u_{\parallel f}$');
        L{2}= horzcat('$u_{\parallel i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar u_{\parallel_f}$ [km/s]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            xlim([UParlim1 UParlim2]);
        end
        hold off
        
    end
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PAR VEL MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if IONVPERPVECflag == 1 % Perpendicular Energy Plots
    
    if (WPerp12plotflag == 1)
        
        % Perpendicular1 Energy Plots
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        if (MomentFilterFlag == 1)
            set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
        else
            set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 2)
        else
            subplot(1, 2, 2)
        end
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                (6.242e18).*(Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr'));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$w_{\perp 1}$ [km/s]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp1limflag == 1)
                caxis([WPerp1lim1 WPerp1lim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, ...
                (6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$w_{\perp 1}$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp1limflag == 1)
                caxis([WPerp1lim1 WPerp1lim2]);
            end
            hold off
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 1)
        else
            subplot(1, 2, 1)
        end
        plot((6.242e18).*real(M2Perp1Tav), rrAB(1, :), 'k', 'LineWidth', linewidth)
        hold on
        plot((6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot((6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar w_{\perp 1}$');
        L{2}= horzcat('$w_{\perp 1 i}$');
        L{3}= horzcat('$w_{\perp 1 f}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$w_{\perp 1}$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WPerp1limflag == 1)
            xlim([WPerp1lim1 WPerp1lim2]);
        end
        hold off

        if (MomentFilterFlag == 1)

            subplot(2, 2, 4)
            pp= pcolor(ndAB/60, rrAB, ...
                (6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp1FiltAvrg));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$w_{\perp 1 f}$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp1limflag == 1)
                caxis([WPerp1lim1 WPerp1lim2]);
            end
            hold off

            clear LH
            subplot(2, 2, 3)
            plot((6.242e18).*real(M2Perp1FiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
            hold on
            plot((6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
            set(gca, 'FontSize', fsize);
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
            LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
            L{1}= horzcat('$\bar w_{\perp 1 f}$');
            L{2}= horzcat('$w_{\perp 1 i}$');
            legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
            xlabel(horzcat('$\bar w_{\perp 1 f}$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp1limflag == 1)
                xlim([WPerp1lim1 WPerp1lim2]);
            end
            hold off

        end

        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP1 ENERGY MOMENT CONTOUR.png')]);
        end
    
        % Perpendicular2 Energy Plots
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        if (MomentFilterFlag == 1)
            set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
        else
            set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 2)
        else
            subplot(1, 2, 2)
        end
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                (6.242e18).*(Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr'));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$w_{\perp 2}$ [km/s]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp2limflag == 1)
                caxis([WPerp2lim1 WPerp2lim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, ...
                (6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$w_{\perp 2}$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp2limflag == 1)
                caxis([WPerp2lim1 WPerp2lim2]);
            end
            hold off
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 1)
        else
            subplot(1, 2, 1)
        end
        plot((6.242e18).*real(M2Perp2Tav), rrAB(1, :), 'k', 'LineWidth', linewidth)
        hold on
        plot((6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot((6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar w_{\perp 2}$');
        L{2}= horzcat('$w_{\perp 2 i}$');
        L{3}= horzcat('$w_{\perp 2 f}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$w_{\perp 2}$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WPerp2limflag == 1)
            xlim([WPerp2lim1 WPerp2lim2]);
        end
        hold off

        if (MomentFilterFlag == 1)

            subplot(2, 2, 4)
            pp= pcolor(ndAB/60, rrAB, ...
                (6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp2FiltAvrg));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$w_{\perp 2 f}$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp2limflag == 1)
                caxis([WPerp2lim1 WPerp2lim2]);
            end
            hold off

            clear LH
            subplot(2, 2, 3)
            plot((6.242e18).*real(M2Perp2FiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
            hold on
            plot((6.242e18).*real(Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
            set(gca, 'FontSize', fsize);
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
            LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
            L{1}= horzcat('$\bar w_{\perp 2 f}$');
            L{2}= horzcat('$w_{\perp 2 i}$');
            legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
            xlabel(horzcat('$\bar w_{\perp 2 f}$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerp2limflag == 1)
                xlim([WPerp2lim1 WPerp2lim2]);
            end
            hold off

        end

        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP2 ENERGY MOMENT CONTOUR.png')]);
        end
        
    end
    
else
    
    if (WPerpplotflag == 1) % Perpendicular Energy Plots
        
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        
        subplot(1, 2, 2)
        if pcolorflag == 0
            if (logWflag == 1)
                contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                    real(log10((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr')));
            end
            if (logWflag == 0)
                contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                    real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr'));
            end
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
            end
            if (logWflag == 1)
                xlabel(cbar, horzcat('log$_{10}(w_\perp)$ [eV]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            if (logWflag == 0)
                xlabel(cbar, horzcat('$w_\perp$ [eV]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerplimflag == 1)
                caxis([WPerplim1 WPerplim2]);
            end
            hold off
        else
            if (logWflag == 1)
                pp= pcolor(ndAB/60, rrAB, real(log10((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr)));
            end
            if (logWflag == 0)
                pp= pcolor(ndAB/60, rrAB, real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr));
            end
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(g)', 'interpreter', 'latex', 'FontSize', 25)
            end
            if (logWflag == 1)
                xlabel(cbar, horzcat('log$_{10}(w_\perp)$ [eV]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            if (logWflag == 0)
                xlabel(cbar, horzcat('$w_\perp$ [eV]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (WPerplimflag == 1)
                caxis([WPerplim1 WPerplim2]);
            end
            hold off
        end
        
        subplot(1, 2, 1)
        if (logWflag == 1)
            plot(real(log10((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpTav)), rrAB, 'k', 'LineWidth', linewidth)
            hold on
            plot(real(log10((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(1, :))), rrAB, 'g', 'LineWidth', linewidth)
            plot(real(log10((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(end, :))), rrAB, 'r', 'LineWidth', linewidth)
        end
        if (logWflag == 0)
            plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpTav), rrAB, 'k', 'LineWidth', linewidth)
            hold on
            plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(1, :)), rrAB, 'g', 'LineWidth', linewidth)
            plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2PerpphRTr(end, :)), rrAB, 'r', 'LineWidth', linewidth)
        end
        set(gca, 'FontSize', fsize);
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        if figlabelflag == 1
            text(textx, texty, '(g)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logWflag == 1)
            xlabel(horzcat('log$_{10}(w_\perp)$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            L{1}= horzcat('log$_{10}(\bar w_\perp)$');
            L{2}= horzcat('log$_{10}(w_{\perp i})$');
            L{3}= horzcat('log$_{10}(w_{\perp f})$');
        end
        if (logWflag == 0)
            xlabel(horzcat('$w_\perp$ [eV]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            L{1}= horzcat('$\bar w_\perp$');
            L{2}= horzcat('$w_{\perp i}$');
            L{3}= horzcat('$w_{\perp f}$');
        end
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WPerplimflag == 1)
            xlim([WPerplim1 WPerplim2]);
        end
        hold off
                
    end
    
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP ENERGY MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if (WParplotflag == 1) % Parallel Energy Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2ParphRTr'));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$w_\parallel$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WParlimflag == 1)
            caxis([WParlim1 WParlim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2ParphRTr));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$w_\parallel$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WParlimflag == 1)
            caxis([WParlim1 WParlim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real((6.242e18).*M2ParTav), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2ParphRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2ParphRTr(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar w_\parallel$');
    L{2}= horzcat('$w_{\parallel i}$');
    L{3}= horzcat('$w_{\parallel f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$w_\parallel$ [eV]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (WParlimflag == 1)
        xlim([WParlim1 WParlim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2ParFiltAvrg));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$w_{\parallel f}$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WParlimflag == 1)
            caxis([WParlim1 WParlim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real((6.242e18).*M2ParFiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real((6.242e18).*Rank(r).Specie(s).FluxTube(f).M2ParphRTr(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar w_{\parallel f}$');
        L{2}= horzcat('$w_{\parallel i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar w_{\parallel f}$ [eV]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (WParlimflag == 1)
            xlim([WParlim1 WParlim2]);
        end
        hold off
        
    end
        
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PAR ENERGY MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if IONVPERPVECflag == 1 % Perpendicular Flux Plots
    
    if (FPerp12plotflag == 1)
        
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
        
        subplot(2, 2, 2)
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                real(log10(Rank(r).Specie(s).FluxTube(f).Perp1Flux')));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('log$_{10}(j_{\perp 1})$ [m$^{-2} \cdot $s$^{-1}$]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (FluxPerp1limflag == 1)
                caxis([FluxPerp1lim1 FluxPerp1lim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, real(log10(Rank(r).Specie(s).FluxTube(f).Perp1Flux)));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            xlabel(cbar, horzcat('log$_{10}(j_{\perp 1})$ [m$^{-2} \cdot $s$^{-1}$]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (FluxPerp1limflag == 1)
                caxis([FluxPerp1lim1 FluxPerp1lim2]);
            end
            hold off
        end
        
        subplot(2, 2, 1)
        plot(real(log10(Perp1FluxTav)), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(real(log10(Rank(r).Specie(s).FluxTube(f).Perp1Flux(1, :))), rrAB, 'g', 'LineWidth', linewidth)
        plot(real(log10(Rank(r).Specie(s).FluxTube(f).Perp1Flux(end, :))), rrAB, 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('log$_{10}(\bar j_{\perp 1})$');
        L{2}= horzcat('log$_{10}(j_{\perp 1 i})$');
        L{3}= horzcat('log$_{10}(j_{\perp 1 f})$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('log$_{10}(j_{\perp 1})$ [m$^{-2} \cdot $s$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (FluxPerp1limflag == 1)
            xlim([FluxPerp1lim1 FluxPerp1lim2]);
        end
        hold off
                        
        subplot(2, 2, 4)
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                real(log10(Rank(r).Specie(s).FluxTube(f).Perp2Flux')));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('log$_{10}(j_{\perp 2})$ [m$^{-2} \cdot $s$^{-1}$]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (FluxPerp2limflag == 1)
                caxis([FluxPerp2lim1 FluxPerp2lim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, real(log10(Rank(r).Specie(s).FluxTube(f).Perp2Flux)));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            xlabel(cbar, horzcat('log$_{10}(j_{\perp 2})$ [m$^{-2} \cdot $s$^{-1}$]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (FluxPerp2limflag == 1)
                caxis([FluxPerp2lim1 FluxPerp2lim2]);
            end
            hold off
        end
        
        subplot(2, 2, 3)
        plot(real(log10(Perp2FluxTav)), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(real(log10(Rank(r).Specie(s).FluxTube(f).Perp2Flux(1, :))), rrAB, 'g', 'LineWidth', linewidth)
        plot(real(log10(Rank(r).Specie(s).FluxTube(f).Perp2Flux(end, :))), rrAB, 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('log$_{10}(\bar j_{\perp 2})$');
        L{2}= horzcat('log$_{10}(j_{\perp 2 i})$');
        L{3}= horzcat('log$_{10}(j_{\perp 2 f})$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('log$_{10}(j_{\perp 2})$ [m$^{-2} \cdot $s$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (FluxPerp2limflag == 1)
            xlim([FluxPerp2lim1 FluxPerp2lim2]);
        end
        hold off
                
        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP12 FLUX MOMENT CONTOUR.png')]);
        end
        
    end
    
else
    
    if (FPerpplotflag == 1)
        
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        
        subplot(2, 2, 2)
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                real(log10(Rank(r).Specie(s).FluxTube(f).PerpFlux')));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('log$_{10}(j_\perp)$ [m$^{-2} \cdot $s$^{-1}$]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (FluxPerplimflag == 1)
                caxis([FluxPerplim1 FluxPerplim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, real(log10(Rank(r).Specie(s).FluxTube(f).PerpFlux)));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            xlabel(cbar, horzcat('log$_{10}(j_\perp)$ [m$^{-2} \cdot $s$^{-1}$]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (FluxPerplimflag == 1)
                caxis([FluxPerplim1 FluxPerplim2]);
            end
            hold off
        end
        
        subplot(2, 2, 1)
        plot(real(log10(PerpFluxTav)), rrAB, 'k', 'LineWidth', linewidth)
        hold on
        plot(real(log10(Rank(r).Specie(s).FluxTube(f).PerpFlux(1, :))), rrAB, 'g', 'LineWidth', linewidth)
        plot(real(log10(Rank(r).Specie(s).FluxTube(f).PerpFlux(end, :))), rrAB, 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('log$_{10}(\bar j_\perp)$');
        L{2}= horzcat('log$_{10}(j_\perp i})$');
        L{3}= horzcat('log$_{10}(j_{\perp f})$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('log$_{10}(j_\perp)$ [m$^{-2} \cdot $s$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (FluxPerplimflag == 1)
            xlim([FluxPerplim1 FluxPerplim2]);
        end
        hold off
                
        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP FLUX MOMENT CONTOUR.png')]);
        end
        
    end
    
end

% -------------------------------------------------------

if (FParplotflag == 1) % Parallel Flux Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    
    subplot(1, 2, 2)
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            -sign(Rank(r).Specie(s).FluxTube(f).ParFlux').* ...
            abs(real(log10(Rank(r).Specie(s).FluxTube(f).ParFlux'))));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('log$_{10}(j_\parallel)$ [m$^{-2} \cdot $s$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (FluxParlimflag == 1)
            caxis([FluxParlim1 FluxParlim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, -sign(Rank(r).Specie(s).FluxTube(f).ParFlux).* ...
            abs(real(log10(Rank(r).Specie(s).FluxTube(f).ParFlux))));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(j)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('log$_{10}(j_\parallel)$ [m$^{-2} \cdot $s$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (FluxParlimflag == 1)
            caxis([FluxParlim1 FluxParlim2]);
        end
        hold off
    end
    
    subplot(1, 2, 1)
    plot(-sign(ParFluxTav).* ...
        abs(real(log10(ParFluxTav))), rrAB, 'k', 'LineWidth', linewidth)
    hold on
    plot(-sign(Rank(r).Specie(s).FluxTube(f).ParFlux(1, :)).*abs(real(log10(Rank(r).Specie(s).FluxTube(f).ParFlux(1, :)))), ...
        rrAB, 'g', 'LineWidth', linewidth)
    plot(-sign(Rank(r).Specie(s).FluxTube(f).ParFlux(end, :)).*abs(real(log10(Rank(r).Specie(s).FluxTube(f).ParFlux(end, :)))), ...
        rrAB, 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(j)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('log$_{10}(\bar j_\parallel)$');
    L{2}= horzcat('log$_{10}(j_{\parallel i})$');
    L{3}= horzcat('log$_{10}(j_{\parallel f})$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('log$_{10}(j_\parallel)$ [m$^{-2} \cdot $s$^{-1}$]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (FluxParlimflag == 1)
        xlim([FluxParlim1 FluxParlim2]);
    end
    hold off
        
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PAR FLUX MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

if IONVPERPVECflag == 1 % Perpendicular Temperature Plots
    
    if (TPerp12plotflag == 1)
        
        % Perpendicular1 Temp Plots
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        if (MomentFilterFlag == 1)
            set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
        else
            set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 2)
        else
            subplot(1, 2, 2)
        end
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                (Rank(r).Specie(s).FluxTube(f).Perp1Temp'));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$T_{\perp 1}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp1limflag == 1)
                caxis([TPerp1lim1 TPerp1lim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, ...
                real(Rank(r).Specie(s).FluxTube(f).Perp1Temp));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$T_{\perp 1}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp1limflag == 1)
                caxis([TPerp1lim1 TPerp1lim2]);
            end
            hold off
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 1)
        else
            subplot(1, 2, 1)
        end
        plot(real(Perp1TempTav), rrAB(1, :), 'k', 'LineWidth', linewidth)
        hold on
        plot(real(Rank(r).Specie(s).FluxTube(f).Perp1Temp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot(real(Rank(r).Specie(s).FluxTube(f).Perp1Temp(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar T_{\perp 1}$');
        L{2}= horzcat('$T_{\perp 1 i}$');
        L{3}= horzcat('$T_{\perp 1 f}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$T_{\perp 1}$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TPerp1limflag == 1)
            xlim([TPerp1lim1 TPerp1lim2]);
        end
        hold off

        if (MomentFilterFlag == 1)

            subplot(2, 2, 4)
            pp= pcolor(ndAB/60, rrAB, ...
                real(Rank(r).Specie(s).FluxTube(f).Perp1TempFiltAvrg));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$T_{\perp 1 f}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp1limflag == 1)
                caxis([TPerp1lim1 TPerp1lim2]);
            end
            hold off

            clear LH
            subplot(2, 2, 3)
            plot(real(Perp1TempFiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
            hold on
            plot(real(Rank(r).Specie(s).FluxTube(f).Perp1Temp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
            set(gca, 'FontSize', fsize);
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
            LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
            L{1}= horzcat('$\bar T_{\perp 1 f}$');
            L{2}= horzcat('$T_{\perp 1 i}$');
            legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
            xlabel(horzcat('$\bar T_{\perp 1 f}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp1limflag == 1)
                xlim([TPerp1lim1 TPerp1lim2]);
            end
            hold off

        end

        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP1 TEMP MOMENT CONTOUR.png')]);
        end
                
        % Perpendicular2 Temp Plots
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        if (MomentFilterFlag == 1)
            set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
        else
            set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 2)
        else
            subplot(1, 2, 2)
        end
        if pcolorflag == 0
            contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                (Rank(r).Specie(s).FluxTube(f).Perp2Temp'));
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$T_{\perp 2}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp2limflag == 1)
                caxis([TPerp2lim1 TPerp2lim2]);
            end
            hold off
        else
            pp= pcolor(ndAB/60, rrAB, ...
                real(Rank(r).Specie(s).FluxTube(f).Perp2Temp));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$T_{\perp 2}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp2limflag == 1)
                caxis([TPerp2lim1 TPerp2lim2]);
            end
            hold off
        end

        if (MomentFilterFlag == 1)
            subplot(2, 2, 1)
        else
            subplot(1, 2, 1)
        end
        plot(real(Perp2TempTav), rrAB(1, :), 'k', 'LineWidth', linewidth)
        hold on
        plot(real(Rank(r).Specie(s).FluxTube(f).Perp2Temp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        plot(real(Rank(r).Specie(s).FluxTube(f).Perp2Temp(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar T_{\perp 2}$');
        L{2}= horzcat('$T_{\perp 2 i}$');
        L{3}= horzcat('$T_{\perp 2 f}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$T_{\perp 2}$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TPerp2limflag == 1)
            xlim([TPerp2lim1 TPerp2lim2]);
        end
        hold off

        if (MomentFilterFlag == 1)

            subplot(2, 2, 4)
            pp= pcolor(ndAB/60, rrAB, ...
                real(Rank(r).Specie(s).FluxTube(f).Perp2TempFiltAvrg));
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            xlabel(cbar, horzcat('$T_{\perp 2 f}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp2limflag == 1)
                caxis([TPerp2lim1 TPerp2lim2]);
            end
            hold off

            clear LH
            subplot(2, 2, 3)
            plot(real(Perp2TempFiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
            hold on
            plot(real(Rank(r).Specie(s).FluxTube(f).Perp2Temp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
            set(gca, 'FontSize', fsize);
            if figlabelflag == 1
                text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
            end
            LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
            LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
            L{1}= horzcat('$\bar T_{\perp 2 f}$');
            L{2}= horzcat('$T_{\perp 2 i}$');
            legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
            xlabel(horzcat('$\bar T_{\perp 2 f}$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'xaxisLocation', 'top')
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerp2limflag == 1)
                xlim([TPerp2lim1 TPerp2lim2]);
            end
            hold off

        end

        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP2 TEMP MOMENT CONTOUR.png')]);
        end
        
    end
    
else
    
    if (TPerpplotflag == 1)
        
        fignum= fignum+ 1; % Assign figure number
        fig(fignum)= figure(fignum);
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
        
        subplot(1, 2, 2)
        if pcolorflag == 0
            if (logTflag == 1)
                contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                    real(log10(Rank(r).Specie(s).FluxTube(f).PerpTemp')));
            end
            if (logTflag == 0)
                contourf(ndAB(:, 1)/60, rrAB(1, :), ...
                    real(Rank(r).Specie(s).FluxTube(f).PerpTemp'));
            end
            cbar= colorbar;
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(i)', 'interpreter', 'latex', 'FontSize', 25)
            end
            if (logTflag == 1)
                xlabel(cbar, horzcat('log$_{10}(T_\perp)$ [K]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            if (logTflag == 0)
                xlabel(cbar, horzcat('$T_\perp$ [K]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerplimflag == 1)
                caxis([TPerplim1 TPerplim2]);
            end
            hold off
        else
            if (logTflag == 1)
                pp= pcolor(ndAB/60, rrAB, real(log10(Rank(r).Specie(s).FluxTube(f).PerpTemp)));
            end
            if (logTflag == 0)
                pp= pcolor(ndAB/60, rrAB, real(Rank(r).Specie(s).FluxTube(f).PerpTemp));
            end
            set(pp, 'EdgeColor', 'none');
            shading interp
            colormap(fignum, cmapset);
            cbar= colorbar(cbarloc);
            set(gca, 'FontSize', fsize);
            xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
            if figlabelflag == 1
                text(textx, texty, '(k)', 'interpreter', 'latex', 'FontSize', 25)
            end
            if (logTflag == 1)
                xlabel(cbar, horzcat('log$_{10}(T_\perp)$ [K]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            if (logTflag == 0)
                xlabel(cbar, horzcat('$T_\perp$ [K]'), ...
                    'interpreter', 'latex', 'FontSize', 25)
            end
            set(gca, 'xaxisLocation', 'top')
            ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
            set(gca, 'yminortick', 'on','xminortick', 'on');
            set(gcf, 'color','white');
            ylim([ylim1 ylim2]);
            if (TPerplimflag == 1)
                caxis([TPerplim1 TPerplim2]);
            end
            hold off
        end
        
        subplot(1, 2, 1)
        if (logTflag == 1)
            plot(real(log10(Rank(r).Specie(s).FluxTube(f).PerpTempTav)), rrAB, 'k', 'LineWidth', linewidth)
            hold on
            plot(real(log10(Rank(r).Specie(s).FluxTube(f).PerpTemp(1, :))), rrAB, 'g', 'LineWidth', linewidth)
            plot(real(log10(Rank(r).Specie(s).FluxTube(f).PerpTemp(end, :))), rrAB, 'r', 'LineWidth', linewidth)
        end
        if (logTflag == 0)
            plot(real(Rank(r).Specie(s).FluxTube(f).PerpTempTav), rrAB, 'k', 'LineWidth', linewidth)
            hold on
            plot(real(Rank(r).Specie(s).FluxTube(f).PerpTemp(1, :)), rrAB, 'g', 'LineWidth', linewidth)
            plot(real(Rank(r).Specie(s).FluxTube(f).PerpTemp(end, :)), rrAB, 'r', 'LineWidth', linewidth)
        end
        set(gca, 'FontSize', fsize);
        LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
        if figlabelflag == 1
            text(textx, texty, '(k)', 'interpreter', 'latex', 'FontSize', 25)
        end
        if (logTflag == 1)
            xlabel(horzcat('log$_{10}(T_\perp)$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            L{1}= horzcat('log$_{10}(\bar T_\perp)$');
            L{2}= horzcat('log$_{10}(T_{\perp i})$');
            L{3}= horzcat('log$_{10}(T_{\perp f})$');
        end
        if (logTflag == 0)
            xlabel(horzcat('$T_\perp$ [K]'), ...
                'interpreter', 'latex', 'FontSize', 25)
            L{1}= horzcat('$\bar T_\perp$');
            L{2}= horzcat('$T_{\perp i}$');
            L{3}= horzcat('$T_{\perp f}$');
        end
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TPerplimflag == 1)
            xlim([TPerplim1 TPerplim2]);
        end
        hold off
                
        if FIGSAVEflag == 1
            saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PERP TEMP MOMENT CONTOUR.png')]);
        end
        
    end
    
end

% -------------------------------------------------------

if (TParplotflag == 1) % Parallel Temperature Plots
    
    fignum= fignum+ 1; % Assign figure number
    fig(fignum)= figure(fignum);
    if (MomentFilterFlag == 1)
        set(fig(fignum), 'Position', [10 10 xfigsize+ dxfigsize yfigsize+ dyfigsize])
    else
        set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 2)
    else
        subplot(1, 2, 2)
    end
    if pcolorflag == 0
        contourf(ndAB(:, 1)/60, rrAB(1, :), ...
            real(Rank(r).Specie(s).FluxTube(f).ParTemp'));
        cbar= colorbar;
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$T_\parallel$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TParlimflag == 1)
            caxis([TParlim1 TParlim2]);
        end
        hold off
    else
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).ParTemp));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$T_\parallel$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TParlimflag == 1)
            caxis([TParlim1 TParlim2]);
        end
        hold off
    end
    
    if (MomentFilterFlag == 1)
        subplot(2, 2, 1)
    else
        subplot(1, 2, 1)
    end
    plot(real(ParTempTav), rrAB(1, :), 'k', 'LineWidth', linewidth)
    hold on
    plot(real(Rank(r).Specie(s).FluxTube(f).ParTemp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
    plot(real(Rank(r).Specie(s).FluxTube(f).ParTemp(end, :)), rrAB(1, :), 'r', 'LineWidth', linewidth)
    set(gca, 'FontSize', fsize);
    if figlabelflag == 1
        text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
    end
    LH(1)= plot(nan, nan, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
    LH(3)= plot(nan, nan, 'r', 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\bar T_\parallel$');
    L{2}= horzcat('$T_{\parallel i}$');
    L{3}= horzcat('$T_{\parallel f}$');
    legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
    xlabel(horzcat('$T_\parallel$ [K]'), ...
        'interpreter', 'latex', 'FontSize', 25)
    ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gcf, 'color','white');
    ylim([ylim1 ylim2]);
    if (TParlimflag == 1)
        xlim([TParlim1 TParlim2]);
    end
    hold off
    
    if (MomentFilterFlag == 1)
        
        subplot(2, 2, 4)
        pp= pcolor(ndAB/60, rrAB, ...
            real(Rank(r).Specie(s).FluxTube(f).ParTempFiltAvrg));
        set(pp, 'EdgeColor', 'none');
        shading interp
        colormap(fignum, cmapset);
        cbar= colorbar(cbarloc);
        set(gca, 'FontSize', fsize);
        xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        xlabel(cbar, horzcat('$T_{\parallel f}$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TParlimflag == 1)
            caxis([TParlim1 TParlim2]);
        end
        hold off
        
        clear LH
        subplot(2, 2, 3)
        plot(real(ParTempFiltAvrgTav), rrAB(1, :), 'b', 'LineWidth', linewidth)
        hold on
        plot(real(Rank(r).Specie(s).FluxTube(f).ParTemp(1, :)), rrAB(1, :), 'g', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        if figlabelflag == 1
            text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        end
        LH(1)= plot(nan, nan, 'b', 'MarkerSize', 10, 'LineWidth', 2);
        LH(2)= plot(nan, nan, 'g', 'MarkerSize', 10, 'LineWidth', 2);
        L{1}= horzcat('$\bar T_{\parallel f}$');
        L{2}= horzcat('$T_{\parallel i}$');
        legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
        xlabel(horzcat('$\bar T_{\parallel f}$ [K]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (TParlimflag == 1)
            xlim([TParlim1 TParlim2]);
        end
        hold off
        
    end
        
    if FIGSAVEflag == 1
        saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PAR TEMP MOMENT CONTOUR.png')]);
    end
    
end

% -------------------------------------------------------

disp('ION MOMENT PLOTS COMPLETE')

%%
tic

clc
close all

% IMPORT ION DISTRIBUTION FUNCTIONS
NoiseFilterflag= 0;
NoiseLimNph= 5;  % 10000;
RedNaNflag= 0; % Set to 1 for zero values to NaN

Nbeg= NNtT+ 1;
Nend= NNtT+ 1;

for Vperp1ind= 1:1:NVperp1G
    Vperp1GCp(Vperp1ind)= ...
        Rank(root).Specie(s).FluxTube(f).QCell(1).Vperp1GCp(Vperp1ind, 1, 1);
end
for Vperp2ind= 1:1:NVperp2G
    Vperp2GCp(Vperp2ind)= ...
        Rank(root).Specie(s).FluxTube(f).QCell(1).Vperp2GCp(1, Vperp2ind, 1);
end
for Vparind= 1:1:NVparG
    VparGCp(Vparind)= ...
        Rank(root).Specie(s).FluxTube(f).QCell(1).VparGCp(1, 1, Vparind);
end

for nn= Nbeg:1:Nend
    for Qind= NqLB:1:NqUB
        
        N2PerpphRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
            num2str(f), '_', num2str(Qind), '_', 'N2PerpphRTfort.bin'));
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTf(:, :, :, nn)= ...
            fread(N2PerpphRTID, NVperp1G*NVperp2G*NVparG, 'real*8');
        fclose(N2PerpphRTID);
                
        Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, :, :, nn)= ...
            reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTf(:, :, :, nn), ...
            NVperp1G, NVperp2G, NVparG); % (Vperp2= Vx, Vperp1= Vy, Vpar= Vz)
                
        Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn)= ...
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, :, :, nn)./ ...
            Rank(r).Specie(s).FluxTube(f).QCell(1).d3vCp(:, :, :);
        
        Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(:, :, :, nn)= ...
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, :, :, nn).* ...
            (Rank(root).Specie(s).FluxTube(f).d3xC(Qind)/ ...
            Rank(root).Specie(s).FluxTube(f).nsnormfacT(1));
                
        Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormERRORRTp(:, :, :, nn)= ...
            (Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(:, :, :, nn)).^(-1e0/2e0);
        
        if (NoiseFilterflag == 1)
            for Vperp1ind= 1:1:NVperp1G
                for Vperp2ind= 1:1:NVperp2G
                    for Vparind= 1:1:NVparG
                        if (Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                                N2PerpphReNormRT(Vperp1ind, Vperp2ind, Vparind, nn) < NoiseLimNph)
                            Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                                N2PerpphReNormRT(Vperp1ind, Vperp2ind, Vparind, nn)= 0e0;
                        end
                    end
                end
            end
                    
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, :, :, nn)= ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(:, :, :, nn)./ ...
                (Rank(root).Specie(s).FluxTube(f).d3xC(Qind)/ ...
                Rank(root).Specie(s).FluxTube(f).nsnormfacT(1));
            
            Rank(r).Specie(s).FluxTube(f).NqRTp(nn, Qind)= ...
            sum(sum(sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, :, :, nn))));
        
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn)= ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, :, :, nn)./ ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).d3vCp(:, :, :);
            
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormERRORRTp(:, :, :, nn)= ...
                (Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(:, :, :, nn)).^(-1e0/2e0);
        end
        
        Rank(root).Specie(s).FluxTube(f).N2PerpphReNormRTp(:, :, :, Qind, nn)= ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(:, :, :, nn);
        
        for Vperp1ind= 1:1:NVperp1G
            for Vperp2ind= 1:1:NVperp2G
                for Vparind= 1:1:NVparG
                    if (isinf(Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                            N2PerpphReNormERRORRTp(Vperp1ind, Vperp2ind, Vparind, nn)) == 1) | ...
                            (isnan(Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                            N2PerpphReNormERRORRTp(Vperp1ind, Vperp2ind, Vparind, nn)) == 1)
                        Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                            N2PerpphReNormERRORRTp(Vperp1ind, Vperp2ind, Vparind, nn)= 0e0;
                    end
                end
            end
        end
        
        if QExchangeFlag == 1
            NphENARTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', ...
                num2str(f), '_', num2str(Qind), '_', 'NphENARTfort.bin'));
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphENARTf(:, :, :, nn)= ...
                fread(NphENARTID, NVpG*NVqG*NVphiG, 'real*8');
            fclose(NphENARTID);
            
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphENARTp(:, :, :, nn)= ...
                reshape(Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphENARTf(:, :, :, nn), ...
                NVpG, NVqG, NVphiG);
            
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).FphENARTp(:, :, :, nn)= ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphENARTp(:, :, :, nn)./ ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).d33vCp(:, :, :);
            
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormENART(:, :, :, nn)= ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphENARTp(:, :, :, nn).* ...
                (Rank(root).Specie(s).FluxTube(f).QCell(Qind).d3xC(1)/ ...
                Rank(root).Specie(s).FluxTube(f).nsnormfacT(1));
            
            Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORENARTp(:, :, :, nn)= ...
                (Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormENART(:, :, :, nn)).^(-1e0/2e0);
            
            for Vpind= 1:1:NVpG
                for Vqind= 1:1:NVqG
                    for Vphiind= 1:1:NVphiG
                        if (isinf(Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                                NphReNormERRORENARTp(Vpind, Vqind, Vphiind, nn)) == 1) | ...
                                (isnan(Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                                NphReNormERRORENARTp(Vpind, Vqind, Vphiind, nn)) == 1)
                            Rank(root).Specie(s).FluxTube(f).QCell(Qind). ...
                                NphReNormERRORENARTp(Vpind, Vqind, Vphiind, nn)= 0e0;
                        end
                    end
                end
            end
            
        end
    end
    
    Rank(root).Specie(s).FluxTube(f).TotalParticles(nn)= ...
        sum(sum(sum(sum(Rank(root).Specie(s).FluxTube(f).N2PerpphReNormRTp(:, :, :, :, nn)))));
end

% Compute noise-filtered moments:
if (NoiseFilterflag == 1)
    for Qind= NqLB:1:NqUB
        for nn= Nbeg:1:Nend
                        
            ggg0(:, :, :, nn)= ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            ggg1Perp1(:, :, :, nn)= ...
                Rank(1).Specie(1).FluxTube(1).QCell(1).Vperp1GCp(:, :, :).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            ggg1Perp2(:, :, :, nn)= ...
                Rank(1).Specie(1).FluxTube(1).QCell(1).Vperp2GCp(:, :, :).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            ggg1Par(:, :, :, nn)= ...
                Rank(1).Specie(1).FluxTube(1).QCell(1).VparGCp(:, :, :).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
                                                                            
            ggg0Int(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg0(:, :, :, nn);
            ggg1Perp1Int(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg1Perp1(:, :, :, nn);
            ggg1Perp2Int(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg1Perp2(:, :, :, nn);
            ggg1ParInt(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg1Par(:, :, :, nn);
                        
            ggg0Sum(nn)= sum(sum(sum(ggg0Int(:, :, :, nn))));
            ggg1Perp1Sum(nn)= sum(sum(sum(ggg1Perp1Int(:, :, :, nn))));
            ggg1Perp2Sum(nn)= sum(sum(sum(ggg1Perp2Int(:, :, :, nn))));
            ggg1ParSum(nn)= sum(sum(sum(ggg1ParInt(:, :, :, nn))));
            
            Rank(root).Specie(s).FluxTube(f).M0phRTr(nn, Qind)= ggg0Sum(nn);
            Rank(root).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)= ggg1Perp1Sum(nn)/ggg0Sum(nn);
            Rank(root).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)= ggg1Perp2Sum(nn)/ggg0Sum(nn);
            Rank(root).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)= ggg1ParSum(nn)/ggg0Sum(nn);
                        
            ggg2(:, :, :, nn)= ...
                (Rank(1).Specie(1).FluxTube(1).QCell(1).Vperp1GCp(:, :, :).^2e0+ ...
                Rank(1).Specie(1).FluxTube(1).QCell(1).Vperp2GCp(:, :, :).^2e0+ ...
                (Rank(1).Specie(1).FluxTube(1).QCell(1).VparGCp(:, :, :)- ...
                Rank(root).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)).^2e0).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            ggg2Perp1(:, :, :, nn)= ...
                (Rank(1).Specie(1).FluxTube(1).QCell(1).Vperp1GCp(:, :, :).^2e0).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            ggg2Perp2(:, :, :, nn)= ...
                (Rank(1).Specie(1).FluxTube(1).QCell(1).Vperp2GCp(:, :, :).^2e0).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            ggg2Par(:, :, :, nn)= ...
                ((Rank(1).Specie(1).FluxTube(1).QCell(1).VparGCp(:, :, :)- ...
                Rank(root).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)).^2e0).* ...
                Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
            
            ggg2Int(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg2(:, :, :, nn);
            ggg2Perp1Int(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg2Perp1(:, :, :, nn);
            ggg2Perp2Int(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg2Perp2(:, :, :, nn);
            ggg2ParInt(:, :, :, nn)= ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp1Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVperp2Cp(:, :, :).* ...
                Rank(r).Specie(s).FluxTube(f).QCell(1).dVparCp(:, :, :).* ...
                ggg2Par(:, :, :, nn);
            
            ggg2Sum(nn)= sum(sum(sum(ggg2Int(:, :, :, nn))));
            ggg2Perp1Sum(nn)= sum(sum(sum(ggg2Perp1Int(:, :, :, nn))));
            ggg2Perp2Sum(nn)= sum(sum(sum(ggg2Perp2Int(:, :, :, nn))));
            ggg2ParSum(nn)= sum(sum(sum(ggg2ParInt(:, :, :, nn))));
            
            Rank(root).Specie(s).FluxTube(f).M2phRTr(nn, Qind)= ...
                (mion/2e0)*ggg2Sum(nn)/ggg0Sum(nn);
            Rank(root).Specie(s).FluxTube(f).M2Perp1phRTr(nn, Qind)= ...
                (mion/2e0)*ggg2Perp1Sum(nn)/ggg0Sum(nn);
            Rank(root).Specie(s).FluxTube(f).M2Perp2phRTr(nn, Qind)= ...
                (mion/2e0)*ggg2Perp2Sum(nn)/ggg0Sum(nn);
            Rank(root).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind)= ...
                (mion/2e0)*ggg2ParSum(nn)/ggg0Sum(nn);
            
            if (Rank(root).Specie(s).FluxTube(f).M0phRTr(nn, Qind) == 0e0) 
                Rank(root).Specie(s).FluxTube(f).M0phRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M2phRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M2Perp1phRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M2Perp2phRTr(nn, Qind)= NaN;
                Rank(root).Specie(s).FluxTube(f).M2ParphRTr(nn, Qind)= NaN;
            end
        end        
        
        M0Tav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M0phRTr(:, Qind))/(NNtT+ 1);
        M1Perp1Tav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M1Perp1phRTr(:, Qind))/(NNtT+ 1);
        M1Perp2Tav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M1Perp2phRTr(:, Qind))/(NNtT+ 1);
        M1ParTav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M1ParphRTr(:, Qind))/(NNtT+ 1);
        M2Tav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M2phRTr(:, Qind))/(NNtT+ 1);
        M2Perp1Tav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M2Perp1phRTr(:, Qind))/(NNtT+ 1);
        M2Perp2Tav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M2Perp2phRTr(:, Qind))/(NNtT+ 1);
        M2ParTav(Qind)= ...
            sum(Rank(root).Specie(s).FluxTube(f).M2ParphRTr(:, Qind))/(NNtT+ 1);
    end
    
    Rank(r).Specie(s).FluxTube(f).Perp1Temp(:, :)= ...
        2e0.*Rank(r).Specie(s).FluxTube(f).M2Perp1phRTr(:, :)./kB; % [K]
    Rank(r).Specie(s).FluxTube(f).Perp2Temp(:, :)= ...
        2e0.*Rank(r).Specie(s).FluxTube(f).M2Perp2phRTr(:, :)./kB; % [K]
    Rank(r).Specie(s).FluxTube(f).PerpTemp(:, :)= ...
        (1e0/2e0).*Rank(r).Specie(s).FluxTube(f).Perp1Temp(:, :)+ ...
        (1e0/2e0).*Rank(r).Specie(s).FluxTube(f).Perp2Temp(:, :);
    Rank(r).Specie(s).FluxTube(f).Perp1Flux(:, :)= ...
        Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
        Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(:, :); % [m^2/s]
    Rank(r).Specie(s).FluxTube(f).Perp2Flux(:, :)= ...
        Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
        Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(:, :); % [m^2/s]
    Rank(r).Specie(s).FluxTube(f).PerpFlux(:, :)= ...
        Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
        sqrt(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(:, :).^2e0+ ...
        Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(:, :).^2e0); % [m^2/s]

    Rank(r).Specie(s).FluxTube(f).ParTemp(:, :)= ...
        2e0.*Rank(r).Specie(s).FluxTube(f).M2ParphRTr(:, :)./kB; % [K]
    Rank(r).Specie(s).FluxTube(f).Temp(:, :)= ...
        (1e0/3e0)*Rank(r).Specie(s).FluxTube(f).ParTemp(:, :)+ ...
        (2e0/3e0)*Rank(r).Specie(s).FluxTube(f).PerpTemp(:, :);
    Rank(r).Specie(s).FluxTube(f).ParFlux(:, :)= ...
        Rank(r).Specie(s).FluxTube(f).M0phRTr(:, :).* ...
        Rank(r).Specie(s).FluxTube(f).M1ParphRTr(:, :); % [m^2/s]
end

% Compute reduced distribution functions:
for nn= Nbeg:1:Nend
    for Qind= NqLB:1:NqUB;
        for Vparind= 1:1:NVparG
            FLogVparRed(Vparind, Qind, nn)= ...
                sum(sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, Vparind, nn)));
        end
        for Vperp1ind= 1:1:NVperp1G
            FLogVperp1Red(Vperp1ind, Qind, nn)= ...
                sum(sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(Vperp1ind, :, :, nn)));
        end
        for Vperp2ind= 1:1:NVperp2G
            FLogVperp2Red(Vperp2ind, Qind, nn)= ...
                sum(sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, Vperp2ind, :, nn)));
        end
        for Vperp1ind= 1:1:NVperp1G
            for Vparind= 1:1:NVparG
                NLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(Vperp1ind, :, Vparind, nn));
                FLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(Vperp1ind, :, Vparind, nn));
                NReNormLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= ...
                    sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(Vperp1ind, :, Vparind, nn));
                NReNormERRORLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= ...
                    (NReNormLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn))^(-1e0/2e0);
                if (RedNaNflag == 1)
                    if (NLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn) == 0e0)
                        NLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= NaN;
                    end
                    if (FLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn) == 0e0)
                        FLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= NaN;
                    end
                    if (NReNormLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn) == 0e0)
                        NReNormLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= NaN;
                    end
                    if (NReNormERRORLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn) == 0e0)
                        NReNormERRORLogVperp1VparRed(Vperp1ind, Vparind, Qind, nn)= NaN;
                    end
                end
            end
        end
        for Vperp2ind= 1:1:NVperp2G
            for Vparind= 1:1:NVparG
                NLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(:, Vperp2ind, Vparind, nn));
                FLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, Vperp2ind, Vparind, nn));
                NReNormLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= ...
                    sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(:, Vperp2ind, Vparind, nn));
                NReNormERRORLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= ...
                    (NReNormLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn))^(-1e0/2e0);
                if (RedNaNflag == 1)
                    if (NLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn) == 0e0)
                        NLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= NaN;
                    end
                    if (FLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn) == 0e0)
                        FLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= NaN;
                    end
                    if (NReNormLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn) == 0e0)
                        NReNormLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= NaN;
                    end
                    if (NReNormERRORLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn) == 0e0)
                        NReNormERRORLogVperp2VparRed(Vperp2ind, Vparind, Qind, nn)= NaN;
                    end
                end
                
            end
        end
        for Vperp1ind= 1:1:NVperp1G
            for Vperp2ind= 1:1:NVperp2G
                NLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).N2PerpphRTp(Vperp1ind, Vperp2ind, :, nn));
                FLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(Vperp1ind, Vperp2ind, :, nn));
                NReNormLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= ...
                    sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).N2PerpphReNormRT(Vperp1ind, Vperp2ind, :, nn));
                NReNormERRORLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= ...
                    (NReNormLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn))^(-1e0/2e0);
                if (RedNaNflag == 1)
                    if (NLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn) == 0e0)
                        NLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= NaN;
                    end
                    if (FLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn) == 0e0)
                        FLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= NaN;
                    end
                    if (NReNormLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn) == 0e0)
                        NReNormLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= NaN;
                    end
                    if (NReNormERRORLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn) == 0e0)
                        NReNormERRORLogVperp1Vperp2Red(Vperp1ind, Vperp2ind, Qind, nn)= NaN;
                    end
                end
            end
        end
        if QExchangeFlag == 1
            for Vpind= 1:1:NVpLog
                for Vphiind= 1:1:NVphiLog
                    NLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphENARTp(Vpind, :, Vphiind, nn));
                    FLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FphENARTp(Vpind, :, Vphiind, nn));
                    NReNormLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= ...
                        sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphReNormENART(Vpind, :, Vphiind, nn));
                    NReNormERRORLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORENARTp(Vpind, :, Vphiind, nn));
                    if (RedNaNflag == 1)
                        if (NLogVpVphiENARed(Vpind, Vphiind, Qind, nn) == 0e0)
                            NLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= NaN;
                        end
                        if (FLogVpVphiENARed(Vpind, Vphiind, Qind, nn) == 0e0)
                            FLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= NaN;
                        end
                        if (NReNormLogVpVphiENARed(Vpind, Vphiind, Qind, nn) == 0e0)
                            NReNormLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= NaN;
                        end
                        if (NReNormERRORLogVpVphiENARed(Vpind, Vphiind, Qind, nn) == 0e0)
                            NReNormERRORLogVpVphiENARed(Vpind, Vphiind, Qind, nn)= NaN;
                        end
                    end
                end
            end
            for Vqind= 1:1:NVqLog
                for Vphiind= 1:1:NVphiLog
                    NLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphENARTp(:, Vqind, Vphiind, nn));
                    FLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FphENARTp(:, Vqind, Vphiind, nn));
                    NReNormLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= ...
                        sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphReNormENART(:, Vqind, Vphiind, nn));
                    NReNormERRORLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORENARTp(:, Vqind, Vphiind, nn));
                    if (RedNaNflag == 1)
                        if (NLogVqVphiENARed(Vqind, Vphiind, Qind, nn) == 0e0)
                            NLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= NaN;
                        end
                        if (FLogVqVphiENARed(Vqind, Vphiind, Qind, nn) == 0e0)
                            FLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= NaN;
                        end
                        if (NReNormLogVqVphiENARed(Vqind, Vphiind, Qind, nn) == 0e0)
                            NReNormLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= NaN;
                        end
                        if (NReNormERRORLogVqVphiENARed(Vqind, Vphiind, Qind, nn) == 0e0)
                            NReNormERRORLogVqVphiENARed(Vqind, Vphiind, Qind, nn)= NaN;
                        end
                    end
                    
                end
            end
            for Vpind= 1:1:NVpLog
                for Vqind= 1:1:NVqLog
                    NLogVpVqENARed(Vpind, Vqind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphENARTp(Vpind, Vqind, :, nn));
                    FLogVpVqENARed(Vpind, Vqind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FphENARTp(Vpind, Vqind, :, nn));
                    NReNormLogVpVqENARed(Vpind, Vqind, Qind, nn)= ...
                        sum(Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphReNormENART(Vpind, Vqind, :, nn));
                    NReNormERRORLogVpVqENARed(Vpind, Vqind, Qind, nn)= ...
                        sum(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORENARTp(Vpind, Vqind, :, nn));
                    if (RedNaNflag == 1)
                        if (NLogVpVqENARed(Vpind, Vqind, Qind, nn) == 0e0)
                            NLogVpVqENARed(Vpind, Vqind, Qind, nn)= NaN;
                        end
                        if (FLogVpVqENARed(Vpind, Vqind, Qind, nn) == 0e0)
                            FLogVpVqENARed(Vpind, Vqind, Qind, nn)= NaN;
                        end
                        if (NReNormLogVpVqENARed(Vpind, Vqind, Qind, nn) == 0e0)
                            NReNormLogVpVqENARed(Vpind, Vqind, Qind, nn)= NaN;
                        end
                        if (NReNormERRORLogVpVqENARed(Vpind, Vqind, Qind, nn) == 0e0)
                            NReNormERRORLogVpVqENARed(Vpind, Vqind, Qind, nn)= NaN;
                        end
                    end
                end
            end
        end
    end
end

toc

disp('ION DISTRIBUTION FUNCTIONS IMPORT COMPLETE')

%%

clc
close all

tic

% COMPUTE PITCH ANGLE DISTRIBUTIONS AND DIFFERENTIAL NUMBER AND ENERGY
% FLUXES:

NEG= 1d2;
NalphaG= 1d2;
NthetaG= 1d2;

clear EnergyG alphaG thetaG Vperp2E Vperp1E VparE
clear EnergyRTp VelRTp VelPerp1RTp VelPerp2RTp VelPerpRTp VelParRTp AlphaRTp ThetaRTp
for Vperp1ind= 1:1:NVperp1G
    for Vperp2ind= 1:1:NVperp2G
        for Vparind= 1:1:NVparG
            % Energy Grid [J]
            EnergyRTp(Vperp1ind, Vperp2ind, Vparind)= ...
                (mion/2e0)*(Vperp1GCp(Vperp1ind)^2e0+ ...
                Vperp2GCp(Vperp2ind)^2e0+ VparGCp(Vparind)^2e0);
            % Velocity Grid [m/s]
            VelRTp(Vperp1ind, Vperp2ind, Vparind)= ...
                sqrt(Vperp1GCp(Vperp1ind)^2e0+ Vperp2GCp(Vperp2ind)^2e0+ ...
                VparGCp(Vparind)^2e0);
            % Perp1 Velocity Grid [m/s]
            VelPerp1RTp(Vperp1ind, Vperp2ind, Vparind)= Vperp1GCp(Vperp1ind);
            % Perp2 Velocity Grid [m/s]
            VelPerp2RTp(Vperp1ind, Vperp2ind, Vparind)= Vperp2GCp(Vperp2ind);
            % Perp Velocity Grid [m/s]
            VelPerpRTp(Vperp1ind, Vperp2ind, Vparind)= ...
                sqrt(Vperp1GCp(Vperp1ind)^2e0+ Vperp2GCp(Vperp2ind)^2e0);
            % Par Velocity Grid [m/s]
            VelParRTp(Vperp1ind, Vperp2ind, Vparind)= VparGCp(Vparind);
            % Pitch Angle Grid [degs]
            AlphaRTp(Vperp1ind, Vperp2ind, Vparind)= ...
                abs(90d0- atand(VelPerpRTp(Vperp1ind, Vperp2ind, Vparind)/ ...
                VelParRTp(Vperp1ind, Vperp2ind, Vparind)));
            % Gyro Angle Grid [degs]
            ThetaRTp(Vperp1ind, Vperp2ind, Vparind)= ...
                acosd(VelPerp2RTp(Vperp1ind, Vperp2ind, Vparind)/ ...
                VelPerpRTp(Vperp1ind, Vperp2ind, Vparind));
        end
    end
end

% Set query points:
EnergyG= linspace(min(min(min(EnergyRTp(:, :, :)))), max(max(max(EnergyRTp(:, :, :)))), NEG);
alphaG= linspace(min(min(min(AlphaRTp(:, :, :)))), max(max(max(AlphaRTp(:, :, :)))), NalphaG);
thetaG= linspace(min(min(min(ThetaRTp(:, :, :)))), max(max(max(ThetaRTp(:, :, :)))), NthetaG);
[EnergyGG alphaGG thetaGG]= ndgrid(EnergyG, alphaG, thetaG);
Vperp2E(:, :, :)= sqrt(2e0.*EnergyGG(:, :, :)./mion).*sind(alphaGG(:, :, :)).*cosd(thetaGG(:, :, :));
Vperp1E(:, :, :)= sqrt(2e0.*EnergyGG(:, :, :)./mion).*sind(alphaGG(:, :, :)).*sind(thetaGG(:, :, :));
VparE(:, :, :)= sqrt(2e0.*EnergyGG(:, :, :)./mion).*cosd(alphaGG(:, :, :));

% Set known points:
[Vperp1GCm1 Vperp2GCm1 VparGCm1]= ndgrid(Vperp1GCp, Vperp2GCp, VparGCp);

clear FinterpIn FinterpOut FEInterp
% Interpolate vel distribs:
for nn= Nend
    for Qind= NqLB:1:NqUB
        FinterpIn(:, :, :)= Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpphRTp(:, :, :, nn);
        FinterpOut(:, :, :)= interpn(Vperp1GCm1, Vperp2GCm1, VparGCm1, FinterpIn, Vperp1E, Vperp2E, VparE);
        FEInterp(:, :, :, Qind, nn)= FinterpOut(:, :, :);
    end
end
for nn= Nend
    for Qind= NqLB:1:NqUB
        for Eind= 1:1:NEG
            for alphaind= 1:1:NalphaG
                for thetaind= 1:1:NthetaG
                    if isnan(FEInterp(Eind, alphaind, thetaind, Qind, nn)) == 1
                        FEInterp(Eind, alphaind, thetaind, Qind, nn)= 0e0;
                    end
                end
            end
        end
    end
end

clear FEnergy DNF DEF
for nn= Nend
    for Qind= NqLB:1:NqUB
        % Pitch Angle Distributions [m^-3 eV^-1 sr^-1]
        FEnergy(:, :, :, Qind, nn)= ...
            (1e0/mion).*(sqrt(2e0.*EnergyGG(:, :, :)./mion)).* ...
            FEInterp(:, :, :, Qind, nn);
        % Differential Number Flux [m^-2 s^-1 eV^-1 sr^-1]
        DNF(:, :, :, Qind, nn)= ...
            (2e0.*EnergyGG(:, :, :)./(mion^2e0)).* ...
            FEInterp(:, :, :, Qind, nn);
        % Differential Energy Flux [eV m^-2 s^-1 eV^-1 sr^-1]
        DEF(:, :, :, Qind, nn)= ...
            (2e0.*(EnergyGG(:, :, :).^2e0)./(mion^2e0)).* ...
            FEInterp(:, :, :, Qind, nn);
    end
end

% Compute Reduced Distributions
clear FEnergyEnergyAlphaRed DNFEnergyAlphaRed DEFEnergyAlphaRed ...
    FEnergyEnergyThetaRed DNFEnergyThetaRed DEFEnergyThetaRed ...
    FEnergyAlphaThetaRed DNFAlphaThetaRed DEFAlphaThetaRed

for nn= Nend
    for Qind= NqLB:1:NqUB
        for Eind= 1:1:NEG
            for alphaind= 1:1:NalphaG
                FEEnergyAlphaRed(Eind, alphaind, Qind, nn)= ...
                    sum(FEInterp(Eind, alphaind, :, Qind, nn));
                FEnergyEnergyAlphaRed(Eind, alphaind, Qind, nn)= ...
                    sum(FEnergy(Eind, alphaind, :, Qind, nn));
                DNFEnergyAlphaRed(Eind, alphaind, Qind, nn)= ...
                    sum(DNF(Eind, alphaind, :, Qind, nn));
                DEFEnergyAlphaRed(Eind, alphaind, Qind, nn)= ...
                    sum(DEF(Eind, alphaind, :, Qind, nn));
            end
        end
        for Eind= 1:1:NEG
            for thetaind= 1:1:NthetaG
                FEEnergyThetaRed(Eind, thetaind, Qind, nn)= ...
                    sum(FEInterp(Eind, :, thetaind, Qind, nn));
                FEnergyEnergyThetaRed(Eind, thetaind, Qind, nn)= ...
                    sum(FEnergy(Eind, :, thetaind, Qind, nn));
                DNFEnergyThetaRed(Eind, thetaind, Qind, nn)= ...
                    sum(DNF(Eind, :, thetaind, Qind, nn));
                DEFEnergyThetaRed(Eind, thetaind, Qind, nn)= ...
                    sum(DEF(Eind, :, thetaind, Qind, nn));
            end
        end
        for alphaind= 1:1:NalphaG
            for thetaind= 1:1:NthetaG
                FEAlphaThetaRed(alphaind, thetaind, Qind, nn)= ...
                    sum(FEInterp(:, alphaind, thetaind, Qind, nn));
                FEnergyAlphaThetaRed(alphaind, thetaind, Qind, nn)= ...
                    sum(FEnergy(:, alphaind, thetaind, Qind, nn));
                DNFAlphaThetaRed(alphaind, thetaind, Qind, nn)= ...
                    sum(DNF(:, alphaind, thetaind, Qind, nn));
                DEFAlphaThetaRed(alphaind, thetaind, Qind, nn)= ...
                    sum(DEF(:, alphaind, thetaind, Qind, nn));
            end
        end
    end
end

toc

disp('PITCH ANGLE DISTRIBUTIONS AND DIFFERENTIAL FLUXES COMPLETE')

%%

clc
close all

tic

% COMPUTE LINEAR INTERPOLATIONS OF DISTRIBUTION FUNCTIONS:

% Set known points:
[VparGCm1 Vperp1GCm1]= meshgrid(VparGCp, Vperp1GCp);
[VparGCm2 Vperp2GCm2]= meshgrid(VparGCp, Vperp2GCp);
[Vperp2GCm3 Vperp1GCm3]= meshgrid(Vperp2GCp, Vperp1GCp);

% Set query points:
NVperp1GCq= 1e2;
NVperp2GCq= 1e2;
NVparGCq= 1e2;
Vperp1GCqq= linspace(Vperp1GCp(1), Vperp1GCp(end), NVperp1GCq);
Vperp2GCqq= linspace(Vperp2GCp(1), Vperp2GCp(end), NVperp2GCq);
VparGCqq= linspace(VparGCp(1), VparGCp(end), NVparGCq);

[VparGCq1 Vperp1GCq1]= meshgrid(VparGCqq, Vperp1GCqq);
[VparGCq2 Vperp2GCq2]= meshgrid(VparGCqq, Vperp2GCqq);
[Vperp2GCq3 Vperp1GCq3]= meshgrid(Vperp2GCqq, Vperp1GCqq);

for nn= 1:1:Nend
    for Qind= NqLB:1:NqUB
        
        % Interpolate Raw Ion Counts:
        FinterpIn(:, :)= NReNormLogVperp1VparRed(:, :, Qind, nn);
        FinterpOut(:, :)= interpn(Vperp1GCm1, VparGCm1, FinterpIn, Vperp1GCq1, VparGCq1);
        NReNormLogVperp1VparRedInterp(:, :, Qind, nn)= FinterpOut(:, :);
        NReNormERRORLogVperp1VparRedInterp(:, :, Qind, nn)= ...
            NReNormLogVperp1VparRedInterp(:, :, Qind, nn).^(-1e0/2e0);
        
        FinterpIn(:, :)= NReNormLogVperp2VparRed(:, :, Qind, nn);
        FinterpOut(:, :)= interpn(Vperp2GCm2, VparGCm2, FinterpIn, Vperp2GCq2, VparGCq2);
        NReNormLogVperp2VparRedInterp(:, :, Qind, nn)= FinterpOut(:, :);
        NReNormERRORLogVperp2VparRedInterp(:, :, Qind, nn)= ...
            NReNormLogVperp2VparRedInterp(:, :, Qind, nn).^(-1e0/2e0);
        
        FinterpIn(:, :)= NReNormLogVperp1Vperp2Red(:, :, Qind, nn);
        FinterpOut(:, :)= interpn(Vperp1GCm3, Vperp2GCm3, FinterpIn, Vperp1GCq3, Vperp2GCq3);
        NReNormLogVperp1Vperp2RedInterp(:, :, Qind, nn)= FinterpOut(:, :);
        NReNormERRORLogVperp1Vperp2RedInterp(:, :, Qind, nn)= ...
            NReNormLogVperp1Vperp2RedInterp(:, :, Qind, nn).^(-1e0/2e0);
        
        % Interpolate Distribution Functions:
        FinterpIn(:, :)= FLogVperp1VparRed(:, :, Qind, nn);
        FinterpOut(:, :)= interpn(Vperp1GCm1, VparGCm1, FinterpIn, Vperp1GCq1, VparGCq1);
        FLogVperp1VparRedInterp(:, :, Qind, nn)= FinterpOut(:, :);
                
        FinterpIn(:, :)= FLogVperp2VparRed(:, :, Qind, nn);
        FinterpOut(:, :)= interpn(Vperp2GCm2, VparGCm2, FinterpIn, Vperp2GCq2, VparGCq2);
        FLogVperp2VparRedInterp(:, :, Qind, nn)= FinterpOut(:, :);
                
        FinterpIn(:, :)= FLogVperp1Vperp2Red(:, :, Qind, nn);
        FinterpOut(:, :)= interpn(Vperp1GCm3, Vperp2GCm3, FinterpIn, Vperp1GCq3, Vperp2GCq3);
        FLogVperp1Vperp2RedInterp(:, :, Qind, nn)= FinterpOut(:, :);
        
    end
end

toc

disp('LINEAR INTERPOLATIONS OF DISTRIBUTION FUNCTIONS COMPLETE')

%%
clc
close all

% PLOT DISTRIBUTION FUNCTIONS:

InterpFlag= 1;

Qind1= NqLB+ 3; 
Qind2= NqUB- 1;

pausevar= 1e-9;
FontSize= 25;

cbarpos= 'SouthOutside';
fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1200 800])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' DISTRIBUTION FUNCTIONS VID.avi'));
open(v);

if InterpFlag == 0
    Vperp12lim1= 4;
    Vperp12lim2= 25;
    Vparlim1= 4;
    Vparlim2= 25;
end
if InterpFlag == 1
    Vperp12lim1= 1; %20;
    Vperp12lim2= 100; %81;
    Vparlim1= 1; %20;
    Vparlim2= 100; %81;
end
LinHeatingflag= 0;
LogHeatingflag= 0;

% nn= Nend
% for Qind= NqLB:1:NqUB
for nn= Nend
    
    Qind= Qind1;
    subplot(2, 3, 1)
    if InterpFlag == 0
        pcolor((FLogVperp1VparRed(:, :, Qind, nn)/(max(max(FLogVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FLogVperp1VparRedInterp(:, :, Qind, nn)/(max(max(FLogVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_{\parallel})$ [s$^3$m$^{-6}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 2)
    if InterpFlag == 0
        pcolor((FLogVperp2VparRed(:, :, Qind, nn)/(max(max(FLogVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FLogVperp2VparRedInterp(:, :, Qind, nn)/(max(max(FLogVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 2}, \: v_{\parallel})$ [s$^3$m$^{-6}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 3)
    if InterpFlag == 0
        pcolor((FLogVperp1Vperp2Red(:, :, Qind, nn)/(max(max(FLogVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FLogVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(FLogVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_{\perp 2})$ [s$^3$m$^{-6}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    Qind= Qind2;
    subplot(2, 3, 4)
    if InterpFlag == 0
        pcolor((FLogVperp1VparRed(:, :, Qind, nn)/(max(max(FLogVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FLogVperp1VparRedInterp(:, :, Qind, nn)/(max(max(FLogVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_{\parallel})$ [s$^3$m$^{-6}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 5)
    if InterpFlag == 0
        pcolor((FLogVperp2VparRed(:, :, Qind, nn)/(max(max(FLogVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FLogVperp2VparRedInterp(:, :, Qind, nn)/(max(max(FLogVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 2}, \: v_{\parallel})$ [s$^3$m$^{-6}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 6)
    if InterpFlag == 0
        pcolor((FLogVperp1Vperp2Red(:, :, Qind, nn)/(max(max(FLogVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FLogVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(FLogVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_{\perp 2})$ [s$^3$m$^{-6}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    pause(pausevar)
    frame= getframe(gcf);
    writeVideo(v, frame);
end

close(v);

toc

disp('DISTRIBUTION FUNCTIONS MOVIE COMPLETE')

%%
clc
close all

% PLOT PITCH ANGLE DISTRIBUTIONS:

pausevar= 1e-9;
FontSize= 25;

cbarpos= 'SouthOutside';
fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1200 800])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' PITCH ANGLE DISTRIBS VID.avi'));
open(v);

LinHeatingflag= 0;
LogHeatingflag= 0;
                    
% nn= Nend
% for Qind= NqLB:1:NqUB
for nn= Nend
    
    Qind= Qind1;
    subplot(2, 3, 1)
    if InterpFlag == 0
        pcolor((FFVperp1VparRed(:, :, Qind, nn)/(max(max(FFVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FFVperp1VparRedInterp(:, :, Qind, nn)/(max(max(FFVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f_\alpha(v_{\perp 1}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 2)
    if InterpFlag == 0
        pcolor((FFVperp2VparRed(:, :, Qind, nn)/(max(max(FFVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FFVperp2VparRedInterp(:, :, Qind, nn)/(max(max(FFVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f_\alpha(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 3)
    if InterpFlag == 0
        pcolor((FFVperp1Vperp2Red(:, :, Qind, nn)/(max(max(FFVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FFVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(FFVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f_\alpha(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    Qind= Qind2;
    subplot(2, 3, 4)
    if InterpFlag == 0
        pcolor((FFVperp1VparRed(:, :, Qind, nn)/(max(max(FFVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FFVperp1VparRedInterp(:, :, Qind, nn)/(max(max(FFVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f_\alpha(v_{\perp 1}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 5)
    if InterpFlag == 0
        pcolor((FFVperp2VparRed(:, :, Qind, nn)/(max(max(FFVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FFVperp2VparRedInterp(:, :, Qind, nn)/(max(max(FFVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f_\alpha(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 6)
    if InterpFlag == 0
        pcolor((FFVperp1Vperp2Red(:, :, Qind, nn)/(max(max(FFVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((FFVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(FFVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f_\alpha(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
        
    pause(pausevar)
    frame= getframe(gcf);
    writeVideo(v, frame);
end

close(v);

toc

disp('PITCH ANGLE MOVIE COMPLETE')

%%
clc
close all

% PLOT DIFFERENTIAL NUMBER FLUXES:

pausevar= 1e-9;
FontSize= 25;

cbarpos= 'SouthOutside';
fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1200 800])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' DNF VID.avi'));
open(v);

LinHeatingflag= 0;
LogHeatingflag= 0;

% nn= Nend
% for Qind= NqLB:1:NqUB
for nn= Nend
    
    Qind= Qind1;
    subplot(2, 3, 1)
    if InterpFlag == 0
        pcolor((DNFVperp1VparRed(:, :, Qind, nn)/(max(max(DNFVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DNFVperp1VparRedInterp(:, :, Qind, nn)/(max(max(DNFVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_N(v_{\perp 1}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 2)
    if InterpFlag == 0
        pcolor((DNFVperp2VparRed(:, :, Qind, nn)/(max(max(DNFVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DNFVperp2VparRedInterp(:, :, Qind, nn)/(max(max(DNFVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_N(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 3)
    if InterpFlag == 0
        pcolor((DNFVperp1Vperp2Red(:, :, Qind, nn)/(max(max(DNFVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DNFVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(DNFVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$phi_N(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    Qind= Qind2;
    subplot(2, 3, 4)
    if InterpFlag == 0
        pcolor((DNFVperp1VparRed(:, :, Qind, nn)/(max(max(DNFVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DNFVperp1VparRedInterp(:, :, Qind, nn)/(max(max(DNFVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_N(v_{\perp 1}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 5)
    if InterpFlag == 0
        pcolor((DNFVperp2VparRed(:, :, Qind, nn)/(max(max(DNFVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DNFVperp2VparRedInterp(:, :, Qind, nn)/(max(max(DNFVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_N(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 6)
    if InterpFlag == 0
        pcolor((DNFVperp1Vperp2Red(:, :, Qind, nn)/(max(max(DNFVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DNFVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(DNFVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_N(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
        
    pause(pausevar)
    frame= getframe(gcf);
    writeVideo(v, frame);
end

close(v);

toc

disp('DNF MOVIE COMPLETE')

%%
clc
close all

% PLOT DIFFERENTIAL ENERGY FLUXES:

pausevar= 1e-9;
FontSize= 25;

cbarpos= 'SouthOutside';
fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1200 800])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' DEF VID.avi'));
open(v);

LinHeatingflag= 0;
LogHeatingflag= 0;

% nn= Nend
% for Qind= NqLB:1:NqUB
for nn= Nend
    
    Qind= Qind1;
    subplot(2, 3, 1)
    if InterpFlag == 0
        pcolor((DEFVperp1VparRed(:, :, Qind, nn)/(max(max(DEFVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DEFVperp1VparRedInterp(:, :, Qind, nn)/(max(max(DEFVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_E(v_{\perp 1}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 2)
    if InterpFlag == 0
        pcolor((DEFVperp2VparRed(:, :, Qind, nn)/(max(max(DEFVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DEFVperp2VparRedInterp(:, :, Qind, nn)/(max(max(DEFVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_E(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 3)
    if InterpFlag == 0
        pcolor((DEFVperp1Vperp2Red(:, :, Qind, nn)/(max(max(DEFVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DEFVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(DEFVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$phi_E(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    Qind= Qind2;
    subplot(2, 3, 4)
    if InterpFlag == 0
        pcolor((DEFVperp1VparRed(:, :, Qind, nn)/(max(max(DEFVperp1VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DEFVperp1VparRedInterp(:, :, Qind, nn)/(max(max(DEFVperp1VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_E(v_{\perp 1}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 5)
    if InterpFlag == 0
        pcolor((DEFVperp2VparRed(:, :, Qind, nn)/(max(max(DEFVperp2VparRed(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DEFVperp2VparRedInterp(:, :, Qind, nn)/(max(max(DEFVperp2VparRedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_E(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat(num2str(round((Rank(r).Specie(s).FluxTube(f).rGL(Qind)- RE)*1e-3)), '[km]', '$\le r \le$', ...
        num2str(round((Rank(r).Specie(s).FluxTube(f).rGH(Qind)- RE)*1e-3)), '[km], ', '$t= $', ...
        num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)/3600e0), '[hr]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
    
    subplot(2, 3, 6)
    if InterpFlag == 0
        pcolor((DEFVperp1Vperp2Red(:, :, Qind, nn)/(max(max(DEFVperp1Vperp2Red(:, :, Qind, nn)))))')
    end
    if InterpFlag == 1
        pcolor((DEFVperp1Vperp2RedInterp(:, :, Qind, nn)/(max(max(DEFVperp1Vperp2RedInterp(:, :, Qind, nn)))))')
    end
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$\phi_E(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    xlim([Vperp12lim1 Vperp12lim2]);
    ylim([Vparlim1 Vparlim2]);
    if InterpFlag == 0
        yticks([4 9 14.5 20 25]);
        yticklabels({num2str(round(-VparGCp(4)*1e-3, 1)), ...
            num2str(round(-VparGCp(9)*1e-3, 1)), ...
            num2str(round(-mean(VparGCp(14:15))*1e-3, 1)), ...
            num2str(round(-VparGCp(20)*1e-3, 1)), ...
            num2str(round(-VparGCp(25)*1e-3, 1))});
        xticks([4 9 14.5 20 25]);
        xticklabels({num2str(round(Vperp1GCp(4)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(9)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCp(14:15))*1e-3, 1)), ...
            num2str(round(Vperp1GCp(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCp(25)*1e-3, 1))});
    end
    if InterpFlag == 1
        yticks([20 35 50.5 66 81]);
        yticklabels({num2str(round(-VparGCqq(20)*1e-3, 1)), ...
            num2str(round(-VparGCqq(35)*1e-3, 1)), ...
            num2str(round(-mean(VparGCqq(50:51))*1e-3, 1)), ...
            num2str(round(-VparGCqq(66)*1e-3, 1)), ...
            num2str(round(-VparGCqq(81)*1e-3, 1))});
        xticks([20 35 50.5 66 81]);
        xticklabels({num2str(round(Vperp1GCqq(20)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(35)*1e-3, 1)), ...
            num2str(round(mean(Vperp1GCqq(50:51))*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(66)*1e-3, 1)), ...
            num2str(round(Vperp1GCqq(81)*1e-3, 1))});
    end
        
    pause(pausevar)
    frame= getframe(gcf);
    writeVideo(v, frame);
end

close(v);

toc

disp('DEF MOVIE COMPLETE')

%%
clc
close all

% PLOT ION DISTRIBUTION FUNCTIONS IN VELOCITY-SPACE:

pausevar= 1e-9;
FontSize= 25;
cbarpos= 'SouthOutside';
fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1200 800])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION DISTRIB FNC VID.avi'));
open(v);

LinHeatingflag= 1;

Qind= 10;
conlim1p1= min(min(min(NLogVperp1VparRed(:, :, Qind, :))));
conlim2p1= max(max(max(NLogVperp1VparRed(:, :, Qind, :))));
conlim1p2= min(min(min(NLogVperp2VparRed(:, :, Qind, :))));
conlim2p2= max(max(max(NLogVperp2VparRed(:, :, Qind, :))));
conlim1= max(conlim1p1, conlim1p2);
conlim2= max(conlim2p1, conlim2p2);
conlim112= min(min(min(NLogVperp1Vperp2Red(:, :, Qind, :))));
conlim212= max(max(max(NLogVperp1Vperp2Red(:, :, Qind, :))));

conlim1p1f= min(min(min(FLogVperp1VparRed(:, :, Qind, :))));
conlim2p1f= max(max(max(FLogVperp1VparRed(:, :, Qind, :))));
conlim1p2f= min(min(min(FLogVperp2VparRed(:, :, Qind, :))));
conlim2p2f= max(max(max(FLogVperp2VparRed(:, :, Qind, :))));
conlim1f= max(conlim1p1f, conlim1p2f);
conlim2f= max(conlim2p1f, conlim2p2f);
conlim112f= min(min(min(FLogVperp1Vperp2Red(:, :, Qind, :))));
conlim212f= max(max(max(FLogVperp1Vperp2Red(:, :, Qind, :))));
for nn= 1:1:Nend
    subplot(2, 3, 1)
    pcolor(NLogVperp1VparRed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    caxis([conlim1 conlim2])
    shading interp
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    if LinHeatingflag == 0
        xlim([Vperp12lim1 Vperp12lim2])
        ylim([Vparlim1 Vparlim2])
        yticks([6 11 16 20.5 25 30 35]);
        yticklabels({num2str(round(-VparGCp(6)*1e-3, 2)), ...
            num2str(round(-VparGCp(11)*1e-3, 2)), ...
            num2str(round(-VparGCp(16)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(25)*1e-3, 2)), ...
            num2str(round(-VparGCp(30)*1e-3, 2)), ...
            num2str(round(-VparGCp(35)*1e-3, 2))});
        xticks([2 4 6 10 15]);
        xticklabels({num2str(round(Vperp1GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    if LinHeatingflag == 1
        yticks([2 10 20.5 31 39]);
        yticklabels({num2str(round(-VparGCp(2)*1e-3, 2)), ...
            num2str(round(-VparGCp(10)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(31)*1e-3, 2)), ...
            num2str(round(-VparGCp(39)*1e-3, 2))});
        xticks([5 10 15]);
        xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    
    subplot(2, 3, 2)
    pcolor(NLogVperp2VparRed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    shading interp
    caxis([conlim1 conlim2])
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat('$r \sim$', num2str(round((Rank(r).Specie(s).FluxTube(f).QCell(Qind).rGC- RE)*1e-3)), ' [km]', ...
        ', $t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)), ' [s]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    if LinHeatingflag == 0
        xlim([Vperp12lim1 Vperp12lim2])
        ylim([Vparlim1 Vparlim2])
        yticks([6 11 16 20.5 25 30 35]);
        yticklabels({num2str(round(-VparGCp(6)*1e-3, 2)), ...
            num2str(round(-VparGCp(11)*1e-3, 2)), ...
            num2str(round(-VparGCp(16)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(25)*1e-3, 2)), ...
            num2str(round(-VparGCp(30)*1e-3, 2)), ...
            num2str(round(-VparGCp(35)*1e-3, 2))});
        xticks([2 4 6 10 15]);
        xticklabels({num2str(round(Vperp2GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
    end
    if LinHeatingflag == 1
        yticks([2 10 20.5 31 39]);
        yticklabels({num2str(round(-VparGCp(2)*1e-3, 2)), ...
            num2str(round(-VparGCp(10)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(31)*1e-3, 2)), ...
            num2str(round(-VparGCp(39)*1e-3, 2))});
        xticks([5 10 15]);
        xticklabels({num2str(round(Vperp2GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
    end
    
    subplot(2, 3, 3)
    pcolor(NLogVperp1Vperp2Red(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    shading interp
    caxis([conlim112 conlim212])
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    if LinHeatingflag == 0
        xlim([Vperp12lim1 Vperp12lim2])
        ylim([Vperp12lim1 Vperp12lim2])
        yticks([2 4 6 10 15]);
        yticklabels({num2str(round(Vperp2GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
        xticks([2 4 6 10 15]);
        xticklabels({num2str(round(Vperp1GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    if LinHeatingflag == 1
        yticks([5 10 15]);
        yticklabels({num2str(round(Vperp2GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
        xticks([5 10 15]);
        xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    
    subplot(2, 3, 4)
    pcolor(FLogVperp1VparRed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    shading interp
    caxis([conlim1f conlim2f])
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    if LinHeatingflag == 0
        xlim([Vperp12lim1 Vperp12lim2])
        ylim([Vparlim1 Vparlim2])
        yticks([6 11 16 20.5 25 30 35]);
        yticklabels({num2str(round(-VparGCp(6)*1e-3, 2)), ...
            num2str(round(-VparGCp(11)*1e-3, 2)), ...
            num2str(round(-VparGCp(16)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(25)*1e-3, 2)), ...
            num2str(round(-VparGCp(30)*1e-3, 2)), ...
            num2str(round(-VparGCp(35)*1e-3, 2))});
        xticks([2 4 6 10 15]);
        xticklabels({num2str(round(Vperp1GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    if LinHeatingflag == 1
        yticks([2 10 20.5 31 39]);
        yticklabels({num2str(round(-VparGCp(2)*1e-3, 2)), ...
            num2str(round(-VparGCp(10)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(31)*1e-3, 2)), ...
            num2str(round(-VparGCp(39)*1e-3, 2))});
        xticks([5 10 15]);
        xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    
    subplot(2, 3, 5)
    pcolor(FLogVperp2VparRed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    shading interp
    caxis([conlim1f conlim2f])
    xlabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\parallel}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 2}, \: v_{\parallel})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    if LinHeatingflag == 0
        xlim([Vperp12lim1 Vperp12lim2])
        ylim([Vparlim1 Vparlim2])
        yticks([6 11 16 20.5 25 30 35]);
        yticklabels({num2str(round(-VparGCp(6)*1e-3, 2)), ...
            num2str(round(-VparGCp(11)*1e-3, 2)), ...
            num2str(round(-VparGCp(16)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(25)*1e-3, 2)), ...
            num2str(round(-VparGCp(30)*1e-3, 2)), ...
            num2str(round(-VparGCp(35)*1e-3, 2))});
        xticks([2 4 6 10 15]);
        xticklabels({num2str(round(Vperp2GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
    end
    if LinHeatingflag == 1
        yticks([2 10 20.5 31 39]);
        yticklabels({num2str(round(-VparGCp(2)*1e-3, 2)), ...
            num2str(round(-VparGCp(10)*1e-3, 2)), ...
            num2str(round(-mean(VparGCp(20:21))*1e-3, 2)), ...
            num2str(round(-VparGCp(31)*1e-3, 2)), ...
            num2str(round(-VparGCp(39)*1e-3, 2))});
        xticks([5 10 15]);
        xticklabels({num2str(round(Vperp2GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
    end
    
    subplot(2, 3, 6)
    pcolor(FLogVperp1Vperp2Red(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    shading interp
    caxis([conlim112f conlim212f])
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_{\perp 2})$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    if LinHeatingflag == 0
        xlim([Vperp12lim1 Vperp12lim2])
        ylim([Vperp12lim1 Vperp12lim2])
        yticks([2 4 6 10 15]);
        yticklabels({num2str(round(Vperp2GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
        xticks([2 4 6 10 15]);
        xticklabels({num2str(round(Vperp1GCp(2)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(4)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(6)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    if LinHeatingflag == 1
        yticks([5 10 15]);
        yticklabels({num2str(round(Vperp2GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp2GCp(15)*1e-3, 2))});
        xticks([5 10 15]);
        xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
            num2str(round(Vperp1GCp(15)*1e-3, 2))});
    end
    
    pause(pausevar)
    frame= getframe(gcf);
end

writeVideo(v, frame);
close(v);

toc

disp('ION DISTRIBUTION FUNCTIONS MOVIE COMPLETE')

%%
clc
close all

% PLOT RAW ENA COUNTS IN VELOCITY-SPACE:

pausevar= 1e-9;
FontSize= 25;

cbarpos= 'SouthOutside';
fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1200 800])

set(gca, 'nextplot', 'replacechildren');
% v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' RAW ENA COUNTS VID.avi'));
% open(v);

Vpphilim1= 1;
Vpphilim2= 7;
Vqlim1= 8;
Vqlim2= 33;
Heatingflag= 0;

Qind= 2;
for nn= 1:1:Nend
    subplot(2, 3, 1)
    pcolor(NReNormLogVpVphiENARed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_p$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_\phi$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N_r(v_p, \: v_\phi)$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    %     if Heatingflag == 0
    %         xlim([Vpphilim1 Vpphilim2])
    %         ylim([Vqlim1 Vqlim2])
    %         yticks([6 11 16 25 30 35]);
    %         yticklabels({num2str(round(-VqGCp(6)*1e-3, 2)), ...
    %             num2str(round(-VqGCp(11)*1e-3, 2)), ...
    %             num2str(round(-VqGCp(16)*1e-3, 2)), ...
    %             num2str(round(-VqGCp(25)*1e-3, 2)), ...
    %             num2str(round(-VqGCp(30)*1e-3, 2)), ...
    %             num2str(round(-VqGCp(35)*1e-3, 2))});
    %         xticks([2 4 6 10 15]);
    %         xticklabels({num2str(round(VpGCp(2)*1e-3, 2)), ...
    %             num2str(round(VpGCp(4)*1e-3, 2)), ...
    %             num2str(round(VpGCp(6)*1e-3, 2)), ...
    %             num2str(round(VpGCp(10)*1e-3, 2)), ...
    %             num2str(round(VpGCp(15)*1e-3, 2))});
    %     end
    %     if Heatingflag == 1
    % %         xlim([Vperp12lim1 Vperp12lim2])
    % %         ylim([Vparlim1 Vparlim2])
    %         yticks([10 20 30 40]);
    %         yticklabels({num2str(round(-VparGCp(10)*1e-3, 2)), ...
    %             num2str(round(-VparGCp(20)*1e-3, 2)), ...
    %             num2str(round(-VparGCp(30)*1e-3, 2)), ...
    %             num2str(round(-VparGCp(40)*1e-3, 2))});
    %         xticks([5 10 15]);
    %         xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
    %             num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
    %             num2str(round(Vperp1GCp(15)*1e-3, 2))});
    %     end
    
    subplot(2, 3, 2)
    pcolor(NReNormLogVqVphiENARed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_q$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_\phi$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N_r(v_q, \: v_\phi)$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    title(horzcat('$r \sim$', num2str(round((Rank(r).Specie(s).FluxTube(f).QCell(Qind).rGC- RE)*1e-3)), ' [km]', ...
        ', $t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)), ' [s]'), ...
        'interpreter', 'latex', 'FontSize', 25);
    %     xlim([Vperp12lim1 Vperp12lim2])
    %     ylim([Vparlim1 Vparlim2])
    %     yticks([5 10 15]);
    %     yticklabels({num2str(round(VqGCp(5)*1e-3, 2)), ...
    %         num2str(round(VqGCp(10)*1e-3, 2)), ...
    %         num2str(round(VqGCp(15)*1e-3, 2))});
    %     xticks([5 10 15]);
    %     xticklabels({num2str(round(VphiGCp(5)*1e-3, 2)), ...
    %         num2str(round(VphiGCp(10)*1e-3, 2)), ...
    %         num2str(round(VphiGCp(15)*1e-3, 2))});
    
    subplot(2, 3, 3)
    pcolor(NReNormLogVpVqENARed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    shading interp
    colormap jet
    xlabel(horzcat('$v_p$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_q$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N_r(v_p, \: v_q)$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    %     xlim([Vperp12lim1 Vperp12lim2])
    %     ylim([Vperp12lim1 Vperp12lim2])
    %     yticks([2 4 6 10 15]);
    %     yticklabels({num2str(round(Vperp2GCp(2)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(4)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(6)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(15)*1e-3, 2))});
    %     xticks([2 4 6 10 15]);
    %     xticklabels({num2str(round(Vperp1GCp(2)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(4)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(6)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(15)*1e-3, 2))});
    
    subplot(2, 3, 4)
    pcolor(NReNormERRORLogVpVphiENARed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    xlabel(horzcat('$v_p$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_\phi$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N_r^{-1/2}(v_p, \: v_\phi)$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    %     yticks([6 11 16 25 30 35]);
    %     yticklabels({num2str(round(-VparGCp(6)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(11)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(16)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(25)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(30)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(35)*1e-3, 2))});
    %     xticks([5 10 15]);
    %     xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(15)*1e-3, 2))});
    
    subplot(2, 3, 5)
    pcolor(NReNormERRORLogVqVphiENARed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    xlabel(horzcat('$v_q$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_\phi$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N_r^{-1/2}(v_q, \: v_\phi)$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'YDir', 'reverse')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    %     yticks([6 11 16 25 30 35]);
    %     yticklabels({num2str(round(-VparGCp(6)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(11)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(16)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(25)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(30)*1e-3, 2)), ...
    %         num2str(round(-VparGCp(35)*1e-3, 2))});
    %     xticks([5 10 15]);
    %     xticklabels({num2str(round(Vperp2GCp(5)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(15)*1e-3, 2))});
    
    subplot(2, 3, 6)
    pcolor(NReNormERRORLogVpVqENARed(:, :, Qind, nn)')
    cbarqq0= colorbar(cbarpos);
    colormap jet
    xlabel(horzcat('$v_p$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_q$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$N_r^{-1/2}(v_p, \: v_q)$'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'xaxisLocation', 'top')
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    %     yticks([5 10 15]);
    %     yticklabels({num2str(round(Vperp2GCp(5)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(10)*1e-3, 2)), ...
    %         num2str(round(Vperp2GCp(15)*1e-3, 2))});
    %     xticks([5 10 15]);
    %     xticklabels({num2str(round(Vperp1GCp(5)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(10)*1e-3, 2)), ...
    %         num2str(round(Vperp1GCp(15)*1e-3, 2))});
    
    pause(pausevar)
    frame= getframe(gcf);
end

% writeVideo(v, frame);
% close(v);
%
toc

disp('RAW ENA COUNTS MOVIE COMPLETE')

%%
tic

clc
close all

% PLOT ION REDUCED VELOCITY DISTRIBUTIONS VIDEO:

qVind= 5; %VISIONS-1 altitude

colorp= {'k.' 'r.' 'b.' 'm.' 'g.' 'y.' 'c.'};
color= repmat(colorp, 1, ranksize);
pausevar= 1e-9;
linewidth= 2; fsize= 20;

fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1600 750])

FontSize= 25;
Vperp1lim1= 0e0;
Vperp1lim2= 50e0;
Vperp2lim1= Vperp1lim1;
Vperp2lim2= Vperp1lim2;
Vparlim1= -5e0;
Vparlim2= 5e0;
Flim1= 0e0;
Flim2= 5e0;
FElim1= 0e0;
FElim2= 5e0;

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION REDUCED PHASE SPACE DISTRIB FNC VID.avi'));
open(v);

for nn= Nend
    
    % ----------------------------------------------------
    %                 Qind= qVind;
    Qind= 2;
    
    % ----------------------------------------------------
    
    %     contourf(NReNormLinVperp1VparRed(:, :, 2, 1))
    %     xlim([1 10])
    %     ylim([70 75])
    %     colorbar
    %     caxis([0 1e3]);
    subplot(2, 3, 1)
    pcolor(squeeze(Vperp1lin3(:, 1, :))*1e-3, ...
        squeeze(Vparlin3(:, 1, :))*1e-3, NReNormLinVperp1VparRed(:, :, Qind, nn))
    %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
    cbarqq0= colorbar('EastOutside');
    %     caxis([0 1e3]);
    colormap jet
    %     xlim([Vperp1lim1 Vperp1lim2])
    %     ylim([Vparlim1 Vparlim2])
    %     zlim([0 1e4])
    xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    ylabel(horzcat('$v_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_\parallel)$ [s$^2 \cdot$ m$^{-5}$]'), ...
        'interpreter', 'latex', 'FontSize', FontSize)
    set(gca, 'yminortick', 'on','xminortick', 'on');
    set(gca, 'FontSize', FontSize);
    set(gcf, 'color','white');
    grid on
    hold off
    
    %     subplot(2, 3, 4)
    %     pcolor(squeeze(Vperp1log(:, 1, :))*1e-3, ...
    %         squeeze(Vparlog(:, 1, :))*1e-3, NReNormLogVperp1VparRed(:, :, Qind, nn))
    %     %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
    %     cbarqq0= colorbar('EastOutside');
    % %     caxis([Flim1 Flim2]);
    %     colormap jet
    %     xlim([Vperp1lim1 Vperp1lim2])
    %     ylim([Vparlim1 Vparlim2])
    % %     zlim([Flim1 Flim2])
    %     xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    %     ylabel(horzcat('$v_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    %     xlabel(cbarqq0, horzcat('$f(v_{\perp 1}, \: v_\parallel)$ [s$^2 \cdot$ m$^{-5}$]'), ...
    %         'interpreter', 'latex', 'FontSize', FontSize)
    %     set(gca, 'yminortick', 'on','xminortick', 'on');
    %     set(gca, 'FontSize', FontSize);
    %     set(gcf, 'color','white');
    %     grid on
    %     hold off
    
    %     subplot(2, 3, 2)
    %     pcolor(squeeze(Vperp1lin(:, 1, :))*1e-3, ...
    %         squeeze(Vparlin(:, 1, :))*1e-3, NReNormLinVperp1VparRed(:, :, Qind, nn))
    %     %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
    %     cbarqq0= colorbar('EastOutside');
    % %     caxis([Flim1 Flim2]);
    %     colormap jet
    % %     xlim([Vperp1lim1 Vperp1lim2])
    % %     ylim([Vparlim1 Vparlim2])
    % %     zlim([Flim1 Flim2])
    %     title(horzcat('$r \sim$', num2str(round((Rank(r).Specie(s).FluxTube(f).QCell(Qind).rGC- RE)*1e-3)), ' [km]', ...
    %         ', $t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)), ' [s]'), ...
    %         'interpreter', 'latex', 'FontSize', 25);
    %     xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    %     ylabel(horzcat('$v_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    %     xlabel(cbarqq0, horzcat('$N_m(v_{\perp 1}, \: v_\parallel)$ [unitless]'), ...
    %         'interpreter', 'latex', 'FontSize', FontSize)
    %     set(gca, 'yminortick', 'on','xminortick', 'on');
    %     set(gca, 'FontSize', FontSize);
    %     set(gcf, 'color','white');
    %     grid on
    %     hold off
    %
    %     subplot(2, 3, 3)
    %     pcolor(squeeze(Vperp1lin(:, 1, :))*1e-3, ...
    %         squeeze(Vparlin(:, 1, :))*1e-3, NReNormERRORLinVperp1VparRed(:, :, Qind, nn))
    %     %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
    %     cbarqq0= colorbar('EastOutside');
    % %     caxis([Flim1 Flim2]);
    %     colormap jet
    % %     xlim([Vperp1lim1 Vperp1lim2])
    % %     ylim([Vparlim1 Vparlim2])
    % %     zlim([Flim1 Flim2])
    %     xlabel(horzcat('$v_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    %     ylabel(horzcat('$v_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', FontSize)
    %     xlabel(cbarqq0, horzcat('$N_{error}(v_{\perp 1}, \: v_\parallel)$ [unitless]'), ...
    %         'interpreter', 'latex', 'FontSize', FontSize)
    %     set(gca, 'yminortick', 'on','xminortick', 'on');
    %     set(gca, 'FontSize', FontSize);
    %     set(gcf, 'color','white');
    %     grid on
    %     hold off
    
    % ----------------------------------------------------
    
    pause(pausevar)
    frame= getframe(gcf);
    
    writeVideo(v, frame);
    
    % ----------------------------------------------------
    
end

close(v);

toc

disp('ION REDUCED VELOCITY DISTRIBUTION MOVIE COMPLETE')

%%
for nn= 1:1:NNtT+ 1
    for Qind= NqLB:1:NqUB
        for Vparind= 1:1:NVparG
            for Vperpind= 1:1:NVperpG
                
                if ((Vparind == 1) | (Vparind == 2))
                    if ((Rank(r).Specie(s).FluxTube(f).QCell(1).VparGLp(1, 1) <= VparlinGC(Vparind, 1)) & ...
                            (VparlinGC(Vparind, 1) <= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, 2)))
                        
                        xLinInterp= VparlinGC(Vparind, 1);
                        xLinInterp1= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGLp(1, 1);
                        xLinInterp2= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, 2);
                        yLinInterp1(Qind, Vperpind, Vparind, nn)= Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphRTp(Vperpind, 2, nn);
                        yLinInterp2(Qind, Vperpind, Vparind, nn)= Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphRTp(Vperpind, 3, nn);
                        
                        yLinInterp(Qind, Vperpind, Vparind, nn)= yLinInterp1(Qind, Vperpind, Vparind, nn)+ ...
                            (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2(Qind, Vperpind, Vparind, nn)- ...
                            yLinInterp1(Qind, Vperpind, Vparind, nn)));
                        
                        NphVparLin(Qind, Vperpind, Vparind, nn)= yLinInterp(Qind, Vperpind, Vparind, nn);
                    end
                end
                
                if ((Vparind ~= 1) & (Vparind ~= 2) & ...
                        (Vparind ~= NVparG))
                    if ((Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, Vparind- 1) < ...
                            VparlinGC(Vparind, 1)) & ...
                            (VparlinGC(Vparind, 1) <= ...
                            Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, Vparind)))
                        
                        xLinInterp= VparlinGC(Vparind, 1);
                        xLinInterp1= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, Vparind- 1);
                        xLinInterp2= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, Vparind);
                        yLinInterp1(Qind, Vperpind, Vparind, nn)= Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphRTp(Vperpind, Vparind- 1, nn);
                        yLinInterp2(Qind, Vperpind, Vparind, nn)= Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphRTp(Vperpind, Vparind+ 1, nn);
                        
                        yLinInterp(Qind, Vperpind, Vparind, nn)= yLinInterp1(Qind, Vperpind, Vparind, nn)+ ...
                            (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2(Qind, Vperpind, Vparind, nn)- ...
                            yLinInterp1(Qind, Vperpind, Vparind, nn)));
                        
                        NphVparLin(Qind, Vperpind, Vparind, nn)= yLinInterp(Qind, Vperpind, Vparind, nn);
                        
                    end
                end
                
                if (Vparind == NVparG)
                    if ((Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, Vparind- 1) < ...
                            VparlinGC(Vparind, 1)) & ...
                            (VparlinGC(Vparind, 1) <= ...
                            Rank(r).Specie(s).FluxTube(f).QCell(1).VparGHp(1, Vparind)))
                        
                        xLinInterp= VparlinGC(Vparind, 1);
                        xLinInterp1= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGCp(1, Vparind- 1);
                        xLinInterp2= Rank(r).Specie(s).FluxTube(f).QCell(1).VparGHp(1, Vparind);
                        yLinInterp1(Qind, Vperpind, Vparind, nn)= Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphRTp ...
                            (Vperpind, NVparG- 1, nn);
                        yLinInterp2(Qind, Vperpind, Vparind, nn)= Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphRTp ...
                            (Vperpind, NVparG, nn);
                        
                        yLinInterp(Qind, Vperpind, Vparind, nn)= yLinInterp1(Qind, Vperpind, Vparind, nn)+ ...
                            (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2(Qind, Vperpind, Vparind, nn)- ...
                            yLinInterp1(Qind, Vperpind, Vparind, nn)));
                        
                        NphVparLin(Qind, Vperpind, Vparind, nn)= yLinInterp(Qind, Vperpind, Vparind, nn);
                        
                    end
                end
            end
        end
    end
end

for nn= 1:1:NNtT+ 1
    for Qind= NqLB:1:NqUB
        for Vparind= 1:1:NVparG
            for Vperpind= 1:1:NVperpG
                
                if ((Vperpind == 1) | (Vperpind == 2))
                    if ((Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGLp(1, 1) <= VperplinGC(1, Vperpind)) & ...
                            (VperplinGC(1, Vperpind) <= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(2, 1)))
                        
                        xLinInterp= VperplinGC(1, Vperpind);
                        xLinInterp1= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGLp(1, 1);
                        xLinInterp2= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(2, 1);
                        yLinInterp1(Qind, Vperpind, Vparind, nn)= NphVparLin(Qind, 2, Vparind, nn);
                        yLinInterp2(Qind, Vperpind, Vparind, nn)= NphVparLin(Qind, 3, Vparind, nn);
                        
                        yLinInterp(Qind, Vperpind, Vparind, nn)= yLinInterp1(Qind, Vperpind, Vparind, nn)+ ...
                            (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2(Qind, Vperpind, Vparind, nn)- ...
                            yLinInterp1(Qind, Vperpind, Vparind, nn)));
                        
                        NphLin(Qind, Vperpind, Vparind, nn)= yLinInterp(Qind, Vperpind, Vparind, nn);
                    end
                end
                
                if ((Vperpind ~= 1) & (Vperpind ~= 2) & (Vperpind ~= NVperpG))
                    if ((Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(Vperpind- 1, 1) < ...
                            VperplinGC(1, Vperpind)) & ...
                            (VperplinGC(1, Vperpind) <= ...
                            Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(Vperpind, 1)))
                        
                        xLinInterp= VperplinGC(1, Vperpind);
                        xLinInterp1= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(Vperpind- 1, 1);
                        xLinInterp2= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(Vperpind, 1);
                        yLinInterp1(Qind, Vperpind, Vparind, nn)= NphVparLin(Qind, Vperpind- 1, Vparind, nn);
                        yLinInterp2(Qind, Vperpind, Vparind, nn)= NphVparLin(Qind, Vperpind+ 1, Vparind, nn);
                        
                        yLinInterp(Qind, Vperpind, Vparind, nn)= yLinInterp1(Qind, Vperpind, Vparind, nn)+ ...
                            (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2(Qind, Vperpind, Vparind, nn)- ...
                            yLinInterp1(Qind, Vperpind, Vparind, nn)));
                        
                        NphLin(Qind, Vperpind, Vparind, nn)= yLinInterp(Qind, Vperpind, Vparind, nn);
                        
                    end
                end
                
                if (Vperpind == NVperpG)
                    if ((Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(Vperpind- 1, 1) < ...
                            VperplinGC(1, Vperpind)) & ...
                            (VperplinGC(1, Vperpind) <= ...
                            Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGHp(Vperpind, 1)))
                        
                        xLinInterp= VperplinGC(1, Vperpind);
                        xLinInterp1= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGCp(Vperpind- 1, 1);
                        xLinInterp2= Rank(r).Specie(s).FluxTube(f).QCell(1).VperpGHp(Vperpind, 1);
                        yLinInterp1(Qind, Vperpind, Vparind, nn)= NphVparLin ...
                            (Qind, NVperpG- 1, Vparind, nn);
                        yLinInterp2(Qind, Vperpind, Vparind, nn)= NphVparLin ...
                            (Qind, NVperpG, Vparind, nn);
                        
                        yLinInterp(Qind, Vperpind, Vparind, nn)= yLinInterp1(Qind, Vperpind, Vparind, nn)+ ...
                            (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2(Qind, Vperpind, Vparind, nn)- ...
                            yLinInterp1(Qind, Vperpind, Vparind, nn)));
                        
                        NphLin(Qind, Vperpind, Vparind, nn)= yLinInterp(Qind, Vperpind, Vparind, nn);
                        
                    end
                end
            end
        end
    end
end

for nn= 1:1:NNtT+ 1
    for Qind= NqLB:1:NqUB
        for Vparind= 1:1:NVparG
            for Vperpind= 1:1:NVperpG
                Rank(r).Specie(s).FluxTube(f).QCell(Qind).NphLinRT(Vperpind, Vparind, nn)= ...
                    NphLin(Qind, Vperpind, Vparind, nn);
            end
        end
    end
end

%% -------------------------------------------------------
clc
close all

figure

subplot(1, 3, 1)
contourf(Rank(root).Specie(s).FluxTube(f).QCell(5).FNorm2log10(:, :, 2)'); colorbar
subplot(1, 3, 2)
contourf(Rank(root).Specie(s).FluxTube(f).QCell(5).NphReNormRT2log10(:, :, 2)'); colorbar
subplot(1, 3, 3)
contourf(Rank(root).Specie(s).FluxTube(f).QCell(5).NphReNormERRORRTp2log10(:, :, 2)'); colorbar

%%
tic

clc
close all

% PLOT ION VELOCITY DISTRIBUTIONS VIDEO:

qVind= 5; %VISIONS-1 altitude

colorp= {'k.' 'r.' 'b.' 'm.' 'g.' 'y.' 'c.'};
color= repmat(colorp, 1, ranksize);
pausevar= 1e-9;
Vgridlogflag= 1;
linewidth= 2; fsize= 20;

Qindlen= ((NqUB- NqLB)+ 1);
for Qind= 1:1:Qindlen;
    rr(Qind)= Rank(r).Specie(s).FluxTube(f).rGCp(Qind);
end
nd(:)= Rank(r).Specie(s).FluxTube(f).TimeT(:);
[rrAB ndAB]= meshgrid(rr, nd);

ylim1= 0; %200;
ylim2= 0.8d3;

Nlim1= 0; Nlim2= 12.5; Nlimflag= 0;
UPerplim1= 0; UPerplim2= 1; UPerplimflag= 1;
UParlim1= -2; UParlim2= 2; UParlimflag= 1;

fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1600 750])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION PHASE SPACE DISTRIB FNC VID.avi'));
open(v);

if IONVPERPVECflag == 1
    
    for nn= NNtT+ 1
        
        % ----------------------------------------------------
        %                 Qind= qVind;
        Qind= 1;
        
        % ----------------------------------------------------
        
        subplot(2, 3, 6)
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            xlim([UParlim1 UParlim2]);
        end
        hold off
        
        subplot(2, 3, 5)
        plot(real(Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_\perp$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerplimflag == 1)
            xlim([UPerplim1 UPerplim2]);
        end
        hold off
        
        subplot(2, 3, 4)
        plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, :)), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        plot(log10(NqRTp(nn, :)), rrAB, 'ro', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('log$_{10}(n)$ [m$^{-3}$]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Nlimflag == 1)
            xlim([Nlim1 Nlim2]);
        end
        hold off
        
        subplot(2, 3, 1)
        pp= isosurface(FPerpNorm2log10p(:, :, :));
        xlim1= 6;
        xlim2= NVperpG+ 30;
        xlim([xlim1 xlim2]);
        xticks([15.5 25.5 35.5 45.5 55.5]);
        if Vgridlogflag == 1
            xticklabels({num2str(round(mean(VperpGCplog(15:16)), 1)), ...
                num2str(round(mean(VperpGCplog(25:26)), 1)), ...
                num2str(round(mean(VperpGCplog(35:36)), 1)), ...
                num2str(round(mean(VperpGCplog(45:46)), 1)), ...
                num2str(round(mean(VperpGCplog(55:56)), 1))});
        else
            xticklabels({num2str(round(mean(VperpGCp(15:16)), 1)), ...
                num2str(round(mean(VperpGCp(25:26)), 1)), ...
                num2str(round(mean(VperpGCp(35:36)), 1)), ...
                num2str(round(mean(VperpGCp(45:46)), 1)), ...
                num2str(round(mean(VperpGCp(55:56)), 1))});
        end
        ylim1= 4;
        ylim2= 29;
        ylim([ylim1 ylim2]);
        yticks([4.5 10.5 18 25.5 32.5]);
        if Vgridlogflag == 1
            yticklabels({num2str(-round(mean(VparGCplog(4:5)), 1)), ...
                num2str(-round(mean(VparGCplog(10:11)), 1)), ...
                num2str(-round(mean(VparGCplog(18)), 1)), ...
                num2str(-round(mean(VparGCplog(25:26)), 1)), ...
                num2str(-round(mean(VparGCplog(32:33)), 1))});
        else
            yticklabels({num2str(-round(mean(VparGCp(4:5)), 1)), ...
                num2str(-round(mean(VparGCp(10:11)), 1)), ...
                num2str(-round(mean(VparGCp(18)), 1)), ...
                num2str(-round(mean(VparGCp(25:26)), 1)), ...
                num2str(-round(mean(VparGCp(32:33)), 1))});
        end
        set(gca, 'Ydir', 'reverse')
        text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        cmin= min(min(min(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FNorm2log10(:, :, :))));
        cmax= max(max(max(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FNorm2log10(:, :, :))));
        caxis([-8 -2]);
        xlabel('log$_{10}(\mathbf{v_\perp})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}(f)$', ' [s$^3 \cdot$ m$^{-6}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        %                 set(gca, 'YScale', 'log')
        %                 set(gca, 'XScale', 'log')
        grid on
        hold off
        
        subplot(2, 3, 2)
        pp= contourf(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormRT2log10(:, :, nn)', 5);
        xlim1= 6;
        xlim2= NVperpG+ 30;
        ylim1= 4;
        ylim2= 29;
        xlim([xlim1 xlim2]);
        ylim([ylim1 ylim2]);
        xticks([15.5 25.5 35.5 45.5 55.5]);
        xticklabels({num2str(round(mean(VperpGCplog(15:16)), 1)), ...
            num2str(round(mean(VperpGCplog(25:26)), 1)), ...
            num2str(round(mean(VperpGCplog(35:36)), 1)), ...
            num2str(round(mean(VperpGCplog(45:46)), 1)), ...
            num2str(round(mean(VperpGCplog(55:56)), 1))});
        yticks([4.5 10.5 18 25.5 32.5]);
        yticklabels({num2str(-round(mean(VparGCplog(4:5)), 1)), ...
            num2str(-round(mean(VparGCplog(10:11)), 1)), ...
            num2str(-round(mean(VparGCplog(18)), 1)), ...
            num2str(-round(mean(VparGCplog(25:26)), 1)), ...
            num2str(-round(mean(VparGCplog(32:33)), 1))});
        set(gca, 'Ydir', 'reverse')
        text(xlim1- 8, ylim1- 5, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        title(horzcat('$r \sim$', num2str(round((Rank(r).Specie(s).FluxTube(f).QCell(Qind).rGC- RE)*1e-3)), ' [km]', ...
            ', $t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)), ' [s]'), ...
            'interpreter', 'latex', 'FontSize', 25);
        cmin= min(min(min(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormRT2log10(:, :, :))));
        cmax= max(max(max(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormRT2log10(:, :, :))));
        caxis([cmin cmax]);
        xlabel('log$_{10}(\mathbf{v_\perp})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}(\mathcal{N^\prime})$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        %                 set(gca, 'YScale', 'log')
        %                 set(gca, 'XScale', 'log')
        grid on
        hold off
        
        subplot(2, 3, 3)
        pp= contourf(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORRTp2log10(:, :, nn)', 5);
        xlim1= 6;
        xlim2= NVperpG+ 30;
        ylim1= 4;
        ylim2= 29;
        xlim([xlim1 xlim2]);
        ylim([ylim1 ylim2]);
        xticks([15.5 25.5 35.5 45.5 55.5]);
        xticklabels({num2str(round(mean(VperpGCplog(15:16)), 1)), ...
            num2str(round(mean(VperpGCplog(25:26)), 1)), ...
            num2str(round(mean(VperpGCplog(35:36)), 1)), ...
            num2str(round(mean(VperpGCplog(45:46)), 1)), ...
            num2str(round(mean(VperpGCplog(55:56)), 1))});
        yticks([4.5 10.5 18 25.5 32.5]);
        yticklabels({num2str(-round(mean(VparGCplog(4:5)), 1)), ...
            num2str(-round(mean(VparGCplog(10:11)), 1)), ...
            num2str(-round(mean(VparGCplog(18)), 1)), ...
            num2str(-round(mean(VparGCplog(25:26)), 1)), ...
            num2str(-round(mean(VparGCplog(32:33)), 1))});
        set(gca, 'Ydir', 'reverse')
        text(xlim1- 8, ylim1- 5, '(c)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        cmin= min(min(min(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORRTp2log10(:, :, :))));
        cmax= max(max(max(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORRTp2log10(:, :, :))));
        caxis([cmin cmax]);
        xlabel('log$_{10}(\mathbf{v_\perp})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}(\mathcal{N^\prime}^{-1/2})$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        %                 set(gca, 'YScale', 'log')
        %                 set(gca, 'XScale', 'log')
        grid on
        hold off
        
        % ----------------------------------------------------
        
        pause(pausevar)
        frame= getframe(gcf);
        
        writeVideo(v, frame);
        
        % ----------------------------------------------------
        
    end
    
else
    
    for nn= NNtT+ 1
        
        % ----------------------------------------------------
        %                 Qind= qVind;
        Qind= 1;
        
        % ----------------------------------------------------
        
        subplot(2, 3, 6)
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            xlim([UParlim1 UParlim2]);
        end
        hold off
        
        subplot(2, 3, 5)
        plot(real(Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(Rank(r).Specie(s).FluxTube(f).M1PerpphRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_\perp$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerplimflag == 1)
            xlim([UPerplim1 UPerplim2]);
        end
        hold off
        
        subplot(2, 3, 4)
        plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, :)), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(log10(Rank(r).Specie(s).FluxTube(f).M0phRTr(nn, Qind)), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        plot(log10(NqRTp(nn, :)), rrAB, 'ro', 'LineWidth', linewidth)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('log$_{10}(n)$ [m$^{-3}$]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Nlimflag == 1)
            xlim([Nlim1 Nlim2]);
        end
        hold off
        
        subplot(2, 3, 1)
        pp= contourf(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FNorm2log10(:, :, nn)', 5);
        xlim1= 6;
        xlim2= NVperpG+ 30;
        xlim([xlim1 xlim2]);
        xticks([15.5 25.5 35.5 45.5 55.5]);
        if Vgridlogflag == 1
            xticklabels({num2str(round(mean(VperpGCplog(15:16)), 1)), ...
                num2str(round(mean(VperpGCplog(25:26)), 1)), ...
                num2str(round(mean(VperpGCplog(35:36)), 1)), ...
                num2str(round(mean(VperpGCplog(45:46)), 1)), ...
                num2str(round(mean(VperpGCplog(55:56)), 1))});
        else
            xticklabels({num2str(round(mean(VperpGCp(15:16)), 1)), ...
                num2str(round(mean(VperpGCp(25:26)), 1)), ...
                num2str(round(mean(VperpGCp(35:36)), 1)), ...
                num2str(round(mean(VperpGCp(45:46)), 1)), ...
                num2str(round(mean(VperpGCp(55:56)), 1))});
        end
        ylim1= 4;
        ylim2= 29;
        ylim([ylim1 ylim2]);
        yticks([4.5 10.5 18 25.5 32.5]);
        if Vgridlogflag == 1
            yticklabels({num2str(-round(mean(VparGCplog(4:5)), 1)), ...
                num2str(-round(mean(VparGCplog(10:11)), 1)), ...
                num2str(-round(mean(VparGCplog(18)), 1)), ...
                num2str(-round(mean(VparGCplog(25:26)), 1)), ...
                num2str(-round(mean(VparGCplog(32:33)), 1))});
        else
            yticklabels({num2str(-round(mean(VparGCp(4:5)), 1)), ...
                num2str(-round(mean(VparGCp(10:11)), 1)), ...
                num2str(-round(mean(VparGCp(18)), 1)), ...
                num2str(-round(mean(VparGCp(25:26)), 1)), ...
                num2str(-round(mean(VparGCp(32:33)), 1))});
        end
        set(gca, 'Ydir', 'reverse')
        text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        cmin= min(min(min(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FNorm2log10(:, :, :))));
        cmax= max(max(max(Rank(root).Specie(s).FluxTube(f).QCell(Qind).FNorm2log10(:, :, :))));
        caxis([-8 -2]);
        xlabel('log$_{10}(\mathbf{v_\perp})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}(f)$', ' [s$^3 \cdot$ m$^{-6}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        %                 set(gca, 'YScale', 'log')
        %                 set(gca, 'XScale', 'log')
        grid on
        hold off
        
        subplot(2, 3, 2)
        pp= contourf(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormRT2log10(:, :, nn)', 5);
        xlim1= 6;
        xlim2= NVperpG+ 30;
        ylim1= 4;
        ylim2= 29;
        xlim([xlim1 xlim2]);
        ylim([ylim1 ylim2]);
        xticks([15.5 25.5 35.5 45.5 55.5]);
        xticklabels({num2str(round(mean(VperpGCplog(15:16)), 1)), ...
            num2str(round(mean(VperpGCplog(25:26)), 1)), ...
            num2str(round(mean(VperpGCplog(35:36)), 1)), ...
            num2str(round(mean(VperpGCplog(45:46)), 1)), ...
            num2str(round(mean(VperpGCplog(55:56)), 1))});
        yticks([4.5 10.5 18 25.5 32.5]);
        yticklabels({num2str(-round(mean(VparGCplog(4:5)), 1)), ...
            num2str(-round(mean(VparGCplog(10:11)), 1)), ...
            num2str(-round(mean(VparGCplog(18)), 1)), ...
            num2str(-round(mean(VparGCplog(25:26)), 1)), ...
            num2str(-round(mean(VparGCplog(32:33)), 1))});
        set(gca, 'Ydir', 'reverse')
        text(xlim1- 8, ylim1- 5, '(b)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        title(horzcat('$r \sim$', num2str(round((Rank(r).Specie(s).FluxTube(f).QCell(Qind).rGC- RE)*1e-3)), ' [km]', ...
            ', $t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn)), ' [s]'), ...
            'interpreter', 'latex', 'FontSize', 25);
        cmin= min(min(min(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormRT2log10(:, :, :))));
        cmax= max(max(max(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormRT2log10(:, :, :))));
        caxis([cmin cmax]);
        xlabel('log$_{10}(\mathbf{v_\perp})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}(\mathcal{N^\prime})$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        %                 set(gca, 'YScale', 'log')
        %                 set(gca, 'XScale', 'log')
        grid on
        hold off
        
        subplot(2, 3, 3)
        pp= contourf(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORRTp2log10(:, :, nn)', 5);
        xlim1= 6;
        xlim2= NVperpG+ 30;
        ylim1= 4;
        ylim2= 29;
        xlim([xlim1 xlim2]);
        ylim([ylim1 ylim2]);
        xticks([15.5 25.5 35.5 45.5 55.5]);
        xticklabels({num2str(round(mean(VperpGCplog(15:16)), 1)), ...
            num2str(round(mean(VperpGCplog(25:26)), 1)), ...
            num2str(round(mean(VperpGCplog(35:36)), 1)), ...
            num2str(round(mean(VperpGCplog(45:46)), 1)), ...
            num2str(round(mean(VperpGCplog(55:56)), 1))});
        yticks([4.5 10.5 18 25.5 32.5]);
        yticklabels({num2str(-round(mean(VparGCplog(4:5)), 1)), ...
            num2str(-round(mean(VparGCplog(10:11)), 1)), ...
            num2str(-round(mean(VparGCplog(18)), 1)), ...
            num2str(-round(mean(VparGCplog(25:26)), 1)), ...
            num2str(-round(mean(VparGCplog(32:33)), 1))});
        set(gca, 'Ydir', 'reverse')
        text(xlim1- 8, ylim1- 5, '(c)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        cmin= min(min(min(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORRTp2log10(:, :, :))));
        cmax= max(max(max(Rank(root).Specie(s).FluxTube(f).QCell(Qind).NphReNormERRORRTp2log10(:, :, :))));
        caxis([cmin cmax]);
        xlabel('log$_{10}(\mathbf{v_\perp})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}(\mathcal{N^\prime}^{-1/2})$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        %                 set(gca, 'YScale', 'log')
        %                 set(gca, 'XScale', 'log')
        grid on
        hold off
        
        % ----------------------------------------------------
        
        pause(pausevar)
        frame= getframe(gcf);
        
        writeVideo(v, frame);
        
        % ----------------------------------------------------
        
    end
    
end

close(v);

toc

disp('ION VELOCITY DISTRIBUTION MOVIE COMPLETE')

%% -------------------------------------------------------
tic
clc
close all

Qindlen= ((NqUB- NqLB)+ 1);
for Qind= 1:1:Qindlen;
    rr(Qind)= Rank(r).Specie(s).FluxTube(f).rGCp(Qind);
end
nd(:)= Rank(r).Specie(s).FluxTube(f).TimeT(:);
[rrAB ndAB]= meshgrid(rr, nd);

cbarloc= 'EastOutside'; logAEflag= 0; logTflag= 1; logWflag= 0; fignum= 1; figlabelflag= 0; FIGSAVEflag= 1;
pcolorflag= 1; cmapset= 'jet'; xfigsize= 750; yfigsize= 350; fsize= 20; linewidth= 2;
dxfigsize= 250; dyfigsize= 250; legloc= 'NorthEast';
ylim1= 0.5d3;
ylim2= 0.7d3;

fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);

set(fig(fignum), 'Position', [10 10 xfigsize yfigsize])

subplot(1, 2, 1)
pp= pcolor(ndAB/60, rrAB, ...
    Rank(r).Specie(s).FluxTube(f).nuIonNeutRT);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(fignum, cmapset);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\nu_{in}h$'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);

subplot(1, 2, 2)
pp= pcolor(ndAB/60, rrAB, ...
    Rank(r).Specie(s).FluxTube(f).sigmaIonNeutRT);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(fignum, cmapset);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\sigma_{in}$'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);

%%
clc
close all

for Qind= NqLB:1:NqUB
    if (Qind == 1)
        %         EAPressureR1(Qind)= (2e0/abs(Rank(r).Specie(s).FluxTube(f).QCell(Qind+ 1).hqC(1)+ ...
        %             Rank(r).Specie(s).FluxTube(f).QCell(Qind).hqC(1)))* ...
        %             (kB*Rank(r).Specie(s).FluxTube(f).TeNT(1)/ ...
        %             (qion*Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind)));
        %         EAPressureR2(Qind)= (abs(Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind+ 1)- ...
        %             Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind)))/ ...
        %             (abs(Rank(r).Specie(s).FluxTube(f).QCell(Qind+ 1).qGC(1)- ...
        %             Rank(r).Specie(s).FluxTube(f).QCell(Qind).qGC(1)));
        
        EAPressureR0(Qind)= (kB*Rank(r).Specie(s).FluxTube(f).TeNT(1)/ ...
            (qion*Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind)));
        EAPressureR1(Qind)= (2e0/abs(Rank(r).Specie(s).FluxTube(f).QCell(Qind+ 1).hqC(1)+ ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).hqC(1)));
        EAPressureR2(Qind)= ((abs(Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind+ 1)- ...
            Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind)))/ ...
            (abs(Rank(r).Specie(s).FluxTube(f).QCell(Qind+ 1).qGC(1)- ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).qGC(1))));
    end
    if (Qind ~= 1)
        %         EAPressureR1(Qind)= (2e0/(Rank(r).Specie(s).FluxTube(f).QCell(Qind- 1).hqC(1)+ ...
        %             Rank(r).Specie(s).FluxTube(f).QCell(Qind).hqC(1)))* ...
        %             (kB*Rank(r).Specie(s).FluxTube(f).TeNT(1)/ ...
        %             (qion*Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind)));
        %         EAPressureR2(Qind)= (Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind- 1)- ...
        %             Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind))/ ...
        %             (abs(Rank(r).Specie(s).FluxTube(f).QCell(Qind- 1).qGC(1)- ...
        %             Rank(r).Specie(s).FluxTube(f).QCell(Qind).qGC(1)));
        
        EAPressureR0(Qind)= (kB*Rank(r).Specie(s).FluxTube(f).TeNT(1)/ ...
            (qion*Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind)));
        EAPressureR1(Qind)= (2e0/(Rank(r).Specie(s).FluxTube(f).QCell(Qind- 1).hqC(1)+ ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).hqC(1)));
        EAPressureR2(Qind)= ((Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind- 1)- ...
            Rank(r).Specie(s).FluxTube(f).M0phRTr(1, Qind))/ ...
            (abs(Rank(r).Specie(s).FluxTube(f).QCell(Qind- 1).qGC(1)- ...
            Rank(r).Specie(s).FluxTube(f).QCell(Qind).qGC(1))));
    end
    EAPressure(Qind)= EAPressureR0(Qind)*EAPressureR1(Qind)*EAPressureR2(Qind);
end

subplot(4, 1, 1)
plot(EAPressureR0(:), rr)
ylim([500 750])
subplot(4, 1, 2)
plot(EAPressureR1(:).*EAPressureR2(:), rr)
ylim([500 750])
subplot(4, 1, 3)
plot(EAPressureR2(:), rr)
ylim([500 750])
subplot(4, 1, 4)
plot(Rank(r).Specie(s).FluxTube(f).EAPressureRTr(end, :), rr, 'kx-')
hold on
plot(EAPressure(:), rr, 'b')
ylim([500 750])

%%
tic

clc
close all

% PLOT ION REDUCED VELOCITY DISTRIBUTIONS VIDEO:

qVind= 5; %VISIONS-1 altitude

colorp= {'k.' 'r.' 'b.' 'm.' 'g.' 'y.' 'c.'};
color= repmat(colorp, 1, ranksize);
pausevar= 1e-9;
Vgridlogflag= 1;
linewidth= 2; fsize= 20;
LOGREDflag= 1;

fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1600 750])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION REDUCED PHASE SPACE DISTRIB FNC VID.avi'));
open(v);

if IONVPERPVECflag == 1
    
    for nn= NNtT+ 1
        
        % ----------------------------------------------------
        %                 Qind= qVind;
        Qind= 2;
        
        % ----------------------------------------------------
        
        subplot(2, 3, 1)
        if LOGREDflag == 1
            pp= contourf(F2PerpVperp1VparRedlog10(:, :, Qind, nn)');
        else
            pp= contourf(F2PerpVperp1VparRed(:, :, Qind, nn)');
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGREDflag == 1
            cmin= min(min(min(min(F2PerpVperp1VparRedlog10(:, :, :, :)))));
            cmax= max(max(max(max(F2PerpVperp1VparRedlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(F2PerpVperp1VparRed(:, :, :, :)))));
            cmax= max(max(max(max(F2PerpVperp1VparRed(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        xlim([1 NVperp1G]);
        xticks([2 4 6]);
        %         xticklabels({round(log10(Vperp1GCp(2))), round(log10(Vperp1GCp(4))), ...
        %             round(log10(Vperp1GCp(6)))});
        xticklabels({round(Vperp1GCp(2)), round(Vperp1GCp(4)), ...
            round(Vperp1GCp(6))});
        %         yticks([5 10 15]);
        %         yticklabels({round(-log10(VparGCp(5))), 0e0, ...
        %             round(log10(VparGCp(15)))});
        xlabel('log$_{10}(v_{\perp 1})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[f(v_{\perp 1}, \: v_\parallel)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        %         set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid off
        
        %         subplot(2, 3, 4)
        %         if LOGREDflag == 1
        %             pp= contourf(F2PerpEVperp1VparRedlog10(:, :, Qind, nn)');
        %         else
        %             pp= contourf(F2PerpEVperp1VparRed(:, :, Qind, nn)');
        %         end
        %         %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        %         cbarqq0= colorbar('SouthOutside');
        %         if LOGREDflag == 1
        %             cmin= min(min(min(min(F2PerpEVperp1VparRedlog10(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpEVperp1VparRedlog10(:, :, :, :)))));
        %         else
        %             cmin= min(min(min(min(F2PerpEVperp1VparRed(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpEVperp1VparRed(:, :, :, :)))));
        %         end
        %         caxis([cmin cmax]);
        %         xlim([1 NVperp1G]);
        %         xticks([2 4 6]);
        %         xticklabels({round(PitchAngleGCp(2)), round(PitchAngleGCp(4)), ...
        %             round(PitchAngleGCp(6))});
        %         yticks([5 10 15]);
        %         yticklabels({round(GyroAngleGCp(5)), round(GyroAngleGCp(10)), ...
        %             round(GyroAngleGCp(15))});
        %         xlabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        %         ylabel('$\theta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        %         xlabel(cbarqq0, horzcat('log$_{10}[f(\alpha, \: \theta)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
        %             'interpreter', 'latex', 'FontSize', 25)
        %         set(gca, 'xaxisLocation', 'top');
        %         set(gca, 'yminortick', 'on','xminortick', 'on');
        %         set(gca, 'FontSize', 25);
        %         set(gcf, 'color','white');
        %         grid on
        %
        %         subplot(2, 3, 2)
        %         if LOGREDflag == 1
        %             pp= contourf(F2PerpVperp2VparRedlog10(:, :, Qind, nn)');
        %         else
        %             pp= contourf(F2PerpVperp2VparRed(:, :, Qind, nn)');
        %         end
        %         %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        %         cbarqq0= colorbar('SouthOutside');
        %         if LOGREDflag == 1
        %             cmin= min(min(min(min(F2PerpVperp2VparRedlog10(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpVperp2VparRedlog10(:, :, :, :)))));
        %         else
        %             cmin= min(min(min(min(F2PerpVperp2VparRed(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpVperp2VparRed(:, :, :, :)))));
        %         end
        %         caxis([cmin cmax]);
        %         xlim([1 NVperp2G]);
        %         xticks([2 4 6]);
        %         xticklabels({round(log10(Vperp2GCp(2))), round(log10(Vperp2GCp(4))), ...
        %             round(log10(Vperp2GCp(6)))});
        %         yticks([5 10 15]);
        %         yticklabels({round(-log10(VparGCp(5))), 0e0, ...
        %             round(log10(VparGCp(15)))});
        %         xlabel('log$_{10}(v_{\perp 2})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        %         ylabel('log$_{10}(\mathbf{v_\parallel})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        %         xlabel(cbarqq0, horzcat('log$_{10}[f(v_{\perp 2}, \: v_\parallel)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
        %             'interpreter', 'latex', 'FontSize', 25)
        %         set(gca, 'xaxisLocation', 'top');
        %         set(gca, 'yminortick', 'on','xminortick', 'on');
        %         set(gca, 'FontSize', 25);
        %         set(gcf, 'color','white');
        %         grid on
        %
        %         subplot(2, 3, 5)
        %         if LOGREDflag == 1
        %             pp= contourf(F2PerpEVperp2VparRedlog10(:, :, Qind, nn)');
        %         else
        %             pp= contourf(F2PerpEVperp2VparRed(:, :, Qind, nn)');
        %         end
        %         %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        %         cbarqq0= colorbar('SouthOutside');
        %         if LOGREDflag == 1
        %             cmin= min(min(min(min(F2PerpEVperp2VparRedlog10(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpEVperp2VparRedlog10(:, :, :, :)))));
        %         else
        %             cmin= min(min(min(min(F2PerpEVperp2VparRed(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpEVperp2VparRed(:, :, :, :)))));
        %         end
        %         caxis([cmin cmax]);
        %         xlim([1 NVperp2G]);
        %         xticks([2 4 6]);
        %         xticklabels({round(EnergyGCp(2)*(6.242e18)), round(EnergyGCp(4)*(6.242e18)), ...
        %             round(EnergyGCp(6)*(6.242e18))});
        %         yticks([5 10 15]);
        %         yticklabels({round(GyroAngleGCp(5)), round(GyroAngleGCp(10)), ...
        %             round(GyroAngleGCp(15))});
        %         xlabel('$E$ [eV]', 'interpreter', 'latex', 'FontSize', 25)
        %         ylabel('$\theta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        %         xlabel(cbarqq0, horzcat('log$_{10}[f(E, \: \theta)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
        %             'interpreter', 'latex', 'FontSize', 25)
        %         set(gca, 'xaxisLocation', 'top');
        %         set(gca, 'yminortick', 'on','xminortick', 'on');
        %         set(gca, 'FontSize', 25);
        %         set(gcf, 'color','white');
        %         grid on
        %
        %         subplot(2, 3, 3)
        %         if LOGREDflag == 1
        %             pp= contourf(F2PerpVperp1Vperp2Redlog10(:, :, Qind, nn)');
        %         else
        %             pp= contourf(F2PerpVperp1Vperp2Red(:, :, Qind, nn)');
        %         end
        %         cbarqq0= colorbar('SouthOutside');
        %         if LOGREDflag == 1
        %             cmin= min(min(min(min(F2PerpVperp1Vperp2Redlog10(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpVperp1Vperp2Redlog10(:, :, :, :)))));
        %         else
        %             cmin= min(min(min(min(F2PerpVperp1Vperp2Red(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpVperp1Vperp2Red(:, :, :, :)))));
        %         end
        %         caxis([cmin cmax]);
        %         xlim([1 NVperp1G]);
        %         xticks([2 4 6]);
        %         xticklabels({round(log10(Vperp1GCp(2))), round(log10(Vperp1GCp(4))), ...
        %             round(log10(Vperp1GCp(6)))});
        %         ylim([1 NVperp2G]);
        %         yticks([2 4 6]);
        %         yticklabels({round(log10(Vperp2GCp(2))), round(log10(Vperp2GCp(4))), ...
        %             round(log10(Vperp2GCp(6)))});
        %         xlabel('log$_{10}(v_{\perp 1})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        %         ylabel('log$_{10}(v_{\perp 2})$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        %         xlabel(cbarqq0, horzcat('log$_{10}[f(v_{\perp 1}, \: v_{\perp 2})]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
        %             'interpreter', 'latex', 'FontSize', 25)
        %         set(gca, 'xaxisLocation', 'top');
        %         set(gca, 'yminortick', 'on','xminortick', 'on');
        %         set(gca, 'FontSize', 25);
        %         set(gcf, 'color','white');
        %         grid on
        %
        %         subplot(2, 3, 6)
        %         if LOGREDflag == 1
        %             pp= contourf(F2PerpEVperp1Vperp2Redlog10(:, :, Qind, nn)');
        %         else
        %             pp= contourf(F2PerpEVperp1Vperp2Red(:, :, Qind, nn)');
        %         end
        %         %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        %         cbarqq0= colorbar('SouthOutside');
        %         if LOGREDflag == 1
        %             cmin= min(min(min(min(F2PerpEVperp1Vperp2Redlog10(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpEVperp1Vperp2Redlog10(:, :, :, :)))));
        %         else
        %             cmin= min(min(min(min(F2PerpEVperp1Vperp2Red(:, :, :, :)))));
        %             cmax= max(max(max(max(F2PerpEVperp1Vperp2Red(:, :, :, :)))));
        %         end
        %         caxis([cmin cmax]);
        %         xlim([1 NVperp1G]);
        %         xticks([2 4 6]);
        %         xticklabels({round(EnergyGCp(2)*(6.242e18)), round(EnergyGCp(4)*(6.242e18)), ...
        %             round(EnergyGCp(6)*(6.242e18))});
        %         ylim([1 NVperp2G]);
        %         yticks([2 4 6]);
        %         yticklabels({round(PitchAngleGCp(2)), round(PitchAngleGCp(4)), ...
        %             round(PitchAngleGCp(6))});
        %         xlabel('log$_{10}(E)$ [eV]', 'interpreter', 'latex', 'FontSize', 25)
        %         ylabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        %         xlabel(cbarqq0, horzcat('log$_{10}[f(E, \: \alpha)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
        %             'interpreter', 'latex', 'FontSize', 25)
        %         set(gca, 'xaxisLocation', 'top');
        %         set(gca, 'yminortick', 'on','xminortick', 'on');
        %         set(gca, 'FontSize', 25);
        %         set(gcf, 'color','white');
        %         grid on
        
        % ----------------------------------------------------
        
        pause(pausevar)
        frame= getframe(gcf);
        
        writeVideo(v, frame);
        
        % ----------------------------------------------------
        
    end
end

close(v);

toc

disp('ION REDUCED VELOCITY DISTRIBUTION MOVIE COMPLETE')

%%
tic

clc
close all

% PLOT ION DIFFERENTIAL FLUXES VIDEO:

qVind= 5; %VISIONS-1 altitude

colorp= {'k.' 'r.' 'b.' 'm.' 'g.' 'y.' 'c.'};
color= repmat(colorp, 1, ranksize);
pausevar= 1e-9;
Vgridlogflag= 1;
linewidth= 2; fsize= 20;
LOGDIFFflag= 1;

fignum= 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1600 750])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION DIFFERENTIAL FLUXES VID.avi'));
open(v);

if IONVPERPVECflag == 1
    
    for nn= NNtT+ 1
        
        % ----------------------------------------------------
        %                 Qind= qVind;
        Qind= 2;
        
        % ----------------------------------------------------
        
        subplot(2, 3, 1)
        if LOGDIFFflag == 1
            pp= contourf(DiffNumberFluxVperp1VparRedlog10(:, :, Qind, nn));
        else
            pp= contourf(DiffNumberFluxVperp1VparRed(:, :, Qind, nn));
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGDIFFflag == 1
            cmin= min(min(min(min(DiffNumberFluxVperp1VparRedlog10(:, :, :, :)))));
            cmax= max(max(max(max(DiffNumberFluxVperp1VparRedlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(DiffNumberFluxVperp1VparRed(:, :, :, :)))));
            cmax= max(max(max(max(DiffNumberFluxVperp1VparRed(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        %         xticks([5 10 15]);
        %         xticklabels({round(log10(EnergyGCp(5))), 0e0, ...
        %             round(log10(EnergyGCp(15)))});
        %         ylim([1 NVperp1G]);
        %         xlim([1 NVperp1G]);
        %         yticks([2 4 6]);
        %         yticklabels({round(GyroAngleGCp(2)), round(GyroAngleGCp(4)), ...
        %             round(GyroAngleGCp(6))});
        xlabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\theta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[\phi_N(\alpha, \: \theta)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        
        subplot(2, 3, 2)
        if LOGREDflag == 1
            pp= contourf(DiffNumberFluxVperp2VparRedlog10(:, :, Qind, nn));
        else
            pp= contourf(DiffNumberFluxVperp2VparRed(:, :, Qind, nn));
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGREDflag == 1
            cmin= min(min(min(min(DiffNumberFluxVperp2VparRedlog10(:, :, :, :)))));
            cmax= max(max(max(max(DiffNumberFluxVperp2VparRedlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(DiffNumberFluxVperp2VparRed(:, :, :, :)))));
            cmax= max(max(max(max(DiffNumberFluxVperp2VparRed(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        xticks([5 10 15]);
        xticklabels({round(log10(EnergyGCp(5))), 0e0, ...
            round(log10(EnergyGCp(15)))});
        ylim([1 NVperp1G]);
        yticks([2 4 6]);
        yticklabels({round(GyroAngleGCp(2)), round(GyroAngleGCp(4)), ...
            round(GyroAngleGCp(6))});
        xlabel('log$_{10}(E)$ [J]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\theta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[\phi_N(E, \: \theta)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        
        subplot(2, 3, 3)
        if LOGREDflag == 1
            pp= contourf(DiffNumberFluxVperp1Vperp2Redlog10(:, :, Qind, nn));
        else
            pp= contourf(DiffNumberFluxVperp1Vperp2Red(:, :, Qind, nn));
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGREDflag == 1
            cmin= min(min(min(min(DiffNumberFluxVperp1Vperp2Redlog10(:, :, :, :)))));
            cmax= max(max(max(max(DiffNumberFluxVperp1Vperp2Redlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(DiffNumberFluxVperp1Vperp2Red(:, :, :, :)))));
            cmax= max(max(max(max(DiffNumberFluxVperp1Vperp2Red(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        xticks([5 10 15]);
        xticklabels({round(log10(EnergyGCp(5))), 0e0, ...
            round(log10(EnergyGCp(15)))});
        ylim([1 NVperp1G]);
        yticks([2 4 6]);
        yticklabels({round(PitchAngleGCp(2)), round(PitchAngleGCp(4)), ...
            round(PitchAngleGCp(6))});
        xlabel('log$_{10}(E)$ [J]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[\phi_N(E, \: \alpha)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        
        subplot(2, 3, 4)
        if LOGDIFFflag == 1
            pp= contourf(DiffEnergyFluxVperp1VparRedlog10(:, :, Qind, nn));
        else
            pp= contourf(DiffEnergyFluxVperp1VparRed(:, :, Qind, nn));
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGDIFFflag == 1
            cmin= min(min(min(min(DiffEnergyFluxVperp1VparRedlog10(:, :, :, :)))));
            cmax= max(max(max(max(DiffEnergyFluxVperp1VparRedlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(DiffEnergyFluxVperp1VparRed(:, :, :, :)))));
            cmax= max(max(max(max(DiffEnergyFluxVperp1VparRed(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        %         xticks([5 10 15]);
        %         xticklabels({round(log10(EnergyGCp(5))), 0e0, ...
        %             round(log10(EnergyGCp(15)))});
        %         ylim([1 NVperp1G]);
        %         xlim([1 NVperp1G]);
        %         yticks([2 4 6]);
        %         yticklabels({round(GyroAngleGCp(2)), round(GyroAngleGCp(4)), ...
        %             round(GyroAngleGCp(6))});
        xlabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\theta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[\phi_E(\alpha, \: \theta)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        
        subplot(2, 3, 5)
        if LOGREDflag == 1
            pp= contourf(DiffEnergyFluxVperp2VparRedlog10(:, :, Qind, nn));
        else
            pp= contourf(DiffEnergyFluxVperp2VparRed(:, :, Qind, nn));
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGREDflag == 1
            cmin= min(min(min(min(DiffEnergyFluxVperp2VparRedlog10(:, :, :, :)))));
            cmax= max(max(max(max(DiffEnergyFluxVperp2VparRedlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(DiffEnergyFluxVperp2VparRed(:, :, :, :)))));
            cmax= max(max(max(max(DiffEnergyFluxVperp2VparRed(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        xticks([5 10 15]);
        xticklabels({round(log10(EnergyGCp(5))), 0e0, ...
            round(log10(EnergyGCp(15)))});
        ylim([1 NVperp1G]);
        yticks([2 4 6]);
        yticklabels({round(GyroAngleGCp(2)), round(GyroAngleGCp(4)), ...
            round(GyroAngleGCp(6))});
        xlabel('log$_{10}(E)$ [J]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\theta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[\phi_E(E, \: \theta)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        
        subplot(2, 3, 6)
        if LOGREDflag == 1
            pp= contourf(DiffEnergyFluxVperp1Vperp2Redlog10(:, :, Qind, nn));
        else
            pp= contourf(DiffEnergyFluxVperp1Vperp2Red(:, :, Qind, nn));
        end
        %         text(xlim1- 8, ylim1- 5, '(a)', 'interpreter', 'latex', 'FontSize', 25)
        cbarqq0= colorbar('SouthOutside');
        if LOGREDflag == 1
            cmin= min(min(min(min(DiffEnergyFluxVperp1Vperp2Redlog10(:, :, :, :)))));
            cmax= max(max(max(max(DiffEnergyFluxVperp1Vperp2Redlog10(:, :, :, :)))));
        else
            cmin= min(min(min(min(DiffEnergyFluxVperp1Vperp2Red(:, :, :, :)))));
            cmax= max(max(max(max(DiffEnergyFluxVperp1Vperp2Red(:, :, :, :)))));
        end
        caxis([cmin cmax]);
        xticks([5 10 15]);
        xticklabels({round(log10(EnergyGCp(5))), 0e0, ...
            round(log10(EnergyGCp(15)))});
        ylim([1 NVperp1G]);
        yticks([2 4 6]);
        yticklabels({round(PitchAngleGCp(2)), round(PitchAngleGCp(4)), ...
            round(PitchAngleGCp(6))});
        xlabel('log$_{10}(E)$ [J]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        xlabel(cbarqq0, horzcat('log$_{10}[\phi_E(E, \: \alpha)]$', ' [s$^2 \cdot$ m$^{-5}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        
        % ----------------------------------------------------
        
        pause(pausevar)
        frame= getframe(gcf);
        
        writeVideo(v, frame);
        
        % ----------------------------------------------------
        
    end
end

close(v);

toc

disp('ION DIFFERENTIAL FLUXES MOVIE COMPLETE')

%%
tic

clc
close all

% PLOT 3D ION VELOCITY DISTRIBUTIONS VIDEO:

clear Vperp2iso Vperp1iso Vpariso PitchAngleiso Energyiso GyroAngleiso

% [Energyiso PitchAngleiso GyroAngleiso]= meshgrid(EnergyGCp(:, :, :).*(6.242e18), ...
%     PitchAngleGCp(:, :, :), GyroAngleGCp(:, :, :));
% [Vperp2iso Vperp1iso Vpariso]= meshgrid(Vperp2GCp(:).*1e-3, ...
%     Vperp1GCp(:).*1e-3, VparGCp(:).*1e-3);

qVind= 5; %VISIONS-1 altitude
Qindlen= ((NqUB- NqLB)+ 1);
for Qind= 1:1:Qindlen;
    rr(Qind)= Rank(r).Specie(s).FluxTube(f).rGCp(Qind);
end
nd(:)= Rank(r).Specie(s).FluxTube(f).TimeT(:);
[rrAB ndAB]= meshgrid(rr, nd);
pausevar= 1e-9;
linewidth= 2; fsize= 20;

ylim1= 0; %200;
ylim2= 0.8d3;

% Plot 3D energy/pitch angle space distribution functions:

fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1600 750])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION 3D ENERGY PITCH ANGLE DISTRIB FNC VID.avi'));
open(v);

for nn= NNtT
    for Qind= 1
        
        clear F2PerpERTplog10p DiffNumberFluxlog10p DiffEnergyFluxlog10p
        
        F2PerpERTplog10p(:, :, :)= Rank(root).Specie(s).FluxTube(f).QCell(Qind).F2PerpERTplog10(:, :, :, nn);
        DiffNumberFluxlog10p(:, :, :)= Rank(root).Specie(s).FluxTube(f).QCell(Qind).DiffNumberFluxlog10(:, :, :, nn);
        DiffEnergyFluxlog10p(:, :, :)= Rank(root).Specie(s).FluxTube(f).QCell(Qind).DiffEnergyFluxlog10(:, :, :, nn);
        
        for Vperp1ind= 1:1:NV
            if (Vperp1ind > NVperp1G)
                F2PerpERTplog10p(Vperp1ind, :, :)= NaN;
                DiffNumberFluxlog10p(Vperp1ind, :, :)= NaN;
                DiffEnergyFluxlog10p(Vperp1ind, :, :)= NaN;
            else
                F2PerpERTplog10p(Vperp1ind, :, :)= F2PerpERTplog10p(Vperp1ind, :, :);
                DiffNumberFluxlog10p(Vperp1ind, :, :)= DiffNumberFluxlog10p(Vperp1ind, :, :);
                DiffEnergyFluxlog10p(Vperp1ind, :, :)= DiffEnergyFluxlog10p(Vperp1ind, :, :);
            end
        end
        for Vperp2ind= 1:1:NV
            if (Vperp2ind > NVperp2G)
                F2PerpERTplog10p(:, Vperp2ind, :)= NaN;
                DiffNumberFluxlog10p(:, Vperp2ind, :)= NaN;
                DiffEnergyFluxlog10p(:, Vperp2ind, :)= NaN;
            else
                F2PerpERTplog10p(:, Vperp2ind, :)= F2PerpERTplog10p(:, Vperp2ind, :);
                DiffNumberFluxlog10p(:, Vperp2ind, :)= DiffNumberFluxlog10p(:, Vperp2ind, :);
                DiffEnergyFluxlog10p(:, Vperp2ind, :)= DiffEnergyFluxlog10p(:, Vperp2ind, :);
            end
        end
        for Vparind= 1:1:NV
            if (Vparind > NVparG)
                F2PerpERTplog10p(:, :, Vparind)= NaN;
                DiffNumberFluxlog10p(:, :, Vparind)= NaN;
                DiffEnergyFluxlog10p(:, :, Vparind)= NaN;
            else
                F2PerpERTplog10p(:, :, Vparind)= F2PerpERTplog10p(:, :, Vparind);
                DiffNumberFluxlog10p(:, :, Vparind)= DiffNumberFluxlog10p(:, :, Vparind);
                DiffEnergyFluxlog10p(:, :, Vparind)= DiffEnergyFluxlog10p(:, :, Vparind);
            end
        end
        
        subplot(2, 3, 1)
        isosurface(EnergyGCp, PitchAngleGCp, GyroAngleGCp, F2PerpERTplog10p);
        % isonormals(FPerpNormlog10p, iso)
        view(3)
        camlight;
        lighting gouraud;
        xlabel('$E$ [eV]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        zlabel('$\vartheta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        title(horzcat('log$_{10}[f(E, \alpha, \vartheta)]$', ' [s$^3 \cdot$ m$^{-6} \cdot$ sR$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        % hold off
        
        subplot(2, 3, 2)
        isosurface(Energyiso, PitchAngleiso, GyroAngleiso, DiffNumberFluxlog10p);
        % isonormals(FPerpNormlog10p, iso)
        view(3)
        camlight;
        lighting gouraud;
        xlabel('$E$ [eV]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        zlabel('$\vartheta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        title(horzcat('log$_{10}(\phi_N)$', ' [s$^2 \cdot$ m$^{-5} \cdot$ sR$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        % hold off
        
        subplot(2, 3, 3)
        isosurface(Energyiso, PitchAngleiso, GyroAngleiso, DiffEnergyFluxlog10p);
        % isonormals(FPerpNormlog10p, iso)
        view(3)
        camlight;
        lighting gouraud;
        xlabel('$E$ [eV]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$\alpha$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        zlabel('$\vartheta$ [degs]', 'interpreter', 'latex', 'FontSize', 25)
        title(horzcat('log$_{10}(\phi_E)$', ' [J$ \cdot$ s$^2 \cdot$ m$^{-5} \cdot$ sR$^{-1}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        % hold off
        
        subplot(2, 3, 4)
        plot(EnergyM(nn, :)*(6.242e18), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(EnergyM(nn, Qind)*(6.242e18), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$E$ [eV]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (Wlimflag == 1)
            xlim([Wlim1 Wlim2]);
        end
        hold off
        
        subplot(2, 3, 5)
        plot(PitchAngleM(nn, :), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(PitchAngleM(nn, Qind), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$\alpha$ [degs]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        %         if (UPerp2limflag == 1)
        %             xlim([UPerp2lim1 UPerp2lim2]);
        %         end
        hold off
        
        subplot(2, 3, 6)
        plot(GyroAngleM(nn, :), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(GyroAngleM(nn, Qind), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$\theta$ [degs]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        %         if (UPerp2limflag == 1)
        %             xlim([UPerp2lim1 UPerp2lim2]);
        %         end
        hold off
        
        pause(pausevar)
        frame= getframe(gcf);
        
        writeVideo(v, frame);
        
    end
end

close(v);

% Plot 3D velocity space distribution functions:

fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1600 750])

set(gca, 'nextplot', 'replacechildren');
v= VideoWriter(horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION 3D PHASE SPACE DISTRIB FNC VID.avi'));
open(v);

for nn= NNtT
    for Qind= 1
        
        FPerpNormlog10p(:, :, :)= Rank(1).Specie(1).FluxTube(1).QCell(Qind).F2PerpNormlog10(:, :, :, nn);
        N2PerpphReNormRTlog10p(:, :, :)= Rank(1).Specie(1).FluxTube(1).QCell(Qind).N2PerpphReNormRTlog10(:, :, :, nn);
        N2PerpphReNormERRORRTplog10p(:, :, :)= Rank(1).Specie(1).FluxTube(1).QCell(Qind).N2PerpphReNormERRORRTplog10(:, :, :, nn);
        
        subplot(2, 3, 1)
        isosurface(Vperp2iso, Vperp1iso, Vpariso, FPerpNormlog10p);
        % isonormals(FPerpNormlog10p, iso)
        view(3)
        camlight;
        lighting gouraud;
        xlabel('$v_{\perp 2}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$v_{\perp 1}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        zlabel('$v_\parallel$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        title(horzcat('log$_{10}[f(v_{\perp 1}, v_{\perp 2}, v_\parallel)]$', ' [s$^3 \cdot$ m$^{-6}$]'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        % hold off
        
        subplot(2, 3, 2)
        isosurface(Vperp2iso, Vperp1iso, Vpariso, N2PerpphReNormRTlog10p);
        % isonormals(FPerpNormlog10p, iso)
        view(3)
        camlight;
        lighting gouraud;
        xlabel('$v_{\perp 2}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$v_{\perp 1}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        zlabel('$v_\parallel$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        title(horzcat('log$_{10}(\mathcal{N^\prime})$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        % hold off
        
        subplot(2, 3, 3)
        isosurface(Vperp2iso, Vperp1iso, Vpariso, N2PerpphReNormERRORRTplog10p);
        % isonormals(FPerpNormlog10p, iso)
        view(3)
        camlight;
        lighting gouraud;
        xlabel('$v_{\perp 2}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        ylabel('$v_{\perp 1}$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        zlabel('$v_\parallel$ [km/s]', 'interpreter', 'latex', 'FontSize', 25)
        title(horzcat('log$_{10}(\mathcal{N^\prime}^{-1/2})$'), ...
            'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top');
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gca, 'FontSize', 25);
        set(gcf, 'color','white');
        grid on
        % hold off
        
        subplot(2, 3, 4)
        plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_{\perp 2}$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp2limflag == 1)
            xlim([UPerp2lim1 UPerp2lim2]);
        end
        hold off
        
        subplot(2, 3, 5)
        plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(e)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_{\perp 1}$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UPerp1limflag == 1)
            xlim([UPerp1lim1 UPerp1lim2]);
        end
        hold off
        
        subplot(2, 3, 6)
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, :)*1e-3), rrAB, 'k', 'LineWidth', linewidth)
        hold on;
        plot(real(-Rank(r).Specie(s).FluxTube(f).M1ParphRTr(nn, Qind)*1e-3), rrAB(1, Qind), 'r*', 'LineWidth', linewidth, 'MarkerSize', 20)
        set(gca, 'FontSize', fsize);
        %     if figlabelflag == 1
        %         text(textx, texty, '(f)', 'interpreter', 'latex', 'FontSize', 25)
        %     end
        xlabel(horzcat('$u_\parallel$ [km/s]'), 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'xaxisLocation', 'top')
        ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
        set(gca, 'yminortick', 'on','xminortick', 'on');
        set(gcf, 'color','white');
        ylim([ylim1 ylim2]);
        if (UParlimflag == 1)
            xlim([UParlim1 UParlim2]);
        end
        hold off
        
        pause(pausevar)
        frame= getframe(gcf);
        
        writeVideo(v, frame);
        
    end
end

close(v);

toc

disp('3D ION DISTRIBUTION FUNCTIONS PLOTS COMPLETE')

%% -------------------------------------------------------

clc
close all

tic

% BOUNDARY CONDITIONS PLOTS:
Qindlen= ((NqUB- NqLB)+ 1);
for Qind= 1:1:Qindlen;
    rr(Qind)= Rank(r).Specie(s).FluxTube(f).rGCp(Qind);
end
nd(:)= Rank(r).Specie(s).FluxTube(f).TimeT(:);
[rrAB ndAB]= meshgrid(rr, nd);

Nend= NNtT+ 1;
for nn= 1:1:Nend
    
    NsnTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'NsnTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).NsnT(nn)= fread(NsnTID, 1, 'integer*8');
    fclose(NsnTID);
    
    NsnRRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'NsnRRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).NsnRRT(nn)= fread(NsnRRTID, 1, 'integer*8');
    fclose(NsnRRTID);
    
    NqLBoutfluxIonRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'NqLBoutfluxIonRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).NqLBoutfluxIonRT(nn)= fread(NqLBoutfluxIonRTID, 1, 'real*8');
    fclose(NqLBoutfluxIonRTID);
    
    LBoutfluxIonRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'LBoutfluxIonRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).LBoutfluxIonRT(nn)= fread(LBoutfluxIonRTID, 1, 'real*8');
    fclose(LBoutfluxIonRTID);
    
    NqUBoutfluxIonRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'NqUBoutfluxIonRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).NqUBoutfluxIonRT(nn)= fread(NqUBoutfluxIonRTID, 1, 'real*8');
    fclose(NqUBoutfluxIonRTID);
    
    UBoutfluxIonRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'UBoutfluxIonRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).UBoutfluxIonRT(nn)= fread(UBoutfluxIonRTID, 1, 'real*8');
    fclose(UBoutfluxIonRTID);
    
    LBNetDensityTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'LBNetDensityTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).LBNetDensityT(nn)= fread(LBNetDensityTID, 1, 'real*8');
    fclose(LBNetDensityTID);
    
    UBNetDensityTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(nn), '_', num2str(s), '_', num2str(f), '_', ...
        'UBNetDensityTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).UBNetDensityT(nn)= fread(UBNetDensityTID, 1, 'real*8');
    fclose(UBNetDensityTID);
        
end
fignum= 0;
fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1500 950])
r= root;
s= 1;
f= 1;

clear dNs1 dNs2 dNs3 NsnT
dNs1(:)= Rank(r).Specie(s).FluxTube(f).NqLBoutfluxIonRT(:)/Rank(r).Specie(s).FluxTube(f).nsnormfacT+ ...
    Rank(r).Specie(s).FluxTube(f).NqUBoutfluxIonRT(:)/Rank(r).Specie(s).FluxTube(f).nsnormfacT;
dNs2(:)= Rank(r).Specie(s).FluxTube(f).LBNetDensityT(:);
dNs3(:)= Rank(r).Specie(s).FluxTube(f).UBNetDensityT(:);
NsnT(:)= Rank(r).Specie(s).FluxTube(f).NsnT(:);

subplot(3, 2, 1)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).NqLBoutfluxIonRT(1:size(ndAB(:, 1)))/Rank(r).Specie(s).FluxTube(f).nsnormfacT, 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
title('Boundary Escape Fluxes', 'interpreter', 'latex', 'FontSize', 25);
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^o_{LB}$ (unitless)'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 2)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).LBoutfluxIonRT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$j^o_{LB}$ [m$^{-2}$ $\cdot$ s$^{-1}$]'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 3)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).NqUBoutfluxIonRT(1:size(ndAB(:, 1)))/Rank(r).Specie(s).FluxTube(f).nsnormfacT, 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^o_{UB}$ (unitless)'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 4)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).UBoutfluxIonRT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$j^o_{UB}$ [m$^{-2}$ $\cdot$ s$^{-1}$]'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 5)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).NqLBoutfluxIonRT(1:size(ndAB(:, 1)))/Rank(r).Specie(s).FluxTube(f).nsnormfacT+ ...
    Rank(r).Specie(s).FluxTube(f).NqUBoutfluxIonRT(1:size(ndAB(:, 1)))/Rank(r).Specie(s).FluxTube(f).nsnormfacT, 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^o_{LB}+ N^o_{UB}$ (unitless)'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 6)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).LBoutfluxIonRT(1:size(ndAB(:, 1)))+ ...
    Rank(r).Specie(s).FluxTube(f).UBoutfluxIonRT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$j^o_{LB}+ j^o_{UB}$ [m$^{-2}$ $\cdot$ s$^{-1}$]'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' BOUNDARY CONDITIONS PLOT 1.png')]);

fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1500 950])

subplot(3, 2, 1)
plot(ndAB(:, 1)/60, ranksize*Rank(r).Specie(s).FluxTube(f).LBNetDensityT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
hold on
%         plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).LBNominalNormDensityT(1)* ...
%             Rank(r).Specie(s).FluxTube(f).d3xCLBp/Rank(r).Specie(s).FluxTube(f).nsnormfacT, ...
%             'r*', 'MarkerSize', 20, 'LineWidth', 2)
grid on
title('Total Cluster Particle Numbers', 'interpreter', 'latex', 'FontSize', 25);
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^i_{LB}$'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 2)
plot(ndAB(:, 1)/60, ranksize*Rank(r).Specie(s).FluxTube(f).UBNetDensityT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
hold on
%         plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).UBNominalNormDensityT(1)* ...
%             Rank(r).Specie(s).FluxTube(f).d3xCUBp/ ...
%             Rank(r).Specie(s).FluxTube(f).nsnormfacT, 'rx', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^i_{UB}$'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 3)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).NqLBoutfluxIonRT(1:size(ndAB(:, 1)))/Rank(r).Specie(s).FluxTube(f).nsnormfacT+ ...
    Rank(r).Specie(s).FluxTube(f).NqUBoutfluxIonRT(1:size(ndAB(:, 1)))/Rank(r).Specie(s).FluxTube(f).nsnormfacT, 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^o_{LB}+ N^o_{UB}$ (unitless)'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 5)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).NsnT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N_s$ (unitless)'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(3, 2, 6)
plot(ndAB(:, 1)/60, Rank(r).Specie(s).FluxTube(f).NsnRRT(1:size(ndAB(:, 1))), 'k', 'MarkerSize', 20, 'LineWidth', 2)
hold on
% plot(ndAB(:, 1)/60, Rank(root).Specie(s).FluxTube(f).TotalParticles(1:size(ndAB(:, 1))), 'k-.', 'MarkerSize', 20, 'LineWidth', 2)
grid on
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$N^C_s$ (unitless)'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' BOUNDARY CONDITIONS PLOT 2.png')]);

toc

disp('BOUNDARY CONDITIONS PLOTS COMPLETE')

%% -------------------------------------------------------

clc
close all

tic

nnind= 1; % NNtT+ 1;

% INITIAL DENSITY PLOTS:

s= 1;
f= 1;
NqIC= abs(NqICB- NqICA)+ 1e0;

fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1500 1500])

for Qind= 1:1:NqIC
    Rank(1).Specie(s).FluxTube(f).rGCICpp(Qind)= (Rank(1).Specie(s).FluxTube(f).QCell(Qind).rGC- RE)/RE;
end

rrIC(:)= Rank(1).Specie(s).FluxTube(f).rGCICpp(:);

for Qind= 1:1:NqIC
    
    Tearg1(Qind)= Rank(root).Specie(s).FluxTube(f).TeNT(nnind);
    
    Temparg1(Qind)= Rank(root).Specie(s).FluxTube(f).Temp(nnind, Qind);
    PerpTemparg1(Qind)= Rank(root).Specie(s).FluxTube(f).PerpTemp(nnind, Qind);
    ParTemparg1(Qind)= Rank(root).Specie(s).FluxTube(f).ParTemp(nnind, Qind);
    
    Temparg2(Qind)= Rank(root).Specie(s).FluxTube(f).Tsp(Qind);
    PerpTemparg2(Qind)= Rank(root).Specie(s).FluxTube(f).TsPerpp(Qind);
    ParTemparg2(Qind)= Rank(root).Specie(s).FluxTube(f).TsParp(Qind);
    
end

subplot(1, 3, 1)
plot(rrIC(:), PerpTemparg2(:), 'gd-', 'MarkerSize', 20, 'LineWidth', 2);
hold on
plot(rrIC(:), PerpTemparg1(:), 'bo-', 'MarkerSize', 20, 'LineWidth', 2)
grid on
title(horzcat('$t= $', num2str(Rank(r).Specie(s).FluxTube(f).TimeT(nnind)), ...
    ' [min]'), 'interpreter', 'latex', 'FontSize', 25);    xlabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$T_\perp$ [K]'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(1, 3, 2)
plot(rrIC(:), ParTemparg2(:), 'gd-', 'MarkerSize', 20, 'LineWidth', 2);
hold on
plot(rrIC(:), ParTemparg1(:), 'bo-', 'MarkerSize', 20, 'LineWidth', 2)
grid on
title(horzcat('$t= $', num2str(Rank(r).Specie(s).FluxTube(f).TimeT(nnind)), ...
    ' [min]'), 'interpreter', 'latex', 'FontSize', 25);    xlabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$T_\parallel$ [K]'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

subplot(1, 3, 3)
plot(rrIC(:), Temparg2(:), 'gd-', 'MarkerSize', 20, 'LineWidth', 2);
hold on
plot(rrIC(:), Temparg1(:), 'bo-', 'MarkerSize', 20, 'LineWidth', 2)
grid on
title(horzcat('$t= $', num2str(Rank(r).Specie(s).FluxTube(f).TimeT(nnind)), ...
    ' [min]'), 'interpreter', 'latex', 'FontSize', 25);    xlabel('$r$ [$R_E$]', 'interpreter', 'latex', 'FontSize', 25);
ylabel(horzcat('$T$ [K]'), 'interpreter', 'latex', 'FontSize', 25);
set(gca, 'FontSize', 25);
set(gcf, 'color','white');
hold off

saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' INITIAL DENSITY PLOT.png')]);

toc

disp('INITIAL DENSITY PLOT COMPLETE')

%% ----------------------------------------------------
clc
close all

nn= 1;
plot(abs(Rank(root).Specie(s).FluxTube(f).PEg(nn, :)), rr, 'kx')
hold on
plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn, :)+ ...
    Rank(root).Specie(s).FluxTube(f).PEa(nn, :)), rr, 'bx')
ylim([500 1d3])

%%

clc
close all

tic

% PLOT ION POTENTIAL WELL PLOTS:

GRAVITYflag= 1;
THERMALflag= 1;
AMBIPOLARflag= 1;
MIRRORflag= 1;
EPARALLELflag= 1;
PEaInertialflag= 1;
elecflag= 1;
Qindlen= NqUB- NqLB+ 1;
Nqend= Qindlen;

ylim1= 500; ylim2= 0.8d3;

if MIRRORflag == 0
    xlim1= 0; xlim2= 0.5e0;
else
    %     xlim1= 0; xlim2= 25e0;
    xlim1= 0; xlim2= 5e0;
end

for Qind= 1:1:Nqend;
    rr2(Qind)= Rank(r).Specie(s).FluxTube(f).rGCp(Qind);
end
nd2(:)= Rank(r).Specie(s).FluxTube(f).TimeT(:);
[rrAB2 ndAB2]= meshgrid(rr2, nd2);

GG= 667e-13; % Universal gravitational constant [N m^2/kg^2]
ME= 598e22; % Earth mass [kg]
qion= 1.602e-19; % Charge [C]
mEarth= 8.06e15; % Earth magnetic moment [T m^3]
mu0= 1.257e-6; % Permeability of Free Space [m kg s^-2 A^-2]
eps0= 8.854e-12; % Permittivity of Free Space [m^-3 kg^-1 s^4 A^2]
mion= (16e0)*(1.67e-27); % Mass [kg]
zns0= RE+ 370.78e3; % For VISIONS-1 Case Study
ns0= 6.678e10;
Qindns0= 3; % Rank(1).Specie(1).Qindns0;
nsnormfac= 6.2e17;
Tfac= 1e0;
r= root;
s= 1; f= 1;

for nn= 1:1:NNtT+ 1
    EPmagRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(s), '_', ...
        num2str(f), '_', 'EPmagRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).EPmagRTr(:)= fread(EPmagRTID, ...
        (NqUB- NqLB)+ 1, 'real*8');
    fclose(EPmagRTID);
    EGmagRTID= fopen(horzcat(datadir, num2str(r- 1), '_', num2str(s), '_', ...
        num2str(f), '_', 'EGmagRTfort.bin'));
    Rank(r).Specie(s).FluxTube(f).EGmagRTr(:)= fread(EGmagRTID, ...
        (NqUB- NqLB)+ 1, 'real*8');
    fclose(EGmagRTID);
end

Rank(r).Specie(s).FluxTube(f).Bmag(:)= ...
    mEarth.*sqrt(Rank(r).Specie(s).FluxTube(f).ellGC(:))./ ...
    (Rank(root).Specie(s).FluxTube(f).rGC(:).^3);

fG(:)= (qion.*Rank(r).Specie(s).FluxTube(f).Bmag(:)./mion)./(2e0*pi);

for nn= 1:1:NNtT+ 1
    for Qind= NqLB:1:NqUB
        rhoG(nn, Qind)= sqrt(Rank(r).Specie(s).FluxTube(f).M1Perp1phRTr(nn, Qind)^2e0+ ...
            Rank(r).Specie(s).FluxTube(f).M1Perp2phRTr(nn, Qind)^2e0)/fG(Qind);
    end
end

fG1= (qion*Rank(r).Specie(s).FluxTube(f).Bmag(1)/mion)/(2e0*pi)
fGend= (qion*Rank(r).Specie(s).FluxTube(f).Bmag(end)/mion)/(2e0*pi)
tauG1= 1e0/fG1;
tauGend= 1e0/fGend;

% Compute gravitational potential energy well
if GRAVITYflag == 1
    for nn= 1:1:NNtT+ 1
        Rank(root).Specie(s).FluxTube(f).PEgA(nn, :)= ...
            Rank(r).Specie(s).FluxTube(f).EGmagRTr(:).*(qion/mion);

        Rank(root).Specie(s).FluxTube(f).PEgArg(nn, :)= ...
            Rank(root).Specie(s).FluxTube(f).PEgA(nn, :)'.* ...
            Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
            Rank(root).Specie(s).FluxTube(f).dqC(:);
        for Qind= 1:1:Nqend
            Rank(root).Specie(s).FluxTube(f).PEgSum(nn, Qind)= ...
                sum(Rank(root).Specie(s).FluxTube(f).PEgArg(nn, 1:Qind));
        end
        PEg(nn, :)= ...
            ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEgSum(nn, :);
    end
end
    
% Compute gravitational potential energy well
% if GRAVITYflag == 1
%     gpar(:)= (2e0*GG*ME).*cos(Rank(root).Specie(s).FluxTube(f).thetaGC(:))./ ...
%         ((Rank(r).Specie(s).FluxTube(f).rGC(:).^2e0).* ...
%         sqrt(Rank(root).Specie(s).FluxTube(f).ellGC(:)));
%     Fg(:)= mion.*gpar(:);
%     FgArg(:)= Fg(:).*Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
%         Rank(root).Specie(s).FluxTube(f).dqC(:);
%     for Qind= 1:1:Nqend
%         FgSum(Qind)= sum(FgArg(1:Qind));
%     end
%     for nn= 1:1:NNtT+ 1
%         PEg(nn, :)= (6.242e18).*FgSum(:);
%     end
% end

for nn= 1:1:NNtT+ 1
    Densp(nn, :)= Rank(1).Specie(s).FluxTube(f).M0phRTr(nn, :);
    if (MIRRORflag == 1) | (EPARALLELflag == 1)
        Tempp(nn, :)= Rank(1).Specie(s).FluxTube(f).Temp(nn, :);
    else
%         Tempp(nn, :)= Rank(root).Specie(s).FluxTube(f).TsParp(:);
        Tempp(nn, :)= Rank(1).Specie(s).FluxTube(f).Temp(nn, :);
    end
end
    
% Compute thermal potential energy well
if THERMALflag == 1
    for nn= 1:1:NNtT+ 1
        for Qind= NqLB:1:NqUB
            if Qind == NqLB            
                dn(nn, Qind)= abs(-25e0*Densp(nn, Qind)+ 48e0*Densp(nn, Qind+ 1) ...
                    - 36e0*Densp(nn, Qind+ 2)+ 16e0*Densp(nn, Qind+ 3) ...
                    - 3e0*Densp(nn, Qind+ 4));
                dq(Qind)= 12e0*abs(Rank(root).Specie(s).FluxTube(f).dqC(Qind));
                dnds(nn, Qind)= (1e0/Rank(root).Specie(s).FluxTube(f).hqC(Qind))* ...
                    (dn(nn, Qind)/dq(Qind));
                
                Ft(nn, Qind)= ((1e0*kB*Tempp(nn, Qind))/ ...
                    (Densp(nn, Qind)))* ...
                    (dnds(nn, Qind));
            end
            if Qind == NqUB         
                dn(nn, Qind)= abs(-25e0*Densp(nn, Qind)+ 48e0*Densp(nn, Qind- 1) ...
                    - 36e0*Densp(nn, Qind- 2)+ 16e0*Densp(nn, Qind- 3) ...
                    - 3e0*Densp(nn, Qind- 4));
                dq(Qind)= -12e0*abs(Rank(root).Specie(s).FluxTube(f).dqC(Qind));
                dnds(nn, Qind)= (1e0/Rank(root).Specie(s).FluxTube(f).hqC(Qind))* ...
                    (dn(nn, Qind)/dq(Qind));
                
                Ft(nn, Qind)= ((1e0*kB*Tempp(nn, Qind))/ ...
                    (Densp(nn, Qind)))* ...
                    (dnds(nn, Qind));
            end
            if (Qind == NqLB+ 1) | (Qind == NqUB- 1)     
                dn(nn, Qind)= abs(Densp(nn, Qind+ 1)- Densp(nn, Qind- 1));
                dq(Qind)= 2e0*abs(Rank(root).Specie(s).FluxTube(f).dqC(Qind));
                dnds(nn, Qind)= (1e0/Rank(root).Specie(s).FluxTube(f).hqC(Qind))* ...
                    (dn(nn, Qind)/dq(Qind));
                
                Ft(nn, Qind)= ((1e0*kB*Tempp(nn, Qind))/ ...
                    (Densp(nn, Qind)))* ...
                    (dnds(nn, Qind));
            end
            if (Qind ~= NqLB) & (Qind ~= NqUB) & (Qind ~= NqLB+ 1) & (Qind ~= NqUB- 1)     
                dn(nn, Qind)= abs(Densp(nn, Qind- 2)- 8e0*Densp(nn, Qind- 1) ...
                    + 8e0*Densp(nn, Qind+ 1)- Densp(nn, Qind+ 2));
                dq(Qind)= 12e0*abs(Rank(root).Specie(s).FluxTube(f).dqC(Qind));
                dnds(nn, Qind)= (1e0/Rank(root).Specie(s).FluxTube(f).hqC(Qind))* ...
                    (dn(nn, Qind)/dq(Qind));
                
                Ft(nn, Qind)= ((1e0*kB*Tempp(nn, Qind))/ ...
                    (Densp(nn, Qind)))* ...
                    (dnds(nn, Qind));
            end
            
            if isnan(Ft(nn, Qind)) == 1 | isinf(Ft(nn, Qind)) == 1
                Ft(nn, Qind)= 0e0;
            end
            
            % Compute thermal potential energy well
            FtArg(nn, Qind)= Ft(nn, Qind)*Rank(root).Specie(s).FluxTube(f).hqC(Qind)* ...
                Rank(root).Specie(s).FluxTube(f).dqC(Qind);
            FtSum(nn, Qind)= sum(FtArg(nn, 1:Qind));
            PEt(nn, Qind)= (6.242e18)*FtSum(nn, Qind);
        
        end
    end
end

% if THERMALflag == 1
%     for nn= 1:1:NNtT+ 1
%         % Compute thermal potential energy well
%         Rank(root).Specie(s).FluxTube(f).PEtA(nn, :)= ...
%             Rank(r).Specie(s).FluxTube(f).EThermalRTr(nn, :).*(qion/mion);
% 
%         Rank(root).Specie(s).FluxTube(f).PEtArg(nn, :)= ...
%             Rank(root).Specie(s).FluxTube(f).PEtA(nn, :)'.* ...
%             Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
%             Rank(root).Specie(s).FluxTube(f).dqC(:);
%         for Qind= 1:1:Nqend
%             Rank(root).Specie(s).FluxTube(f).PEtSum(nn, Qind)= ...
%                 sum(Rank(root).Specie(s).FluxTube(f).PEtArg(nn, 1:Qind));
%         end
%         PEt(nn, :)= ...
%             ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEtSum(nn, :);
%     end
% end

for nn= 1:1:NNtT+ 1
    % Compute parallel electric potential energy well
    if EPARALLELflag == 1
        Rank(root).Specie(s).FluxTube(f).PEpA(nn, :)= ...
            Rank(r).Specie(s).FluxTube(f).EPmagRTr(:).*(qion/mion);

        Rank(root).Specie(s).FluxTube(f).PEpArg(nn, :)= ...
            Rank(root).Specie(s).FluxTube(f).PEpA(nn, :)'.* ...
            Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
            Rank(root).Specie(s).FluxTube(f).dqC(:);
        for Qind= 1:1:Nqend
            Rank(root).Specie(s).FluxTube(f).PEpSum(nn, Qind)= ...
                sum(Rank(root).Specie(s).FluxTube(f).PEpArg(nn, 1:Qind));
        end
        PEp(nn, :)= ...
            ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEpSum(nn, :);
    end
            
    % Compute mirror force potential energy well
    if MIRRORflag == 1
        if IONVPERPVECflag == 1
            Rank(root).Specie(s).FluxTube(f).mu(nn, :)= ...
                (mion.*(Rank(root).Specie(s).FluxTube(f).M1Perp1phRTr(nn, :).^2e0+ ...
                Rank(root).Specie(s).FluxTube(f).M1Perp2phRTr(nn, :).^2e0))./ ...
                (2e0.*Rank(root).Specie(s).FluxTube(f).Bmag(:)');
        else
            Rank(root).Specie(s).FluxTube(f).mu(nn, :)= ...
                (mion.*Rank(root).Specie(s).FluxTube(f).M1PerpphRTr(nn, :).^2e0)./ ...
                (2e0.*Rank(root).Specie(s).FluxTube(f).Bmag(:));
        end
        Rank(root).Specie(s).FluxTube(f).dBds(nn, :)= ...
            (((-6e0*mEarth).*cos(Rank(root).Specie(s).FluxTube(f).thetaGC(:)))./ ...
            (Rank(root).Specie(s).FluxTube(f).rGC(Qind).^4e0))- ...
            ((3e0*mEarth).*cos(Rank(root).Specie(s).FluxTube(f).thetaGC(:)).* ...
            (sin(Rank(root).Specie(s).FluxTube(f).thetaGC(:)).^2e0)./ ...
            (Rank(root).Specie(s).FluxTube(f).ellGC(:).* ...
            Rank(root).Specie(s).FluxTube(f).rGC(:).^4e0));
        Rank(root).Specie(s).FluxTube(f).PEmA(nn, :)= ...
            abs((-Rank(root).Specie(s).FluxTube(f).mu(nn, :)./mion).* ...
            Rank(root).Specie(s).FluxTube(f).dBds(nn, :));
        Rank(root).Specie(s).FluxTube(f).PEmArg(nn, :)= ...
            Rank(root).Specie(s).FluxTube(f).PEmA(nn, :).* ...
            Rank(root).Specie(s).FluxTube(f).hqC(:)'.* ...
            Rank(root).Specie(s).FluxTube(f).dqC(:)';
        for Qind= 1:1:Nqend
            Rank(root).Specie(s).FluxTube(f).PEmSum(nn, Qind)= ...
                sum(Rank(root).Specie(s).FluxTube(f).PEmArg(nn, 1:Qind));
        end
        PEm(nn, :)= ...
            ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEmSum(nn, :);
    end
        
    % Compute ambipolar electric potential energy well
    if AMBIPOLARflag == 1
        Rank(root).Specie(s).FluxTube(f).PEaAPressure(nn, :)= ...
            Rank(r).Specie(s).FluxTube(f).EAPressureRTr(nn, :).*(qion/mion);
        if PEaInertialflag == 1
            Rank(root).Specie(s).FluxTube(f).PEaAInertial(nn, :)= ...
                Rank(r).Specie(s).FluxTube(f).EAInertialRTr(nn, :).*(qion/mion);
        end
        Rank(root).Specie(s).FluxTube(f).PEaA(nn, :)= ...
            Rank(r).Specie(s).FluxTube(f).EAMagRTr(nn, :).*(qion/mion);

        Rank(root).Specie(s).FluxTube(f).PEaArgPressure(nn, :)= ...
            Rank(root).Specie(s).FluxTube(f).PEaAPressure(nn, :)'.* ...
            Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
            Rank(root).Specie(s).FluxTube(f).dqC(:);
        if PEaInertialflag == 1
            Rank(root).Specie(s).FluxTube(f).PEaArgInertial(nn, :)= ...
                Rank(root).Specie(s).FluxTube(f).PEaAInertial(nn, :)'.* ...
                Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
                Rank(root).Specie(s).FluxTube(f).dqC(:);
        end
        Rank(root).Specie(s).FluxTube(f).PEaArg(nn, :)= ...
            Rank(root).Specie(s).FluxTube(f).PEaA(nn, :)'.* ...
            Rank(root).Specie(s).FluxTube(f).hqC(:).* ...
            Rank(root).Specie(s).FluxTube(f).dqC(:);
        for Qind= 1:1:Nqend
            Rank(root).Specie(s).FluxTube(f).PEaSumPressure(nn, Qind)= ...
                sum(Rank(root).Specie(s).FluxTube(f).PEaArgPressure(nn, 1:Qind));
            if PEaInertialflag == 1
                Rank(root).Specie(s).FluxTube(f).PEaSumInertial(nn, Qind)= ...
                    sum(Rank(root).Specie(s).FluxTube(f).PEaArgInertial(nn, 1:Qind));
            end
            Rank(root).Specie(s).FluxTube(f).PEaSum(nn, Qind)= ...
                sum(Rank(root).Specie(s).FluxTube(f).PEaArg(nn, 1:Qind));
        end
        Rank(root).Specie(s).FluxTube(f).PEaPressure(nn, :)= ...
            ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEaSumPressure(nn, :);
        if PEaInertialflag == 1
            Rank(root).Specie(s).FluxTube(f).PEaInertial(nn, :)= ...
                ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEaSumInertial(nn, :);
        end
        PEa(nn, :)= ...
            ((6.242e18)*mion).*Rank(root).Specie(s).FluxTube(f).PEaSum(nn, :);
    end
        
        %         if Densp(nn, Qind) == 0
        %             if GRAVITYflag == 1
        %                 Rank(root).Specie(s).FluxTube(f).PEg(nn, Qind)= NaN;
        %             end
        %             if THERMALflag == 1
        %                 Rank(root).Specie(s).FluxTube(f).PEt(nn, Qind)= NaN;
        %             end
        %             if MIRRORflag == 1
        %                 Rank(root).Specie(s).FluxTube(f).PEm(nn, Qind)= NaN;
        %             end
        %             if EPARALLELflag == 1
        %                 Rank(root).Specie(s).FluxTube(f).PEp(nn, Qind)= NaN;
        %             end
        %             if AMBIPOLARflag == 1
        %                 Rank(root).Specie(s).FluxTube(f).PEaPressure(nn, Qind)= NaN;
        %                 if PEaInertialflag == 1
        %                     Rank(root).Specie(s).FluxTube(f).PEaInertial(nn, Qind)= NaN;
        %                 end
        %                 Rank(root).Specie(s).FluxTube(f).PEa(nn, Qind)= NaN;
        %             end
        %         end
end
%%
clc
close all
Nend= NNtT+ 1;

figure(1)
nn= Nend;
plot(abs(PEg(nn, 1:Qindlen)), rrAB2(1, :), 'k', 'MarkerSize', 10, 'LineWidth', 2);
hold on
plot(abs(PEt(nn, 1:Qindlen)), rrAB2(1, :), 'r.-', 'MarkerSize', 10, 'LineWidth', 2);
plot(abs(PEt(nn, 1:Qindlen))+ abs(PEa(nn, 1:Qindlen)), rrAB2(1, :), 'g.-', 'MarkerSize', 10, 'LineWidth', 2);
plot(abs(PEa(nn, 1:Qindlen)), rrAB2(1, :), 'b.-', 'MarkerSize', 10, 'LineWidth', 2);
plot(abs(PEp(nn, 1:Qindlen)), rrAB2(1, :), 'b*', 'MarkerSize', 10, 'LineWidth', 2);

% ylim([500 1200])
% plot(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(1, 1:Qindlen)), rrAB2, 'bx-', 'MarkerSize', 10, 'LineWidth', 2);
% hold on
% plot(abs(Rank(r).Specie(s).FluxTube(f).EAMagRTr(4, 1:Qindlen)), rrAB2, 'r.-', 'MarkerSize', 10, 'LineWidth', 2);

% plot(abs(PEt(1, 1:Qindlen)), rrAB2, 'b.-', 'MarkerSize', 10, 'LineWidth', 2);
% hold on
% plot(abs(PEt(end, 1:Qindlen)), rrAB2, 'r.-', 'MarkerSize', 10, 'LineWidth', 2);

%%
r= root; s= 1; f= 1;
fsize= 20; textx= -4.5; texty= 7; textx1= 0; texty1= 0.5; PEglim1= 0; PEglim2= 4;

legloc= 'SouthEast'; nn1= 1; nn2= round(NNtT/4); nn3= round(NNtT/2); nn4= NNtT;
xfigsize= 750; yfigsize= 350;

char1= 'ko-'; char2= 'bx-'; char3= 'rs-'; char4= 'g-'; char5= 'm-'; char6= 'ms-'; char7= 'r.-';

for Qind= 1:1:Qindlen;
    if AMBIPOLARflag == 1
        PEaTAv(Qind)= sum(Rank(root).Specie(s).FluxTube(f).PEa(2:NNtT+ 1, Qind))/(NNtT);
    end
    if THERMALflag == 1
        PEtTAv(Qind)= sum(Rank(root).Specie(s).FluxTube(f).PEt(2:NNtT+ 1, Qind))/(NNtT);
    end
    %     if (AMBIPOLARflag == 1) & (THERMALflag == 1)
    %         PEatTAv(Qind)= PEaTAv(Qind)+ PEtTAv(Qind);
    %         if isnan(PEatTAv(Qind)) == 1
    %             Rank(root).Specie(s).FluxTube(f).PEt(:, Qind)= NaN;
    %             Rank(root).Specie(s).FluxTube(f).PEt(:, Qind)= NaN;
    %             Rank(root).Specie(s).FluxTube(f).PEaPressure(:, Qind)= NaN;
    %             Rank(root).Specie(s).FluxTube(f).PEaInertial(:, Qind)= NaN;
    %             Rank(root).Specie(s).FluxTube(f).PEa(:, Qind)= NaN;
    %             if GRAVITYflag == 1
    %                 Rank(root).Specie(s).FluxTube(f).PEg(:, Qind)= NaN;
    %             end
    %         end
    %     end
end

clear L LH
fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 xfigsize+ 500 yfigsize+ 450])

subplot(2, 2, 1)
if GRAVITYflag == 1
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEg(nn1, 1:Qindlen)), rrAB2, char1, 'MarkerSize', 10, 'LineWidth', 2);
end
hold on
if (elecflag == 1)
    if THERMALflag == 1
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn1, 1:Qindlen)), rrAB2, char3, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn1, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn1, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn1, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn1, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn1, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
            if PEaInertialflag == 1
                plot(abs(Rank(root).Specie(s).FluxTube(f).PEaInertial(nn1, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
            end
        end
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEa(nn1, 1:Qindlen)), rrAB2, char4, 'MarkerSize', 10, 'LineWidth', 2);
        if (PEaInertialflag == 0)
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn1, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn1, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        else
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn1, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEa(nn1, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
else
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn1, 1:Qindlen)), rrAB2, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
end
title(horzcat('$t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn1)/60), ' [min]'), 'interpreter', 'latex', 'FontSize', 25)
xlabel('[eV]', 'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
if GRAVITYflag == 1
    LH(1)= plot(nan, nan, char(char1), 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\phi_g$');
end
if (elecflag == 1)
    if THERMALflag == 1
        LH(3)= plot(nan, nan, char(char3), 'MarkerSize', 10, 'LineWidth', 2);
        L{3}= horzcat('$\phi_t$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_p$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
        LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
        L{6}= horzcat('$\phi_p$');
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
            L{5}= horzcat('$\phi_a^p$');
            if PEaInertialflag == 1
                LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
                L{6}= horzcat('$\phi_a^i$');
            end
        end
        LH(4)= plot(nan, nan, char(char4), 'MarkerSize', 10, 'LineWidth', 2);
        L{4}= horzcat('$\phi_a$');
    end
    LH(2)= plot(nan, nan, char(char2), 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t+ \phi_a$');
else
    LH(2)= plot(nan, nan, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t$');
end
legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
xlim([xlim1 xlim2;]);
set(gca, 'FontSize', 25);
hold off

subplot(2, 2, 2)
if GRAVITYflag == 1
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEg(nn2, 1:Qindlen)), rrAB2, char1, 'MarkerSize', 10, 'LineWidth', 2);
end
hold on
if (elecflag == 1)
    if THERMALflag == 1
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn2, 1:Qindlen)), rrAB2, char3, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn2, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn2, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn2, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn2, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn2, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
            if PEaInertialflag == 1
                plot(abs(Rank(root).Specie(s).FluxTube(f).PEaInertial(nn2, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
            end
        end
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEa(nn2, 1:Qindlen)), rrAB2, char4, 'MarkerSize', 10, 'LineWidth', 2);
        if (PEaInertialflag == 0)
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn2, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn2, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        else
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn2, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEa(nn2, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
else
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn2, 1:Qindlen)), rrAB2, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
end
title(horzcat('$t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn2)/60), ' [min]'), 'interpreter', 'latex', 'FontSize', 25)
xlabel('[eV]', 'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
if GRAVITYflag == 1
    LH(1)= plot(nan, nan, char(char1), 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\phi_g$');
end
if (elecflag == 1)
    if THERMALflag == 1
        LH(3)= plot(nan, nan, char(char3), 'MarkerSize', 10, 'LineWidth', 2);
        L{3}= horzcat('$\phi_t$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_p$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
        LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
        L{6}= horzcat('$\phi_p$');
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
            L{5}= horzcat('$\phi_a^p$');
            if PEaInertialflag == 1
                LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
                L{6}= horzcat('$\phi_a^i$');
            end
        end
        LH(4)= plot(nan, nan, char(char4), 'MarkerSize', 10, 'LineWidth', 2);
        L{4}= horzcat('$\phi_a$');
    end
    LH(2)= plot(nan, nan, char(char2), 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t+ \phi_a$');
else
    LH(2)= plot(nan, nan, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t$');
end
legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
xlim([xlim1 xlim2;]);
set(gca, 'FontSize', 25);
hold off

subplot(2, 2, 3)
if GRAVITYflag == 1
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEg(nn3, 1:Qindlen)), rrAB2, char1, 'MarkerSize', 10, 'LineWidth', 2);
end
hold on
if (elecflag == 1)
    if THERMALflag == 1
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn3, 1:Qindlen)), rrAB2, char3, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn3, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn3, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn3, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn3, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn3, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
            if PEaInertialflag == 1
                plot(abs(Rank(root).Specie(s).FluxTube(f).PEaInertial(nn3, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
            end
        end
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEa(nn3, 1:Qindlen)), rrAB2, char4, 'MarkerSize', 10, 'LineWidth', 2);
        if (PEaInertialflag == 0)
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn3, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn3, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        else
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn3, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEa(nn3, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
else
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn3, 1:Qindlen)), rrAB2, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
end
title(horzcat('$t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn3)/60), ' [min]'), 'interpreter', 'latex', 'FontSize', 25)
xlabel('[eV]', 'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
if GRAVITYflag == 1
    LH(1)= plot(nan, nan, char(char1), 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\phi_g$');
end
if (elecflag == 1)
    if THERMALflag == 1
        LH(3)= plot(nan, nan, char(char3), 'MarkerSize', 10, 'LineWidth', 2);
        L{3}= horzcat('$\phi_t$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_p$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
        LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
        L{6}= horzcat('$\phi_p$');
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
            L{5}= horzcat('$\phi_a^p$');
            if PEaInertialflag == 1
                LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
                L{6}= horzcat('$\phi_a^i$');
            end
        end
        LH(4)= plot(nan, nan, char(char4), 'MarkerSize', 10, 'LineWidth', 2);
        L{4}= horzcat('$\phi_a$');
    end
    LH(2)= plot(nan, nan, char(char2), 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t+ \phi_a$');
else
    LH(2)= plot(nan, nan, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t$');
end
legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
xlim([xlim1 xlim2;]);
set(gca, 'FontSize', 25);
hold off

subplot(2, 2, 4)
if GRAVITYflag == 1
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEg(nn4, 1:Qindlen)), rrAB2, char1, 'MarkerSize', 10, 'LineWidth', 2);
end
hold on
if (elecflag == 1)
    if THERMALflag == 1
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn4, 1:Qindlen)), rrAB2, char3, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn4, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn4, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEm(nn4, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEp(nn4, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn4, 1:Qindlen)), rrAB2, char5, 'MarkerSize', 10, 'LineWidth', 2);
            if PEaInertialflag == 1
                plot(abs(Rank(root).Specie(s).FluxTube(f).PEaInertial(nn4, 1:Qindlen)), rrAB2, char6, 'MarkerSize', 10, 'LineWidth', 2);
            end
        end
        plot(abs(Rank(root).Specie(s).FluxTube(f).PEa(nn4, 1:Qindlen)), rrAB2, char4, 'MarkerSize', 10, 'LineWidth', 2);
        if (PEaInertialflag == 0)
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn4, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEaPressure(nn4, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        else
            plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn4, 1:Qindlen))+ ...
                abs(Rank(root).Specie(s).FluxTube(f).PEa(nn4, 1:Qindlen)), rrAB2, char2, 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
else
    plot(abs(Rank(root).Specie(s).FluxTube(f).PEt(nn4, 1:Qindlen)), rrAB2, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
end
title(horzcat('$t= $', num2str(Rank(1).Specie(1).FluxTube(1).TimeT(nn4)/60), ' [min]'), 'interpreter', 'latex', 'FontSize', 25)
xlabel('[eV]', 'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
if GRAVITYflag == 1
    LH(1)= plot(nan, nan, char(char1), 'MarkerSize', 10, 'LineWidth', 2);
    L{1}= horzcat('$\phi_g$');
end
if (elecflag == 1)
    if THERMALflag == 1
        LH(3)= plot(nan, nan, char(char3), 'MarkerSize', 10, 'LineWidth', 2);
        L{3}= horzcat('$\phi_t$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 0)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
    end
    if (MIRRORflag == 0) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_p$');
    end
    if (MIRRORflag == 1) & (EPARALLELflag == 1)
        LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
        L{5}= horzcat('$\phi_m$');
        LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
        L{6}= horzcat('$\phi_p$');
    end
    if AMBIPOLARflag == 1
        if MIRRORflag == 0 & EPARALLELflag == 0
            LH(5)= plot(nan, nan, char(char5), 'MarkerSize', 10, 'LineWidth', 2);
            L{5}= horzcat('$\phi_a^p$');
            if PEaInertialflag == 1
                LH(6)= plot(nan, nan, char(char6), 'MarkerSize', 10, 'LineWidth', 2);
                L{6}= horzcat('$\phi_a^i$');
            end
        end
        LH(4)= plot(nan, nan, char(char4), 'MarkerSize', 10, 'LineWidth', 2);
        L{4}= horzcat('$\phi_a$');
    end
    LH(2)= plot(nan, nan, char(char2), 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t+ \phi_a$');
else
    LH(2)= plot(nan, nan, 'b-', 'MarkerSize', 10, 'LineWidth', 2);
    L{2}= horzcat('$\phi_t$');
end
legend(LH, L, 'interpreter', 'latex', 'FontSize', 15, 'location', legloc)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
xlim([xlim1 xlim2;]);
set(gca, 'FontSize', 25);
hold off

saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION POTENTIAL WELL PLOTS.png')]);

toc

disp('ION POTENTIAL WELL PLOTS COMPLETE')

%% ----------------------------------------------------

clc
close all

tic

% PLOT ION POTENTIAL WELL CONTOURS:

fignum= fignum+ 1; % Assign figure number
fig(fignum)= figure(fignum);
set(fig(fignum), 'Position', [10 10 1050 1500])

clim1= 0e0;
clim2= 1e0;

if MIRRORflag == 0
    PEm(:, :)= 0e0;
end
if EPARALLELflag == 0
    PEp(:, :)= 0e0;
end

PEup(:, :)= PEt(:, :)+ ...
    PEm(:, :)+ PEa(:, :);

PEdown(:, :)= PEg(:, :)+ PEp(:, :);

sub1= subplot(2, 3, 1);
pp= pcolor(ndAB2/60, rrAB2, abs(PEg));
set(pp, 'EdgeColor', 'none');
shading interp
colormap(sub1, jet);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
% text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\phi_g$ [eV]'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
if MIRRORflag == 0
%     caxis([clim1 clim2]);
else
%     caxis([0 1]);
end
hold off

sub1= subplot(2, 3, 2);
pp= pcolor(ndAB2/60, rrAB2, PEt);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(sub1, jet);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
% text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\phi_t$ [eV]'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
if MIRRORflag == 0
%     caxis([-clim2 clim1]);
else
%     caxis([-1 0]);
end
hold off

sub1= subplot(2, 3, 3);
pp= pcolor(ndAB2/60, rrAB2, PEa);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(sub1, jet);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
% text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\phi_a$ [eV]'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
if MIRRORflag == 0
%     caxis([-clim2 clim1]);
else
%     caxis([-1 0]);
end
hold off

sub1= subplot(2, 3, 4);
pp= pcolor(ndAB2/60, rrAB2, PEp);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(sub1, jet);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
% text(textx, texty, '(b)', 'interpreter', 'latex', 'FontSize', 25)
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\phi_p$ [eV]'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
if MIRRORflag == 0
%     caxis([clim1 clim2]);
else
%     caxis([1 50]);
    %     caxis([1 10]);
end
hold off

sub1= subplot(2, 3, 5);
pp= pcolor(ndAB2/60, rrAB2, PEdown);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(sub1, jet);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
% text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\phi_{\downarrow}$ [eV]'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
if MIRRORflag == 0
%     caxis([-clim2 clim1]);
else
%     caxis([-1 0]);
end
hold off

sub1= subplot(2, 3, 6);
pp= pcolor(ndAB2/60, rrAB2, PEup);
set(pp, 'EdgeColor', 'none');
shading interp
colormap(sub1, jet);
cbar= colorbar('SouthOutside');
set(gca, 'FontSize', fsize);
% text(textx, texty, '(a)', 'interpreter', 'latex', 'FontSize', 25)
xlabel('$t$ [min]', 'interpreter', 'latex', 'FontSize', 25)
xlabel(cbar, horzcat('$\phi_{\uparrow}$ [eV]'), ...
    'interpreter', 'latex', 'FontSize', 25)
ylabel('$r$ [km]', 'interpreter', 'latex', 'FontSize', 25)
set(gca, 'xaxisLocation', 'top')
set(gca, 'yminortick', 'on','xminortick', 'on');
set(gcf, 'color','white');
ylim([ylim1 ylim2]);
if MIRRORflag == 0
%     caxis([clim1 clim2]);
else
    caxis([0 0.1]);
    %     caxis([0 10]);
end
hold off

saveas(figure(fignum), [horzcat(dataexpdir, 'Sim ', num2str(simnum), ' ION POTENTIAL WELL CONTOURS.png')]);

toc

disp('ION POTENTIAL WELL CONTOURS COMPLETE')
