module KineticMainParams

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	MAIN PARAMETERS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use ModelInput

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

include 'mpif.h' ! Include file for MPI commands

! ----------------------------------------------------

! Define date and time variables
integer(4), dimension(8) :: dateint
character(10), dimension(3)  :: datechar
character(50) :: monthstring, daystring, yearstring, hourstring, minutestring, secondstring

! Define all MPI variables
integer(4) status(MPI_STATUS_SIZE)
integer(4) :: len, ierr, rank
character(10) :: rankstring

! Define all running indices
integer(kind= dp) :: s, f, Qind, FAindIC, Vperp1ind, Vperp2ind, Vperpind, Vparind, Vpind, Vqind, Vphiind, &
	j, n, nn, jj, jnn, rr, MAfilterQind1, MAfilterQind2
character(50) :: Nsstring, sstring, fstring, Qindstring, Vperp1indstring, Vperp2indstring, Vperpindstring, Vparindstring, &
	Vpindstring, Vqindstring, Vphiindstring, nnstring, jstring, expstring, &
	paramstring, pdriftmaxstring, pdriftmeanstring, Timestring

! Define all simulation real-time values and strings
real(kind= dp) :: TotEnd, S1End, S2End, S3End, S4End, S5End, S6End, S7End, S811End, &
	S812End, S813End, S81End, S83End, S84End, S861End, S862End, S863End, &
 S864End, S865End, S866End, S87End, S88End, S8End, SDE2End, SDE3End, KS0End

character(50) :: Totstring, S1string, S2string, S3string, S4string, S5string, &
	S6string, S7string, S811string, S812string, S813string, S81string, S82string, &
	S83string, S84string, S861string, S862string, S863string, S864string, S865string, &
	S866string, S87string, S88string, S8string, SDE2string, SDE3string, KS0string

! ----------------- 1 SIMULATION PARAMETERIZATION -----------------

! Define number of particle species and flux tubes
integer(kind= dp), parameter :: Stot= 1d0 ! Number of particle species
integer(kind= dp), dimension(1) :: Nf ! Number of flux tubes per particle species
integer(kind= dp), dimension(:), allocatable :: Nfp

! Define imported ion and ENA config-space grid dimensions from MATLAB
integer(kind= dp), dimension(:, :), allocatable :: NqGp
integer(kind= dp), dimension(1) :: NqG ! Number of config-space grid cells
real(kind= dp), dimension(:, :), allocatable :: HiCratioMeanp, Tep
real(kind= dp), dimension(1) :: HiCratioMean, Te

! Define imported ion phase-space grid dimensions from MATLAB
integer(kind= dp), dimension(:, :, :), allocatable :: NVperp1Gp, NVperp2Gp, NVperpGp, NVparGp
integer(kind= dp), dimension(1) :: NVperp1G, NVperp2G, NVperpG, NVparG ! Number of ion vel-space grid cells
real(kind= dp) :: Vperp12NlinRange, VparNlinRange, Vperp12Gridlinspace, Vperp12sigma, Vperp12sigmaFac, &
	VparGridlinspace, Vparsigma, VparsigmaFac, VpphiNlinRange, VqNlinRange, &
	VpphiGridlinspace, Vpphisigma, VpphisigmaFac, VqGridlinspace, Vqsigma, VqsigmaFac

! Define imported reference density parameters
integer(kind= dp) :: nsnormCLBInput, nsnormCUBInput
integer(kind= dp) :: LBREPLENISHflagInput, UBREPLENISHflagInput
real(kind= dp), dimension(:), allocatable :: DensityInput, TemperatureInput, &
	EAInertialInput, EAPressureInput, EAmagInput

! Define imported ENA phase-space grid dimensions from MATLAB
integer(kind= dp), dimension(:, :, :), allocatable :: NVpGp, NVqGp, NVphiGp
integer(kind= dp), dimension(1) :: NVpG, NVqG, NVphiG ! Number of ENA vel-space grid cells

! Define all simulation time parameters
real(kind= dp), dimension(1) :: h
integer(kind= dp), dimension(1) :: NNt, ndatfac, Q0NNt, Q0ndatfac
real(kind= dp), dimension(1) :: IonNoiseLimit, ENANoiseLimit
integer(kind= dp), parameter :: dt= 1d0

! Define imported particle masses and charges from MATLAB
real(kind= dp), dimension(1) :: ms, qs
integer(kind= dp), dimension(1) :: Qindns0
real(kind= dp), dimension(:), allocatable :: msp, qsp
integer(kind= dp), dimension(:), allocatable :: Qindns0p

! Define imported initial temperatures and configuration-space grid parameters from MATLAB
real(kind= dp), dimension(:, :, :), allocatable :: hq, dq, qGCp, hqCp, dpCp, dqCp, dphiCp, &
	TsPerpp, TsParp, rGCp, phiGCp, thetaGCp, ellGCp, qGLp, qGHp, pGCp, d3xCp
real(kind= dp), dimension(1) :: qGC, hqC, dpC, dqC, dphiC, TsPerp, TsPar, &
	rGC, phiGC, thetaGC, ellGC, qGL, qGH, pGC, d3xC

! Define imported ion velocity-space grid parameters from MATLAB
real(kind= dp), dimension(:), allocatable :: VperpGLpp, VperpGHpp, VperpGCpp, dVperpGpp, &
	VparGLpp, VparGHpp, VparGCpp, dVparGpp, d3vCpp
real(kind= dp), dimension(:, :), allocatable :: VperpGLp, VperpGHp, VperpGCp, dVperpGp, &
	VparGLp, VparGHp, VparGCp, dVparGp, d3vCp
real(kind= dp), dimension(:), allocatable :: Vperp1GLpp, Vperp1GHpp, Vperp1GCpp, dVperp1Gpp, &
	Vperp2GLpp, Vperp2GHpp, Vperp2GCpp, dVperp2Gpp, VparGL2pp, VparGH2pp, VparGC2pp, dVparG2pp, d3vC2pp
real(kind= dp), dimension(:, :, :), allocatable :: Vperp1GLp, Vperp1GHp, Vperp1GCp, dVperp1Gp, &
	Vperp2GLp, Vperp2GHp, Vperp2GCp, dVperp2Gp, VparGL2p, VparGH2p, VparGC2p, dVparG2p, d3vC2p
real(kind= dp), dimension(1) :: dVperpG, VperpGL, VperpGH, VperpGC, Vperp1GL, Vperp1GH, Vperp1GC, dVperp1G, &
	Vperp2GL, Vperp2GH, Vperp2GC, dVperp2G, VparGL, VparGH, VparGC, dVparG, d3vC

! Define imported ENA velocity-space grid parameters from MATLAB
real(kind= dp), dimension(:), allocatable :: VpGLpp, VpGHpp, VpGCpp, VqGLpp, VqGHpp, VqGCpp, &
	VphiGLpp, VphiGHpp, VphiGCpp, hVpCpp, hVqCpp, hVphiCpp, dVpCpp, dVqCpp, dVphiCpp, d33vCpp
real(kind= dp), dimension(:, :, :), allocatable :: VpGLp, VpGHp, VpGCp, VqGLp, VqGHp, VqGCp, &
	VphiGLp, VphiGHp, VphiGCp, hVpCp, hVqCp, hVphiCp, dVpCp, dVqCp, dVphiCp, d33vCp
real(kind= dp), dimension(1) :: VpGL, VpGH, VpGC, VqGL, VqGH, VqGC, VphiGL, VphiGH, VphiGC, &
	hVpC, hVqC, hVphiC, dVpC, dVqC, dVphiC, d33vC

! Define range of configuration-space field-aligned grid cells
integer(kind= dp), dimension(1) :: NqIC, NqLB, NqUB

! Define drift limits over each time-step and on entire simulation duration
real(kind= dp), dimension(1) :: pdriftLim, rdriftLim, thetadriftLim, phidriftLim

! Define reference parallel potential drop parameters
real(kind= dp), dimension(1) :: EPar0p

! Define ICR heating parameters
real(kind= dp), dimension(1) :: OmegaG0p

integer(kind= dp), dimension(1) :: ranksize
integer(4) :: ranksizep

! ----------------- 2 DENSITY PROFILE A -----------------

! Define ion density normalization constant and ion scale-height density profile variables
real(kind= dp), dimension(:), allocatable :: ICbbp
real(kind= dp), dimension(1) :: ICbb, IC0bb, gC, argC, nsC, nsnormC
real(kind= dp), dimension(1) :: gCNeut, argCNeut, nsCNeut, nsnormCNeut
integer(kind= dp), dimension(1) :: NsFARR, NsFARRpTSum
integer(kind= dp), dimension(:), allocatable :: NsFARRp, NsFARRpTMP
integer(kind= dp), dimension(1) :: NsRR
real(kind= dp), dimension(1) :: Ns, NsR
integer(kind= dp), dimension(1) :: NsRRTot, nsnormCLB, nsnormCUB
integer(kind= dp), dimension(:), allocatable :: nsnormCNeut0
integer(kind= dp), dimension(:), allocatable :: NsRRTotpp
integer(kind= dp), dimension(1) :: NsTqind, NsTFAind
real(kind= dp), dimension(1) :: NsFAp, NsFARp
logical(kind= dp), dimension(:), allocatable :: NsFARRmsk

! ----------------- 3_1 GAUSSIAN RANDOM NUMBER GENERATOR -----------------

real(kind= dp), dimension(1) :: UniformRN1, UniformRN2, GaussianRN

! ----------------- 3_2 POISSON RANDOM NUMBER GENERATOR -----------------

real(kind= dp), dimension(1) :: MeanPoisson, ExpPoisson, UniformProd, RandUni
integer(4), dimension(1) :: RandPoisson

! ----------------- 5 DIPOLE POLYNOMIAL SOLVER -----------------

! Define all quartic dipole polynomial solver variables
real(kind= dp), dimension(1) :: ApIC, BpIC, CpIC, DpIC, &
	AbIC, BbIC, CbIC, D3IC, D4IC, AbbIC, BbbIC, DeltaIC, EpsilonIC, sigmaIC, &
	muIC, piiIC, r3IC, theta3IC, qtest3IC, ptest3IC, rfinalIC, thetafinalIC, &
	phifinalIC, xfinalIC, yfinalIC, zfinalIC, qfinalIC, pfinalIC, ellfinalIC
complex(kind= dp), dimension(1) :: ThetapIC, nuIC, gamma3IC, sigma23IC

! ----------------- 6 VELOCITY DISTRIBUTION -----------------

! Define all initial velocity distribution variables
real(kind= dp), dimension(1) :: Vxrandn1, Vxrandn2, Vyrandn1, Vyrandn2, Vzrandn1, Vzrandn2
real(kind= dp), dimension(1) :: thetaMB, phiMB, ellMB, TsMB

! ----------------- 6_1 MB VELOCITY DISTRIBUTION -----------------

integer(kind= dp), dimension(1) :: ICflag, LBflag
real(kind= dp), dimension(1) :: Vx0pp, Vy0pp, Vz0pp, &
	sigmaVx0, sigmaVy0, sigmaVz0, Vx0p, Vy0p, Vz0p, Vp0p, Vq0p, Vphi0p, &
	Vperp1MB, Vperp2MB, VperpMB, VparMBp, VparMB, VxMB, VyMB, VzMB
real(kind= dp), parameter :: V0BMoffset= 1d-50

! ----------------- 8 3D KINETIC SOLVER -----------------

! Define all kinetic solver variables
real(kind= dp), parameter :: ICRBMoffset= 1d-50
real(kind= dp), parameter :: EXPParticleSize= 7d0 ! Number of test particles for exporting
integer(kind= dp), dimension(1) :: dexpparticle, nnFlag
integer(kind= dp), dimension(1) :: NsTK, dNsTK1Ion, dNsTK1ENA, dNsTK1, dNsTK2, dNsTK3, dNsTK
real(kind= dp), dimension(:), allocatable :: NsTKp, NsTKRp
real(kind= dp), dimension(1) :: Time, TimeN
real(kind= dp), dimension(:), allocatable :: AEAmag, AGmag, AEPmag, AEAmagN, AGmagN, AEPmagN, x, y, z, &
	Vperp1, Vperp2, Vperp, Vpar, Vx, Vy, Vz
real(kind= dp), dimension(:), allocatable :: AEAmagTMP, AGmagTMP, AEPmagTMP, AEAmagNTMP, AGmagNTMP, AEPmagNTMP, &
	xTMP, yTMP, zTMP, Vperp1TMP, Vperp2TMP, VperpTMP, VparTMP, VxTMP, VyTMP, VzTMP
logical(kind= dp), dimension(:), allocatable :: LBoutfluxIonMsk, UBoutfluxIonMsk, &
	LBreplenishIonMsk, UBreplenishIonMsk, LBoutfluxENAMsk, UBoutfluxENAMsk
logical(kind= dp), dimension(:), allocatable :: LBoutfluxIonMskTMP, UBoutfluxIonMskTMP, &
 	LBreplenishIonMskTMP, UBreplenishIonMskTMP, LBoutfluxENAMskTMP, UBoutfluxENAMskTMP
real(kind= dp), dimension(:), allocatable :: BCreal
integer(kind= dp), dimension(:), allocatable :: BCinteger
logical(kind= dp), dimension(:), allocatable :: BClogical

! ----------------- 8_1 KINETIC RK4 UPDATE -----------------

real(kind= dp), dimension(:), allocatable :: xN, yN, zN, VxN, VyN, VzN, &
	xNTMP, yNTMP, zNTMP, VxNTMP, VyNTMP, VzNTMP
real(kind= dp), dimension(1) :: thetak1, phik1, ellk1
real(kind= dp), dimension(1) :: k1Vx, k1Vy, k1Vz, xk2, yk2, zk2, &
	k2Vx, k2Vy, k2Vz, k3Vx, k3Vy, k3Vz, xk4, yk4, zk4, k4Vx, k4Vy, k4Vz, &
	k1x, k1y, k1z, Vxk2, Vyk2, Vzk2, k2x, k2y, k2z, k3x, k3y, &
	k3z, Vxk4, Vyk4, Vzk4, k4x, k4y, k4z, xNp, yNp, zNp, rNp, thetaNp, &
	rN, thetaN, phiN, ellN, VxNp, VyNp, VzNp, VparNp, VpNp, VphiNp

! ----------------- 8_1_1 KINETIC UPDATE A-----------------

real(kind= dp), dimension(:), allocatable :: Vperp1N, Vperp2N, VperpN, Vperp1NTMP, Vperp2NTMP, VperpNTMP
real(kind= dp), dimension(:), allocatable :: Vphi, Vp, VphiTMP, VpTMP
real(kind= dp), dimension(1) :: rk1, pk1, qk1, Rperpk1, &
	Bmagk1, OmegaGk1, dBdsk1, muk1, AMxk1p, AMyk1p, AMzk1p, AMxk1, AMyk1, AMzk1, AMpark1, AMpk1, AMphik1, &
	AGxk1p, AGyk1p, AGzk1p, AGxk1, AGyk1, AGzk1, AGpark1, AGpk1, AGphik1, &
	AEAxk1p, AEAyk1p, AEAzk1p, AEAxk1, AEAyk1, AEAzk1, AEApark1, AEApk1, AEAphik1, AEAmagSk1, AGmagSk1, &
	AEPxk1p, AEPyk1p, AEPzk1p, AEPxk1, AEPyk1, AEPzk1, AEPpark1, AEPpk1, AEPphik1, AEPmagSk1, &
	Axk1, Ayk1, Azk1, Vr, Vtheta, EpsilonPar, &
	EpsilonPerp, alpha, VxR, VyR, VzR, lambdaPerp, EtaLH, XiPerp1, XiPerp2, S0, OmegaG0, ChiPerp1, &
	ChiPerp2, GammaPerp1, GammaPerp2, SigmaPerp, DPerp1, DPerp2, DVPerp1, DVPerp2, &
	epsPerp, tauPerp
integer(kind= dp), dimension(:), allocatable :: Qindk1, Qindk0, Vparindk1, Vperpindk1, Vperp1indk1, Vperp2indk1, &
	Vparindk0, Vperpindk0, Vperp1indk0, Vperp2indk0, Vpindk1, Vqindk1, Vphiindk1
integer(kind= dp), dimension(:), allocatable :: Qindk1TMP, Qindk0TMP, Vparindk1TMP, Vperp1indk1TMP, Vperp2indk1TMP, &
	Vperpindk1TMP, Vparindk0TMP, Vperp1indk0TMP, Vperp2indk0TMP, Vperpindk0TMP, Vpindk1TMP, Vqindk1TMP, Vphiindk1TMP
real(kind= dp), dimension(1) :: U1GammaPerp1randn, U2GammaPerp1randn, &
	U1GammaPerp2randn, U2GammaPerp2randn

! ----------------- 8_1_2 KINETIC UPDATE B-----------------

real(kind= dp), dimension(1) :: rk2, thetak2, phik2, qk2, pk2, &
	Rperpk2, ellk2, Bmagk2, OmegaGk2, dBdsk2, muk2, &
	AMxk2p, AMyk2p, AMzk2p, AMxk2, AMyk2, AMzk2, AMpark2, AMpk2, AMphik2, &
	AGxk2p, AGyk2p, AGzk2p, AGxk2, AGyk2, AGzk2, AGpark2, AGpk2, AGphik2, &
	AEAxk2p, AEAyk2p, AEAzk2p, AEAxk2, AEAyk2, AEAzk2, AEApark2, AEApk2, AEAphik2, AEAmagSk2, AGmagSk2, &
	AEPxk2p, AEPyk2p, AEPzk2p, AEPxk2, AEPyk2, AEPzk2, AEPpark2, AEPpk2, AEPphik2, AEPmagSk2, &
	Axk2, Ayk2, Azk2

! ----------------- 8_1_3 KINETIC UPDATE C-----------------

real(kind= dp), dimension(:), allocatable :: pk4, rk4, thetak4, phik4, ellk4, qk4, &
	pk4TMP, rk4TMP, thetak4TMP, phik4TMP, ellk4TMP, qk4TMP
real(kind= dp), dimension(1) ::  Rperpk4, Bmagk4, OmegaGk4, dBdsk4, muk4, &
	AMxk4p, AMyk4p, AMzk4p, AMxk4, AMyk4, AMzk4, AMpark4, AMpk4, AMphik4, &
	AGxk4p, AGyk4p, AGzk4p, AGxk4, AGyk4, AGzk4, AGpark4, AGpk4, AGphik4, &
	AEAxk4p, AEAyk4p, AEAzk4p, AEAxk4, AEAyk4, AEAzk4, AEApark4, AEApk4, AEAphik4, AEAmagSk4, AGmagSk4, &
	AEPxk4p, AEPyk4p, AEPzk4p, AEPxk4, AEPyk4, AEPzk4, AEPpark4, AEPpk4, AEPphik4, AEPmagSk4, &
	Axk4, Ayk4, Azk4

! ----------------- 8_1_4 RK4 DIPOLE POLYNOMIAL SOLVER -----------------

! Define all RK4 quartic dipole polynomial solver variables
real(kind= dp), dimension(1) :: qNp, pNp, ApRK4, BpRK4, CpRK4, DpRK4, &
	AbRK4, BbRK4, CbRK4, D3RK4, D4RK4, AbbRK4, BbbRK4, DeltaRK4, EpsilonRK4, sigmaRK4, &
	muRK4, piiRK4, r3RK4, theta3RK4, qtest3RK4, ptest3RK4, rfinalRK4, thetafinalRK4, &
	phifinalRK4, xfinalRK4, yfinalRK4, zfinalRK4, qfinalRK4, pfinalRK4
complex(kind= dp), dimension(1) :: ThetapRK4, nuRK4, gamma3RK4, sigma23RK4

! ----------------- 8_1_1_1 ION-NEUTRAL COLLISIONS -----------------

integer(kind= dp), dimension(1) :: ENAcount
real(kind= dp), dimension(1) :: Vpp, Vphip, Vparp
logical(kind= dp), dimension(:), allocatable :: ENAflag, ENAflagTMP
integer(kind= dp), dimension(:), allocatable :: ENAflagN0ind, ENAflagN1ind, ENAflagN0indTMP, ENAflagN1indTMP

! ----------------- 8_3 PARTICLE COUNTS -----------------

! Define all ion and ENA config space statistical particle counts variables:

real(kind= dp), dimension(1) :: NqLBoutfluxIon, NqLBoutfluxIonR, NqUBoutfluxIon, NqUBoutfluxIonR, &
	NqLBoutfluxENA, NqLBoutfluxENAR, NqUBoutfluxENA, NqUBoutfluxENAR
real(kind= dp), dimension(1) :: Nq, NqENA, NqR, NqENAR
real(kind= dp) :: jcount

! Define all ion phase space statistical particle counts variables:

real(kind= dp), dimension(1) :: Nph, NphR, N2Perpph, N2perpphR

! Define all ion limits:

real(kind= dp), dimension(:), allocatable :: pdriftion, qdriftion, phidriftion, &
	rdriftENA, thetadriftENA, phidriftENA, pdriftionTMP, qdriftionTMP, phidriftionTMP, &
	rdriftENATMP, thetadriftENATMP, phidriftENATMP
real(kind= dp), dimension(1) :: pdriftionMax, pdriftionMean, qdriftionMax, phidriftionMax, pdriftionMaxR, &
 	pdriftionMeanR, qdriftionMaxR, phidriftionMaxR, rdriftENAMax, thetadriftENAMax, phidriftENAMax, &
	rdriftENAMaxR, thetadriftENAMaxR, phidriftENAMaxR
real(kind= dp), dimension(1) :: VperpMin, Vperp1Min, Vperp2Min, VperpMax, Vperp1Max, Vperp2Max, &
	VparMin, VparMax, VperpMinR, Vperp1MinR, Vperp2MinR, VperpMaxR, Vperp1MaxR, Vperp2MaxR, &
	VparMinR, VparMaxR

! Define all ENA phase space statistical particle counts variables:

real(kind= dp), dimension(1) :: NphENA, NphENAR

! Define all ion and ENA fluid reference variables:

real(kind= dp), dimension(:), allocatable :: Vperp1REF, Vperp2REF, VperpREF, VparREF, &
	Vperp1REFsig, Vperp2REFsig, VperpREFsig, VparREFsig
real(kind= dp), dimension(1) :: VperpREFp, Vperp1REFp, Vperp2REFp, VparREFp, &
	VperpREFsigp, Vperp1REFsigp, Vperp2REFsigp, VparREFsigp, &
	VperpREFR, Vperp1REFR, Vperp2REFR, VparREFR, &
	VperpREFsigR, Vperp1REFsigR, Vperp2REFsigR, VparREFsigR
real(kind= dp), dimension(:), allocatable :: VpENAREF, VqENAREF, VphiENAREF, &
	VpENAREFsig, VqENAREFsig, VphiENAREFsig
real(kind= dp), dimension(1) :: VpENAREFp, VqENAREFp, VphiENAREFp, &
 	VpENAREFsigp, VqENAREFsigp, VphiENAREFsigp, &
	VpENAREFR, VqENAREFR, VphiENAREFR, &
	VpENAREFsigR, VqENAREFsigR, VphiENAREFsigR

! ----------------- 8_5 2D VELOCITY SPACE INTEGRATOR -----------------

real(kind= dp), dimension(:, :, :), allocatable :: ggg, Sum1
real(kind= dp), dimension(:, :, :, :), allocatable :: ggg2Perp, Sum12Perp
real(kind= dp), dimension(:), allocatable :: II, MM

! ----------------- 8_6_6_1 ION-NEUTRAL COLLISION FREQUENCY -----------------

real(kind= dp), dimension(1) :: sigmaIonNeut, VxNeutrandn1, VxNeutrandn2, &
	VyNeutrandn1, VyNeutrandn2, VzNeutrandn1, VzNeutrandn2, VxNeutp, VyNeutp, VzNeutp, &
	VxNeut, VyNeut, VzNeut, sigmaVxNeut, sigmaVyNeut, sigmaVzNeut, &
	VelIonNeut, MionVx, MionVy, MionVz, nuIonNeutM

! ----------------- 8_7 AMBIPOLAR ELECTRIC FIELD -----------------

real(kind= dp), dimension(1) :: EGravR, EAInertialR1, EAInertialR2, EAInertialR, &
	EAPressureR1, EAPressureR2, EAPressureR, &
	xLinInterp, xLinInterp1, xLinInterp2, yLinInterp, yLinInterp1, yLinInterp2, EAmagInterp, EGmagInterp, &
	M0FiltSum1, M0FiltSum2, M1ParFiltSum1, M1ParFiltSum2

! ----------------- 8_7_1 MOMENT FILTER -----------------

integer(kind= dp), dimension(1) :: MAfilterPt
real(kind= dp), dimension(:, :), allocatable :: MomentFiltIn, MomentFiltOut
real(kind= dp), dimension(1) :: FiltSum1, FiltSum2, FiltAvrg1, FiltAvrg2

! ----------------- 8_7_2 GRAVITATIONAL FIELD -----------------

real(kind= dp), dimension(:), allocatable :: ICbbpG
real(kind= dp), dimension(1) :: ICbbG, IC0bbG, gCG

! ----------------- 8_8 PARALLEL ELECTRIC FIELD -----------------

real(kind= dp), dimension(1) :: EParR, EPmagInterp

! ----------------- 8_9 3D VELOCITY SPACE INTEGRATOR -----------------

real(kind= dp), dimension(:, :, :, :), allocatable :: gggENA, Sum1ENA
real(kind= dp), dimension(:), allocatable :: IIENA, MMENA

! ----------------- 9 MPI REDUCE SUM -----------------

real(kind= dp), dimension(1) :: MPIRedSumIn, MPIRedSumOut

! ----------------- 0 DATA EXPORT -----------------

real(kind= dp), dimension(100000) :: RNtest1, RNtest2

integer(4) :: expint ! Data export file unit
! Define all exported data file names:
character(50) :: RNtest1file, RNtest2file, NsnTfile, NsnRRTfile, NqLBoutfluxIonRTfile, LBoutfluxIonRTfile, &
	NqUBoutfluxIonRTfile, UBoutfluxIonRTfile, NqLBoutfluxENARTfile, LBoutfluxENARTfile, &
	NqUBoutfluxENARTfile, UBoutfluxENARTfile, LBNetDensityTfile, UBNetDensityTfile, &
	nsnormfacTfile, TimeTfile, TeNTfile, N2PerpphRTfile, NphRTfile, NphENARTfile, M0phRTfile, M0phENARTfile, &
	M1PerpphRTfile, M1Perp1phRTfile, M1Perp2phRTfile, M1ParphRTfile, sigmaIonNeutRTfile, nuIonNeutRTfile, &
	M2phRTfile, M2PerpphRTfile, M2Perp1phRTfile, M2Perp2phRTfile, M2ParphRTfile, M1PphENARTfile, &
	M1QphENARTfile, M1PHIphENARTfile, M2phENARTfile, M2PphENARTfile, M2QphENARTfile, M2PHIphENARTfile, &
	PhiParRTfile, EAInertialRTfile, EAPressureRTfile, EAmagRTfile, EGmagRTfile, EPmagRTfile, &
	M0FiltAvrgRTfile, M1Perp1FiltAvrgRTfile, M1Perp2FiltAvrgRTfile, M1ParFiltAvrgRTfile, M2ParFiltAvrgRTfile, &
	nsnormCLBTfile, nsnormCUBTfile, DensityOutputRTfile, TemperatureOutputRTfile, &
	EAInertialOutputRTfile, EAPressureOutputRTfile, EAmagOutputRTfile, &
	LBREPLENISHflagTfile, UBREPLENISHflagTfile

! ----------------- ALL DERIVED DATA TYPES -----------------

type V3Celltype ! Vel-space grid cell derived data type: rank 3, indices [s, f, Qind, (Vpind, Vqind, Vphiind)]
	! Simulation Parameterization:
	real(kind= dp), dimension(1) :: VpGLT, VpGHT, VpGCT, VqGLT, VqGHT, VqGCT, &
		VphiGLT, VphiGHT, VphiGCT, hVpCT, hVqCT, hVphiCT, dVpCT, dVqCT, dVphiCT, d33vCT
	! Particle Counts:
	real(kind= dp), dimension(:), allocatable :: NphENART
	! ENA Distribution Functions:
	real(kind= dp), dimension(:), allocatable :: FphENART
	! ENA Moments:
	real(kind= dp), dimension(:), allocatable :: g0phENART

end type V3Celltype

type V2PerpCelltype ! Vel-space grid cell derived data type: rank 3, indices [s, f, Qind, (Vperp1ind, Vperp2ind, Vparind)]
	! Simulation Parameterization:
	real(kind= dp), dimension(1) :: Vperp1GLT, Vperp1GHT, Vperp1GCT, dVperp1GT, &
		Vperp2GLT, Vperp2GHT, Vperp2GCT, dVperp2GT, VparGLT, VparGHT, VparGCT, dVparGT, d3vCT
	! Particle Counts:
	real(kind= dp), dimension(:), allocatable :: N2PerpphRT
	! Ion Distribution Functions:
	real(kind= dp), dimension(:), allocatable :: F2PerpphRT
	! Ion Moments:
	real(kind= dp), dimension(:), allocatable :: g02PerpphRT

end type V2PerpCelltype

type VCelltype ! Vel-space grid cell derived data type: rank 2, indices [s, f, Qind, (Vperpind, Vparind)]
  ! Simulation Parameterization:
	real(kind= dp), dimension(1) :: dVperpGT, dVparGT, VperpGLT, VperpGHT, VperpGCT, &
		VparGLT, VparGHT, VparGCT, d3vCT
	! Particle Counts:
	real(kind= dp), dimension(:), allocatable :: NphRT
	! Ion Distribution Functions:
	real(kind= dp), dimension(:), allocatable :: FphRT
	! Ion Moments:
	real(kind= dp), dimension(:), allocatable :: g0phRT

end type VCelltype

type QCellICtype ! Config-space initialization grid cell derived data type: rank 1, indices [s, f, QindIC]
	! Density Profile A:
	integer(kind= dp), dimension(1) :: NsFARRT, NsFAT
	! Density Profile B:
	real(kind= dp), dimension(:), allocatable :: qniICT, pniICT
	! Dipole Polynomial Solver:
	real(kind= dp), dimension(:), allocatable :: rfinalICT, thetafinalICT, phifinalICT, &
		xfinalICT, yfinalICT, zfinalICT, qfinalICT, pfinalICT, ellfinalICT
	! Velocity Distribution:
	real(kind= dp), dimension(:), allocatable :: Vperp1ICT, Vperp2ICT, VperpICT, VparICT, VxICT, &
		VyICT, VzICT

end type QCellICtype

type QCell0type ! Config-space initialization grid cell derived data type: rank 1, indices [s, f, Qind0]
	! Simulation Parameterization:
	real(kind= dp), dimension(1) :: qGC0T, hqC0T, dpC0T, dqC0T, dphiC0T, Ts0T, TsPerp0T, TsPar0T, &
		rGC0T, phiGC0T, thetaGC0T, ellGC0T, qGL0T, qGH0T, pGC0T, d3xC0T
	! Density Profile A:
	real(kind= dp), dimension(1) :: nsnormCT
	real(kind= dp), dimension(1) :: nsnormCNeut0T

end type QCell0type

type QCelltype ! Config-space grid cell derived data type: rank 1, indices [s, f, Qind]
  ! Simulation Parameterization:
	real(kind= dp), dimension(1) :: DensityInputT, TemperatureInputT, &
		EAInertialInputT, EAPressureInputT, EAmagInputT
	integer(kind= dp), dimension(1) :: NVperp1GT, NVperp2GT, NVperpGT, NVparGT
	integer(kind= dp), dimension(1) :: NVpGT, NVqGT, NVphiGT
	real(kind= dp), dimension(1) :: qGCT, hqCT, dpCT, dqCT, dphiCT, TsT, TsPerpT, TsParT, rGCT, &
		phiGCT, thetaGCT, ellGCT, qGLT, qGHT, pGCT, d3xCT
	real(kind= dp), dimension(1) :: lambdaPerppT, EtaLHpT, XiPerp1pT, XiPerp2pT, S0pT, &
		OmegaG0pT, ChiPerp1pT, ChiPerp2pT
	! Density Profile A:
	real(kind= dp), dimension(1) :: nsnormCNeutT
	! Particle Counts:
	real(kind= dp), dimension(:), allocatable :: NqT, NqENAT, NqRT, NqENART
	real(kind= dp), dimension(:, :, :), allocatable :: NphReNormT, NphReNormRT
	real(kind= dp), dimension(:, :, :, :), allocatable :: N2PerpphReNormT, N2PerpphReNormRT
	real(kind= dp), dimension(:, :, :, :), allocatable :: NphReNormENAT, NphReNormENART
	real(kind= dp), dimension(:, :, :), allocatable :: NphRTp
	real(kind= dp), dimension(:, :, :, :), allocatable :: N2PerpphRTp
	real(kind= dp), dimension(:, :, :, :), allocatable :: NphENARTp
	! Ion Distribution Functions:
	real(kind= dp), dimension(:, :, :), allocatable :: FphRTp
	real(kind= dp), dimension(:, :, :, :), allocatable :: F2PerpphRTp
	! ENA Distribution Functions:
	real(kind= dp), dimension(:, :, :, :), allocatable :: FphENARTp

	type(VCelltype), dimension(:, :), allocatable :: VCellT
	type(V2PerpCelltype), dimension(:, :, :), allocatable :: V2PerpCellT
	type(V3Celltype), dimension(:, :, :), allocatable :: V3CellT
end type QCelltype

type FluxTubetype ! Flux tube number derived data type: rank 1, indices [s, f]
  ! Simulation Parameterization:
	integer(kind= dp), dimension(1) :: nsnormCLBInputT, nsnormCUBInputT
	real(kind= dp), dimension(:), allocatable :: ellqCT
	real(kind= dp), dimension(1) :: SUMellqCT
	integer(kind= dp), dimension(1) :: NqGT, NqG0T
	real(kind= dp), dimension(:), allocatable :: TeNT
	real(kind= dp), dimension(1) :: HiCratioMeanT, TeT, LBNominalDensityT, UBNominalDensityT, &
		d3xCLBT, d3xCUBT, sigmaLBT, sigmaUBT
	real(kind= dp), dimension(1) :: AT, BT, hT
	integer(kind= dp), dimension(1) :: NtT, NNtT, ndatfacT, Q0NNtT, Q0ndatfacT
	real(kind= dp), dimension(1) :: IonNoiseLimitT, ENANoiseLimitT
	integer(kind= dp), dimension(1) :: IONNOISEflagT, FLUIDIONEXPORTflagT, FLUIDENAEXPORTflagT, &
		LBCONDITIONflagT, UBCONDITIONflagT, LBREPLENISHflagT, UBREPLENISHflagT, &
		DENSITYPROFILEflagT, STATICINJECTIONflagT, DENSITYOUTPUTflagT, DENSITYINPUTflagT, &
		ICRflagT, MIRRORflagT, GRAVflagT, EAMBSELFCONSISTflagT, EAMBSIGNflagT, EAMBflagT, &
		EAINERTIALflagT, EAPRESSUREflagT, MOMENTFILTERflagT, EPARflagT, &
		ION2VPERPflagT, PHASEIONDISTRIBflagT, PHASEDENSITYIONMOMENTflagT, &
		PHASEVELPERPIONMOMENTflagT, PHASEVELPARIONMOMENTflagT, PHASEENERGYIONMOMENTflagT, &
		PHASEENERGYPERPIONMOMENTflagT, PHASEENERGYPARIONMOMENTflagT, PHASEENERGYPERPIONMOMENTCENTERflagT, &
		PHASEENERGYPARIONMOMENTCENTERflagT, FLUIDIONREFflagT, &
		ENANOISEflagT, QEXCHANGEflagT
	integer(kind= dp), dimension(1) :: NqICAT, NqICBT, NqICT
	real(kind= dp), dimension(1) :: zns0T, ns0T, nsnormfacT, zns0NeutT, ns0NeutT
	real(kind= dp), dimension(1) :: pdriftLimT, rdriftLimT, thetadriftLimT, &
		phidriftLimT
	real(kind= dp), dimension(:), allocatable :: rGCTp, hqCTp, dqCTp
	real(kind= dp), dimension(1) :: PhiPar0BT
	! Density Profile A:
	integer(kind= dp), dimension(1) :: nsnormCLBT, nsnormCUBT
	integer(kind= dp), dimension(:), allocatable :: NsFARRpT, NsFApT, NsFARpT
	integer(kind= dp), dimension(1) :: NsRRT, NsRT, NsT
	integer(kind= dp), dimension(:), allocatable :: NsnT, NsnRRT
	! Initial Conditions:
	real(kind= dp), dimension(:), allocatable :: x0T, y0T, z0T, r0T, theta0T, phi0T, q0T, &
		p0T, Vperp10T, Vperp20T, Vperp0T, Vpar0T, Vx0T, Vy0T, Vz0T
	! Kinetic Solver:
	integer(kind= dp), dimension(1) :: expparticle1T, expparticle2T, expparticle3T, &
		expparticle4T, expparticle5T, expparticle6T, expparticle7T
	real(kind= dp), dimension(:), allocatable :: TimeT
	! Particle Counts:
	real(kind= dp), dimension(:), allocatable :: NqReNormLBoutfluxIonT, NqReNormUBoutfluxIonT, &
		NqReNormLBreplenishIonT, NqReNormUBreplenishIonT, &
		NqReNormLBoutfluxENAT, NqReNormUBoutfluxENAT
	real(kind= dp), dimension(:), allocatable :: NqLBoutfluxIonRT, NqUBoutfluxIonRT, &
		LBoutfluxIonRT, UBoutfluxIonRT, NqLBoutfluxENART, NqUBoutfluxENART, &
		LBoutfluxENART, UBoutfluxENART, LBNetDensityT, UBNetDensityT
	real(kind= dp), dimension(:, :), allocatable :: NqReNormT, NqReNormENAT
	real(kind= dp), dimension(:, :), allocatable :: VperpREFpT, Vperp1REFpT, Vperp2REFpT, VparREFpT, &
		VperpREFsigpT, Vperp1REFsigpT, Vperp2REFsigpT, VparREFsigpT, &
		VperpREFRT, Vperp1REFRT, Vperp2REFRT, VparREFRT, &
		VperpREFsigRT, Vperp1REFsigRT, Vperp2REFsigRT, VparREFsigRT, &
		VpENAREFpT, VqENAREFpT, VphiENAREFpT, VpENAREFsigpT, VqENAREFsigpT, VphiENAREFsigpT, &
		VpENAREFRT, VqENAREFRT, VphiENAREFRT, VpENAREFsigRT, VqENAREFsigRT, VphiENAREFsigRT
	real(kind= dp), dimension(1) :: pdriftionMaxRT, pdriftionMeanRT, qdriftionMaxRT, phidriftionMaxRT, &
		rdriftENAMaxRT, thetadriftENAMaxRT, phidriftENAMaxRT
	real(kind= dp), dimension(1) :: VperpMinT, Vperp1MinT, Vperp2MinT, VperpMaxT, Vperp1MaxT, Vperp2MaxT, &
		VparMinT, VparMaxT, VperpMinRT, Vperp1MinRT, Vperp2MinRT, VperpMaxRT, Vperp1MaxRT, Vperp2MaxRT, VparMinRT, VparMaxRT
	real(kind= dp), dimension(:, :), allocatable :: NqTp, NqENATp, NqRTp, NqENARTp
	! Ion Moments:
	real(kind= dp), dimension(:, :), allocatable :: M0phRT, M1PerpphRT, M1Perp1phRT, M1Perp2phRT, &
		M1ParphRT, M2phRT, M2ParphRT, M2PerpphRT, M2Perp1phRT, M2Perp2phRT
	real(kind= dp), dimension(:), allocatable :: DensityOutputRT, TemperatureOutputRT, &
		EAInertialOutputRT, EAPressureOutputRT, EAmagOutputRT
	! Ion Neutral Collision Frequency:
	real(kind= dp), dimension(:, :), allocatable :: sigmaIonNeutRT, nuIonNeutRT, nuIonNeutRankSumRT, nuIonNeutPoiT
	! Ambipolar Electric Field:
	real(kind= dp), dimension(:, :), allocatable :: LambdaDRT, EAInertialRT, EAPressureRT, EAmagRT
	! Gravitational Field:
	real(kind= dp), dimension(:), allocatable :: EGmagRT
	! Moment Filter:
	real(kind= dp), dimension(:, :), allocatable :: MomentFiltInT, MomentFiltOutT, &
		M0FiltAvrgRT, M1Perp1FiltAvrgRT, M1Perp2FiltAvrgRT, M1ParFiltAvrgRT, M2ParFiltAvrgRT
	! Parallel Electric Field:
	real(kind= dp), dimension(:), allocatable :: PhiParRT
	real(kind= dp), dimension(:), allocatable :: EPmagRT
	! ENA Moments:
	real(kind= dp), dimension(:, :), allocatable :: M0phENART, M1PphENART, &
		M1QphENART, M1PHIphENART, M2phENART, M2PphENART, M2QphENART, M2PHIphENART

	type(QCellICtype), dimension(:), allocatable :: QCellICT
	type(QCelltype), dimension(:), allocatable :: QCellT
	type(QCell0type), dimension(:), allocatable :: QCell0T
end type FluxTubetype

type Specietype ! Particle species number derived data type: rank 1, index [s]
 	! Simulation Parameterization:
    integer(kind= dp), dimension(1) :: NfT
    real(kind= dp), dimension(1) :: msT, qsT
		integer(kind= dp), dimension(1) :: Qindns0T
	! Density Profile A:
	integer(kind= dp), dimension(1) :: NsRRTotpST
	integer(kind= dp), dimension(:), allocatable :: NsRRTotpT
	type(FluxTubetype), dimension(:), allocatable :: FluxTubeT
end type Specietype
type(Specietype), dimension(:), allocatable :: SpecieT

type Ranktype ! Rank number derived data type: rank 1, index [rr+ 1]
	! Density Profile A:
	integer(kind= dp), dimension(1) :: NsTqindT, NsTFAindT
end type Ranktype
type(Ranktype), dimension(:), allocatable :: RankT

! ----------------------------------------------------

contains

! ----------------------------------------------------

! DEFINE ALL PARTICLE TRANSLATIONAL ENERGIZATION TERMS:

! ----------------------------------------------------

! DEFINE ALL PARTICLE TRANSLATIONAL COORDINATES:

subroutine rsub(r, x, y, z)
	implicit none
	real(kind= dp), intent(out) :: r
	real(kind= dp), intent(in) :: x, y, z
		r= sqrt(x**2d0+ y**2d0+ z**2d0)
end subroutine rsub

subroutine thetasub(theta, z, r)
	implicit none
	real(kind= dp), intent(out) :: theta
	real(kind= dp), intent(in) :: z, r
		theta= acos(z/r)
end subroutine thetasub

subroutine phisub(phi, x, y)
	implicit none
	real(kind= dp), intent(out) :: phi
	real(kind= dp), intent(in) :: x, y
		phi= atan2(y, x)
end subroutine phisub

subroutine qsub(q, r, theta)
	implicit none
	real(kind= dp), intent(out) :: q
	real(kind= dp), intent(in) :: r, theta
		q= (RE**2d0)*cos(theta)/r**2d0
end subroutine qsub

subroutine psub(p, r, theta)
	implicit none
	real(kind= dp), intent(out) :: p
	real(kind= dp), intent(in) :: r, theta
		p= r/(RE*(sin(theta))**2d0)
end subroutine psub

subroutine Rperpsub(Rperp, p, phi)
	implicit none
	real(kind= dp), intent(out) :: Rperp
	real(kind= dp), intent(in) :: p, phi
		Rperp= sqrt(p**2d0+ phi**2d0)
end subroutine Rperpsub

subroutine ellsub(ell, theta)
	implicit none
	real(kind= dp), intent(out) :: ell
	real(kind= dp), intent(in) :: theta
		ell= 1d0+ 3d0*(cos(theta))**2d0
end subroutine ellsub

! ----------------------------------------------------

! DEFINE MIRROR FORCE PARTICLE ACCELERATION TERMS:

subroutine Bmagsub(Bmag, r, ell)
	implicit none
	real(kind= dp), intent(out) :: Bmag
	real(kind= dp), intent(in) :: r, ell
		Bmag= m*sqrt(ell)/r**3d0
end subroutine Bmagsub

subroutine OmegaGsub(OmegaG, qs, ms, Bmag)
	implicit none
	real(kind= dp), intent(out) :: OmegaG
	real(kind= dp), intent(in) :: qs, ms, Bmag
		OmegaG= qs*Bmag/ms ! [rads/s]
end subroutine OmegaGsub

subroutine dBdssub(dBds, r, theta, ell)
	implicit none
	real(kind= dp), intent(out) :: dBds
	real(kind= dp), intent(in) :: r, theta, ell
		dBds= (-6d0*m*cos(theta)/r**4d0)- (3d0*m*cos(theta)*(sin(theta)**2d0)/(ell*r**4d0))
end subroutine dBdssub

subroutine musub(mu, ms, Bmag, Vperp)
	implicit none
	real(kind= dp), intent(out) :: mu
	real(kind= dp), intent(in) :: ms, Bmag, Vperp
		mu= ms*(Vperp**2d0)/(2d0*Bmag)
end subroutine musub

subroutine AMxsub(AMx, ms, mu, dBds, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AMx
	real(kind= dp), intent(in) :: ms, mu, dBds, theta, phi, ell
		AMx= ((-mu*dBds)*3d0*cos(phi)*cos(theta)*sin(theta)/(sqrt(ell)))/ms
end subroutine AMxsub

subroutine AMysub(AMy, ms, mu, dBds, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AMy
	real(kind= dp), intent(in) :: ms, mu, dBds, theta, phi, ell
		AMy= ((-mu*dBds)*3d0*sin(phi)*cos(theta)*sin(theta)/(sqrt(ell)))/ms
end subroutine AMysub

subroutine AMzsub(AMz, ms, mu, dBds, theta, ell)
	implicit none
	real(kind= dp), intent(out) :: AMz
	real(kind= dp), intent(in) :: ms, mu, dBds, theta, ell
		AMz= ((-mu*dBds)*(3d0*(cos(theta))**2d0- 1d0)/(sqrt(ell)))/ms
end subroutine AMzsub

! ----------------------------------------------------

! DEFINE GRAVITATIONAL FORCE PARTICLE ACCELERATION TERMS:

subroutine AGxsub(AGx, AGmagS, AGmag, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AGx
	real(kind= dp), intent(in) :: AGmagS, AGmag, theta, phi, ell
		AGx= AGmagS*abs(AGmag)*(3d0*cos(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
end subroutine AGxsub

subroutine AGysub(AGy, AGmagS, AGmag, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AGy
	real(kind= dp), intent(in) :: AGmagS, AGmag, theta, phi, ell
		AGy= AGmagS*abs(AGmag)*(3d0*sin(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
end subroutine AGysub

subroutine AGzsub(AGz, AGmagS, AGmag, theta, ell)
	implicit none
	real(kind= dp), intent(out) :: AGz
	real(kind= dp), intent(in) :: AGmagS, AGmag, theta, ell
		AGz= AGmagS*abs(AGmag)*((3d0*(cos(theta))**2d0- 1d0)/(sqrt(ell)))
end subroutine AGzsub



subroutine AGENAxsub(AGx, mNeut, r, theta, phi)
	implicit none
	real(kind= dp), intent(out) :: AGx
	real(kind= dp), intent(in) :: mNeut, r, theta, phi
		AGx= ((-GG*ME*mNeut/(r**2d0))*sin(theta)*cos(phi))/mNeut
end subroutine AGENAxsub

subroutine AGENAysub(AGy, mNeut, r, theta, phi)
	implicit none
	real(kind= dp), intent(out) :: AGy
	real(kind= dp), intent(in) :: mNeut, r, theta, phi
		AGy= ((-GG*ME*mNeut/(r**2d0))*sin(theta)*sin(phi))/mNeut
end subroutine AGENAysub

subroutine AGENAzsub(AGz, mNeut, r, theta)
	implicit none
	real(kind= dp), intent(out) :: AGz
	real(kind= dp), intent(in) :: mNeut, r, theta
		AGz= ((-GG*ME*mNeut/(r**2d0))*cos(theta))/mNeut
end subroutine AGENAzsub

! ----------------------------------------------------

! DEFINE AMBIPOLAR ELECTRIC FIELD PARTICLE ACCELERATION TERMS:

subroutine AEAxsub(AEAx, AEAmagS, AEAmag, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AEAx
	real(kind= dp), intent(in) :: AEAmagS, AEAmag, theta, phi, ell
		if (EAMBSIGNflag == 0) then
			AEAx= AEAmagS*abs(AEAmag)*(3d0*cos(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
		end if
		if (EAMBSIGNflag == 1) then
			AEAx= AEAmagS*AEAmag*(3d0*cos(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
		end if
end subroutine AEAxsub

subroutine AEAysub(AEAy, AEAmagS, AEAmag, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AEAy
	real(kind= dp), intent(in) :: AEAmagS, AEAmag, theta, phi, ell
		if (EAMBSIGNflag == 0) then
			AEAy= AEAmagS*abs(AEAmag)*(3d0*sin(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
		end if
		if (EAMBSIGNflag == 1) then
			AEAy= AEAmagS*AEAmag*(3d0*sin(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
		end if
end subroutine AEAysub

subroutine AEAzsub(AEAz, AEAmagS, AEAmag, theta, ell)
	implicit none
	real(kind= dp), intent(out) :: AEAz
	real(kind= dp), intent(in) :: AEAmagS, AEAmag, theta, ell
		if (EAMBSIGNflag == 0) then
			AEAz= AEAmagS*abs(AEAmag)*((3d0*(cos(theta))**2d0- 1d0)/(sqrt(ell)))
		end if
		if (EAMBSIGNflag == 1) then
			AEAz= AEAmagS*AEAmag*((3d0*(cos(theta))**2d0- 1d0)/(sqrt(ell)))
		end if
end subroutine AEAzsub

! ----------------------------------------------------

! DEFINE PARALLEL ELECTRIC FIELD PARTICLE ACCELERATION TERMS:

subroutine AEPxsub(AEPx, AEPmagS, AEPmag, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AEPx
	real(kind= dp), intent(in) :: AEPmagS, AEPmag, theta, phi, ell
		AEPx= AEPmagS*abs(AEPmag)*(3d0*cos(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
end subroutine AEPxsub

subroutine AEPysub(AEPy, AEPmagS, AEPmag, theta, phi, ell)
	implicit none
	real(kind= dp), intent(out) :: AEPy
	real(kind= dp), intent(in) :: AEPmagS, AEPmag, theta, phi, ell
		AEPy= AEPmagS*abs(AEPmag)*(3d0*sin(phi)*cos(theta)*sin(theta)/(sqrt(ell)))
end subroutine AEPysub

subroutine AEPzsub(AEPz, AEPmagS, AEPmag, theta, ell)
	implicit none
	real(kind= dp), intent(out) :: AEPz
	real(kind= dp), intent(in) :: AEPmagS, AEPmag, theta, ell
		AEPz= AEPmagS*abs(AEPmag)*((3d0*(cos(theta))**2d0- 1d0)/(sqrt(ell)))
end subroutine AEPzsub

! ----------------------------------------------------

! DEFINE CARTESIAN COORDINATE NET PARTICLE FORCE TERMS:

subroutine Axsub(Ax, AMx, AGx, AEAx, AEPx)
	implicit none
	real(kind= dp), intent(out) :: Ax
	real(kind= dp), intent(in) :: AMx, AGx, AEAx, AEPx
		Ax= AMx+ AGx+ AEAx+ AEPx
end subroutine Axsub

subroutine Aysub(Ay, AMy, AGy, AEAy, AEPy)
	implicit none
	real(kind= dp), intent(out) :: Ay
	real(kind= dp), intent(in) :: AMy, AGy, AEAy, AEPy
		Ay= AMy+ AGy+ AEAy+ AEPy
end subroutine Aysub

subroutine Azsub(Az, AMz, AGz, AEAz, AEPz)
	implicit none
	real(kind= dp), intent(out) :: Az
	real(kind= dp), intent(in) :: AMz, AGz, AEAz, AEPz
		Az= AMz+ AGz+ AEAz+ AEPz
end subroutine Azsub

! ----------------------------------------------------

! DEFINE LINEAR INTERPOLATION SUBROUTINE:

subroutine LinInterpsub(yLinInterp, xLinInterp, yLinInterp1, yLinInterp2, xLinInterp1, xLinInterp2)
	implicit none
	real(kind= dp), intent(out) :: yLinInterp
	real(kind= dp), intent(in) :: xLinInterp, yLinInterp1, yLinInterp2, xLinInterp1, xLinInterp2
		yLinInterp= yLinInterp1+ (((xLinInterp- xLinInterp1)/(xLinInterp2- xLinInterp1))*(yLinInterp2- yLinInterp1))
end subroutine LinInterpsub

! ----------------------------------------------------

! TERMINAL OUTPUT FONT COLORS:

! [90m=dark grey           [30m=black
! [91m=peach               [31m=red
! [92m=light green         [32m=green
! [93m=light yellow        [33m=yellow
! [94m=light blue          [34m=blue
! [95m=pink                [35m=purple
! [96m=light aqua          [36m=aqua
! [97m=pearl white

! ----------------------------------------------------

end module KineticMainParams
