import json
import numpy as np
import matplotlib.pyplot as plt

NNt= 13
Stot= 1
Nf= 1
NqG= 25
NVparG= 28
NVperp1G= 28
NVperp2G= 28
RE= 6380e3

Qindval= 7- 1
print('Qindval= ', Qindval)

if ((NVparG != NVperp1G) or (NVparG != NVperp2G) or (NVperp1G != NVperp2G)):
	print('ERROR: UNEQUAL VALUES OF NVparG= ', NVparG, ', NVperp1G= ', NVperp1G, ', NVperp2G= ', NVperp2G)

ranksize= 1
SpinupFlag= 1

if (SpinupFlag == 1):
	#datadir= '/Users/robertalbarran/Library/Mobile Documents/com~apple~CloudDocs/Desktop/ACADEMICS/DISSERTATION/VPCFinalSims/W0F2ncscOUTPUTa0/KMIODataW0F2ncscOUTPUTa0/'
	#datadir= '/Users/robertalbarran/Library/Mobile Documents/com~apple~CloudDocs/Desktop/ACADEMICS/DISSERTATION/VPCFinalSims/W0F1ncscOUTPUTa0/KMIODataW0F1ncscOUTPUTa0/'
	datadir= '/Users/robertalbarran/Desktop/KAOS_M1/KAOSDataSpinup/'
if (SpinupFlag == 0):
	datadir= '/Users/robertalbarran/Desktop/KAOS_M1/KAOSDataOutput/'

EAmag= np.zeros(Stot*Nf*NNt*NqG)
EAmag= EAmag.reshape(Stot, Nf, NNt, NqG)
M0= np.zeros(Stot*Nf*NNt*NqG)
M0= M0.reshape(Stot, Nf, NNt, NqG)
M1Perp1= np.zeros(Stot*Nf*NNt*NqG)
M1Perp1= M1Perp1.reshape(Stot, Nf, NNt, NqG)
M1Perp2= np.zeros(Stot*Nf*NNt*NqG)
M1Perp2= M1Perp2.reshape(Stot, Nf, NNt, NqG)
M1Par= np.zeros(Stot*Nf*NNt*NqG)
M1Par= M1Par.reshape(Stot, Nf, NNt, NqG)
Nph= np.zeros(Stot*Nf*NNt*NqG*NVperp1G*NVperp2G*NVparG)
Nph= Nph.reshape(Stot, Nf, NNt, NqG, NVperp1G*NVperp2G*NVparG)
NphVperp1Vpar= np.zeros(Stot*Nf*NNt*NqG*NVperp1G*NVparG)
NphVperp1Vpar= NphVperp1Vpar.reshape(Stot, Nf, NVperp1G, NVparG, NqG, NNt)
NphVperp2Vpar= np.zeros(Stot*Nf*NNt*NqG*NVperp2G*NVparG)
NphVperp2Vpar= NphVperp2Vpar.reshape(Stot, Nf, NVperp2G, NVparG, NqG, NNt)
NphVperp1Vperp2= np.zeros(Stot*Nf*NNt*NqG*NVperp1G*NVperp2G)
NphVperp1Vperp2= NphVperp1Vperp2.reshape(Stot, Nf, NVperp1G, NVperp2G, NqG, NNt)
Time= np.zeros(Stot*Nf*NNt)
Time= Time.reshape(Stot, Nf, NNt)
rGCGT= np.zeros(Stot*Nf*NNt*NqG)
rGCGT= rGCGT.reshape(Stot, Nf, NNt, NqG)
pGCGT= np.zeros(Stot*Nf*NNt*NqG)
pGCGT= pGCGT.reshape(Stot, Nf, NNt, NqG)
phiGCGT= np.zeros(Stot*Nf*NNt*NqG)
phiGCGT= phiGCGT.reshape(Stot, Nf, NNt, NqG)

s= 0
f= 0

for nn in range(0, NNt- 1):
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_TimeTfort.bin'
	with open(path, 'rb') as f:
		Time[0, 0, nn]= np.fromfile(f)
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_rGCGTfort.bin'
	with open(path, 'rb') as f:
		rGCGT[0, 0, nn, :]= np.fromfile(f)
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_pGCGTfort.bin'
	with open(path, 'rb') as f:
		pGCGT[0, 0, nn, :]= np.fromfile(f)
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_phiGCGTfort.bin'
	with open(path, 'rb') as f:
		phiGCGT[0, 0, nn, :]= np.fromfile(f)

	for Qind in range(0, NqG- 1):
		path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_' + str(Qind+ 1) + '_N2PerpphRTfort.bin'
		with open(path, 'rb') as f:
			Nph[0, 0, nn, Qind, :]= np.fromfile(f)

Nph= np.transpose(Nph)
Nph= Nph.reshape(Stot, Nf, NVperp1G, NVperp2G, NVparG, NqG, NNt)
rGCGT= rGCGT- RE
for nn in range(0, NNt- 1):
	for Qind in range(0, NqG- 1):
		for Vperp1ind in range(0, NVperp1G- 1):
			for Vparind in range(0, NVparG- 1):
				NphVperp1Vpar[0, 0, Vperp1ind, Vparind, Qind, nn]= np.sum(Nph[0, 0, Vperp1ind, :, Vparind, Qind, nn])

		for Vperp2ind in range(0, NVperp2G- 1):
			for Vparind in range(0, NVparG- 1):
				NphVperp2Vpar[0, 0, Vperp2ind, Vparind, Qind, nn]= np.sum(Nph[0, 0, :, Vperp2ind, Vparind, Qind, nn])

		for Vperp1ind in range(0, NVperp1G- 1):
			for Vperp2ind in range(0, NVperp2G- 1):
				NphVperp1Vperp2[0, 0, Vperp1ind, Vperp2ind, Qind, nn]= np.sum(Nph[0, 0, Vperp1ind, Vperp2ind, :, Qind, nn])

for nn in range(0, NNt- 1):
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_EAmagRTfort.bin'
	with open(path, 'rb') as f:
		EAmag[0, 0, nn, :]= np.fromfile(f)
	
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M0phRTfort.bin'
	with open(path, 'rb') as f:
		M0[0, 0, nn, :]= np.fromfile(f)
	
	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M1Perp1phRTfort.bin'
	with open(path, 'rb') as f:
		M1Perp1[0, 0, nn, :]= np.fromfile(f)

	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M1Perp2phRTfort.bin'
	with open(path, 'rb') as f:
		M1Perp2[0, 0, nn, :]= np.fromfile(f)

	path= datadir + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M1ParphRTfort.bin'
	with open(path, 'rb') as f:
		M1Par[0, 0, nn, :]= np.fromfile(f)

#---------------------------------
# PLOT DATA:

alt= np.zeros(NqG)
alt[0]= 1
for r in range(0, NqG- 1):
        alt[r+ 1]= alt[r]+ 1

nn= NNt- 2
s= 1- 1
f= 1- 1
print('rGCGT= ', 1e-3*rGCGT[s, f, nn, :])
TTr, rrT= np.meshgrid(Time[s, f, :], 1e-3*rGCGT[s, f, nn, :])

fig,axarr= plt.subplots(1, 2)
figA= 0
EAplt= axarr[figA].pcolormesh(Time[s, f, :]/3600e0, 1e-3*rGCGT[s, f, nn, :], np.transpose(EAmag[s, f, :, :]), shading= 'gouraud')
fig.colorbar(EAplt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1])
axarr[figA].set_title(Title)
axarr[figA].set_xlabel('Time [hr]')
axarr[figA].set_ylabel('$r$ [km]')
figB= 1
axarr[figB].plot(EAmag[s, f, nn, :], 1e-3*rGCGT[s, f, nn, :])
axarr[figB].set_xlabel('$E_A$ [V/m]')
axarr[figB].set_ylabel('$r$ [km]')
plt.show()

fig,axarr= plt.subplots(1, 3)
figA= 0
M0plt= axarr[figA].pcolormesh(Time[s, f, :]/3600e0, 1e-3*rGCGT[s, f, nn, :], np.transpose(np.log(M0[s, f, :, :])), shading= 'gouraud')
fig.colorbar(M0plt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) 
axarr[figA].set_title(Title)
axarr[figA].set_xlabel('Time [hr]')
axarr[figA].set_ylabel('$r$ [km]')
figB= 1
axarr[figB].plot(M0[s, f, 1, :], 1e-3*rGCGT[s, f, nn, :], 'k-', label= 'initial time')
axarr[figB].plot(M0[s, f, nn, :], 1e-3*rGCGT[s, f, nn, :], 'k.',label= 'current time')
axarr[figB].set_xlabel('$n$ [m$^{-3}$]')
axarr[figB].set_ylabel('$r$ [km]')
axarr[figB].legend()
figC= 2
axarr[figC].plot(Time[s, f, :]/3600e0, M0[s, f, :, Qindval])
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) + ', r [km]= ' + str(round(1e-3*rGCGT[s, f, nn, Qindval]))
axarr[figC].set_title(Title)
axarr[figC].set_xlabel('Time [hr]')
axarr[figC].set_ylabel('$n$ [m$^{-3}$]')
plt.show()

fig,axarr= plt.subplots(1, 2)
figA= 0
M1Perp1plt= axarr[figA].pcolormesh(Time[s, f, :]/3600e0, 1e-3*rGCGT[s, f, nn, :], np.transpose(M1Perp1[s, f, :, :]), shading= 'gouraud')
fig.colorbar(M1Perp1plt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) 
axarr[figA].set_title(Title)
axarr[figA].set_xlabel('Time [hr]')
axarr[figA].set_ylabel('$r$ [km]')
figB= 1
axarr[figB].plot(M1Perp1[s, f, nn, :], 1e-3*rGCGT[s, f, nn, :])
axarr[figB].set_xlabel('$u_{\perp 1}$ [m/s]')
axarr[figB].set_ylabel('$r$ [km]')
plt.show()

fig,axarr= plt.subplots(1, 2)
figA= 0
M1Perp2plt= axarr[figA].pcolormesh(Time[s, f, :]/3600e0, 1e-3*rGCGT[s, f, nn, :], np.transpose(M1Perp2[s, f, :, :]), shading= 'gouraud')
fig.colorbar(M1Perp2plt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) 
axarr[figA].set_title(Title)
axarr[figA].set_xlabel('Time [hr]')
axarr[figA].set_ylabel('$r$ [km]')
figB= 1
axarr[figB].plot(M1Perp2[s, f, nn, :], 1e-3*rGCGT[s, f, nn, :])
axarr[figB].set_xlabel('$u_{\perp 2}$ [m/s]')
axarr[figB].set_ylabel('$r$ [km]')
plt.show()

fig,axarr= plt.subplots(1, 3)
figA= 0
M1Parplt= axarr[figA].pcolormesh(Time[s, f, :]/3600e0, 1e-3*rGCGT[s, f, nn, :], np.transpose(M1Par[s, f, :, :]), shading= 'gouraud')
fig.colorbar(M1Parplt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1])
axarr[figA].set_title(Title)
axarr[figA].set_xlabel('Time [hr]')
axarr[figA].set_ylabel('$r$ [km]')
figB= 1
axarr[figB].plot(M1Par[s, f, nn, :], 1e-3*rGCGT[s, f, nn, :])
axarr[figB].set_xlabel('$u_\parallel$ [m/s]')
axarr[figB].set_ylabel('$r$ [km]')
figC= 2
axarr[figC].plot(Time[s, f, :]/3600e0, M1Par[s, f, :, Qindval])
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) + ', r [km]= ' + str(round(1e-3*rGCGT[s, f, nn, Qindval]))
axarr[figC].set_title(Title)
axarr[figC].set_xlabel('Time [hr]')
axarr[figC].set_ylabel('$u_\parallel$ [m/s]')
plt.show()

#for nn in range(1, NNt- 1):
nn= 1
fig,axarr= plt.subplots(1, 3)
NphVperp1Vparplt= axarr[0].pcolormesh(NphVperp1Vpar[s, f, :, :, Qindval, nn], shading= 'gouraud')
fig.colorbar(NphVperp1Vparplt)
axarr[0].set_xlabel('$N(v_{\perp 1},\:v_{\parallel})$')
#plt.title('Time [s]= ', Time[s, f, nn])

NphVperp2Vparplt= axarr[1].pcolormesh(NphVperp1Vperp2[s, f, :, :, Qindval, nn], shading= 'gouraud')
fig.colorbar(NphVperp2Vparplt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) + ', r [km]= ' + str(round(1e-3*rGCGT[s, f, nn, Qindval]))
axarr[1].set_title(Title)
axarr[1].set_xlabel('$N(v_{\perp 2},\:v_{\parallel})$')

NphVperp1Vperp2plt= axarr[2].pcolormesh(NphVperp2Vpar[s, f, :, :, Qindval, nn], shading= 'gouraud')
fig.colorbar(NphVperp1Vperp2plt)
axarr[2].set_xlabel('$N(v_{\perp 1},\:v_{\perp 2})$')
plt.show()

nn= NNt- 2
fig,axarr= plt.subplots(1, 3)
NphVperp1Vparplt= axarr[0].pcolormesh(NphVperp1Vpar[s, f, :, :, Qindval, nn], shading= 'gouraud')
fig.colorbar(NphVperp1Vparplt)
axarr[0].set_xlabel('$N(v_{\perp 1},\:v_{\parallel})$')
#plt.title('Time [s]= ', Time[s, f, nn])

NphVperp2Vparplt= axarr[1].pcolormesh(NphVperp1Vperp2[s, f, :, :, Qindval, nn], shading= 'gouraud')
fig.colorbar(NphVperp2Vparplt)
Title= 'L [RE]= ' + str(pGCGT[s, f, nn, 1]) + ', r [km]= ' + str(round(1e-3*rGCGT[s, f, nn, Qindval])) \
	+ ', Time [hr]= ' + str(Time[0, 0, nn]/3600e0)
axarr[1].set_title(Title)
axarr[1].set_xlabel('$N(v_{\perp 2},\:v_{\parallel})$')

NphVperp1Vperp2plt= axarr[2].pcolormesh(NphVperp2Vpar[s, f, :, :, Qindval, nn], shading= 'gouraud')
fig.colorbar(NphVperp1Vperp2plt)
axarr[2].set_xlabel('$N(v_{\perp 1},\:v_{\perp 2})$')
plt.show()

print('End Sim Time [hr]= ', (25e3/60e0)/60e0)
print('Current Time [hr]= ', (Time[s, f, NNt- 2]/3600e0)) 
