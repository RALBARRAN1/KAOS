import json
import numpy as np
import matplotlib.pyplot as plt

NNt= 9
Stot= 1
Nf= 1
NqG= 12
NVparG= 28
NVperp1G= 28
NVperp2G= 28

Qindval= 6

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

EAmag= np.zeros(NNt*NqG)
EAmag= EAmag.reshape(NNt, NqG)
M0= np.zeros(NNt*NqG)
M0= M0.reshape(NNt, NqG)
M1Perp1= np.zeros(NNt*NqG)
M1Perp1= M1Perp1.reshape(NNt, NqG)
M1Perp2= np.zeros(NNt*NqG)
M1Perp2= M1Perp2.reshape(NNt, NqG)
M1Par= np.zeros(NNt*NqG)
M1Par= M1Par.reshape(NNt, NqG)
Nph= np.zeros(NNt*NqG*NVperp1G*NVperp2G*NVparG)
Nph= Nph.reshape(NNt, NqG, NVperp1G*NVperp2G*NVparG)
NphVperp1Vpar= np.zeros(NNt*NqG*NVperp1G*NVparG)
NphVperp1Vpar= NphVperp1Vpar.reshape(NVperp1G, NVparG, NqG, NNt)
NphVperp2Vpar= np.zeros(NNt*NqG*NVperp2G*NVparG)
NphVperp2Vpar= NphVperp2Vpar.reshape(NVperp2G, NVparG, NqG, NNt)
NphVperp1Vperp2= np.zeros(NNt*NqG*NVperp1G*NVperp2G)
NphVperp1Vperp2= NphVperp1Vperp2.reshape(NVperp1G, NVperp2G, NqG, NNt)
Time= np.zeros(NNt)

for nn in range(0, NNt- 1):
	path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_TimeTfort.bin'
	with open(path, 'rb') as f:
		Time[nn]= np.fromfile(f)

for nn in range(0, NNt- 1):
	for Qind in range(0, NqG- 1):
		path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_' + str(Qind+ 1) + '_N2PerpphRTfort.bin'
		with open(path, 'rb') as f:
			Nph[nn, Qind, :]= np.fromfile(f)

Nph= np.transpose(Nph)
print(Nph.shape)

Nph= Nph.reshape(NVperp1G, NVperp2G, NVparG, NqG, NNt)
print(Nph.shape)
for nn in range(0, NNt- 1):
	for Qind in range(0, NqG- 1):
		for Vperp1ind in range(0, NVperp1G- 1):
			for Vparind in range(0, NVparG- 1):
				NphVperp1Vpar[Vperp1ind, Vparind, Qind, nn]= np.sum(Nph[Vperp1ind, :, Vparind, Qind, nn])

for nn in range(0, NNt- 1):
        for Qind in range(0, NqG- 1):
                for Vperp2ind in range(0, NVperp2G- 1):
                        for Vparind in range(0, NVparG- 1):
                                NphVperp2Vpar[Vperp2ind, Vparind, Qind, nn]= np.sum(Nph[:, Vperp2ind, Vparind, Qind, nn])

for nn in range(0, NNt- 1):
        for Qind in range(0, NqG- 1):
                for Vperp1ind in range(0, NVperp1G- 1):
                        for Vperp2ind in range(0, NVperp2G- 1):
                                NphVperp1Vperp2[Vperp1ind, Vperp2ind, Qind, nn]= np.sum(Nph[Vperp1ind, Vperp2ind, :, Qind, nn])

for nn in range(0, NNt- 1):
	path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_EAmagRTfort.bin'
	with open(path, 'rb') as f:
		EAmag[nn, :]= np.fromfile(f)
	
	path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M0phRTfort.bin'
	with open(path, 'rb') as f:
		M0[nn, :]= np.fromfile(f)
	
	path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M1Perp1phRTfort.bin'
	with open(path, 'rb') as f:
		M1Perp1[nn, :]= np.fromfile(f)

	path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M1Perp2phRTfort.bin'
	with open(path, 'rb') as f:
		M1Perp2[nn, :]= np.fromfile(f)

	path= datadir + str(0) + '_' + str(nn+ 1) + '_' + str(1) + '_' + str(1) + '_M1ParphRTfort.bin'
	with open(path, 'rb') as f:
		M1Par[nn, :]= np.fromfile(f)

#---------------------------------
# PLOT DATA:

alt= np.zeros(NqG)
alt[0]= 1
for r in range(0, NqG- 1):
        alt[r+ 1]= alt[r]+ 1

nn= NNt- 2

fig,axarr= plt.subplots(1, 2)
figA= 0
EAplt= axarr[figA].pcolormesh(np.transpose(EAmag[:, :]))
fig.colorbar(EAplt)
figB= 1
axarr[figB].plot(EAmag[nn, :], alt[:])
axarr[figB].set_xlabel('$E_A$ [V/m]')
plt.show()

fig,axarr= plt.subplots(1, 3)
figA= 0
M0plt= axarr[figA].pcolormesh(np.transpose(M0[:, :]))
fig.colorbar(M0plt)
figB= 1
axarr[figB].plot(M0[1, :], alt[:], 'k-', label= 'initial time')
axarr[figB].plot(M0[nn, :], alt[:], 'k.',label= 'current time')
axarr[figB].set_xlabel('$n$ [m$^{-3}$]')
axarr[figB].legend()
figC= 2
axarr[figC].plot(Time[:], M0[:, Qindval])
axarr[figC].set_xlabel('Time [s]')
axarr[figC].set_ylabel('$n$ [m$^{-3}$]')
plt.show()

fig,axarr= plt.subplots(1, 2)
figA= 0
M1Perp1plt= axarr[figA].pcolormesh(np.transpose(M1Perp1[:, :]))
fig.colorbar(M1Perp1plt)
figB= 1
axarr[figB].plot(M1Perp1[nn, :], alt[:])
axarr[figB].set_xlabel('$u_{\perp 1}$ [m/s]')
plt.show()

fig,axarr= plt.subplots(1, 2)
figA= 0
M1Perp2plt= axarr[figA].pcolormesh(np.transpose(M1Perp2[:, :]))
fig.colorbar(M1Perp2plt)
figB= 1
axarr[figB].plot(M1Perp2[nn, :], alt[:])
axarr[figB].set_xlabel('$u_{\perp 2}$ [m/s]')
plt.show()

fig,axarr= plt.subplots(1, 3)
figA= 0
M1Parplt= axarr[figA].pcolormesh(np.transpose(M1Par[:, :]))
fig.colorbar(M1Parplt)
figB= 1
axarr[figB].plot(M1Par[nn, :], alt[:])
axarr[figB].set_xlabel('$u_\parallel$ [m/s]')
figC= 2
axarr[figC].plot(Time[:], M1Par[:, Qindval])
axarr[figC].set_xlabel('Time [s]')
axarr[figC].set_ylabel('$u_\parallel$ [m/s]')
plt.show()

for nn in range(NNt- 5, NNt- 1):

	fig,axarr= plt.subplots(1, 3)
	NphVperp1Vparplt= axarr[0].pcolormesh(NphVperp1Vpar[:, :, Qindval, nn])
	fig.colorbar(NphVperp1Vparplt)
	axarr[0].set_xlabel('$N(v_{\perp 1},\:v_{\parallel})$')
	#plt.title('Time [s]= ', Time[nn])

	NphVperp2Vparplt= axarr[1].pcolormesh(NphVperp1Vperp2[:, :, Qindval, nn])
	fig.colorbar(NphVperp2Vparplt)
	axarr[1].set_xlabel('$N(v_{\perp 2},\:v_{\parallel})$')

	NphVperp1Vperp2plt= axarr[2].pcolormesh(NphVperp2Vpar[:, :, Qindval, nn])
	fig.colorbar(NphVperp1Vperp2plt)
	axarr[2].set_xlabel('$N(v_{\perp 1},\:v_{\perp 2})$')
	plt.show()

print('End Sim Time [hr]= ', (0.18e5/60e0)/60e0)
print('Current Time [hr]= ', (Time[NNt- 2]/60e0)/60e0) 
