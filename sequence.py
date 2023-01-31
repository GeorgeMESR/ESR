import numpy as np
import scipy.constants as const
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import scipy.fftpack as fftp
import wavefun_print

def Ro0(WF,T):
	ro0=np.zeros((WF.N,WF.N),dtype='complex')
	for i in range(WF.N):
		ro0[i,i]=np.exp(-WF.Energy[i]*1e6/T/const.physical_constants['kelvin-hertz relationship'][0])
	ro0=ro0/np.trace(ro0)
	return ro0

R_opperator= lambda t, matr : np.matrix(linalg.expm(-1j*2.0*np.pi*matr*t))

def NormM(Mx,My,Mz):
	Mx-=(np.max(Mx)+np.min(Mx))/2.0
	My-=(np.max(My)+np.min(My))/2.0
	return Mx,My,Mz

def MakeExp2Decay(Mx,My,Mz):
	N=len(Mz)
	x=np.linspace(0,N,N)
	exp2=np.exp(-(x**2)*4.0/N/N)
	Mx*=exp2
	My*=exp2
	return Mx,My,Mz

def AverageByGauss(fun,dH):
	N=32
	db=np.linspace(-dH,dH,N)
	gauss=lambda x, dH : np.exp(-x*x/2.0/dH/dH)
	Mx,My,Mz=fun(db[0])
	nc=gauss(db[0],dH)
	Mx*=nc
	My*=nc
	Mz*=nc
	for i in range(1,N):
		Mx1, My1, Mz1 = fun(db[i])
		nc = gauss(db[i], dH)
		Mx += Mx1 * nc
		My += My1 * nc
		Mz += Mz1 * nc
	return (Mx,My,Mz)

def AverageByGaussFFT(fun,dH, dt):
	Mx,My,Mz=fun(0)
	gauss=lambda x, dH : np.exp(-x*x/2.0/dH/dH/2.8/2.8)

	yFFT = fftp.fft(Mx + 1j * My)
	yA = yFFT
	xFFT = fftp.fftfreq(len(Mx), d=dt)
	d = np.array([xFFT, yA])
	d1 = d.transpose()
	# d1 = d1[d1[:, 0].argsort()]

	N=256
	db=np.linspace(-dH*3,dH*3,N)
	k1=(db[1]-db[0])/(d1[1,0]-d1[0,0])
	FM1=np.zeros(len(Mx),dtype='complex')
	for i in range(len(Mx)):
		for j in range(len(db)):
			gs = gauss(db[j], dH)
			j1=np.argmin(np.abs(d1[:,0]-d1[i,0]-db[j]))
			FM1[j1]+=d1[i,1]*gs
	kgs=sum(gauss(db,dH))
	FM1=FM1/complex(kgs)
	# plt.plot(d1[:,0],np.abs(d1[:,1]),d1[:,0],np.abs(FM1))
	# plt.show()
	Y=fftp.ifft(FM1)
	Mx=Y.real
	My=Y.imag
	return Mx,My,Mz


def OneDim(fun, t):
	Mx=np.zeros(len(t))
	My=np.zeros(len(t))
	Mz=np.zeros(len(t))
	for i in range(len(t)):
		(Mx[i],My[i],Mz[i])=fun(t[i])
	(Mx, My, Mz) = NormM(Mx, My, Mz)
	return Mx,My,Mz

def TwoDim(fun, t1, t2):
	Mx=np.zeros((len(t1),len(t2)))
	My=np.zeros((len(t1),len(t2)))
	Mz=np.zeros((len(t1),len(t2)))
	for i in range(len(t1)):
		for j in range(len(t2)):
			Mx[i,j],My[i,j],Mz[i,j]=fun(t1[i],t2[j])
	(Mx, My, Mz) = NormM(Mx, My, Mz)
	return Mx,My,Mz

def FID(WF, T, matrH0, matrHF, matrHP, dB,  tpulse, tadc):
	# first pi/2 pulse
	ro0=Ro0(WF,T)
	R_FirstPulse = R_opperator(tpulse, matrH0+matrHF*dB+matrHP)
	ro1=R_FirstPulse*ro0*np.transpose(np.conjugate(R_FirstPulse))
	if(type(tadc)=='float'):
		tadc=[tadc]
	Mx=np.zeros(len(tadc))
	My=np.zeros(len(tadc))
	Mz=np.zeros(len(tadc))
	for i in range(len(tadc)):
		R_FID_time = R_opperator(tadc[i], matrH0+matrHF*dB)
		ro2=R_FID_time*ro1*np.transpose(np.conjugate(R_FID_time))
		matrMx = ro2 * WF.Sx[0]
		matrMy = ro2 * WF.Sy[0]
		matrMz = ro2 * WF.Sz[0]
		# wavefun_print.printHM_HTML('ro2', WF,'R_FirstPulse',R_FirstPulse, 'ro1', ro1, 'ro2', ro2, 'matrMx', matrMx, 'matrMy',matrMy, 'matrMz', matrMz)
		Mx[i] = ((matrMx).trace()).real
		My[i] = ((matrMy).trace()).real
		Mz[i] = ((matrMz).trace()).real
	if(len(tadc)==1):
		Mx = Mx[0]
		My = My[0]
		Mz = Mz[0]
	return Mx,My,Mz



def Hane_Echo(WF, T, matrH0,matrHF, matrHP,  dB, tau, tpulse, tadc):
	# first pi/2 pulse
	ro0=Ro0(WF,T)
	R_FirstPulse = R_opperator(tpulse, matrH0+matrHF*dB+matrHP)
	R_FirstTau = R_opperator(tau-tpulse, matrH0+matrHF*dB)
	R_SecPulse = R_opperator(tpulse*2, matrH0+matrHF*dB+matrHP)
	ro1=R_SecPulse*R_FirstTau*R_FirstPulse*ro0*np.transpose(np.conjugate(R_FirstPulse))*np.transpose(np.conjugate(R_FirstTau))*np.transpose(np.conjugate(R_SecPulse))
	# wavefun_print.printHM_HTML('ro', WF, 'matrHF',matrHF*dB, 'ro0', ro0, 'R_FirstPulse',R_FirstPulse, 'ro1',ro1, 'R_FirstTau', R_FirstTau, 'R_SecPulse', R_SecPulse )
	if(type(tadc)==float):
		tadc=[tadc]
		if(tadc[0]==0):
			tadc[0]=tau - tpulse*0.375 # if 0 then center of echo
	Mx=np.zeros(len(tadc))
	My=np.zeros(len(tadc))
	Mz=np.zeros(len(tadc))
	for i in range(len(tadc)):
		R_SecTau = R_opperator(tadc[i], matrH0+matrHF*dB)
		ro2=R_SecTau*ro1*np.transpose(np.conjugate(R_SecTau))
		matrMx = ro2 * WF.Sx[0]
		matrMy = ro2 * WF.Sy[0]
		matrMz = ro2 * WF.Sz[0]
		# wavefun_print.printHM_HTML('ro2', WF, 'ro2', ro2, 'matrMx', matrMx, 'matrMy',matrMy, 'matrMz', matrMz)
		Mx[i] = ((matrMx).trace()).real
		My[i] = ((matrMy).trace()).real
		Mz[i] = ((matrMz).trace()).real
	if(len(tadc)==1):
		Mx = Mx[0]
		My = My[0]
		Mz = Mz[0]
	return Mx,My,Mz



def MakePlotSeq(t,Mx,My):
	yFFT = fftp.fft(Mx + 1j * My)
	yA = np.abs(yFFT)
	xFFT = fftp.fftfreq(len(t), d=((t[1] - t[0])))
	d = np.array([xFFT, yA])
	d1 = d.transpose()
	d1 = d1[d1[:, 0].argsort()]

	gridsize = (2, 1)
	axis1 = plt.subplot2grid(gridsize, (0, 0))
	axis2 = plt.subplot2grid(gridsize, (1, 0))
	axis1.set_xlabel('$t, \mu s$', fontsize=20)
	axis1.set_ylabel('', fontsize=20)
	axis2.set_xlabel('f, MHz', fontsize=20)
	axis2.set_ylabel('', fontsize=20)


	axis1.plot(t, Mx, 'g', t, My, 'y',t, np.abs(Mx + 1j * My), 'r')
	axis2.plot(d1[:, 0] , d1[:, 1], 'k')

	plt.savefig("fft.png")
	plt.show()

