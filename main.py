import numpy as np
import wavefun
import SpinOpperators
import matplotlib.pyplot as plt
import wavefun_print
import wavefun_set
import wavefun
import sequence


cmMHz=2.997/1000.
nu=93.301 *1e3# MHz
g = 1.9899
B=3350
D=2800
E=-100
Azz=2
Axx=2
Ayy=2
gbetaDh=13.99624604872811*g
gammaN=1.07197/350

S=1/2
I=1/2

tpi2=40/1000. #us
tpi=tpi2*2
nu1=0.25/tpi2
B1=nu1/gbetaDh
tpiRF=82
nuRF=0.5/tpiRF
B1RF=nuRF/gammaN
print('B1=%g BRF=%g' % (B1,B1RF))

T=1
tau=0.24 #us

wS=wavefun_set.WaveFunBasic(S)
wI=wavefun_set.WaveFunBasic(I)
w=wavefun_set.WaveFunCombine(wS,wI)
w_rot=wavefun.RotationWF1(0,0,0,w)




wavefun_print.printHM_HTML('S',w_rot, 'Sz', w_rot.Sz[0], 'Sx', w_rot.Sx[0], 'Sy', w_rot.Sy[0], 'Sp', w_rot.Sp[0], 'Sm', w_rot.Sm[0])
wavefun_print.printHM_HTML('I',w_rot, 'Sz', w_rot.Sz[1], 'Sx', w_rot.Sx[1], 'Sy', w_rot.Sy[1], 'Sp', w_rot.Sp[1], 'Sm', w_rot.Sm[1])

Ham_fine=D*(w_rot.Sz[0]*w_rot.Sz[0]-w_rot.JJ1[0]*w_rot.E[0]/3.0)# +E*(w_rot.Sx[0]*w_rot.Sx[0] - w_rot.Sy[0]*w_rot.Sy[0])
Ham_SFS=Azz*w_rot.Sz[0]*w_rot.Sz[1] +  Axx*w_rot.Sx[0]*w_rot.Sx[1] + Ayy*w_rot.Sy[0]*w_rot.Sy[1]
Ham_EZeem=gbetaDh*w.Sz[0]
Ham_NZeem=gammaN*w.Sz[1]

wavefun_print.printHM_HTML('Ham_fine',w, 'WFfun Matrix', w.WF_matrix, 'Fine', Ham_fine, 'SFS', Ham_SFS, 'Ham_EZeem', Ham_EZeem*B)

w1=wavefun.ApllyHamiltonian(w, Ham_fine+Ham_SFS+Ham_EZeem*B+Ham_NZeem*B)
# w1=wavefun.ApllyHamiltonian(w, Ham_EZeem*B+Ham_NZeem*B)
# w2=wavefun.ExtructSubMatrix(w1,[0,1,2,3,4,5])
wavefun_print.printHM_HTML('Ham_fine1',w1, 'WFfun Matrix', w1.WF_matrix, 'Sz', w1.Sz[0], 'Sx', w1.Sx[0], 'E', wavefun.MakeHamFromEnergy(w1.Energy))

# wavefun_print.plotEnergy(w1)

matrH0=Ham_fine+Ham_SFS+gammaN*w.Sz[1]*B
matrHF=gbetaDh*w.Sz[0]
matrHP=gbetaDh*B1*w.Sy[0]
matrHRF=gammaN*B1RF*w.Sx[1]

N=1024
t=np.linspace(0,(N-1)*0.004,N)
dB0=-0.0
# FIDAver=lambda dB :sequence.FID(w,T,matrH0,matrHF,matrHP,dB0+dB,tpi2,t)
# Mx,My,Mz=FIDAver(0)
# # Mx,My,Mz=sequence.MakeExp2Decay(Mx,My,Mz)
# # Mx,My,Mz=sequence.AverageByGauss(FIDAver,0.05)
# sequence.MakePlotSeq(t,Mx,My)

# HanEchoAver=lambda dB :sequence.Hane_Echo(w1,T,matrH0,matrHF,matrHP,dB0+dB,tau,tpi2,t)
# Mx,My,Mz=HanEchoAver(0)
# # Mx,My,Mz=sequence.AverageByGauss(HanEchoAver,0.1)
# sequence.MakePlotSeq(t,Mx,My)


HanEchoESEEM=lambda t :sequence.Hane_Echo(w1,T,matrH0,matrHF,matrHP,dB0,tau+t,tpi2,0.0)
N=8192
tau1=np.linspace(0,(N-1)*0.004,N)
Mx,My,Mz=sequence.OneDim(HanEchoESEEM,tau1)

sequence.MakePlotSeq(tau1,Mx,My)
