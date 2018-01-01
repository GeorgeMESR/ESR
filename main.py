import wavefun
import SpinOpperators
import numpy as np
import pylab

cmMHz=2.997/1000.
g = 1.9899
B=300
B02=-914.8493/3*cmMHz
B04=-24.4023/60*cmMHz
B06=0.1905/1260*cmMHz
B44=134.2854/60*cmMHz
gbetaDh=0.01399624604872811*g
S=[]
S.append(7/2)
z=SpinOpperators.MakeZeeman(S)
f=SpinOpperators.MakeFineStructure(S)
ham=wavefun.HamiltonianMatrix(S)
Z=z.Sz.M*gbetaDh
F=f.O_02.M*B02+f.O_04.M*B04+f.O_06.M*B06+f.O_44.M*B44
NSp=1024
BMax=1400
Sp={'B':np.zeros(NSp),'E':np.zeros((NSp,ham.wavefunctions.N))}
Sp['B']=np.zeros(NSp)
Sp['E']=np.zeros((NSp,ham.wavefunctions.N))
for i in range(NSp):
    Sp['B'][i]=float(i)/float(NSp-1)*BMax
    ham.M=Z* Sp['B'][i]+F
    w, v=np.linalg.eig(ham.M)
    Sp['E'][i]=w

pylab.plot(Sp['B'],Sp['E'])
pylab.show()
#print (ham.printHM())





