import numpy as np
import scipy
import math
import wavefun_print
import wavefun_set
from scipy.spatial.transform import Rotation as R

def MakeHamFromEnergy(energy :np.ndarray ):
	ham=np.zeros((len(energy),len(energy)),dtype='complex')
	for i in range(len(energy)):
		ham[i,i]=energy[i]
	return ham


def CompareWaveFunctionsLeft(wf1, wf2):
	n2=len(wf2.WFarray)
	wf1a=wf1.WFarray[0:n2]
	return np.array_equal(wf1a,wf2.WFarray)

def CompareWaveFunctionsRight(wf1, wf2):
	n1=len(wf1.WFarray)
	n2=len(wf2.WFarray)
	wf1a=wf1.WFarray[n1-n2:n1]
	return np.array_equal(wf1a,wf2.WFarray)

def RotationMatrice(Q,F,G):
	r=R.from_euler('zyx', [F, Q, G], degrees=True)
	return r.as_matrix()

def RotationWF1(Q,F,G,WF):
	M=RotationMatrice(Q,F,G)
	# print(M)
	WF1=WF.copy()
	tp=str(type(WF))
	if(tp.find('WaveFunBasic')!=-1):
		WF1.Sx=WF.Sx*M[0,0]+WF.Sy*M[1,0]+WF.Sz*M[2,0]
		WF1.Sy=WF.Sx*M[0,1]+WF.Sy*M[1,1]+WF.Sz*M[2,1]
		WF1.Sz=WF.Sx*M[0,2]+WF.Sy*M[1,2]+WF.Sz*M[2,2]
		WF1.Sp=WF1.Sx+1j*WF1.Sy
		WF1.Sm=WF1.Sx-1j*WF1.Sy
	if(tp.find('WaveFunCombine')!=-1):
		for i in range(len(WF.WFs)):
			WF1.Sx[i]=WF.Sx[i]*M[0,0]+WF.Sy[i]*M[1,0]+WF.Sz[i]*M[2,0]
			WF1.Sy[i]=WF.Sx[i]*M[0,1]+WF.Sy[i]*M[1,1]+WF.Sz[i]*M[2,1]
			WF1.Sz[i]=WF.Sx[i]*M[0,2]+WF.Sy[i]*M[1,2]+WF.Sz[i]*M[2,2]
			WF1.Sp[i]=WF1.Sx[i]+1j*WF1.Sy[i]
			WF1.Sm[i]=WF1.Sx[i]-1j*WF1.Sy[i]
	return WF1

def RotationWFZ(Q,F,WF):
	cosQ=math.cos(Q*math.pi/180.0)
	sinQ=math.sin(Q*math.pi/180.0)
	cosF = math.cos(F*math.pi/180.0)
	sinF = math.sin(F*math.pi/180.0)
	WF1=WF.copy()
	WF1.Sz=WF.Sz*cosQ-WF.Sx*sinQ
	WF1.Sx=WF.Sz*(cosF*sinQ)+WF.Sx*(cosF*cosQ)-WF.Sy*sinF
	WF1.Sy=WF.Sz*(sinF*sinQ)+WF.Sx*(sinF*cosQ)+WF.Sy*cosF
	WF1.Sp=WF1.Sx+1j*WF1.Sy
	WF1.Sm=WF1.Sx-1j*WF1.Sy
	return WF1

def ApllyHamiltonian(oldWF, hamilton):
	newWF=oldWF.copy()
	E=MakeHamFromEnergy(oldWF.Energy)+hamilton
	(E1, M) = np.linalg.eig(E)
	for i in range(newWF.N):
		A=M[i,:]
		j=np.argmax(M[i,:])
		if(j!=i):
			t=M[:,j].copy()
			M[:,j]=M[:,i].copy()
			M[:, i]=t.copy()
			t=E1[j]
			E1[j]=E1[i]
			E1[i]=t
	newWF.Energy=E1
	newWF.WF_matrix=M*oldWF.WF_matrix
	for i in range(len(oldWF.WFs)):
		newWF.Sz[i]=M.transpose()*oldWF.Sz[i]*M
		newWF.Sp[i]=M.transpose()*oldWF.Sp[i]*M
		newWF.Sm[i]=M.transpose()*oldWF.Sm[i]*M
		newWF.Sx[i] = (newWF.Sp[i] + newWF.Sm[i]) * 0.5
		newWF.Sy[i] = (newWF.Sp[i] - newWF.Sm[i]) * (0.5 * (1j))
	for i in range(oldWF.N):
		newWF.WFpr[i]=''
		for j in range(oldWF.N):
			s1=wavefun_print.printcomplex(M[j,i])
			if(s1=='1.000'):
				newWF.WFpr[i]+='|'+oldWF.WFprBase[j]+'>+\n'
				break
			if(s1!='0'):
				newWF.WFpr[i]+=s1+'*|'+oldWF.WFprBase[j]+'>+\n'
		newWF.WFpr[i]=newWF.WFpr[i][0:-2]
	return newWF

def ExtructSubMatrix(WF, number_rows):
	newWF=WF.copy()
	newWF.N=len(number_rows)
	newWF.J= newWF.N/2.0-1.0
	newWF.JJ1=newWF.J*(newWF.J+1)
	newWF.E = np.identity(newWF.N, dtype='complex')
	newWF.WF_matrix = newWF.E
	newWF.Energy=np.zeros(newWF.N, dtype='complex')
	newWF.Sz = [np.matrix(np.zeros((newWF.N, newWF.N)), complex) for k in range(len(newWF.WFs))]
	newWF.Sp = [np.matrix(np.zeros((newWF.N, newWF.N)), complex) for k in range(len(newWF.WFs))]
	newWF.Sm = [np.matrix(np.zeros((newWF.N, newWF.N)), complex) for k in range(len(newWF.WFs))]
	newWF.Sx = [np.matrix(np.zeros((newWF.N, newWF.N)), complex) for k in range(len(newWF.WFs))]
	newWF.Sy = [np.matrix(np.zeros((newWF.N, newWF.N)), complex) for k in range(len(newWF.WFs))]
	for i in range(newWF.N):
		newWF.Energy[i]=WF.Energy[number_rows[i]]
		for j in range(newWF.N):
			newWF.WF_matrix[i,j]=WF.WF_matrix[number_rows[i],number_rows[j]]
			for k in range(len(newWF.WFs)):
				newWF.Sz[k][i,j]=WF.Sz[k][number_rows[i],number_rows[j]]
				newWF.Sp[k][i,j]=WF.Sp[k][number_rows[i],number_rows[j]]
				newWF.Sm[k][i,j]=WF.Sm[k][number_rows[i],number_rows[j]]
				newWF.Sx[k][i,j]=WF.Sx[k][number_rows[i],number_rows[j]]
				newWF.Sy[k][i,j]=WF.Sy[k][number_rows[i],number_rows[j]]
	return newWF



def CalcB(x1,x2,y1,y2,z1,z2,hv):
	if((y1+y2-z1-z2)!=0):
		return (hv*(x2-x1)-y1*x2+y2*x1+z1*x2-z2*x1)/(-y1+y2+z1-z2)
	else:
		return None

def IsCalcBX(x:float,x1:float,x2:float):
	if(x1>x2):
		return (x<x1 and x>x2)
	else:
		return (x<x2 and x>x1)

def FindLine(NSp,Sp,WF,hv):
	Br1=[]
	m1=[]
	m2=[]
	for i in range(1,NSp):
		for j in range(WF.N):
			for k in range(WF.N):
				if(j==k):
					continue
				B1=CalcB(Sp['B'][i-1],Sp['B'][i], Sp['E'][i-1][j], Sp['E'][i][j], Sp['E'][i-1][k], Sp['E'][i][k],hv)
				if(B1!=None and IsCalcBX(B1, Sp['B'][i-1],Sp['B'][i])):
					Br1.append(B1)
					m1.append(Sp['E'][i-1][j])
					m2.append(Sp['E'][i-1][k])
	NBres=Br1.__len__()
	return NBres, Br1, m1, m2

