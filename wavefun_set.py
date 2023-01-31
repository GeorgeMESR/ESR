import numpy as np
import scipy
import math
from scipy.spatial.transform import Rotation as R
import wavefun_print

def MakeMSz(res, WF):
	for i in range(WF.N):
		res[i, i] = WF.WF[i]

def MakeMSPlus(res, WF):
	for i in range(WF.N - 1):
		res[i, i + 1] = np.sqrt(WF.JJ1 - WF.WF[i] * WF.WF[i + 1])

def MakeMSMinus(res, WF):
	for i in range(WF.N - 1):
		res[i + 1, i] = np.sqrt(WF.JJ1 - WF.WF[i] * WF.WF[i + 1])

class WaveFunBasic:
	def __init__(self):
		self.J = int(0)
		self.JJ1 = int(0)
		self.WF=[]
		self.Energy=[]
		self.WFpr=[]
		self.WFprBase = []
		self.N=int(2)
		self.E =  np.identity(self.N, dtype='complex')
		self.WF_matrix = self.E
		self.Sz = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sp = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sm = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sx = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sy = np.matrix(np.zeros((self.N, self.N)), complex)
	def __init__(self,J):
		self.J = J
		self.JJ1 = J*(J+1)
		self.N = int(2.0 * J + 1.0)
		self.Energy = np.zeros(self.N)
		self.WF = []
		self.WFpr = []
		self.E =  np.eye(self.N, dtype='complex')
		self.WF_matrix = self.E
		self.Sz = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sp = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sm = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sx = np.matrix(np.zeros((self.N, self.N)), complex)
		self.Sy = np.matrix(np.zeros((self.N, self.N)), complex)
		self.MakeWaveFun(J)
	def copy(self):
		retWF=WaveFunBasic(self.J)
		retWF.J=self.J
		retWF.JJ1=self.JJ1
		retWF.N=self.N
		retWF.Energy=self.Energy
		retWF.WF = [w for w in self.WF]
		retWF.WFpr = [w for w in self.WFpr]
		retWF.WFprBase = [ w for w in self.WFprBase]
		retWF.WF_matrix=self.WF_matrix
		retWF.E=self.E
		retWF.Sz=self.Sz
		retWF.Sp=self.Sp
		retWF.Sm=self.Sm
		retWF.Sx=self.Sx
		retWF.Sy=self.Sy
		return retWF

	def MakeWaveFun(self, J):
		self.J = J
		self.JJ1 = J * (J + 1)
		self.N = int(2.0 * J + 1.0)
		self.Energy=np.zeros(self.N)
		mJ = - self.J
		for i in range( self.N):
			self.WF.append(mJ)
			self.WFpr.append(wavefun_print.PrintSpin(mJ))
			mJ+=1
		self.WFprBase = [ m for m in self.WFpr]
		MakeMSz(self.Sz,self)
		MakeMSPlus(self.Sp,self)
		MakeMSMinus(self.Sm,self)
		self.Sx=(self.Sp+self.Sm)*0.5
		self.Sy=(self.Sp-self.Sm)*(0.5*(1j))

	def sprintWaveFun(self):
		retstr = ''
		for i in range(self.N):
			retstr+=self.WFpr[i]
			if (i < self.N - 1):
				retstr+=(',')
		return retstr

class WaveFunCombine:
	def __init__(self, *listWF):
		self.WFs=[f for f in listWF]
		if(len(self.WFs)==0):
			self.N=0
			return
		self.J=[f.J for f in self.WFs]
		self.N=1
		for w in self.WFs:
			self.N*=w.N
		self.JJ1 = [f.JJ1 for f in self.WFs]
		self.Energy=np.zeros(self.N)
		self.WF=[[] for i in range(self.N)]
		self.WFpr=['' for i in range(self.N)]
		nrepeat1 = int(self.N / self.WFs[0].N)
		nrepeat2 =1
		for i in range( len(self.WFs)):
			m=0
			for j in range(nrepeat2):
				for k in range(self.WFs[i].N):
					for l in range(nrepeat1):
						self.WF[m].append(self.WFs[i].WF[k])
						self.WFpr[m]+=self.WFs[i].WFpr[k]+' '
						m+=1
			if(i<len(self.WFs)-1):
				nrepeat1=int(nrepeat1/self.WFs[i+1].N)
				nrepeat2=int(self.N/self.WFs[i+1].N/nrepeat1)

		for s in range(len(self.WFpr)):
			self.WFpr[s]=self.WFpr[s].strip()
		self.WFprBase = [ m for m in self.WFpr]
		self.E=[]
		for i in range( len(self.WFs)):
			E=np.ones((1,1),complex)
			for j in range(len(self.WFs)):
				if(i!=j):
					E=np.kron(E,self.WFs[j].E)
				else:
					E = np.kron(E, self.WFs[j].E)
			self.E.append(E)
		self.WF_matrix = np.eye(self.N, dtype='complex')
		self.Sz = []
		for i in range( len(self.WFs)):
			E=np.ones((1,1),complex)
			for j in range(len(self.WFs)):
				if(i!=j):
					E=np.kron(E,self.WFs[j].E)
				else:
					E = np.kron(E, self.WFs[j].Sz)
			self.Sz.append(E)
		self.Sp = []
		for i in range( len(self.WFs)):
			E=np.ones((1,1),complex)
			for j in range(len(self.WFs)):
				if(i!=j):
					E=np.kron(E,self.WFs[j].E)
				else:
					E = np.kron(E, self.WFs[j].Sp)
			self.Sp.append(E)
		self.Sm = []
		for i in range( len(self.WFs)):
			E=np.ones((1,1),complex)
			for j in range(len(self.WFs)):
				if(i!=j):
					E=np.kron(E,self.WFs[j].E)
				else:
					E = np.kron(E, self.WFs[j].Sm)
			self.Sm.append(E)
		self.Sx = []
		self.Sy = []
		for i in range(len(self.WFs)):
			self.Sx.append((self.Sp[i]+self.Sm[i])*0.5)
			self.Sy.append((self.Sp[i]-self.Sm[i])*(0.5*(1j)))
	def copy(self):
		retWF=WaveFunCombine()
		retWF.WFs=[f for f in self.WFs]
		retWF.J=[f.J for f in self.WFs]
		retWF.JJ1=[f.JJ1 for f in self.WFs]
		retWF.N=self.N
		retWF.Energy=self.Energy
		retWF.WF = [w for w in self.WF]
		retWF.WFpr = [w for w in self.WFpr]
		retWF.WFprBase = [ m for m in self.WFprBase]
		retWF.E=self.E
		retWF.Sz=self.Sz
		retWF.Sp=self.Sp
		retWF.Sm=self.Sm
		retWF.Sx=self.Sx
		retWF.Sy=self.Sy
		return retWF

