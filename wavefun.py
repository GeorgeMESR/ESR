import numpy as np

def PrintSpin(mJ):
    mJi = int(mJ * 2)
    if(mJ<0):
        if (mJi % 2 != 0):
            retstr = (str(mJi) + '½')
        else:
            retstr = (str(mJi / 2))
    else:
        if (mJi % 2 != 0):
            retstr = ('+'+str(mJi) + '½')
        else:
            retstr = ('+'+str(mJi / 2))
    return retstr


class WaveFun:
    def __init__(self):
        self.J = np.zeros(0)
        self.WFarray=[]
    def Up_mJ(self, mJ, k):
        if (mJ[k] < self.J[k]):
            mJ[k] += 1.0
        else:
            mJ[k] = -self.J[k]
            if (k+1 < self.J.size):
                self.Up_mJ(mJ, k + 1)
    def MakeWaveFun(self, Jarray):
        self.WFarray=[]
        self.J = np.array(Jarray)
        self.N = 1
        for j in  self.J:
            nJlocal = 2 * j + 1
            self.N *= int(nJlocal)
        mJ = - self.J
        for i in range( self.N):
            self.WFarray.append(np.array(mJ))
            self.Up_mJ(mJ, 0)
    def sprintWaveFun(self, k):
        retstr = ''
        nWF = len(self.WFarray[k])
        for i in range(nWF):
            retstr+=PrintSpin(self.WFarray[k][i])
            if (i < nWF - 1):
                retstr+=(',')
        return retstr

def AppendWaveFunctions(wf1, wf2):
    retWF=WaveFun()
    Jarray=wf1.J
    np.append(Jarray,ham2.wavefunctions.J)
    retWF.MakeWaveFun(Jarray)
    return WF

def CompareWaveFunctionsLeft(wf1, wf2):
    n2=len(wf2.WFarray)
    wf1a=wf1.WFarray[0:n2]
    return np.array_equal(wf1a,wf2.WFarray)

def CompareWaveFunctionsRight(wf1, wf2):
    n1=len(wf1.WFarray)
    n2=len(wf2.WFarray)
    wf1a=wf1.WFarray[n1-n2:n1]
    return np.array_equal(wf1a,wf2.WFarray)

class HamiltonianMatrix:
    def __init__(self,Jarray):
        self.wavefunctions = WaveFun()
        self.wavefunctions.MakeWaveFun(Jarray)
        self.M=np.zeros((self.wavefunctions.N,self.wavefunctions.N))
    def printHM(self):
        retstr=''
        s=self.wavefunctions.sprintWaveFun(0)
        nx=len(s)+2
        if(nx<7):
            nx=7
        s=''
        s=s.rjust(nx,' ')
        s += ' | '
        for i in range(self.wavefunctions.N):
            s+= self.wavefunctions.sprintWaveFun(i).rjust(nx,' ')
            s += ' | '
        retstr+=s+'\n'
        for i in range(self.wavefunctions.N):
            s=self.wavefunctions.sprintWaveFun(i).rjust(nx,' ')
            s += ' | '
            for j in range(self.wavefunctions.N):
                s1='%3.3f' % self.M[i][j]
                s+=str(s1).rjust(nx,' ')
                s+= ' | '
            retstr+=s+'\n'
        return retstr



def AppendHamiltonianMatrix(ham1,ham2):
    Jarray=np.hstack((ham1.wavefunctions.J,ham2.wavefunctions.J))
    retHam=HamiltonianMatrix(Jarray)
    for i in range(ham2.wavefunctions.N):
        for j in range(ham1.wavefunctions.N):
            for k in range(ham1.wavefunctions.N):
                retHam.M[j+i*ham1.wavefunctions.N][k+i*ham1.wavefunctions.N]+=ham1.M[j][k]
    for i in range(ham1.wavefunctions.N):
        for j in range(ham2.wavefunctions.N):
            for k in range(ham2.wavefunctions.N):
                retHam.M[j*ham1.wavefunctions.N+i][k*ham1.wavefunctions.N+i]+=ham2.M[j][k]
    return retHam

