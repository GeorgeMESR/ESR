import wavefun
import numpy as np


def MakeSz(hammat):
    hammat.M=np.zeros((hammat.wavefunctions.N,hammat.wavefunctions.N))
    for i in range(hammat.wavefunctions.N):
        hammat.M[i][i]=hammat.wavefunctions.WFarray[i][0]

def MakeSPlus(hammat):
    hammat.M=np.zeros((hammat.wavefunctions.N,hammat.wavefunctions.N))
    S=hammat.wavefunctions.J[0]
    for i in range(hammat.wavefunctions.N-1):
        hammat.M[i][i+1]=np.sqrt(S*(S+1)-hammat.wavefunctions.WFarray[i+1][0]*hammat.wavefunctions.WFarray[i][0])

def MakeSMinus(hammat):
    hammat.M=np.zeros((hammat.wavefunctions.N,hammat.wavefunctions.N))
    S=hammat.wavefunctions.J[0]
    for i in range(hammat.wavefunctions.N-1):
        hammat.M[i+1][i]=np.sqrt(S*(S+1)-hammat.wavefunctions.WFarray[i+1][0]*hammat.wavefunctions.WFarray[i][0])

def HamOk0(N,M):
    N2=int(N/2)
    for i in range(N2):
        M[N-1-i][N-1-i]=M[i][i]

def HamOkq(N,M):
    for j in range(N):
        for i in range(N-1-j):
            M[N-1-j][N-1-i]=M[i][j]
    for i in range(N):
        for j in range(i+1,N):
            M[j][i]=M[i][j]

def Ham_O_02(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if(hammat.wavefunctions.J==1):
        hammat.M[0][0]=1
        hammat.M[1][1]=-2
    if(hammat.wavefunctions.J==3/2):
        hammat.M[0][0]=3
        hammat.M[1][1]=-3
    if(hammat.wavefunctions.J==5/2):
        hammat.M[0][0]=10
        hammat.M[1][1]=-2
        hammat.M[2][2]=-8
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][0]=21
        hammat.M[1][1]=3
        hammat.M[2][2]=-9
        hammat.M[3][3]=-15
    HamOk0(hammat.wavefunctions.N,hammat.M)



def Ham_O_04(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if(hammat.wavefunctions.J==5/2):
        hammat.M[0][0]=60.0
        hammat.M[1][1]=-180.0
        hammat.M[2][2]=120.0
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][0]=420
        hammat.M[1][1]=-780
        hammat.M[2][2]=-180
        hammat.M[3][3]=540
    HamOk0(hammat.wavefunctions.N,hammat.M)


def Ham_O_06(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][0]=1260
        hammat.M[1][1]=-6300
        hammat.M[2][2]=11340
        hammat.M[3][3]=-6300
    HamOk0(hammat.wavefunctions.N,hammat.M)

def Ham_O_22(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if(hammat.wavefunctions.J==1):
        hammat.M[0][2]=1
    if(hammat.wavefunctions.J==3/2):
        hammat.M[0][2]=1.7320508075688772935274463415059
    if(hammat.wavefunctions.J==5/2):
        hammat.M[0][2]=3.1622776601683793319988935444327
        hammat.M[1][3]=4.2426406871192851464050661726291
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][2]=21
        hammat.M[1][3]=3
        hammat.M[2][4]=-9
    HamOkq(hammat.wavefunctions.N,hammat.M)

def Ham_O_24(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if(hammat.wavefunctions.J==5/2):
        hammat.M[0][2]=28.460498941515413987990041899894
        hammat.M[1][3]=-21.213203435596425732025330863145
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][2]=137.47727084867520019764141581184
        hammat.M[1][3]=13.416407864998738178455042012388
        hammat.M[2][4]=-92.951600308978005244302369594778
    HamOkq(hammat.wavefunctions.N,hammat.M)

def Ham_O_26(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][2]=549.90908339470080079056566324736
        hammat.M[1][3]=-1126.9782606598940069902235290406
        hammat.M[2][4]=650.66120216284603671011658716344
    HamOkq(hammat.wavefunctions.N,hammat.M)

def Ham_O_34(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if(hammat.wavefunctions.J==5/2):
        hammat.M[0][3]=9.4868329805051379959966806332982
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][3]=35.49647869859769625540396974937
        hammat.M[1][4]=26.832815729997476356910084024775
    HamOkq(hammat.wavefunctions.N,hammat.M)

def Ham_O_36(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][3]=425.95774438317235506484763699244
        hammat.M[1][4]=-563.48913032994700349511176452028
    HamOkq(hammat.wavefunctions.N,hammat.M)


def Ham_O_44(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if(hammat.wavefunctions.J==5/2):
        hammat.M[0][4]=26.832815729997476356910084024775
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][4]=70.992957397195392510807939498739
        hammat.M[1][5]=103.92304845413263761164678049035
    HamOkq(hammat.wavefunctions.N,hammat.M)

def Ham_O_46(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][4]=1064.8943609579308876621190924811
        hammat.M[1][5]=-727.46133917892846328152746343247
    HamOkq(hammat.wavefunctions.N,hammat.M)

def Ham_O_66(hammat):
    hammat.M = np.zeros((hammat.wavefunctions.N, hammat.wavefunctions.N))
    if (hammat.wavefunctions.J == 7/2):
        hammat.M[0][5]=952.47047198325261258058167131013
    HamOkq(hammat.wavefunctions.N,hammat.M)

class MakeZeeman:
    def __init__(self, J):
        self.Sz=wavefun.HamiltonianMatrix(J)
        MakeSz(self.Sz)
        self.Sp=wavefun.HamiltonianMatrix(J)
        MakeSPlus(self.Sp)
        self.Sm=wavefun.HamiltonianMatrix(J)
        MakeSMinus(self.Sm)
        self.Sx=wavefun.HamiltonianMatrix(J)
        self.Sx.M=0.5*(self.Sp.M+self.Sm.M)
        self.Sy=wavefun.HamiltonianMatrix(J)
        self.Sy.M=np.complex(0,0.5)*(self.Sp.M-self.Sm.M)

class MakeFineStructure:
    def __init__(self, J):
        if (J[0] == 1 or J[0]==3/2):
            self.O_02=wavefun.HamiltonianMatrix(J)
            Ham_O_02(self.O_02)
            self.O_22 = wavefun.HamiltonianMatrix(J)
            Ham_O_22(self.O_22)
        if (J[0] == 2 or J[0]==5/2):
            self.O_02=wavefun.HamiltonianMatrix(J)
            Ham_O_02(self.O_02)
            self.O_22 = wavefun.HamiltonianMatrix(J)
            Ham_O_22(self.O_22)
            self.O_04=wavefun.HamiltonianMatrix(J)
            Ham_O_04(self.O_04)
            self.O_24 = wavefun.HamiltonianMatrix(J)
            Ham_O_24(self.O_24)
            self.O_34 = wavefun.HamiltonianMatrix(J)
            Ham_O_34(self.O_34)
            self.O_44 = wavefun.HamiltonianMatrix(J)
            Ham_O_44(self.O_44)
        if (J[0] == 3 or J[0]==7/2):
            self.O_02=wavefun.HamiltonianMatrix(J)
            Ham_O_02(self.O_02)
            self.O_22 = wavefun.HamiltonianMatrix(J)
            Ham_O_22(self.O_22)
            self.O_04=wavefun.HamiltonianMatrix(J)
            Ham_O_04(self.O_04)
            self.O_24 = wavefun.HamiltonianMatrix(J)
            Ham_O_24(self.O_24)
            self.O_34 = wavefun.HamiltonianMatrix(J)
            Ham_O_34(self.O_34)
            self.O_44 = wavefun.HamiltonianMatrix(J)
            Ham_O_44(self.O_44)
            self.O_06 = wavefun.HamiltonianMatrix(J)
            Ham_O_06(self.O_06)
            self.O_26 = wavefun.HamiltonianMatrix(J)
            Ham_O_26(self.O_26)
            self.O_36 = wavefun.HamiltonianMatrix(J)
            Ham_O_36(self.O_36)
            self.O_46 = wavefun.HamiltonianMatrix(J)
            Ham_O_46(self.O_46)
            self.O_66 = wavefun.HamiltonianMatrix(J)
            Ham_O_66(self.O_66)
