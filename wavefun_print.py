import numpy as np
import scipy
import math
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt

#Function print spin number: half-integer spin printed as N*1/2 ;
def PrintSpin(mJ):
	minus=0
	if(mJ<0):
		mJ=-mJ
		minus=1
	mJi = int(mJ * 2)
	if (mJi % 2 != 0):
		if(mJi!=1):
			retstr = (str(mJi) + '½')
		else:
			retstr = '½'
	else:
		retstr = (str(mJi / 2))
	if(minus):
		retstr='-'+retstr
	else:
		retstr = '+' + retstr
	return retstr

def printcomplex(val:complex):
	s1=''
	if(abs(val)<0.01):
		return '0'
	else:
		if (abs(val.imag) < 0.01):
			s1 = '%3.3f' % val.real
		else:
			if (abs(val.real) < 0.01):
				if (val.imag > 0):
					s1 = 'j%3.3f' % val.imag
				else:
					s1 = '-j%3.3f' % abs(val.imag)
			else:
				if (val.real > 0):
					s1 = '%1.1f' % val.real
				else:
					s1 = '-%1.1f' % abs(val.real)
				if (val.imag > 0):
					s1 += '+j%1.1f' % val.imag
				else:
					s1 += '-j%1.1f' % abs(val.imag)
	return s1

def printcomplexg(val:complex):
	s1=''
	if(abs(val)<0.00001):
		return '0'
	else:
		if (abs(val.imag) < 0.00001):
			s1 = '%g' % val.real
		else:
			if (abs(val.real) < 0.00001):
				if (val.imag > 0):
					s1 = 'j%g' % val.imag
				else:
					s1 = '-j%g' % abs(val.imag)
			else:
				if (val.real > 0):
					s1 = '%g' % val.real
				else:
					s1 = '-%g' % abs(val.real)
				if (val.imag > 0):
					s1 += '+j%g' % val.imag
				else:
					s1 += '-j%g' % abs(val.imag)
	return s1


def printHM(WF, M : np.matrix):
	retstr=''
	s=WF.WFpr[0]
	nx=len(s)+2
	if(nx<12):
		nx=12
	s=''
	s=s.rjust(nx,' ')
	s += ' | '
	for i in range(WF.N):
		s+= WF.WFpr[i].rjust(nx,' ')
		s += ' | '
	retstr+=s+'\n'
	for i in range(WF.N):
		s=WF.WFpr[i].rjust(nx,' ')
		s += ' | '
		for j in range(WF.N):
			s1=printcomplex(M[i,j])
			s+=str(s1).rjust(nx,' ')
			s+= ' | '
		retstr+=s+'\n'
	return retstr

def printHM_HTML(filename:str, WF, *aM):
	mlist=[]
	madd=['', 0]
	for m in aM:
		if(str(type(m)).find('str')!=-1):
			madd[0]=m
		else:
			madd[1]=m
			mlist.append(madd)
			madd=['', 0]
	with open(filename+'.html','w') as f:
		f.write('<html> <meta http-equiv="Content-Type" content="text/html; charset=utf-8">          <head><style> body, b, center, nobr, br,font {font-family:times new roman,sans serif;font-size:10pt;text-indent: 35px;}           .EditTable {text-align: center; text-indent: 0px;padding: 2px;}            H1 {page-break-before: always;}                    span.right {                  -webkit-transform: rotate(-90deg); /* Chrome y Safari */                  -moz-transform: rotate(-90deg); /* Firefox */                  filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=3); /* Internet Explorer */                  -o-transform: rotate(-90deg); /* Opera */                  display: inline-block; padding-left: 8px; vertical-align: text-top; text-indent: 0px;              }  		@page  		{  			mso-page-orientation: landscape;   			size:29.7cm 21cm;  			margin:2cm 2cm 2cm 2cm;  		}			           @page WordSection1 {}  						            div.WordSection1 {page:WordSection1;}           </style>          </head>          <body> <div class="WordSection1">\n')
		for mo in mlist:
			f.write('<BR>'+mo[0]+'\n')
			M=mo[1]
			f.write('<table cellspacing="-1" border="1" width="100%;">\n')
			f.write('<tr><td>&nbsp;</td>')
			for i in range(WF.N):
				f.write('<td>' + WF.WFpr[i].replace('\n','<BR>')+ '</td>')
			f.write('</tr>\n')
			for i in range(WF.N):
				f.write('<tr>')
				f.write('<td>' + WF.WFpr[i].replace('\n','<BR>') + '</td>')
				for j in range(WF.N):
					s1=printcomplexg(M[i,j])
					if(i==j):
						f.write('<td style="background-color:#A0A0FF">' + s1 + '</td>')
					elif(abs(M[i,j])>0.01):
						f.write('<td style="background-color:#FFA0A0">' + s1 + '</td>')
					else:
						f.write('<td style="background-color:#FFFFFF">' + s1 + '</td>')
				f.write('</tr>\n')
			f.write('</table>\n')
		f.write('</div></body></html>')
		f.close()

def plotEnergy(WF):
	B0=np.array([0,1])
	for i in range(len(WF.Energy)):
		E=np.array([np.real(WF.Energy[i])/1000., np.real(WF.Energy[i])/1000.])
		plt.plot(B0,E)
		plt.text(1,np.real(WF.Energy[i])/1000.,WF.WFpr[i].replace('\n',''))
	k=0
	for i in range(len(WF.Energy)):
		for j in range(i+1,len(WF.Energy)):
			if(np.abs(WF.Sx[0][i,j])>0.1):
				plt.vlines(k/10.+0.1,WF.Energy[i]/1000.,WF.Energy[j]/1000.,colors=[1.0,1.0,0.5])
				plt.text(k/10.+0.1,(WF.Energy[i] + WF.Energy[j]) / 2000., printcomplex((WF.Energy[i] - WF.Energy[j]) / 1000.))
				k+=1

		pass
	plt.show()