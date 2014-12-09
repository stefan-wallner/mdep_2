import numpy as np
from sys import argv

def ampl_to_SDM(ampl):
	"""
	Converts a list of amplitudes to th corresponding spin-density matrix elemtns with one anchor wave.
	Assumes a real parameter list
	"""
	SDM = [ampl[0]**2+ampl[1]**2]
	anc = 1.*ampl[0]+ 1.j+ ampl[1]	
	if not len(ampl)%2 == 0:
		print "Error, odd number of input parameters"
		raise IndexError
	for i in range(1,len(ampl)/2):
		amp = 1.*ampl[2*i] - 1.j*ampl[2*i+1] # the -1.j since complex conjugation is needed
		inter = anc*amp
		SDM.append(inter.real)
		SDM.append(inter.imag)
	return SDM

def get_coma_MC(mean,coma,trafo,n_SAMPLE = 10000):	
	"""
	Transforms the covariance matrix according to trafo
	"""
	new_mean = trafo(mean)
	new_coma = []
	new_dim = len(new_mean)
	outfile = open('outfile.txt','w')
	for i in range(new_dim):
		new_coma.append([0. for i in range(new_dim)])
	for i in range(n_SAMPLE):
		sample = np.random.multivariate_normal(mean,coma)
		new_sample = trafo(sample)
		for i in range(new_dim):
			new_sample[i]-=new_mean[i]
		outfile.write(str(sample[0])+'\n')
		for i in range(new_dim):
			for j in range(new_dim):
				new_coma[i][j]+=new_sample[i]*new_sample[j]
	for i in range(new_dim):
		for j in range(new_dim):
			new_coma[i][j]/=n_SAMPLE
	outfile.close()
	return new_coma

def get_coma(mean,coma):
	"""
	Transforms the coma analytically from amplitudes to SDM entries
	"""
	jac = []
	old_dim = len(mean)
	amp_R = mean[0]
	amp_I = mean[1]
	jac_line_anc = [0. for i in range(old_dim)]
	jac_line_anc[0] = 2*amp_R
	jac_line_anc[1] = 2*amp_I
	jac.append(jac_line_anc)
	for i in range(1,old_dim/2):
		jac_line1 = [0. for j in range(old_dim)]
		jac_line2 = [0. for j in range(old_dim)]
		jac_line1[0] = mean[2*i  ]
		jac_line1[1] = mean[2*i+1]
		jac_line2[0] =-mean[2*i+1]
		jac_line2[1] = mean[2*i  ]

		jac_line1[2*i  ] = amp_R
		jac_line1[2*i+1] = amp_I
		jac_line2[2*i  ] = amp_I
		jac_line2[2*i+1] =-amp_R

		jac.append(jac_line1)
		jac.append(jac_line2)

	new_coma = []
	for i in range(old_dim-1):
		new_coma.append([0. for i in range(old_dim-1)])
	for i in range(old_dim-1):
		for j in range(old_dim-1):
			for k in range(old_dim):
				for l in range(old_dim):
					new_coma[i][j]+=jac[i][k]*jac[j][l]*coma[k][l]
	return new_coma

def dtv(var, n=1000):
	ppp = 0.
	for i in range(n):
		ppp+=(np.random.normal(0.,var))**2.
	ppp/=n
	print var,ppp	

def square(xx):
	return [x**2 for x in xx]


if __name__ == "__main__":
	datFile = argv[1]
	comaFile = argv[2]
	data = []
	with open(datFile,'r') as ininin:
		for line in ininin.readlines():
			linevals = [float(val) for val in line.split()]
			data.append(linevals)
	nPoints = len(data[0])	
	coma = []
	with open(comaFile,'r') as ininin:
		for line in ininin.readlines():
			comaLine = [float(val) for val in line.split()]
			oneComa = []
			for i in range(nPoints):
				lin = []
				for j in range(nPoints):
					lin.append(comaLine[i*nPoints+j])
				oneComa.append(lin)
			coma.append(oneComa)
	with open(datFile+'_SDM','w') as outData:
		with open(comaFile+'_SDM','w') as outComa:
			for i in range(len(data)):
				points = ampl_to_SDM(data[i])
				for point in points:
					outData.write(str(point)+'  ')
				outData.write('\n')
				newComa = get_coma(data[i],coma[i])
				for lin in newComa:
					for val in lin:
						outComa.write(str(val)+'  ')
				outComa.write('\n')







