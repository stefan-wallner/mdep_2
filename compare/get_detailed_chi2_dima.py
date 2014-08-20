from sys import exit

def get_detailed_chi2_dima(inFile_name ='./chi2.dat', nWaves = 13, nTbin = 11):
	inFile = open(inFile_name,'r')
	
	deltas={}
	delchi={}
	
	nPoints = 2*nWaves -1
	lines = inFile.readlines()
	
	for i in range(len(lines)):
		chunks = lines[i].split()
		if chunks[0] == 'mass':
			if len(chunks) == 8: #Deltas
				chunks2 = lines[i+1].split()
				complex_bin = int(chunks[3])+1j*int(chunks[7])
#				if int(chunks[3]) in [55,75,80]:
#					print "bin ",chunks[3],"has mass:",chunks[1]
				if not deltas.has_key(complex_bin):
					deltas[complex_bin]={}
				deltas[complex_bin][int(chunks[5])] = [float(chunks2[2]),float(chunks2[5]),float(chunks2[7])]
			if len(chunks) == 10: #Delchis
				chunks2 = lines[i+2].split()
				chunks3 = lines[i+3].split()
				complex_bin = int(chunks[3]) +1j*int(chunks[9])
				if not delchi.has_key(complex_bin):
					delchi[complex_bin]={}
				complex_point = int(chunks[5]) + 1j*int(chunks[7])
				delchi[complex_bin][complex_point] = [float(chunks2[5]),float(chunks3[4])]
	inFile.close()
	return deltas,delchi

def cross_chi2(data,theory,coma):
	chi2=0.
	pointcount = 0
	for complex_bin in coma.iterkeys():
			for complex_point in coma[complex_bin].iterkeys():
				i = complex_point.real
				j = complex_point.imag
				try:
					data[complex_bin][i]
					data[complex_bin][j]
					theory[complex_bin][i]
					theory[complex_bin][j]
				except KeyError:
					continue
				contr = (theory[complex_bin][i][1]-data[complex_bin][i][0])*(theory[complex_bin][j][1]-data[complex_bin][j][0])*coma[complex_bin][complex_point][0]
				if i==j:
					chi2+=contr
				else:
					chi2+=2*contr
				if not contr == 0.:
					pointcount+=1
	return chi2,pointcount

def pointdiff_cc(a,b):
	a_not_b=[]
	for complex_bin in a.iterkeys():
		for complex_point in a[complex_bin].iterkeys():
			try:
				b[complex_bin][complex_point]
			except KeyError:
				a_not_b.append([complex_bin,complex_point])
	b_not_a=[]
	for complex_bin in b.iterkeys():
		for complex_point in b[complex_bin].iterkeys():
			try:
				a[complex_bin][complex_point]
			except KeyError:
				b_not_a.append([complex_bin,complex_point])
	return a_not_b,b_not_a

def list_chi2(points, dat_the, coma):
	chi2=0.
	for point in points:	
		complex_bin = point[0]
		complex_point=point[1]
		i = complex_point.real
		j = complex_point.imag
		contr = (dat_the[complex_bin][i][0] - dat_the[complex_bin][i][1])*(dat_the[complex_bin][j][0] - dat_the[complex_bin][j][1])*coma[complex_bin][complex_point][0]
		if i==j:
			chi2+=contr
		else:
			chi2+=2*contr
	return chi2

def get_single_bins(points):
	first=[]
	second=[]
	for point in points:
		if not point[0] in first:
			first.append(point[0])
		if not point[1] in second:
			second.append(point[1])
	wav=[]
	for sec in second:
		if not sec.real in wav:
			wav.append(sec.real)
		if not sec.imag in wav:
			wav.append(sec.imag)
	wav.sort()
	return first,second,wav

deltas_d,delchi_d=get_detailed_chi2_dima()
deltas_f,delchi_f=get_detailed_chi2_dima('./chi2_f.dat')
BUILD_FILES=False
if BUILD_FILES:
	for tbin in range(1,12): # 11 t'-bins
		datafile = open('./data_coma/data_'+str(tbin)+".dat",'w')
		comafile = open('./data_coma/coma_'+str(tbin)+".dat",'w')
		for mbin in range(0,100): # 100 m-bins
			complex_bin = mbin+1j*tbin
			data = [0. for i in range(25)]
			coma = [[0. for i in range(25)] for i in range(25)]
			for i in range(1,26):
				try:
					data[i-1] = deltas_d[complex_bin][i][0]
				except KeyError:
					pass
				for j in range(1,26):
					complex_point=i+1.j*j
					try:
						coma[i-1][j-1] = delchi_d[complex_bin][complex_point][0]
					except KeyError:
						pass
			for i in range(25):
				for j in range(25):
					if coma[i][j] == 0.: 
						coma[i][j] = coma[j][i]
			for i in range(25):
				for j in range(25):
					comafile.write(str(coma[i][j])+'   ')
			comafile.write('\n')
			datstrings = [str(dat) for dat in data]
			outstring = '   '.join(datstrings)
			outstring+='\n'
			datafile.write(outstring)
		comafile.close()
		datafile.close()

CHECK_DATA_INPUT = False
if CHECK_DATA_INPUT:
	ndiff =0
	for complex_bin in delchi_f.iterkeys():
		for complex_point in delchi_f[complex_bin].iterkeys():
			try:
				delchi_d[complex_bin][complex_point]
			except KeyError:
				if not delchi_f[complex_bin][complex_point][1] == 0.:
					print "nonzero element not appering in dimas file",complex_bin,complex_point,delchi_f[complex_bin][complex_point]
					ndiff+=1
				continue
			delta  = delchi_f[complex_bin][complex_point][0]-delchi_d[complex_bin][complex_point][0]
			delta2 = delchi_f[complex_bin][complex_point][0]+delchi_d[complex_bin][complex_point][0]
			if min(delta*delta,delta2*delta2) > 1.e-10:
				print "large difference in coma entries",min(delta*delta,delta2*delta2)**.5,"at",complex_bin,complex_point,delchi_f[complex_bin][complex_point][0],delchi_d[complex_bin][complex_point][0]
				ndiff+=1
			i = complex_point.real
			delta = deltas_f[complex_bin][i][0] - deltas_d[complex_bin][i][0]
			delta2= deltas_f[complex_bin][i][0] + deltas_d[complex_bin][i][0]
			if min(delta*delta,delta2*delta2) > 1.e-10:
				print "large difference in data entries",min(delta*delta,delta2*delta2)**.5,"at",complex_bin,i,deltas_f[complex_bin][i][0],deltas_d[complex_bin][i][0]
				ndiff+=1
	if ndiff==0:
		print "All matches"

d_not_f,f_not_d = pointdiff_cc(delchi_d,delchi_f)
print "evaluated by d, not by f",len(d_not_f)
print "evaluated by f, not by d",len(f_not_d)
print "contrib to d-chi2, by non f-points",list_chi2(d_not_f,deltas_d,delchi_d)
print "contrib to f-chi2, by non d-points",list_chi2(f_not_d,deltas_f,delchi_f)

print 'ddd',cross_chi2(deltas_d,deltas_d,delchi_d)
#print 'ddf',cross_chi2(deltas_d,deltas_d,delchi_f)
#print 'dfd',cross_chi2(deltas_d,deltas_f,delchi_d)
#print 'fdd',cross_chi2(deltas_f,deltas_d,delchi_d)
#print 'dff',cross_chi2(deltas_d,deltas_f,delchi_f)
#print 'fdf',cross_chi2(deltas_f,deltas_d,delchi_f)
#print 'ffd',cross_chi2(deltas_f,deltas_f,delchi_d)
print 'fff',cross_chi2(deltas_f,deltas_f,delchi_f)

print get_single_bins(d_not_f)
