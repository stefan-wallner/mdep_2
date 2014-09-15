#!/nfs/hicran/home/fkrinner/private/bin/python
import sys
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
from sys import argv
from convertTextOutput import getRelevantData
from convertTextOutput import getIntegralAverage
import yaml
import os
import datetime
import numpy as np
from numpy import linalg as la
import shutil

f0_list=['f0_0278_0320',  'f0_0320_0360',  'f0_0360_0400',  'f0_0400_0440',  'f0_0440_0480',  'f0_0480_0520',  'f0_0520_0560',  'f0_0560_0600',  'f0_0600_0640',  'f0_0640_0680',  'f0_0680_0720',  'f0_0720_0760',  'f0_0760_0800',  'f0_0800_0840',  'f0_0840_0880',  'f0_0880_0920',  'f0_0920_0930',  'f0_0930_0940',  'f0_0940_0950',  'f0_0950_0960',  'f0_0960_0970',  'f0_0970_0980',  'f0_0980_0990',  'f0_0990_1000',  'f0_1000_1010',  'f0_1010_1020',  'f0_1020_1030',  'f0_1030_1040',  'f0_1040_1050',  'f0_1050_1060',  'f0_1060_1070',  'f0_1070_1080',  'f0_1080_1120',  'f0_1120_1160',  'f0_1160_1200',  'f0_1200_1240',  'f0_1240_1280',  'f0_1280_1320',  'f0_1320_1360',  'f0_1360_1400',  'f0_1400_1440',  'f0_1440_1480',  'f0_1480_1520',  'f0_1520_1560',  'f0_1560_1600',  'f0_1600_1640',  'f0_1640_1680',  'f0_1680_1720',  'f0_1720_1760',  'f0_1760_1800',  'f0_1800_1840',  'f0_1840_1880',  'f0_1880_1920',  'f0_1920_1960',  'f0_1960_2000',  'f0_2000_2040',  'f0_2040_2080',  'f0_2080_2120',  'f0_2120_2160',  'f0_2160_2200',  'f0_2200_2240',  'f0_2240_2280']

def get_data_anchor(waves,up,low,direct,DIVIDE_PHASE_SPACE,int_dir,ACC_CORRECTED):
	"""Gets data points and covariance matrix for anchor wave fitting"""
	nWaves = len(waves)
	if len(waves) != len(up) or len(waves) != len(low):
		print "Lenghts of waves and limits do not match"
		raise IndexError
	raw_data = getRelevantData(waves,direct)
	if DIVIDE_PHASE_SPACE:
		raw_data = normalize_to_integrals(raw_data,waves,int_dir,ACC_CORRECTED)
	data_points=[]
	final_comas_inv=[]
	for point in raw_data:
		mass = (point[0]+point[1])/2
		wave_on=[]						# Check, which waves are active in the current bin
		for i in range(nWaves):					# Set all other elements to zero
			if mass > low[i] and mass < up[i]:		# So that the wrong correlations won't be taken into account
				wave_on.append(1.)
			else:
				wave_on.append(0.)
		for i in range(nWaves):					# Set inactive elements to zero
			point[2][2*i  ]*=wave_on[i]
			point[2][2*i+1]*=wave_on[i]	
			for j in range(nWaves):
				point[3][2*i  ][2*j  ]*=wave_on[i]*wave_on[j]			
				point[3][2*i+1][2*j  ]*=wave_on[i]*wave_on[j]			
				point[3][2*i  ][2*j+1]*=wave_on[i]*wave_on[j]			
				point[3][2*i+1][2*j+1]*=wave_on[i]*wave_on[j]		
		jacobian = []
		data_point=[]
		raw_coma=np.matrix(point[3])
		iijjs=[[0,0]]
		for iiii in range(1,nWaves):
			iijjs.append([iiii,0])
			iijjs.append([0,iiii])
		for iijj in iijjs:
			jacLine =[0.]*2*nWaves
			ii = iijj[0]
			jj = iijj[1]
			if ii == jj: # Intensities
				addpoint = point[2][2*ii]**2 + point[2][2*ii+1]**2
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  = 2*point[2][2*ii]   # dInt_i / dRe_i = 2*Re_i
					jacLine[2*ii+1]= 2*point[2][2*ii+1] # dInt_i / dIm_i = 2*Im_i
			elif ii > jj: # Real Part
				addpoint = point[2][2*ii]*point[2][2*jj] + point[2][2*ii+1]*point[2][2*jj+1]
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  = point[2][2*jj]	    # dRe_ij / dRe_i = Re_j			
					jacLine[2*jj]  = point[2][2*ii]	    # dRe_ij / dRe_j = Re_i
					jacLine[2*ii+1]= point[2][2*jj+1]   # dRe_ij / dIm_i = Im_j
					jacLine[2*jj+1]= point[2][2*ii+1]   # dRe_ij / dIm_j = Im_i
			else: # Imaginary Part		# POTENTIAL MINUS SIGN HERE ... INVESTIGATE
			#			R_j      *   I_i           -     R_i       *     I_j
				addpoint = point[2][2*jj]*point[2][2*ii+1] - point[2][2*ii]*point[2][2*jj+1]
				data_point.append(addpoint)
				if not addpoint ==0.:
					jacLine[2*ii]  =-point[2][2*jj+1]  # dIm_ij / dRe_i =-Im_j				
					jacLine[2*jj]  = point[2][2*ii+1]  # dIm_ij / dRe_j = Im_i
					jacLine[2*ii+1]= point[2][2*jj]	   # dIm_ij / dIm_i = Re_j
					jacLine[2*jj+1]=-point[2][2*ii]    # dIm_ij / dIm_j =-Re_i
			jacobian.append(jacLine)
		jacobian_M = np.matrix(jacobian)
		final_coma = jacobian_M*raw_coma*jacobian_M.T
		final_coma_inv = la.pinv(final_coma)
		final_comas_inv.append(final_coma_inv.tolist())
		data_points.append(data_point)
	return data_points,final_comas_inv


def normalize_to_integrals(raw_data,waves,int_dir,ACC_CORRECTED=True):
	"""Divides the data points in raw_data by their integral"""
	nWaves = len(waves)
	wave_60 = [make_60(wave) for wave in waves]
	for point in raw_data:
		m_low= point[0]		
		m_up = point[1]
		integral_matrix = getIntegralAverage(m_low,m_up,int_dir,acceptanceCorrected=ACC_CORRECTED,normalizeToDiag=False) # Normalization to diag would spoil everything ;) => False
		int_factors=[integral_matrix[wave][0]**.5 for wave in wave_60]
		for i in range(nWaves):
			try:
				point[2][2*i  ]/=int_factors[i]
				point[2][2*i+1]/=int_factors[i]
			except ZeroDivisionError:
				point[2][2*i  ] = 0.
				point[2][2*i+1] = 0.				
			for j in range(nWaves):
				try:
					point[3][2*i  ][2*j  ]/=(int_factors[i]*int_factors[j])
					point[3][2*i  ][2*j+1]/=(int_factors[i]*int_factors[j])
					point[3][2*i+1][2*j  ]/=(int_factors[i]*int_factors[j])
					point[3][2*i+1][2*j+1]/=(int_factors[i]*int_factors[j])
				except ZeroDivisionError:
					point[3][2*i  ][2*j  ] = 0.
					point[3][2*i  ][2*j+1] = 0.
					point[3][2*i+1][2*j  ] = 0.
					point[3][2*i+1][2*j+1] = 0.
	return raw_data					


def de_isobarred_list(name):
	"""Generates a list """
	type_isob = ''
	identifier = name.split()[1]
	if identifier in ['f0','f0(980)','[pipi]_S']:
		type_isob='f0'
	if identifier in ['f2']:
		type_isob='f2'
	if identifier in ['rho']:
		type_isob='rho'
	waves=[]
	if type_isob=='f0':
		for step in f0_list:
			repwave =name.replace(identifier,step)
			if len(repwave) < 60:
				waves.append(repwave)
			else:
				waves.append(repwave[:60])
	else:
		raise Exception
	return waves


def islist(obj):
	"""Checks if 'obj' is a list"""
	try:
		obj+=[]
	except TypeError:
		return False
	return True


def make_60(wave):
	"""Appends ' ' to wave until len(wave) = 60"""
	wave_60 = wave
	while len(wave_60) < 60:
		wave_60+=' '
	return wave_60


def write_fit(name,wave_names,lower_lims,upper_lims,indep_fit,DIVIDE_PHASE_SPACE=False,int_dir='',ACC_CORRECTED=True):
	"""Writes data and coma files"""
	datas,comas=get_data_anchor(wave_names,upper_lims,lower_lims,indep_fit,DIVIDE_PHASE_SPACE,int_dir,ACC_CORRECTED)
	dataname = './data_files/'+name+'_data.dat'
	comaname = './data_files/'+name+'_coma.dat'
	dataFile = open(dataname,'w')
	comaFile= open(comaname,'w')
	for bin in range(len(datas)):
		data = datas[bin]
#		print "len(data)",len(data)
		for ptt in data:
			dataFile.write(str(ptt)+' ')
		coma = comas[bin]
#		print "len(coma)",len(coma)
#		print "len(coma[0])",len(coma[0])
		for line in coma:
			for ptt in line:
				comaFile.write(str(ptt)+' ')
		comaFile.write('\n')
		dataFile.write('\n')
	dataFile.close()
	comaFile.close()
	return [os.path.abspath(dataname),os.path.abspath(comaname)]


def process_config(config_file, date_time_flag=False):
	"""Creates data and coma files for a config file and wirtes the according paths to the config_file"""
	try:
		config_stream=open(config_file,'r')
	except IOError:
		print "Error: Config file ("+config_file+") does not exist"
		return 1
	config = yaml.load(config_stream)
	waves = config["waves"]
	wave_names = []
	upper_lims=[]
	lower_lims=[]
	Mmin = waves[0]["mmin"]
	Mmax = waves[0]["mmax"]
	for wave in waves:	
		isdeis = False
		if len(wave["parametrizations"][0]) == 2 and islist(wave["parametrizations"][0]): # islist is needed, because len(...) == 2 can be true for strings
			isdeis= True
			print wave,"is de-isobarred. Automatically append corresponding waves" 
			deiso_waves = de_isobarred_list(wave["name"]) 
			wave_names+=deiso_waves
			nSteps = len(deiso_waves)
			for i in range(nSteps):
				mmin = wave["mmin"]
				mmax = wave["mmax"]
				if mmin < Mmin:
					raise ValueError("Lower mass limit smaller than anchor mass limit")
				if mmax > Mmax:
					raise ValueError("Upper mass limit bigger than anchor mass limit")
				upper_lims.append(mmax)
				lower_lims.append(mmin)
		if not isdeis:
			wave_names.append(wave["name"])
			mmin = wave["mmin"]
			mmax = wave["mmax"]
			if mmin < Mmin:
				raise ValueError("Lower mass limit smaller than anchor mass limit")
			if mmax > Mmax:
				raise ValueError("Upper mass limit bigger than anchor mass limit")
			upper_lims.append(mmax)
			lower_lims.append(mmin)
	if len(wave_names)!=len(set(wave_names)):
		print "Error: Some wave(s) appear multiple times in 'wave_names'."
		return 1
	add_name_to_file = False
	try:
		name = config["fit_name"]
	except KeyError:
		add_name_to_file=True
		name = raw_input("Fit_name: ")
		config["fit_name"] = name
	apath=[]
	tbinning = config["t_binning"]
	norm_to_int = False
	acc_corr = True
	try:
		norm_to_int = config["divide_integrals"]
		try:
			acc_corr = config["acc_corrected"]
		except KeyError:
			pass
	except KeyError:
		pass

	for tbin in range(len(tbinning)-1):
		t_name = name+"_"+str(tbinning[tbin])+'_'+str(tbinning[tbin+1])
		int_dir=''
		if norm_to_int:
			try:
				int_dir = config["integral_dir"][tbin]
			except (KeyError, IndexError):
				print "No integral directory given"
				return 1
		apath.append(write_fit(t_name,wave_names,lower_lims,upper_lims,config["mass_independent_fit"][tbin],norm_to_int,int_dir,acc_corr))
	config["data_files"]=[]
	config["coma_files"]=[]
	for path in apath:
		config["data_files"].append(path[0])
		config["coma_files"].append(path[1])
	config["data_files_created"]=str(datetime.datetime.now())
	if date_time_flag: # Make the writing a bit more complicated that just yaml.dump(...), to keep comments in the old file.
		out_card_name = (os.path.splitext(config_file)[0]+'_'+str(datetime.datetime.now())+os.path.splitext(config_file)[1]).replace(' ','_').replace(':','.')
		shutil.copyfile(os.path.abspath(config_file),"/"+"/".join(os.path.abspath(config_file).split('/')[:-1])+'/'+out_card_name)
	else:
		out_card_name = config_file
	out_card = open(out_card_name,'a')
	if add_name_to_file:
		out_card.write("fit_name: "+name+'\n')
	out_card.write("data_files_created: "+config["data_files_created"]+'\n')
	out_card.write("data_files:\n")
	for data_file in config["data_files"]:
		out_card.write("- "+data_file+'\n')
	out_card.write("coma_files:\n")
	for coma_file in config["coma_files"]:
		out_card.write("- "+coma_file+'\n')
	out_card.close()
	return 0,out_card_name


if __name__ == "__main__":
	try:
		config_file = argv[1]
	except IndexError:
		config_file = raw_input("Config file:")
	process_config(config_file,True)


