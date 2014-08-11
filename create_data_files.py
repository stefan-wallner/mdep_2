#!/nfs/hicran/home/fkrinner/private/bin/python
import sys
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
from sys import argv
from convertTextOutput import getComaData
import yaml
import os
import datetime

def to_bool(string):
	if string.lower() in ['j','y','ja','yes','true']:
		return True
	else:
		return False

def writeFit(	name,			#-explenatory
		waves,			# wave names to be used
		limits,			# mass limits for the waves
		direct, 		# directory of the fit
		CONJUGATE=True ):	# Conjugate the fit result.
	"""Write text files to be read by the chi2 class"""
	nWaves = len(waves)
	upLims=[]
	loLims=[]
	for i in range(nWaves):
		upLims.append(limits[i][1])
		loLims.append(limits[i][0])
	data = getComaData(waves,upLims,loLims,direct, flagg='ANCHOR_FIRST', eps=1.E-3, CONJUGATE=CONJUGATE)
	pointsSorted=[]
	indices=[0]
	for i in range(1,nWaves):
		indices.append(i*nWaves) 
		indices.append(i)

	reduced_dat = []
	reduced_coma= []
	for bin in range(len(data[0])):
		dats = data[0][bin]
		dats_red=[]
		for i in indices:
			dats_red.append(dats[i])
		reduced_dat.append(dats_red)
		coma = data[1][bin]
		coma_red = []
		for i in indices:
			coma_line = []
			for j in indices:				
				coma_line.append(coma[i][j])
			coma_red.append(coma_line)
		reduced_coma.append(coma_red)
	dataname = name+'_data.dat'
	comaname = name+'_coma.dat'
	datFile = open(dataname,'w')
	comaFile= open(comaname,'w')
	for bin in range(len(reduced_dat)):
		dat = reduced_dat[bin]
		for ptt in dat:
			datFile.write(str(ptt)+'   ')
		coma = reduced_coma[bin]
		for line in coma:
			for ptt in line:
				comaFile.write(str(ptt)+'   ')
		comaFile.write('\n')
		datFile.write('\n')
	comaFile.close()
	datFile.close()
	return [os.path.abspath(dataname),os.path.abspath(comaname)]

try:
	config_file = argv[1]
except IndexError:
	config_file = raw_input("Config file: ")
try:
	config_stream=open(config_file,'r')
except IOError:
	print "Error: Config file ("+config_file+") does not exist"
	sys.exit(1)
config = yaml.load(config_stream)
config_stream.close()
try:
	config["data_files_created"]
	print "Data files already created. Exit"
	sys.exit(0)
except KeyError:
	pass

if not len(config["t_binning"]) -1 == len(config["mass_independent_fit"]):
	print "Error: Number if mass independent fits does not match t' binning"
	sys.exit(1)

try:
	if len(config["t_binning"])-1 ==  len(config["data_files"]) and len(config["t_binning"])-1 ==  len(config["coma_files"]):
		nnf = 0
		for dfile in config["data_files"]:
			if not os.path.isfile(dfile):
				nnf+=1
		for cfile in config["coma_files"]:
			if not os.path.isfile(cfile):
				nnf+=1
		if nnf ==0:
			proceed = raw_input("Already enough data- and coma files given, which all exist. Proceed? (y/n):")
			if not to_bool(proceed):
				sys.exit(0)
except KeyError:
	pass

waves = config["waves"]
wave_names = []
mass_limits=[]
Mmin = waves[0]["mmin"]
Mmax = waves[0]["mmax"]
for wave in waves:
	wave_names.append(wave["name"])
	mmin = wave["mmin"]
	mmax = wave["mmax"]
	if mmin < Mmin:
		raise ValueError("Lower mass limit smaller than anchor mass limit")
	if mmax > Mmax:
		raise ValueError("Upper mass limit bigger than anchor mass limit")
	mass_limits.append([mmin,mmax])
try:
	name = config["fit_name"]
except KeyError:
	name = raw_input("Fit_name: ")
	config["fit_name"] = name
try:
	conj = config["conjugate_fit_result"]
except KeyError: 
	conj = to_bool(raw_input("Conjugate fit result (default false)?(y/n):"))
	config["conjugate_fit_result"] = conj
apath=[]
tbinning = config["t_binning"]
for tbin in range(len(tbinning)-1):
	t_name = name+"_"+str(tbinning[tbin])+'_'+str(tbinning[tbin+1])
	apath.append(writeFit(t_name,wave_names,mass_limits,config["mass_independent_fit"][tbin],conj))

config["data_files"]=[]
config["coma_files"]=[]
for path in apath:
	config["data_files"].append(path[0])
	config["coma_files"].append(path[1])
config["data_files_created"]=str(datetime.datetime.now())
config_out =  yaml.dump(config, default_flow_style=False)

out=open(config_file,'w')
out.write(config_out)
out.close()







