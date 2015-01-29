import rootpy
import ROOT
rootpy.log.basic_config_colorized()
from rootpy.io import root_open
from rootpy.plotting import Canvas, Hist, Legend
from ROOT import TH1D

def get_old_method_files(data_file, waveset, target,sfit):
	"""Writes the data from .root files to be used by the full coma approach of fitmd_2"""
	numberNdigit = 3
	waveNumbers = get_wave_numbers(waveset, sfit)
	nWaves = len(waveset)
	histNames = [None]*nWaves**2
	for i in range(nWaves):
		w1 = getNdigitInt(waveNumbers[i],numberNdigit)
		for j in range(nWaves):
			w2 = getNdigitInt(waveNumbers[j],numberNdigit)
			ii = i*nWaves+j
			if i==j:
#				histName = ['h1'+w1]
				histName = ['h'+str(waveNumbers[i])]
			if i>j:
				histName = ['h2'+w1+w2,'h2'+w2+w1]
			if i<j:
				histName = ['h1'+w1+w2,'h1'+w2+w1]
			histNames[ii] = histName
	isChanged = [1]*nWaves**2
	with root_open(data_file,"READ") as inROOT:
		hists = [None]*nWaves**2
		for ii in range(nWaves**2):
			try:
				hists[ii] = inROOT.get(histNames[ii][0])	
			except rootpy.io.DoesNotExist:
				isChanged[ii] = -1
				hists[ii] = inROOT.get(histNames[ii][1])	
		nBins = hists[0].GetNbinsX()
		with open(target+"_data",'w') as outoutout:
			for bin in range(nBins):
				for i in range(len(hists)):
					hist=hists[i]
					fakk=1
					if histNames[i][0][1] == '2':
						fakk = isChanged[i]
					outoutout.write(str(fakk*hist.GetBinContent(bin+1))+' ')
				outoutout.write("\n")
		with open(target+"_coma",'w') as outoutout:
			for bin in range(nBins):
				for i in range(nWaves**2):
					for j in range(nWaves**2):
						if i == j:
							try:
								err = 1/hists[i].GetBinError(bin+1)
							except ZeroDivisionError:
								err = 0.
							err **= 2
							outoutout.write(str(err)+' ')
						else:
							outoutout.write("0.0 ")
				outoutout.write('\n')


def getNdigitInt(i,ndigit):
	"""Returns a string from i, which has at least length ndigit. Appends '0' until it does"""
	string = str(i)
	while len(string) < ndigit:
		string = '0'+string
	return string


def get_wave_numbers(waveset,sfit):
	"""Returns a list of wave numbers for a given waveset from an sfit file"""
	nWaves = len(waveset)
	waveNumbers = {}
	for wave in waveset:
		waveNumbers[wave] = None
	with open(sfit,'r') as ininin:
		for line in ininin.readlines():
			for wave in waveset:
				if wave in line:
					waveNumbers[wave] = int(line.split()[0])	
	for wave in waveset:
		if not waveNumbers[wave]:
			print wave
			raise KeyError # Wave not found
	numbers = []
	for wave in waveset:
		numbers.append(waveNumbers[wave])
	return numbers







if __name__ == "__main__":
	waveset = [	"1-(1++)0+ rho pi S",
			"1-(1-+)1+ rho pi P",
			"1-(1++)0+ rho pi D",
			"1-(0-+)0+ f0(980) pi S",
			"1-(2++)1+ f2 pi P",
			"1-(2++)1+ rho pi D",
			"1-(2++)2+ rho pi D",
			"1-(2-+)0+ f2 pi D",
			"1-(2-+)0+ f2 pi S",
			"1-(2-+)0+ rho pi F",
			"1-(2-+)1+ f2 pi S",
			"1-(4++)1+ f2 pi F",
			"1-(4++)1+ rho pi G"]
	data_files = [	"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.100000-0.112853.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.112853-0.127471.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.127471-0.144385.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.144385-0.164401.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.164401-0.188816.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.188816-0.219907.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.219907-0.262177.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.262177-0.326380.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.326380-0.448588.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.448588-0.724294.root",
			"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/hfit_0.724294-1.000000.root"]
	target = 	"/nfs/mds/user/fkrinner/data_mdep_fit/13w_exotic_old_method/"
	sfit   = 	"/nfs/hicran/project/compass/analysis/swallner/data_massindependentfits/amplitudes/florian/sfit.dat"
	for i in range(len(data_files)):
		get_old_method_files(data_files[i],waveset, target+str(i),sfit)
		print data_files[i]



