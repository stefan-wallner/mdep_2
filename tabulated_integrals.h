// This file was created by a python script:
// /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/chi_squared/get_tabulated_integrals.py
#ifndef TABULATED_INTEGRALS_TOP_LOL
#define TABULATED_INTEGRALS_TOP_LOL

//1-(0-+)0+ f0(980) pi S	0
//1-(1++)0+ f0(980) pi P	3
//1-(1++)0+ rho pi D	7
//1-(1++)0+ rho pi S	2
//1-(2++)1+ f2 pi P	8
//1-(2++)1+ rho pi D	5
//1-(2++)2+ rho pi D	9
//1-(2-+)0+ f2 pi D	12
//1-(2-+)0+ f2 pi S	4
//1-(2-+)0+ rho pi F	11
//1-(2-+)1+ f2 pi S	10
//1-(4++)1+ f2 pi F	6
//1-(4++)1+ rho pi G	1

double tabulated_integrals(double m3Pi,int wave_number){
	int nnnMax = 13;
	if(wave_number >= nnnMax){
		std::cerr << "tabulated_integrals.h: Error: Wave not tabulated: wave_number = " << wave_number << std::endl;
		return 0.;
	};
/*
	double ms[] = {
0.500, 0.510, 0.520, 0.530, 0.540, 0.550, 0.560, 0.570, 0.580, 0.590, 0.600, 0.610, 0.620, 0.630, 0.640, 0.650, 0.660, 0.670, 0.680, 0.690, 0.700, 0.710, 0.720, 0.730, 0.740, 0.750, 0.760, 0.770, 0.780, 0.790, 0.800, 0.810, 0.820, 0.830, 0.840, 0.850, 0.860, 0.870, 0.880, 0.890, 0.900, 0.910, 0.920, 0.930, 0.940, 0.950, 0.960, 0.970, 0.980, 0.990, 1.000, 1.010, 1.020, 1.030, 1.040, 1.050, 1.060, 1.070, 1.080, 1.090, 1.100, 1.110, 1.120, 1.130, 1.140, 1.150, 1.160, 1.170, 1.180, 1.190, 1.200, 1.210, 1.220, 1.230, 1.240, 1.250, 1.260, 1.270, 1.280, 1.290, 1.300, 1.310, 1.320, 1.330, 1.340, 1.350, 1.360, 1.370, 1.380, 1.390, 1.400, 1.410, 1.420, 1.430, 1.440, 1.450, 1.460, 1.470, 1.480, 1.490, 1.500, 1.510, 1.520, 1.530, 1.540, 1.550, 1.560, 1.570, 1.580, 1.590, 1.600, 1.610, 1.620, 1.630, 1.640, 1.650, 1.660, 1.670, 1.680, 1.690, 1.700, 1.710, 1.720, 1.730, 1.740, 1.750, 1.760, 1.770, 1.780, 1.790, 1.800, 1.810, 1.820, 1.830, 1.840, 1.850, 1.860, 1.870, 1.880, 1.890, 1.900, 1.910, 1.920, 1.930, 1.940, 1.950, 1.960, 1.970, 1.980, 1.990, 2.000, 2.010, 2.020, 2.030, 2.040, 2.050, 2.060, 2.070, 2.080, 2.090, 2.100, 2.110, 2.120, 2.130, 2.140, 2.150, 2.160, 2.170, 2.180, 2.190, 2.200, 2.210, 2.220, 2.230, 2.240, 2.250, 2.260, 2.270, 2.280, 2.290, 2.300, 2.310, 2.320, 2.330, 2.340, 2.350, 2.360, 2.370, 2.380, 2.390, 2.400, 2.410, 2.420, 2.430, 2.440, 2.450, 2.460, 2.470, 2.480, 2.490, 2.500};
*/
	 double ints[13][201] = {
{0.112469E-03, 0.147230E-03, 0.188097E-03, 0.235890E-03, 0.290607E-03, 0.353213E-03, 0.424108E-03, 0.503395E-03, 0.593660E-03, 0.692959E-03, 0.803440E-03, 0.927628E-03, 0.106328E-02, 0.121202E-02, 0.137522E-02, 0.155546E-02, 0.175209E-02, 0.196675E-02, 0.219894E-02, 0.245404E-02, 0.272821E-02, 0.302826E-02, 0.335015E-02, 0.369830E-02, 0.407548E-02, 0.448216E-02, 0.492544E-02, 0.539985E-02, 0.590645E-02, 0.646152E-02, 0.705708E-02, 0.768956E-02, 0.837682E-02, 0.912649E-02, 0.991705E-02, 0.107653E-01, 0.116972E-01, 0.126818E-01, 0.137529E-01, 0.149092E-01, 0.161572E-01, 0.174969E-01, 0.189583E-01, 0.205056E-01, 0.221914E-01, 0.240360E-01, 0.260472E-01, 0.281760E-01, 0.305918E-01, 0.331613E-01, 0.360164E-01, 0.391186E-01, 0.425227E-01, 0.463326E-01, 0.506445E-01, 0.554973E-01, 0.608944E-01, 0.672357E-01, 0.742680E-01, 0.828596E-01, 0.932822E-01, 0.105946E+00, 0.121562E+00, 0.141773E+00, 0.160597E+00, 0.178501E+00, 0.196053E+00, 0.212423E+00, 0.228534E+00, 0.245117E+00, 0.261175E+00, 0.276396E+00, 0.292681E+00, 0.307476E+00, 0.322772E+00, 0.338889E+00, 0.353679E+00, 0.368563E+00, 0.384178E+00, 0.397552E+00, 0.413178E+00, 0.427542E+00, 0.441636E+00, 0.457274E+00, 0.472687E+00, 0.485522E+00, 0.499730E+00, 0.512956E+00, 0.527904E+00, 0.542552E+00, 0.555804E+00, 0.569851E+00, 0.584913E+00, 0.600244E+00, 0.611973E+00, 0.624031E+00, 0.638775E+00, 0.651966E+00, 0.666751E+00, 0.680328E+00, 0.693510E+00, 0.706613E+00, 0.718325E+00, 0.737286E+00, 0.746259E+00, 0.759391E+00, 0.770353E+00, 0.786590E+00, 0.796180E+00, 0.813102E+00, 0.825158E+00, 0.840019E+00, 0.853552E+00, 0.868300E+00, 0.876879E+00, 0.890810E+00, 0.905203E+00, 0.916840E+00, 0.925035E+00, 0.940658E+00, 0.948965E+00, 0.967031E+00, 0.975175E+00, 0.986826E+00, 0.998754E+00, 0.101510E+01, 0.102396E+01, 0.103259E+01, 0.105170E+01, 0.106337E+01, 0.107187E+01, 0.108994E+01, 0.110220E+01, 0.111049E+01, 0.112253E+01, 0.113590E+01, 0.114890E+01, 0.115945E+01, 0.117437E+01, 0.118830E+01, 0.119587E+01, 0.120523E+01, 0.122045E+01, 0.123199E+01, 0.124167E+01, 0.125889E+01, 0.126362E+01, 0.127963E+01, 0.128874E+01, 0.130034E+01, 0.131073E+01, 0.132450E+01, 0.133806E+01, 0.134455E+01, 0.136003E+01, 0.137529E+01, 0.138276E+01, 0.139501E+01, 0.140155E+01, 0.140993E+01, 0.143385E+01, 0.143789E+01, 0.145120E+01, 0.146202E+01, 0.147718E+01, 0.147816E+01, 0.149653E+01, 0.150445E+01, 0.152128E+01, 0.152650E+01, 0.153313E+01, 0.154890E+01, 0.156121E+01, 0.156795E+01, 0.158587E+01, 0.159797E+01, 0.159986E+01, 0.161151E+01, 0.161823E+01, 0.163515E+01, 0.165198E+01, 0.166116E+01, 0.166740E+01, 0.167740E+01, 0.168985E+01, 0.170359E+01, 0.170847E+01, 0.173157E+01, 0.173576E+01, 0.174650E+01, 0.175545E+01, 0.176053E+01, 0.176865E+01, 0.178300E+01, 0.180645E+01, 0.180966E+01, 0.182484E+01, 0.182113E+01, 0.184870E+01, 0.184538E+01, 0.186351E+01},
{0.350498E-17, 0.839381E-17, 0.185180E-16, 0.380120E-16, 0.738135E-16, 0.136209E-15, 0.240888E-15, 0.411131E-15, 0.682125E-15, 0.110028E-14, 0.172201E-14, 0.265968E-14, 0.400186E-14, 0.591511E-14, 0.861891E-14, 0.123596E-13, 0.174730E-13, 0.243816E-13, 0.337534E-13, 0.461997E-13, 0.623874E-13, 0.836771E-13, 0.111079E-12, 0.146468E-12, 0.191591E-12, 0.248811E-12, 0.320103E-12, 0.410474E-12, 0.522496E-12, 0.660264E-12, 0.834207E-12, 0.104498E-11, 0.130520E-11, 0.162042E-11, 0.200840E-11, 0.247145E-11, 0.303289E-11, 0.371570E-11, 0.451678E-11, 0.548556E-11, 0.662860E-11, 0.802412E-11, 0.966139E-11, 0.116018E-10, 0.138816E-10, 0.166184E-10, 0.197631E-10, 0.235473E-10, 0.279273E-10, 0.332064E-10, 0.392375E-10, 0.463572E-10, 0.545754E-10, 0.640945E-10, 0.753174E-10, 0.882240E-10, 0.103007E-09, 0.120315E-09, 0.140396E-09, 0.162355E-09, 0.188792E-09, 0.218468E-09, 0.252638E-09, 0.290772E-09, 0.334598E-09, 0.383467E-09, 0.438570E-09, 0.501468E-09, 0.571171E-09, 0.647788E-09, 0.736689E-09, 0.833400E-09, 0.938834E-09, 0.105713E-08, 0.119024E-08, 0.133389E-08, 0.149641E-08, 0.166542E-08, 0.185852E-08, 0.207051E-08, 0.228824E-08, 0.253854E-08, 0.280132E-08, 0.309350E-08, 0.339728E-08, 0.373244E-08, 0.409310E-08, 0.448210E-08, 0.488341E-08, 0.533521E-08, 0.578727E-08, 0.629132E-08, 0.682474E-08, 0.739635E-08, 0.799343E-08, 0.860528E-08, 0.925320E-08, 0.996484E-08, 0.107126E-07, 0.114261E-07, 0.122325E-07, 0.131181E-07, 0.139388E-07, 0.148723E-07, 0.158008E-07, 0.168412E-07, 0.179146E-07, 0.189158E-07, 0.200272E-07, 0.211298E-07, 0.223526E-07, 0.235693E-07, 0.248549E-07, 0.262245E-07, 0.274031E-07, 0.288299E-07, 0.300646E-07, 0.316602E-07, 0.329614E-07, 0.345423E-07, 0.361641E-07, 0.377454E-07, 0.390932E-07, 0.408252E-07, 0.423967E-07, 0.441733E-07, 0.458885E-07, 0.475557E-07, 0.493250E-07, 0.512519E-07, 0.529927E-07, 0.547761E-07, 0.567120E-07, 0.585268E-07, 0.605406E-07, 0.619956E-07, 0.642234E-07, 0.660332E-07, 0.680479E-07, 0.701185E-07, 0.720289E-07, 0.739296E-07, 0.759264E-07, 0.783565E-07, 0.801772E-07, 0.822487E-07, 0.844286E-07, 0.864871E-07, 0.885523E-07, 0.905315E-07, 0.930177E-07, 0.943048E-07, 0.969357E-07, 0.991674E-07, 0.101435E-06, 0.103660E-06, 0.105154E-06, 0.107665E-06, 0.109647E-06, 0.112432E-06, 0.113837E-06, 0.116205E-06, 0.118040E-06, 0.120086E-06, 0.122861E-06, 0.125024E-06, 0.126691E-06, 0.129073E-06, 0.130742E-06, 0.132991E-06, 0.135408E-06, 0.138034E-06, 0.139691E-06, 0.141612E-06, 0.143975E-06, 0.146460E-06, 0.148801E-06, 0.150566E-06, 0.152593E-06, 0.154789E-06, 0.157113E-06, 0.158500E-06, 0.161430E-06, 0.163174E-06, 0.165319E-06, 0.166616E-06, 0.169795E-06, 0.171650E-06, 0.174492E-06, 0.175620E-06, 0.177205E-06, 0.178961E-06, 0.182370E-06, 0.183882E-06, 0.185885E-06, 0.188089E-06, 0.190134E-06, 0.192040E-06, 0.193986E-06, 0.195972E-06, 0.197865E-06},
{0.828747E-04, 0.119396E-03, 0.165853E-03, 0.224729E-03, 0.296099E-03, 0.383103E-03, 0.487569E-03, 0.610453E-03, 0.756431E-03, 0.923993E-03, 0.112120E-02, 0.134816E-02, 0.160801E-02, 0.190522E-02, 0.224284E-02, 0.263028E-02, 0.307019E-02, 0.356965E-02, 0.412414E-02, 0.476003E-02, 0.547541E-02, 0.628289E-02, 0.718413E-02, 0.820120E-02, 0.935630E-02, 0.106511E-01, 0.121368E-01, 0.137954E-01, 0.156606E-01, 0.178261E-01, 0.202342E-01, 0.229549E-01, 0.260820E-01, 0.296827E-01, 0.336758E-01, 0.382648E-01, 0.434646E-01, 0.492877E-01, 0.559089E-01, 0.630713E-01, 0.710498E-01, 0.792864E-01, 0.880909E-01, 0.969494E-01, 0.106020E+00, 0.115247E+00, 0.124601E+00, 0.133377E+00, 0.142536E+00, 0.151407E+00, 0.160071E+00, 0.168917E+00, 0.177064E+00, 0.185186E+00, 0.193391E+00, 0.201437E+00, 0.208964E+00, 0.217439E+00, 0.225282E+00, 0.232292E+00, 0.240987E+00, 0.247715E+00, 0.254889E+00, 0.262530E+00, 0.269667E+00, 0.277170E+00, 0.283574E+00, 0.290709E+00, 0.298019E+00, 0.304677E+00, 0.311738E+00, 0.318371E+00, 0.325063E+00, 0.331531E+00, 0.338294E+00, 0.345978E+00, 0.353070E+00, 0.358339E+00, 0.364979E+00, 0.372033E+00, 0.376611E+00, 0.384356E+00, 0.390421E+00, 0.396900E+00, 0.403182E+00, 0.409602E+00, 0.416282E+00, 0.421523E+00, 0.427239E+00, 0.434219E+00, 0.439839E+00, 0.445770E+00, 0.452945E+00, 0.459134E+00, 0.465763E+00, 0.470993E+00, 0.476618E+00, 0.483107E+00, 0.489391E+00, 0.494240E+00, 0.499855E+00, 0.506271E+00, 0.511038E+00, 0.517987E+00, 0.522800E+00, 0.530359E+00, 0.537300E+00, 0.541003E+00, 0.547891E+00, 0.552555E+00, 0.559174E+00, 0.564877E+00, 0.571321E+00, 0.577815E+00, 0.581838E+00, 0.588778E+00, 0.591773E+00, 0.600050E+00, 0.603530E+00, 0.610342E+00, 0.617068E+00, 0.622847E+00, 0.626272E+00, 0.632749E+00, 0.637191E+00, 0.644256E+00, 0.650033E+00, 0.654849E+00, 0.660537E+00, 0.668026E+00, 0.672787E+00, 0.677600E+00, 0.683930E+00, 0.689204E+00, 0.695729E+00, 0.697804E+00, 0.705888E+00, 0.710389E+00, 0.715565E+00, 0.722059E+00, 0.727126E+00, 0.731416E+00, 0.737851E+00, 0.745303E+00, 0.748977E+00, 0.754735E+00, 0.761442E+00, 0.765953E+00, 0.770909E+00, 0.776332E+00, 0.784028E+00, 0.784562E+00, 0.792602E+00, 0.799097E+00, 0.804892E+00, 0.810655E+00, 0.813391E+00, 0.820551E+00, 0.824490E+00, 0.833435E+00, 0.834578E+00, 0.841072E+00, 0.846205E+00, 0.850769E+00, 0.859027E+00, 0.862763E+00, 0.867024E+00, 0.873751E+00, 0.875992E+00, 0.881733E+00, 0.888273E+00, 0.895848E+00, 0.899949E+00, 0.903797E+00, 0.910254E+00, 0.917487E+00, 0.922205E+00, 0.926612E+00, 0.931849E+00, 0.937090E+00, 0.943267E+00, 0.945483E+00, 0.953770E+00, 0.957121E+00, 0.963570E+00, 0.964004E+00, 0.975243E+00, 0.978871E+00, 0.988536E+00, 0.988980E+00, 0.991696E+00, 0.995679E+00, 0.100744E+01, 0.101043E+01, 0.101271E+01, 0.101861E+01, 0.102485E+01, 0.103015E+01, 0.103572E+01, 0.103964E+01, 0.104450E+01},
{0.610048E-07, 0.882564E-07, 0.123137E-06, 0.167610E-06, 0.221598E-06, 0.287536E-06, 0.367395E-06, 0.461926E-06, 0.573368E-06, 0.702010E-06, 0.854686E-06, 0.102712E-05, 0.122501E-05, 0.145331E-05, 0.170856E-05, 0.200290E-05, 0.233305E-05, 0.270791E-05, 0.311246E-05, 0.357763E-05, 0.409285E-05, 0.466532E-05, 0.529140E-05, 0.598484E-05, 0.676376E-05, 0.762250E-05, 0.856852E-05, 0.961938E-05, 0.107380E-04, 0.120343E-04, 0.133614E-04, 0.148728E-04, 0.165104E-04, 0.183436E-04, 0.202958E-04, 0.224335E-04, 0.248301E-04, 0.273287E-04, 0.301512E-04, 0.332945E-04, 0.367239E-04, 0.403622E-04, 0.444770E-04, 0.489327E-04, 0.538504E-04, 0.590409E-04, 0.651506E-04, 0.715313E-04, 0.787123E-04, 0.866043E-04, 0.955494E-04, 0.105252E-03, 0.116029E-03, 0.128210E-03, 0.141899E-03, 0.157634E-03, 0.175010E-03, 0.195578E-03, 0.218624E-03, 0.246280E-03, 0.279161E-03, 0.317982E-03, 0.365569E-03, 0.427768E-03, 0.503105E-03, 0.591549E-03, 0.689146E-03, 0.793130E-03, 0.904487E-03, 0.102516E-02, 0.114996E-02, 0.127614E-02, 0.141391E-02, 0.154720E-02, 0.168567E-02, 0.183450E-02, 0.197594E-02, 0.212340E-02, 0.227493E-02, 0.241701E-02, 0.257471E-02, 0.272406E-02, 0.287608E-02, 0.303727E-02, 0.320088E-02, 0.334356E-02, 0.349609E-02, 0.364564E-02, 0.380926E-02, 0.397228E-02, 0.412412E-02, 0.428209E-02, 0.444233E-02, 0.461649E-02, 0.474834E-02, 0.489933E-02, 0.505649E-02, 0.521695E-02, 0.537894E-02, 0.552821E-02, 0.568773E-02, 0.583339E-02, 0.598849E-02, 0.618149E-02, 0.629757E-02, 0.645561E-02, 0.658822E-02, 0.677486E-02, 0.688646E-02, 0.706828E-02, 0.722285E-02, 0.738318E-02, 0.754931E-02, 0.768834E-02, 0.782300E-02, 0.797848E-02, 0.812968E-02, 0.828597E-02, 0.839881E-02, 0.857157E-02, 0.868645E-02, 0.887227E-02, 0.900010E-02, 0.914099E-02, 0.928911E-02, 0.945217E-02, 0.958678E-02, 0.970710E-02, 0.990327E-02, 0.100342E-01, 0.101525E-01, 0.103503E-01, 0.104776E-01, 0.105997E-01, 0.107011E-01, 0.108884E-01, 0.110456E-01, 0.111675E-01, 0.113507E-01, 0.115032E-01, 0.115762E-01, 0.117183E-01, 0.118899E-01, 0.120218E-01, 0.121723E-01, 0.123628E-01, 0.124234E-01, 0.125757E-01, 0.127186E-01, 0.128429E-01, 0.129819E-01, 0.131118E-01, 0.132785E-01, 0.134195E-01, 0.135766E-01, 0.137186E-01, 0.138434E-01, 0.140176E-01, 0.140701E-01, 0.142017E-01, 0.144211E-01, 0.144975E-01, 0.146415E-01, 0.147597E-01, 0.149310E-01, 0.150185E-01, 0.152169E-01, 0.153208E-01, 0.154525E-01, 0.155361E-01, 0.156958E-01, 0.158479E-01, 0.159775E-01, 0.160817E-01, 0.162626E-01, 0.164043E-01, 0.164791E-01, 0.166221E-01, 0.166838E-01, 0.168701E-01, 0.170719E-01, 0.172036E-01, 0.172795E-01, 0.174054E-01, 0.175659E-01, 0.176787E-01, 0.177754E-01, 0.180173E-01, 0.181143E-01, 0.182310E-01, 0.183366E-01, 0.184272E-01, 0.185614E-01, 0.186903E-01, 0.189239E-01, 0.189339E-01, 0.191353E-01, 0.191659E-01, 0.194479E-01, 0.194618E-01, 0.196555E-01},
{0.593278E-10, 0.970416E-10, 0.150977E-09, 0.226788E-09, 0.328167E-09, 0.462846E-09, 0.637948E-09, 0.859807E-09, 0.114025E-08, 0.148270E-08, 0.190999E-08, 0.242116E-08, 0.303981E-08, 0.377479E-08, 0.464166E-08, 0.566477E-08, 0.686285E-08, 0.825859E-08, 0.983478E-08, 0.116679E-07, 0.137778E-07, 0.161602E-07, 0.188506E-07, 0.218846E-07, 0.253250E-07, 0.291512E-07, 0.335085E-07, 0.383122E-07, 0.436356E-07, 0.496797E-07, 0.561551E-07, 0.633452E-07, 0.713194E-07, 0.801985E-07, 0.896625E-07, 0.100167E-06, 0.111739E-06, 0.124122E-06, 0.138003E-06, 0.152948E-06, 0.169564E-06, 0.186975E-06, 0.206314E-06, 0.226941E-06, 0.249411E-06, 0.273607E-06, 0.300681E-06, 0.328428E-06, 0.359883E-06, 0.392526E-06, 0.428810E-06, 0.467068E-06, 0.508290E-06, 0.553220E-06, 0.601910E-06, 0.654943E-06, 0.710453E-06, 0.772575E-06, 0.835933E-06, 0.907404E-06, 0.984971E-06, 0.106577E-05, 0.115045E-05, 0.124725E-05, 0.134776E-05, 0.146097E-05, 0.158047E-05, 0.170688E-05, 0.184769E-05, 0.200145E-05, 0.216288E-05, 0.233916E-05, 0.254058E-05, 0.275127E-05, 0.297514E-05, 0.322824E-05, 0.350221E-05, 0.380127E-05, 0.412334E-05, 0.448816E-05, 0.488556E-05, 0.531540E-05, 0.582147E-05, 0.635030E-05, 0.696437E-05, 0.761360E-05, 0.832097E-05, 0.914010E-05, 0.100257E-04, 0.109993E-04, 0.120674E-04, 0.131798E-04, 0.143841E-04, 0.156673E-04, 0.169615E-04, 0.183375E-04, 0.197084E-04, 0.211371E-04, 0.225144E-04, 0.239917E-04, 0.254514E-04, 0.268264E-04, 0.282939E-04, 0.296751E-04, 0.311240E-04, 0.324638E-04, 0.338793E-04, 0.352007E-04, 0.364780E-04, 0.379074E-04, 0.391925E-04, 0.405589E-04, 0.419043E-04, 0.431451E-04, 0.444176E-04, 0.457928E-04, 0.470493E-04, 0.482956E-04, 0.496016E-04, 0.508002E-04, 0.520946E-04, 0.532427E-04, 0.545219E-04, 0.558081E-04, 0.570305E-04, 0.583158E-04, 0.593981E-04, 0.604209E-04, 0.618086E-04, 0.629085E-04, 0.641162E-04, 0.654750E-04, 0.663210E-04, 0.676885E-04, 0.685690E-04, 0.700584E-04, 0.712939E-04, 0.724517E-04, 0.735608E-04, 0.746348E-04, 0.755313E-04, 0.769342E-04, 0.776433E-04, 0.789611E-04, 0.802502E-04, 0.813600E-04, 0.826199E-04, 0.839624E-04, 0.847976E-04, 0.858896E-04, 0.868131E-04, 0.879745E-04, 0.891388E-04, 0.902189E-04, 0.912531E-04, 0.922725E-04, 0.933789E-04, 0.945925E-04, 0.953311E-04, 0.966491E-04, 0.977256E-04, 0.990696E-04, 0.100183E-03, 0.101009E-03, 0.101883E-03, 0.103247E-03, 0.104090E-03, 0.105223E-03, 0.106653E-03, 0.107269E-03, 0.108383E-03, 0.109763E-03, 0.110750E-03, 0.111847E-03, 0.112606E-03, 0.113680E-03, 0.114977E-03, 0.116027E-03, 0.116758E-03, 0.118219E-03, 0.118963E-03, 0.120021E-03, 0.120868E-03, 0.122446E-03, 0.122728E-03, 0.124410E-03, 0.125046E-03, 0.126542E-03, 0.127041E-03, 0.128136E-03, 0.129208E-03, 0.130329E-03, 0.131196E-03, 0.132251E-03, 0.133968E-03, 0.134992E-03, 0.135547E-03, 0.136600E-03, 0.137598E-03, 0.138743E-03, 0.139729E-03},
{0.895361E-10, 0.164290E-09, 0.284361E-09, 0.468130E-09, 0.739769E-09, 0.112748E-08, 0.166787E-08, 0.240441E-08, 0.339657E-08, 0.470758E-08, 0.637608E-08, 0.858324E-08, 0.113008E-07, 0.147238E-07, 0.189635E-07, 0.242066E-07, 0.305615E-07, 0.382242E-07, 0.475739E-07, 0.588197E-07, 0.718385E-07, 0.875941E-07, 0.105801E-06, 0.127347E-06, 0.152462E-06, 0.181991E-06, 0.215769E-06, 0.255283E-06, 0.299680E-06, 0.351481E-06, 0.412712E-06, 0.481230E-06, 0.560352E-06, 0.651311E-06, 0.756543E-06, 0.873080E-06, 0.100986E-05, 0.116570E-05, 0.133938E-05, 0.154392E-05, 0.176732E-05, 0.203305E-05, 0.233242E-05, 0.266637E-05, 0.304529E-05, 0.347439E-05, 0.394457E-05, 0.447720E-05, 0.506923E-05, 0.574190E-05, 0.646482E-05, 0.726217E-05, 0.810685E-05, 0.901342E-05, 0.100228E-04, 0.110983E-04, 0.122254E-04, 0.134507E-04, 0.147612E-04, 0.160688E-04, 0.176020E-04, 0.190400E-04, 0.206282E-04, 0.222795E-04, 0.240004E-04, 0.257475E-04, 0.275228E-04, 0.294368E-04, 0.313572E-04, 0.332601E-04, 0.353490E-04, 0.374289E-04, 0.395218E-04, 0.415740E-04, 0.438212E-04, 0.461371E-04, 0.485172E-04, 0.505058E-04, 0.528814E-04, 0.554542E-04, 0.573381E-04, 0.599320E-04, 0.623185E-04, 0.647758E-04, 0.670677E-04, 0.694703E-04, 0.720178E-04, 0.744542E-04, 0.767702E-04, 0.794181E-04, 0.815966E-04, 0.841527E-04, 0.868030E-04, 0.893380E-04, 0.919178E-04, 0.940919E-04, 0.964974E-04, 0.991240E-04, 0.101811E-03, 0.103822E-03, 0.106226E-03, 0.109100E-03, 0.110907E-03, 0.113757E-03, 0.115856E-03, 0.118770E-03, 0.121550E-03, 0.123346E-03, 0.125969E-03, 0.128063E-03, 0.130828E-03, 0.133213E-03, 0.135759E-03, 0.138666E-03, 0.140271E-03, 0.142992E-03, 0.144442E-03, 0.147737E-03, 0.149106E-03, 0.152003E-03, 0.154621E-03, 0.157060E-03, 0.158249E-03, 0.161146E-03, 0.162990E-03, 0.165923E-03, 0.168044E-03, 0.170021E-03, 0.172447E-03, 0.175483E-03, 0.177206E-03, 0.179346E-03, 0.182009E-03, 0.183985E-03, 0.186670E-03, 0.187539E-03, 0.190711E-03, 0.192400E-03, 0.194753E-03, 0.197199E-03, 0.199135E-03, 0.200941E-03, 0.203020E-03, 0.206237E-03, 0.207881E-03, 0.210108E-03, 0.212492E-03, 0.214537E-03, 0.216613E-03, 0.218484E-03, 0.221548E-03, 0.221575E-03, 0.225051E-03, 0.227325E-03, 0.229743E-03, 0.232137E-03, 0.232823E-03, 0.235568E-03, 0.237281E-03, 0.240668E-03, 0.241152E-03, 0.243720E-03, 0.245138E-03, 0.246852E-03, 0.250234E-03, 0.252205E-03, 0.253325E-03, 0.255851E-03, 0.257018E-03, 0.259154E-03, 0.261575E-03, 0.264496E-03, 0.265607E-03, 0.267166E-03, 0.269612E-03, 0.272191E-03, 0.274403E-03, 0.275747E-03, 0.277492E-03, 0.279686E-03, 0.281916E-03, 0.282647E-03, 0.285925E-03, 0.287288E-03, 0.289241E-03, 0.289835E-03, 0.293543E-03, 0.295163E-03, 0.298195E-03, 0.298506E-03, 0.299592E-03, 0.300938E-03, 0.305048E-03, 0.305980E-03, 0.307956E-03, 0.310069E-03, 0.311887E-03, 0.313464E-03, 0.315334E-03, 0.316955E-03, 0.318658E-03},
{0.134244E-19, 0.318768E-19, 0.697598E-19, 0.142259E-18, 0.273843E-18, 0.501339E-18, 0.880166E-18, 0.149185E-17, 0.245122E-17, 0.392324E-17, 0.609352E-17, 0.933262E-17, 0.139056E-16, 0.204030E-16, 0.294071E-16, 0.418669E-16, 0.586794E-16, 0.811354E-16, 0.111190E-15, 0.150836E-15, 0.201332E-15, 0.267589E-15, 0.350809E-15, 0.457473E-15, 0.591571E-15, 0.760764E-15, 0.968344E-15, 0.122692E-14, 0.153801E-14, 0.192137E-14, 0.239714E-14, 0.296185E-14, 0.364328E-14, 0.446285E-14, 0.545020E-14, 0.659588E-14, 0.797077E-14, 0.959657E-14, 0.114667E-13, 0.137007E-13, 0.162255E-13, 0.192768E-13, 0.227755E-13, 0.267890E-13, 0.313985E-13, 0.367162E-13, 0.426707E-13, 0.496410E-13, 0.575148E-13, 0.667658E-13, 0.769470E-13, 0.885643E-13, 0.101504E-12, 0.116021E-12, 0.132591E-12, 0.151169E-12, 0.171707E-12, 0.195091E-12, 0.221361E-12, 0.249131E-12, 0.282453E-12, 0.317467E-12, 0.357902E-12, 0.401345E-12, 0.450314E-12, 0.503815E-12, 0.562362E-12, 0.628302E-12, 0.699751E-12, 0.776798E-12, 0.865400E-12, 0.959876E-12, 0.106192E-11, 0.117449E-11, 0.130123E-11, 0.143710E-11, 0.159048E-11, 0.174723E-11, 0.192759E-11, 0.211984E-11, 0.232502E-11, 0.255690E-11, 0.280176E-11, 0.307851E-11, 0.337059E-11, 0.369406E-11, 0.404526E-11, 0.442026E-11, 0.482718E-11, 0.528257E-11, 0.576140E-11, 0.630213E-11, 0.689650E-11, 0.752074E-11, 0.821215E-11, 0.894361E-11, 0.976631E-11, 0.106461E-10, 0.116396E-10, 0.126436E-10, 0.138115E-10, 0.150726E-10, 0.163889E-10, 0.179538E-10, 0.195073E-10, 0.213043E-10, 0.231628E-10, 0.252274E-10, 0.273938E-10, 0.298012E-10, 0.324091E-10, 0.352668E-10, 0.383767E-10, 0.415647E-10, 0.449125E-10, 0.487400E-10, 0.526346E-10, 0.568020E-10, 0.611244E-10, 0.659651E-10, 0.710162E-10, 0.763204E-10, 0.816937E-10, 0.875693E-10, 0.937412E-10, 0.100509E-09, 0.106927E-09, 0.113805E-09, 0.121548E-09, 0.129346E-09, 0.136719E-09, 0.145556E-09, 0.153745E-09, 0.162920E-09, 0.171214E-09, 0.181014E-09, 0.190636E-09, 0.200822E-09, 0.210751E-09, 0.220464E-09, 0.230720E-09, 0.241332E-09, 0.251751E-09, 0.263264E-09, 0.275173E-09, 0.286305E-09, 0.298312E-09, 0.311420E-09, 0.322419E-09, 0.334317E-09, 0.345547E-09, 0.357698E-09, 0.371438E-09, 0.384148E-09, 0.396363E-09, 0.407890E-09, 0.421661E-09, 0.434615E-09, 0.445979E-09, 0.461114E-09, 0.473851E-09, 0.489111E-09, 0.502719E-09, 0.514011E-09, 0.527900E-09, 0.539406E-09, 0.550896E-09, 0.567986E-09, 0.584061E-09, 0.593117E-09, 0.605990E-09, 0.622550E-09, 0.636040E-09, 0.648522E-09, 0.660688E-09, 0.675149E-09, 0.688309E-09, 0.702878E-09, 0.712830E-09, 0.729335E-09, 0.740855E-09, 0.753610E-09, 0.766164E-09, 0.784888E-09, 0.790491E-09, 0.808534E-09, 0.820641E-09, 0.835191E-09, 0.845560E-09, 0.858372E-09, 0.870802E-09, 0.885375E-09, 0.897990E-09, 0.909716E-09, 0.930597E-09, 0.944538E-09, 0.952220E-09, 0.963281E-09, 0.974590E-09, 0.991359E-09, 0.100148E-08},
{0.600725E-10, 0.110384E-09, 0.191055E-09, 0.314733E-09, 0.497639E-09, 0.758959E-09, 0.112407E-08, 0.162256E-08, 0.229394E-08, 0.318021E-08, 0.431863E-08, 0.580251E-08, 0.765780E-08, 0.998333E-08, 0.128858E-07, 0.164365E-07, 0.207787E-07, 0.260306E-07, 0.323941E-07, 0.400312E-07, 0.491175E-07, 0.598409E-07, 0.724052E-07, 0.872652E-07, 0.104624E-06, 0.124895E-06, 0.148402E-06, 0.175960E-06, 0.207428E-06, 0.243644E-06, 0.286167E-06, 0.334407E-06, 0.391370E-06, 0.454705E-06, 0.529830E-06, 0.613982E-06, 0.711845E-06, 0.824492E-06, 0.953811E-06, 0.110347E-05, 0.127226E-05, 0.146881E-05, 0.169303E-05, 0.194952E-05, 0.224135E-05, 0.257962E-05, 0.295480E-05, 0.337657E-05, 0.384764E-05, 0.437742E-05, 0.496572E-05, 0.560618E-05, 0.631280E-05, 0.707393E-05, 0.791481E-05, 0.880129E-05, 0.975705E-05, 0.107885E-04, 0.119199E-04, 0.130349E-04, 0.142675E-04, 0.155613E-04, 0.169283E-04, 0.183500E-04, 0.198011E-04, 0.213590E-04, 0.228934E-04, 0.245396E-04, 0.261863E-04, 0.279090E-04, 0.296989E-04, 0.315602E-04, 0.333279E-04, 0.352593E-04, 0.371457E-04, 0.391325E-04, 0.411553E-04, 0.430734E-04, 0.451483E-04, 0.473612E-04, 0.492415E-04, 0.514486E-04, 0.535933E-04, 0.556978E-04, 0.578270E-04, 0.600430E-04, 0.622587E-04, 0.643502E-04, 0.664526E-04, 0.688256E-04, 0.709245E-04, 0.730220E-04, 0.753081E-04, 0.776275E-04, 0.799978E-04, 0.821992E-04, 0.842516E-04, 0.865856E-04, 0.887716E-04, 0.909097E-04, 0.930657E-04, 0.952853E-04, 0.975147E-04, 0.995808E-04, 0.101790E-03, 0.104178E-03, 0.106885E-03, 0.108461E-03, 0.111072E-03, 0.113019E-03, 0.115394E-03, 0.117317E-03, 0.119564E-03, 0.121932E-03, 0.123859E-03, 0.126129E-03, 0.127921E-03, 0.130494E-03, 0.132355E-03, 0.134694E-03, 0.136974E-03, 0.139072E-03, 0.140844E-03, 0.143083E-03, 0.144956E-03, 0.147432E-03, 0.149570E-03, 0.151624E-03, 0.153524E-03, 0.156090E-03, 0.158029E-03, 0.159448E-03, 0.162050E-03, 0.164039E-03, 0.166727E-03, 0.167722E-03, 0.170185E-03, 0.172239E-03, 0.173889E-03, 0.176245E-03, 0.178502E-03, 0.179971E-03, 0.182601E-03, 0.184840E-03, 0.186195E-03, 0.188141E-03, 0.190848E-03, 0.192199E-03, 0.194171E-03, 0.196275E-03, 0.198842E-03, 0.199746E-03, 0.202216E-03, 0.204507E-03, 0.206104E-03, 0.208534E-03, 0.210162E-03, 0.212067E-03, 0.213750E-03, 0.216431E-03, 0.217462E-03, 0.219606E-03, 0.221777E-03, 0.223443E-03, 0.225940E-03, 0.226971E-03, 0.229214E-03, 0.231400E-03, 0.232467E-03, 0.234616E-03, 0.236465E-03, 0.239020E-03, 0.240782E-03, 0.242296E-03, 0.244521E-03, 0.246947E-03, 0.248020E-03, 0.250091E-03, 0.252189E-03, 0.253841E-03, 0.255928E-03, 0.257188E-03, 0.259523E-03, 0.260843E-03, 0.263260E-03, 0.263742E-03, 0.267186E-03, 0.268626E-03, 0.271791E-03, 0.272319E-03, 0.273437E-03, 0.275019E-03, 0.278790E-03, 0.280130E-03, 0.280314E-03, 0.282280E-03, 0.284932E-03, 0.286749E-03, 0.288937E-03, 0.290246E-03, 0.292013E-03},
{0.337420E-13, 0.614652E-13, 0.105568E-12, 0.172601E-12, 0.270793E-12, 0.410058E-12, 0.602773E-12, 0.862662E-12, 0.121115E-11, 0.166581E-11, 0.224235E-11, 0.299739E-11, 0.392255E-11, 0.507382E-11, 0.649148E-11, 0.822743E-11, 0.103142E-10, 0.128066E-10, 0.158099E-10, 0.193913E-10, 0.234983E-10, 0.284157E-10, 0.340421E-10, 0.405906E-10, 0.481258E-10, 0.568432E-10, 0.667107E-10, 0.779923E-10, 0.905503E-10, 0.104963E-09, 0.121588E-09, 0.139751E-09, 0.160334E-09, 0.183543E-09, 0.209431E-09, 0.237542E-09, 0.269666E-09, 0.305117E-09, 0.343507E-09, 0.387011E-09, 0.433346E-09, 0.486285E-09, 0.544230E-09, 0.606204E-09, 0.674271E-09, 0.750913E-09, 0.831495E-09, 0.920824E-09, 0.101970E-08, 0.112937E-08, 0.124485E-08, 0.137290E-08, 0.150813E-08, 0.165602E-08, 0.182336E-08, 0.200015E-08, 0.219165E-08, 0.240352E-08, 0.262828E-08, 0.287085E-08, 0.314779E-08, 0.343555E-08, 0.375003E-08, 0.409472E-08, 0.446625E-08, 0.487216E-08, 0.531403E-08, 0.579263E-08, 0.631253E-08, 0.687244E-08, 0.749968E-08, 0.817385E-08, 0.891271E-08, 0.972744E-08, 0.106167E-07, 0.115893E-07, 0.126951E-07, 0.138688E-07, 0.151868E-07, 0.166173E-07, 0.182288E-07, 0.199894E-07, 0.220133E-07, 0.243004E-07, 0.268207E-07, 0.296359E-07, 0.327480E-07, 0.363273E-07, 0.402928E-07, 0.449021E-07, 0.499043E-07, 0.556354E-07, 0.619799E-07, 0.689854E-07, 0.765830E-07, 0.851097E-07, 0.941103E-07, 0.103821E-06, 0.114168E-06, 0.124909E-06, 0.137352E-06, 0.149065E-06, 0.162095E-06, 0.175144E-06, 0.188999E-06, 0.202504E-06, 0.217387E-06, 0.231121E-06, 0.245521E-06, 0.260982E-06, 0.275915E-06, 0.291840E-06, 0.308383E-06, 0.323519E-06, 0.338525E-06, 0.356371E-06, 0.371073E-06, 0.387900E-06, 0.404090E-06, 0.420040E-06, 0.437985E-06, 0.452973E-06, 0.469493E-06, 0.487038E-06, 0.503525E-06, 0.521804E-06, 0.536571E-06, 0.551640E-06, 0.569362E-06, 0.586254E-06, 0.602967E-06, 0.622348E-06, 0.636384E-06, 0.652252E-06, 0.667316E-06, 0.685770E-06, 0.703557E-06, 0.720383E-06, 0.736590E-06, 0.752398E-06, 0.767344E-06, 0.784274E-06, 0.799145E-06, 0.815553E-06, 0.834815E-06, 0.851899E-06, 0.866462E-06, 0.886244E-06, 0.898278E-06, 0.914406E-06, 0.930368E-06, 0.944678E-06, 0.962294E-06, 0.977770E-06, 0.996436E-06, 0.100996E-05, 0.102695E-05, 0.104226E-05, 0.105619E-05, 0.107452E-05, 0.109002E-05, 0.110751E-05, 0.112484E-05, 0.113652E-05, 0.115338E-05, 0.116860E-05, 0.118326E-05, 0.120034E-05, 0.121975E-05, 0.122939E-05, 0.124433E-05, 0.126408E-05, 0.127989E-05, 0.129611E-05, 0.130837E-05, 0.132470E-05, 0.133945E-05, 0.135487E-05, 0.136435E-05, 0.138637E-05, 0.140146E-05, 0.141437E-05, 0.142713E-05, 0.145076E-05, 0.145356E-05, 0.147792E-05, 0.148862E-05, 0.150925E-05, 0.151714E-05, 0.153321E-05, 0.154851E-05, 0.156120E-05, 0.157899E-05, 0.159251E-05, 0.161942E-05, 0.163638E-05, 0.164000E-05, 0.165699E-05, 0.166916E-05, 0.168775E-05, 0.169985E-05},
{0.895361E-10, 0.164290E-09, 0.284361E-09, 0.468130E-09, 0.739769E-09, 0.112748E-08, 0.166787E-08, 0.240441E-08, 0.339657E-08, 0.470758E-08, 0.637608E-08, 0.858324E-08, 0.113008E-07, 0.147238E-07, 0.189635E-07, 0.242066E-07, 0.305615E-07, 0.382242E-07, 0.475739E-07, 0.588197E-07, 0.718385E-07, 0.875941E-07, 0.105801E-06, 0.127347E-06, 0.152462E-06, 0.181991E-06, 0.215769E-06, 0.255283E-06, 0.299680E-06, 0.351481E-06, 0.412712E-06, 0.481230E-06, 0.560352E-06, 0.651311E-06, 0.756543E-06, 0.873080E-06, 0.100986E-05, 0.116570E-05, 0.133938E-05, 0.154392E-05, 0.176732E-05, 0.203305E-05, 0.233242E-05, 0.266637E-05, 0.304529E-05, 0.347439E-05, 0.394457E-05, 0.447720E-05, 0.506923E-05, 0.574190E-05, 0.646482E-05, 0.726217E-05, 0.810685E-05, 0.901342E-05, 0.100228E-04, 0.110983E-04, 0.122254E-04, 0.134507E-04, 0.147612E-04, 0.160688E-04, 0.176020E-04, 0.190400E-04, 0.206282E-04, 0.222795E-04, 0.240004E-04, 0.257475E-04, 0.275228E-04, 0.294368E-04, 0.313572E-04, 0.332601E-04, 0.353490E-04, 0.374289E-04, 0.395218E-04, 0.415740E-04, 0.438212E-04, 0.461371E-04, 0.485172E-04, 0.505058E-04, 0.528814E-04, 0.554542E-04, 0.573381E-04, 0.599320E-04, 0.623185E-04, 0.647758E-04, 0.670677E-04, 0.694703E-04, 0.720178E-04, 0.744542E-04, 0.767702E-04, 0.794181E-04, 0.815966E-04, 0.841527E-04, 0.868030E-04, 0.893380E-04, 0.919178E-04, 0.940919E-04, 0.964974E-04, 0.991240E-04, 0.101811E-03, 0.103822E-03, 0.106226E-03, 0.109100E-03, 0.110907E-03, 0.113757E-03, 0.115856E-03, 0.118770E-03, 0.121550E-03, 0.123346E-03, 0.125969E-03, 0.128063E-03, 0.130828E-03, 0.133213E-03, 0.135759E-03, 0.138666E-03, 0.140271E-03, 0.142992E-03, 0.144442E-03, 0.147737E-03, 0.149106E-03, 0.152003E-03, 0.154621E-03, 0.157060E-03, 0.158249E-03, 0.161146E-03, 0.162990E-03, 0.165923E-03, 0.168044E-03, 0.170021E-03, 0.172447E-03, 0.175483E-03, 0.177206E-03, 0.179346E-03, 0.182009E-03, 0.183985E-03, 0.186670E-03, 0.187539E-03, 0.190711E-03, 0.192400E-03, 0.194753E-03, 0.197199E-03, 0.199135E-03, 0.200941E-03, 0.203020E-03, 0.206237E-03, 0.207881E-03, 0.210108E-03, 0.212492E-03, 0.214537E-03, 0.216613E-03, 0.218484E-03, 0.221548E-03, 0.221575E-03, 0.225051E-03, 0.227325E-03, 0.229743E-03, 0.232137E-03, 0.232823E-03, 0.235568E-03, 0.237281E-03, 0.240668E-03, 0.241152E-03, 0.243720E-03, 0.245138E-03, 0.246852E-03, 0.250234E-03, 0.252205E-03, 0.253325E-03, 0.255851E-03, 0.257018E-03, 0.259154E-03, 0.261575E-03, 0.264496E-03, 0.265607E-03, 0.267166E-03, 0.269612E-03, 0.272191E-03, 0.274403E-03, 0.275747E-03, 0.277492E-03, 0.279686E-03, 0.281916E-03, 0.282647E-03, 0.285925E-03, 0.287288E-03, 0.289241E-03, 0.289835E-03, 0.293543E-03, 0.295163E-03, 0.298195E-03, 0.298506E-03, 0.299592E-03, 0.300938E-03, 0.305048E-03, 0.305980E-03, 0.307956E-03, 0.310069E-03, 0.311887E-03, 0.313464E-03, 0.315334E-03, 0.316955E-03, 0.318658E-03},
{0.593278E-10, 0.970416E-10, 0.150977E-09, 0.226788E-09, 0.328167E-09, 0.462846E-09, 0.637948E-09, 0.859807E-09, 0.114025E-08, 0.148270E-08, 0.190999E-08, 0.242116E-08, 0.303981E-08, 0.377479E-08, 0.464166E-08, 0.566477E-08, 0.686285E-08, 0.825859E-08, 0.983478E-08, 0.116679E-07, 0.137778E-07, 0.161602E-07, 0.188506E-07, 0.218846E-07, 0.253250E-07, 0.291512E-07, 0.335085E-07, 0.383122E-07, 0.436356E-07, 0.496797E-07, 0.561551E-07, 0.633452E-07, 0.713194E-07, 0.801985E-07, 0.896625E-07, 0.100167E-06, 0.111739E-06, 0.124122E-06, 0.138003E-06, 0.152948E-06, 0.169564E-06, 0.186975E-06, 0.206314E-06, 0.226941E-06, 0.249411E-06, 0.273607E-06, 0.300681E-06, 0.328428E-06, 0.359883E-06, 0.392526E-06, 0.428810E-06, 0.467068E-06, 0.508290E-06, 0.553220E-06, 0.601910E-06, 0.654943E-06, 0.710453E-06, 0.772575E-06, 0.835933E-06, 0.907404E-06, 0.984971E-06, 0.106577E-05, 0.115045E-05, 0.124725E-05, 0.134776E-05, 0.146097E-05, 0.158047E-05, 0.170688E-05, 0.184769E-05, 0.200145E-05, 0.216288E-05, 0.233916E-05, 0.254058E-05, 0.275127E-05, 0.297514E-05, 0.322824E-05, 0.350221E-05, 0.380127E-05, 0.412334E-05, 0.448816E-05, 0.488556E-05, 0.531540E-05, 0.582147E-05, 0.635030E-05, 0.696437E-05, 0.761360E-05, 0.832097E-05, 0.914010E-05, 0.100257E-04, 0.109993E-04, 0.120674E-04, 0.131798E-04, 0.143841E-04, 0.156673E-04, 0.169615E-04, 0.183375E-04, 0.197084E-04, 0.211371E-04, 0.225144E-04, 0.239917E-04, 0.254514E-04, 0.268264E-04, 0.282939E-04, 0.296751E-04, 0.311240E-04, 0.324638E-04, 0.338793E-04, 0.352007E-04, 0.364780E-04, 0.379074E-04, 0.391925E-04, 0.405589E-04, 0.419043E-04, 0.431451E-04, 0.444176E-04, 0.457928E-04, 0.470493E-04, 0.482956E-04, 0.496016E-04, 0.508002E-04, 0.520946E-04, 0.532427E-04, 0.545219E-04, 0.558081E-04, 0.570305E-04, 0.583158E-04, 0.593981E-04, 0.604209E-04, 0.618086E-04, 0.629085E-04, 0.641162E-04, 0.654750E-04, 0.663210E-04, 0.676885E-04, 0.685690E-04, 0.700584E-04, 0.712939E-04, 0.724517E-04, 0.735608E-04, 0.746348E-04, 0.755313E-04, 0.769342E-04, 0.776433E-04, 0.789611E-04, 0.802502E-04, 0.813600E-04, 0.826199E-04, 0.839624E-04, 0.847976E-04, 0.858896E-04, 0.868131E-04, 0.879745E-04, 0.891388E-04, 0.902189E-04, 0.912531E-04, 0.922725E-04, 0.933789E-04, 0.945925E-04, 0.953311E-04, 0.966491E-04, 0.977256E-04, 0.990696E-04, 0.100183E-03, 0.101009E-03, 0.101883E-03, 0.103247E-03, 0.104090E-03, 0.105223E-03, 0.106653E-03, 0.107269E-03, 0.108383E-03, 0.109763E-03, 0.110750E-03, 0.111847E-03, 0.112606E-03, 0.113680E-03, 0.114977E-03, 0.116027E-03, 0.116758E-03, 0.118219E-03, 0.118963E-03, 0.120021E-03, 0.120868E-03, 0.122446E-03, 0.122728E-03, 0.124410E-03, 0.125046E-03, 0.126542E-03, 0.127041E-03, 0.128136E-03, 0.129208E-03, 0.130329E-03, 0.131196E-03, 0.132251E-03, 0.133968E-03, 0.134992E-03, 0.135547E-03, 0.136600E-03, 0.137598E-03, 0.138743E-03, 0.139729E-03},
{0.169071E-13, 0.354944E-13, 0.692783E-13, 0.127021E-12, 0.222742E-12, 0.374151E-12, 0.606034E-12, 0.951517E-12, 0.145330E-11, 0.217148E-11, 0.317569E-11, 0.455051E-11, 0.643909E-11, 0.894579E-11, 0.122416E-10, 0.165904E-10, 0.221698E-10, 0.293016E-10, 0.384616E-10, 0.498681E-10, 0.640703E-10, 0.819653E-10, 0.103905E-09, 0.131126E-09, 0.164035E-09, 0.204011E-09, 0.252537E-09, 0.310140E-09, 0.380293E-09, 0.463147E-09, 0.563069E-09, 0.680467E-09, 0.818784E-09, 0.983448E-09, 0.117577E-08, 0.140404E-08, 0.167039E-08, 0.198186E-08, 0.234435E-08, 0.277019E-08, 0.326865E-08, 0.384230E-08, 0.452547E-08, 0.530654E-08, 0.622107E-08, 0.727220E-08, 0.849427E-08, 0.992843E-08, 0.115663E-07, 0.134535E-07, 0.156487E-07, 0.181544E-07, 0.209893E-07, 0.242386E-07, 0.279046E-07, 0.320956E-07, 0.367981E-07, 0.420292E-07, 0.480229E-07, 0.546347E-07, 0.620640E-07, 0.700990E-07, 0.790009E-07, 0.891014E-07, 0.998586E-07, 0.111732E-06, 0.124588E-06, 0.138658E-06, 0.153552E-06, 0.170123E-06, 0.187568E-06, 0.206702E-06, 0.226593E-06, 0.247921E-06, 0.271125E-06, 0.295500E-06, 0.321577E-06, 0.348043E-06, 0.376683E-06, 0.406664E-06, 0.437325E-06, 0.470538E-06, 0.505022E-06, 0.540262E-06, 0.577577E-06, 0.617277E-06, 0.657367E-06, 0.696827E-06, 0.738660E-06, 0.783105E-06, 0.828162E-06, 0.872904E-06, 0.922009E-06, 0.971540E-06, 0.102263E-05, 0.107313E-05, 0.112495E-05, 0.117700E-05, 0.123066E-05, 0.128411E-05, 0.133948E-05, 0.139504E-05, 0.145387E-05, 0.151024E-05, 0.156943E-05, 0.163127E-05, 0.169627E-05, 0.175130E-05, 0.181780E-05, 0.187605E-05, 0.193955E-05, 0.200157E-05, 0.206609E-05, 0.213003E-05, 0.219207E-05, 0.225842E-05, 0.231718E-05, 0.238781E-05, 0.245349E-05, 0.251833E-05, 0.258856E-05, 0.265380E-05, 0.271604E-05, 0.278414E-05, 0.284711E-05, 0.291642E-05, 0.298733E-05, 0.305641E-05, 0.311692E-05, 0.319004E-05, 0.325987E-05, 0.331350E-05, 0.338520E-05, 0.345689E-05, 0.353544E-05, 0.357933E-05, 0.365233E-05, 0.372255E-05, 0.377854E-05, 0.385136E-05, 0.392548E-05, 0.398344E-05, 0.406118E-05, 0.413274E-05, 0.418553E-05, 0.424728E-05, 0.433038E-05, 0.438770E-05, 0.444938E-05, 0.452000E-05, 0.459645E-05, 0.464191E-05, 0.471318E-05, 0.479120E-05, 0.484803E-05, 0.492338E-05, 0.498943E-05, 0.504336E-05, 0.510947E-05, 0.518885E-05, 0.523324E-05, 0.530713E-05, 0.536982E-05, 0.543327E-05, 0.549908E-05, 0.555329E-05, 0.562063E-05, 0.569310E-05, 0.573479E-05, 0.581097E-05, 0.586811E-05, 0.595750E-05, 0.600643E-05, 0.606100E-05, 0.613368E-05, 0.620615E-05, 0.625006E-05, 0.631292E-05, 0.638386E-05, 0.644000E-05, 0.651300E-05, 0.655707E-05, 0.663009E-05, 0.668607E-05, 0.675794E-05, 0.677886E-05, 0.687369E-05, 0.692797E-05, 0.702162E-05, 0.704769E-05, 0.708771E-05, 0.714437E-05, 0.727460E-05, 0.730114E-05, 0.732280E-05, 0.737679E-05, 0.747082E-05, 0.752269E-05, 0.758817E-05, 0.764741E-05, 0.770027E-05},
{0.591066E-16, 0.123297E-15, 0.239109E-15, 0.437500E-15, 0.760393E-15, 0.126731E-14, 0.203928E-14, 0.317949E-14, 0.482320E-14, 0.715160E-14, 0.103559E-13, 0.147799E-13, 0.206286E-13, 0.284218E-13, 0.384984E-13, 0.517031E-13, 0.684285E-13, 0.896078E-13, 0.116182E-12, 0.149435E-12, 0.189649E-12, 0.239864E-12, 0.299850E-12, 0.372594E-12, 0.459844E-12, 0.565253E-12, 0.689818E-12, 0.837158E-12, 0.100559E-11, 0.120724E-11, 0.144394E-11, 0.171528E-11, 0.202893E-11, 0.239340E-11, 0.281413E-11, 0.328098E-11, 0.383079E-11, 0.444270E-11, 0.513252E-11, 0.593729E-11, 0.680256E-11, 0.780863E-11, 0.894165E-11, 0.101885E-10, 0.115918E-10, 0.131330E-10, 0.148438E-10, 0.167545E-10, 0.188823E-10, 0.212820E-10, 0.238939E-10, 0.267710E-10, 0.298731E-10, 0.332909E-10, 0.370683E-10, 0.412353E-10, 0.457468E-10, 0.507163E-10, 0.561999E-10, 0.619926E-10, 0.687648E-10, 0.754944E-10, 0.832687E-10, 0.916081E-10, 0.100656E-09, 0.110526E-09, 0.120978E-09, 0.132639E-09, 0.145109E-09, 0.158449E-09, 0.173223E-09, 0.189094E-09, 0.206231E-09, 0.224366E-09, 0.244640E-09, 0.267011E-09, 0.290860E-09, 0.315637E-09, 0.343495E-09, 0.372740E-09, 0.404866E-09, 0.440085E-09, 0.477794E-09, 0.519204E-09, 0.564162E-09, 0.612885E-09, 0.665273E-09, 0.721865E-09, 0.783957E-09, 0.852290E-09, 0.926262E-09, 0.100704E-08, 0.109701E-08, 0.119186E-08, 0.129659E-08, 0.140966E-08, 0.153349E-08, 0.166775E-08, 0.181549E-08, 0.196904E-08, 0.213791E-08, 0.232126E-08, 0.251356E-08, 0.273181E-08, 0.294747E-08, 0.318840E-08, 0.343208E-08, 0.369874E-08, 0.396384E-08, 0.425966E-08, 0.456685E-08, 0.489451E-08, 0.523020E-08, 0.556193E-08, 0.591111E-08, 0.629109E-08, 0.666175E-08, 0.704470E-08, 0.743439E-08, 0.784240E-08, 0.825771E-08, 0.866611E-08, 0.910038E-08, 0.954771E-08, 0.997965E-08, 0.104653E-07, 0.108841E-07, 0.113169E-07, 0.118306E-07, 0.123082E-07, 0.127420E-07, 0.132729E-07, 0.137016E-07, 0.142422E-07, 0.146489E-07, 0.151915E-07, 0.156794E-07, 0.161964E-07, 0.166753E-07, 0.171022E-07, 0.176117E-07, 0.181028E-07, 0.185093E-07, 0.190292E-07, 0.195601E-07, 0.200080E-07, 0.205511E-07, 0.211184E-07, 0.215429E-07, 0.219985E-07, 0.224032E-07, 0.229194E-07, 0.234651E-07, 0.239454E-07, 0.243668E-07, 0.247823E-07, 0.253124E-07, 0.257857E-07, 0.261490E-07, 0.267162E-07, 0.271848E-07, 0.278201E-07, 0.283059E-07, 0.286620E-07, 0.290631E-07, 0.294490E-07, 0.297914E-07, 0.304711E-07, 0.310556E-07, 0.313303E-07, 0.317139E-07, 0.323124E-07, 0.327883E-07, 0.331653E-07, 0.335563E-07, 0.340130E-07, 0.344544E-07, 0.349893E-07, 0.352391E-07, 0.357908E-07, 0.361015E-07, 0.365365E-07, 0.369451E-07, 0.376271E-07, 0.376833E-07, 0.382993E-07, 0.387288E-07, 0.391586E-07, 0.394640E-07, 0.398661E-07, 0.402736E-07, 0.407951E-07, 0.411340E-07, 0.414629E-07, 0.422172E-07, 0.426527E-07, 0.428572E-07, 0.431840E-07, 0.434991E-07, 0.441069E-07, 0.443674E-07}};
/*
	if(m3Pi< ms[0]){
		std::cerr << "tabulated_integrals.h: Error: Phase space not defined at m3Pi = " << m3Pi << std::endl;
		return 0.;
	};
	int length = sizeof(ms)/sizeof(*ms);
	if(m3Pi >= ms[length-1]){
		std::cerr << "tabulated_integrals.h: Error: Phase space not defined at m3Pi = " << m3Pi << std::endl;
		return 0.;
	};
	int pos;
	double mmax;
	double mmin;
	for (int i =1; i< length; i++){
		if(ms[i] > m3Pi){ // Position (*)
			pos = i;
			mmin=ms[i-1];
			mmax=ms[i];
			break;
		};
	};
	double frac = (m3Pi - mmin)/(mmax-mmin);
*/ //Directly calculate the position rather than looping through the list
	double mMin = 0.500;
	double mMax = 2.500;
	double step = 0.010;
	if (m3Pi < mMin or m3Pi > mMax){
		std::cerr << "tabulated_integrals.h: Error: Phase space not defined at m3Pi = " << m3Pi << std::endl;
		return 0.;
	};
	int pos = (m3Pi - mMin)/step + 1.; // The +1. is due to historical reasons at position (*)
	double frac = (m3Pi - mMin-(pos-1)*step)/step;


	return pow(ints[wave_number][pos-1]*(1-frac) + ints[wave_number][pos]*frac,.5);
};
#endif//TABULATED_INTEGRALS_TOP_LOL	