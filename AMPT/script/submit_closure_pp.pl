#!/usr/bin/perl -w

use Cwd;

$package = "pp_closure";
$maindir = getcwd();

$groupnum = 25;

$rundir = "${maindir}/running_pp_grp${groupnum}";
mkdir $rundir;

for ($irun=8000; $irun<9000; $irun++){

	sleep 1;

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Error\n";
	print FILE "Requirements = CPU_Speed>=1\n";
	print FILE "Rank = CPU_Speed\n";
	print FILE "Getenv = true\n";
	print FILE "Priority = +20\n";
#	print FILE "request_memory = 4G\n";
	print FILE "Executable = jobscript\n";
	print FILE "Log = jobscript.log\n";
	print FILE "Output = jobscript.out\n";
	print FILE "Error = jobscript.err\n";
	print FILE "Notify_user = shlim\@rcf.rhic.bnl.gov\n";
#	print FILE "+Experiment = \"phenix\"\n";
#	print FILE "+Job_Type = \"cas\"\n";
	print FILE "Queue\n";
	close(FILE);

	$seednum = int(rand(100000000));

	open(FILE, ">input.ampt");
	print FILE "200\n"; # EFRM (sqrt(S_NN) in GeV if FRAME is CMS)
	print FILE "CMS\n"; # FRAME
	print FILE "P\n"; # PROJ
	print FILE "P\n"; # TARG
	print FILE "1\n"; # IAP (projectile A number)
	print FILE "1\n"; # IZP (projectile Z number)
	print FILE "1\n"; # IAT (target A number) Au: 197 Pb: 208
	print FILE "1\n"; # IZT (target Z number) Au: 79 Pb: 82
	print FILE "20000\n"; # NEVNT (total number of events
#	print FILE "100\n"; # NEVNT (total number of events
	print FILE "0.0\n"; # BMIN (mininum impact parameter in fm) 
	print FILE "15.0\n"; # BMAX (maximum impact parameter in fm, also see below)
	print FILE "4\n"; # ISOFT (D=1): select Default AMPT or String Melting(see below)
	print FILE "150\n"; # NTMAX: number of timesteps (D=150), (D=3 off cascade) 
#	print FILE "3\n"; # NTMAX: number of timesteps (D=150), (D=3 off cascade) 
	print FILE "0.2\n"; # DT: timestep in fm (hadron cascade time= DT*NTMAX) (D=0.2)
	print FILE "0.55\n"; # PARJ(41): parameter a in Lund symmetric splitting function   d+Au: 2.2  p+Pb: 0.5
	print FILE "0.15\n"; # PARJ(42): parameter b in Lund symmetric splitting function   d+Au: 0.5  p+Pb: 0.9
#	print FILE "2.2\n"; # PARJ(41): parameter a in Lund symmetric splitting function   d+Au: 2.2  p+Pb: 0.5
#	print FILE "0.5\n"; # PARJ(42): parameter b in Lund symmetric splitting function   d+Au: 0.5  p+Pb: 0.9
	print FILE "1\n"; # (D=1,yes;0,no) flag for popcorn mechanism(netbaryon stopping)
	print FILE "1.0\n"; # PARJ(5) to control BMBbar vs BBbar in popcorn (D=1.0)
	print FILE "1\n"; # shadowing flag (Default=1,yes; 0,no)
	print FILE "0\n"; # quenching flag (D=0,no; 1,yes)
	print FILE "1.0\n"; # quenching parameter -dE/dx (GeV/fm) in case quenching flag=1
	print FILE "2.0\n"; # p0 cutoff in HIJING for minijet productions (D=2.0)
	print FILE "3.2264d0\n"; # parton screening mass fm^(-1) (D=2.265d0) 0mb:5588.2820d0  0.75mb:6.4528d0  1.5mb:4.5628d0  3mb:3.2264d0   p+Pb: 3.2d0
#	print FILE "6.4528d0\n"; # parton screening mass fm^(-1) (D=2.265d0) 0mb:5588.2820d0  0.75mb:6.4528d0  1.5mb:4.5628d0  3mb:3.2264d0   p+Pb: 3.2d0
#	print FILE "5588.2820d0\n"; # parton screening mass fm^(-1) (D=2.265d0) 0mb:5588.2820d0  0.75mb:6.4528d0  1.5mb:4.5628d0  3mb:3.2264d0   p+Pb: 3.2d0
	print FILE "0\n"; # IZPC: (D=0 forward-angle parton scatterings; 100,isotropic)
	print FILE "0.47140452d0\n"; # alpha in parton cascade (D=0.33d0), see parton screening mass  d+Au: 0.47d0  p+Pb: 0.33d0
#	print FILE "0.0d0\n"; # alpha in parton cascade (D=0.33d0), see parton screening mass  d+Au: 0.47d0  p+Pb: 0.33d0
	print FILE "1d6\n"; # dpcoal in GeV
	print FILE "1d6\n"; # drcoal in fm
	print FILE "11\n"; # ihjsed: take HIJING seed from below (D=0)or at runtime(11)
	print FILE "${seednum}\n"; # random seed for HIJING
	print FILE "${seednum}\n"; # random seed for parton cascade
	print FILE "0\n"; # flag for K0s weak decays (D=0,no; 1,yes)
	print FILE "1\n"; # flag for phi decays at end of hadron cascade (D=1,yes; 0,no)
	print FILE "0\n"; # flag for pi0 decays at end of hadron cascade (D=0,no; 1,yes)
	print FILE "3\n"; # optional OSCAR output (D=0,no; 1,yes; 2&3,more parton info)
	print FILE "0\n"; # flag for perturbative deuteron calculation (D=0,no; 1or2,yes)
	print FILE "1\n"; # integer factor for perturbative deuterons(>=1 & <=10000)
	print FILE "1\n"; # choice of cross section assumptions for deuteron reactions
	print FILE "-7.\n"; # Pt in GeV: generate events with >=1 minijet above this value
	print FILE "1000\n"; # maxmiss (D=1000): maximum # of tries to repeat a HIJING event
	print FILE "3\n"; # flag to turn off initial and final state radiation (D=3)
	print FILE "1\n"; # flag to turn off Kt kick (D=1)
	print FILE "0\n"; # flag to turn on quark pair embedding (D=0,no; 1,yes)
	print FILE "7., 0.\n"; # Initial Px and Py values (GeV) of the embedded quark (u or d)
	print FILE "0., 0.\n"; # Initial x & y values (fm) of the embedded back-to-back q/qbar
	print FILE "1, 5., 0.\n"; # nsembd(D=0), psembd (in GeV),tmaxembd (in radian).
	print FILE "0\n"; # Flag to enable users to modify shadowing (D=0,no; 1,yes)
	print FILE "1.d0\n"; # Factor used to modify nuclear shadowing
	print FILE "1\n"; # Flag for random orientation of reaction plane (D=0,no; 1,yes)
#	print FILE "1\n"; # Flag to implement black-disk initial conditions w/o hard processes (D=0, no; 1, yes) 
	close(FILE);

#	$grpdir = sprintf("%s/%s_grp%03d",$maindir,$package,$groupnum);
	$grpdir = sprintf("%s/%s_grp%03d","/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/shlim/simulation/ClosureSample",$package,$groupnum);

	$randnum = int(rand(300));

	open(FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
	print FILE "sleep ${randnum}; source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/setup.csh\n\n\n";

	print FILE "mkdir -p $grpdir\n\n";

	$condor_wrk_dir = "grp${groupnum}_run${irun}";
	print FILE "mkdir -p \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";
	print FILE "cd \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";

#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/jdok/AMPT_GeoScan_0618/ampt_v2.26t5_dau200/* .\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/01.AMPT/source/Ampt-v1.26t5-v2.26t5_all_systems/* .\n";
	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/01.AMPT/source/Ampt-v1.26t9b-v2.26t9b/* .\n";
	print FILE "rm -rf ana/*\n";
	print FILE "cp -av ${wrkdir}/input.ampt .\n";
	print FILE "./exec ${seednum}\n";
#	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/13.hadron_Rcp/AMPT/FillHisto.C .\n";
	print FILE "cp -av /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/01.AMPT/source/AnaCode/Run.C .\n";
#	print FILE "root -l -b -q 'FillHisto.C(\"ana/ampt.dat\")'\n";
	print FILE "root -l -b -q 'Run.C(\"ana/ampt.dat\")'\n";

	$outputfile = sprintf("%s/AMPT_outfile_%s_grp%03d_%05d.root",$grpdir,$package,$groupnum,$irun);
	print FILE "mv outfile.root $outputfile\n";

	print FILE "rm -rf *\n\n";

	close(FILE);
	chmod 0755, "jobscript";
	system "condor_submit condor";
}

