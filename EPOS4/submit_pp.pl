#!/usr/bin/perl -w

use Cwd;

$package = "pp13TeV";
$maindir = getcwd();

$groupnum = 119;

$rundir = "${maindir}/run_${package}_grp${groupnum}";
mkdir $rundir;

for ($irun=9000; $irun<10000; $irun++){

#	sleep 2;

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Never\n";
#	print FILE "Requirements = CPU_Speed>=1\n";
#	print FILE "Rank = CPU_Speed\n";
	print FILE "Getenv = true\n";
	print FILE "Priority = +20\n";
	print FILE "Executable = jobscript\n";
	print FILE "Log = jobscript.log\n";
	print FILE "Output = jobscript.out\n";
	print FILE "Error = jobscript.err\n";
	print FILE "request_memory = 4G\n";
#	print FILE "Notify_user = shlim\@rcf.rhic.bnl.gov\n";
#	print FILE "+Experiment = \"phenix\"\n";
#	print FILE "+Job_Type = \"cas\"\n";
	print FILE "Queue\n";
	close(FILE);

	$seednum1 = int(rand(100000000));
	$seednum2 = int(rand(100000000));

	$grpdir = sprintf("%s/run_%s_grp%03d",$maindir,$package,$groupnum);
	$outputdir1 = sprintf("%s/out_%s_grp%03d",$maindir,$package,$groupnum);
	$outputdir2 = sprintf("%s/epos4_%s_grp%03d",$maindir,$package,$groupnum);

	$randnum = int(rand(60));

	open(FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
	print FILE "sleep ${randnum}; source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/sphenix_mdc/setup.csh\n\n\n";
	print FILE "source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/06.EPOS4/setup.csh\n";

	print FILE "mkdir -p $grpdir\n\n";
	print FILE "mkdir -p $outputdir1\n\n";
	print FILE "mkdir -p $outputdir2\n\n";

#	$condor_wrk_dir = "grp${groupnum}_run${irun}";
	$condor_wrk_dir = "grp${groupnum}_run${irun}";
	print FILE "mkdir -p \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";
	print FILE "cd \$_CONDOR_SCRATCH_DIR/${condor_wrk_dir}\n";

	print FILE "cp -v /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/Model/EPOS4/ScanEPOS02.C .\n";
	print FILE "cp -v /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/Model/EPOS4/${package}.optns .\n";

	print FILE "\$EPO/scripts/epos -root ${package} ${seednum1} ${seednum2}\n";

	print FILE "ls *.root > file.lst\n";
	print FILE "root -l -q ScanEPOS02.C\n";

#	print FILE "root -l -b -q 'phpythia8_sQCD.C(1500,$seednum)'\n";
#	print FILE "root -l -b -q 'phpythia8_sQCD.C(500000,$seednum)'\n";
#	print FILE "root -l -b -q 'Run_fill_ana_tree.C(\"pythia_mn.root\")'\n";

	$outputfile = sprintf("%s/outfile_%s_grp%03d_%05d.root",$outputdir1,$package,$groupnum,$irun);
	print FILE "mv outfile_epos4.root $outputfile\n";

#	$outputfile = sprintf("%s/epos4_%s_grp%03d_%05d.root",$outputdir2,$package,$groupnum,$irun);
#	print FILE "mv *.root $outputfile\n";

	print FILE "rm -rf *\n\n";

	close(FILE);
	chmod 0755, "jobscript";
	system "condor_submit condor";
}

