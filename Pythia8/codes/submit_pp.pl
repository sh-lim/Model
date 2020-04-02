#!/usr/bin/perl -w

use Cwd;

$package = "pp13TeV_set00";
$maindir = getcwd();
$srcdir = "/alice/home/shlim/work/Pythia/codes";
$datadir = "/alice/data/shlim";

$groupnum = 1;
$progname = "mainEx00";

$rundir = "${maindir}/run_${package}_grp${groupnum}";
mkdir $rundir;

$jobname = sprintf("${package}_grp%03d",${groupnum});

for ($irun=0; $irun<5000; $irun++){

#	sleep 1;

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
	print FILE "JobBatchName = ${jobname}\n";
	print FILE "Log = jobscript.log\n";
	print FILE "Output = jobscript.out\n";
	print FILE "Error = jobscript.err\n";
#	print FILE "Notify_user = shlim\@sdfarm.kr\n";
	print FILE "Queue\n";
	close(FILE);

	$seed = $irun+1;

	open(FILE, ">${progname}.cfg");
	print FILE "Main:numberOfEvents = 50000\n";
	print FILE "Main:timesAllowErrors = 3\n";
	print FILE "Random:setSeed = on\n";
	print FILE "Random:seed = ${seed}\n";
	print FILE "Init:showChangedSettings = on\n";
	print FILE "Init:showChangedParticleData = on\n";
	print FILE "Next:numberCount = 1000\n";
	print FILE "Next:numberShowInfo = 1\n";
	print FILE "Next:numberShowProcess = 1\n";
	print FILE "Next:numberShowEvent = 1\n";
	print FILE "Beams:idA = 2212\n";
	print FILE "Beams:idB = 2212\n";
	print FILE "Beams:eCM = 14000.\n";
	print FILE "SoftQCD:nonDiffractive = on\n";
#	print FILE "SoftQCD:inelastic = on\n";
#	print FILE "HardQCD:all = on\n";
#	print FILE "PhaseSpace:pTHatMin = 10.\n";
	print FILE "Tune:pp = 21\n";
#	print FILE "PartonLevel:MPI = off\n";
#	print FILE "PartonLevel:ISR = off\n";
#	print FILE "PartonLevel:FSR = off\n";
	print FILE "HadronLevel:Decay = off\n";

	print FILE "Fragmentation:setVertices = on\n";
	print FILE "PartonVertex:setVertex = on\n";
#	print FILE "PartonVertex:protonRadius = 0.7\n";
#	print FILE "PartonVertex:emissionWidth = 0.1\n";
#	print FILE "ParticleDecays:limitTau0 = on\n";
#	print FILE "ParticleDecays:tau0Max = 10\n";
#
	print FILE "Ropewalk:RopeHadronization = on\n";
	print FILE "Ropewalk:doShoving = on\n";
	print FILE "Ropewalk:doFlavour = on\n";
	close(FILE);

#	$seednum = int(rand(100000000));

	$grpdir = sprintf("%s/%s_grp%03d",$datadir,$package,$groupnum);

	$randnum = int(rand(60));

	open(FILE, ">jobscript");
	print FILE "#!/bin/bash\n";
	print FILE "sleep ${randnum}; source /pool/kiafenv\n\n\n";

	print FILE "mkdir -p $grpdir\n\n";

	print FILE "cp -av ${srcdir}/${progname} .\n";

	print FILE "./${progname}\n";

	$outputfile = sprintf("%s/outfile_%s_grp%03d_%05d.root",$grpdir,$package,$groupnum,$irun);
	print FILE "mv Pythia8_event.root $outputfile\n";

	print FILE "rm -rf ${progname}\n\n";

	close(FILE);
	chmod 0755, "jobscript";
	system "condor_submit -append \"Accounting_Group=group_alice\" condor";
}

