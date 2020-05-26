#!/usr/bin/perl -w

use Cwd;

$package = "pp13TeV_set02";
$maindir = getcwd();

$srcdir = "/alice/home/shlim/work/Pythia/codes";
$datadir = "/alice/data/shlim";

$groupnum = 1;
$try = 100;

$filedir = sprintf("${datadir}/${package}_grp%03d",${groupnum});

$files_per_job = 50; #pp13TeV pythia (string shoving), 3500 files

$rundir = sprintf("${maindir}/runhist_${package}_grp%03d_try%03d",${groupnum},${try});
mkdir $rundir;

$jobname = sprintf("hist_${package}_grp%03d_try%03d",${groupnum},${try});

for ($irun=0; $irun<100; $irun++){

	$wrkdir = "${rundir}/wrk_${irun}";
	mkdir $wrkdir;

	chdir $wrkdir;
	open(FILE, ">condor");
	print FILE "Universe = vanilla\n";
	print FILE "Notification = Never\n";
#	print FILE "Requirements = CPU_Speed>=2\n";
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

	open(FILE, ">jobscript");
	print FILE "#!/bin/bash\n";
	print FILE "source /pool/kiafenv\n\n\n";

#	print FILE "cp -av ${srcdir}/make_hist_pp13TeV_02.C .\n";
#	print FILE "root -l -b -q 'make_hist_pp13TeV_02.C+g(\"file.lst\")'\n";

	print FILE "cp -av ${srcdir}/make_hist_pp13TeV_03.C .\n";
	print FILE "root -l -b -q 'make_hist_pp13TeV_03.C+g(\"file.lst\")'\n";

	$filename = sprintf("%s/outfile_%s_grp%03d_%05d.root",$rundir,$package,$groupnum,$irun);
	print FILE "mv -v outfile_hist.root $filename\n\n";

	print FILE "rm -rf *.C\n\n";

	close(FILE);
	chmod 0755, "jobscript";

	open(FILE,">file.lst");
	for ($iseg=0; $iseg<$files_per_job; $iseg++){
		$filename = sprintf("%s/outfile_%s_grp%03d_%05d.root",$filedir,$package,$groupnum,$files_per_job*$irun+$iseg);
		if ( -e $filename ){
			print FILE $filename."\n";
		}
	}
	close(FILE);

#	system "condor_submit -append \"Accounting_Group=group_alice\" condor";

}
