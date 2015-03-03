#!/usr/bin/perl -w
use Cwd 'abs_path';
use Term::ReadKey;

if($#ARGV < 0) {
    print"ERROR: First argument has to be path of ggH folder!
Script stopped...\n";
    die;
}
$path = $ARGV[0];

unless(-d $path){
    print"ERROR: folder \"$path\" of ggH output does not exist.
Script stopped...\n";
    die;
}

$current_path = Cwd::cwd();

# FeynHiggs:
if($#ARGV > 1) {
    $fhpath = $ARGV[2];
}
else{
    system("curl \"http://wwwth.mpp.mpg.de/members/heinemey/feynhiggs/cFeynHiggs.html\" -o \"cFeynHiggs.html\" 2> $current_path/FH_curl.log");
    if($? == 0) {
	$line = `grep '\\.tar\\.gz' cFeynHiggs.html | head -1`;
	chomp ($line);
	system("rm cFeynHiggs.html");
	$start=index($line,"newversion/FeynHiggs");
	$end=index($line,".tar.gz\">FeynHiggs");
	$filename=substr($line,$start+11,$end-4-$start);
	$fhversion=substr($line,$start+11,$end-11-$start);
	print"\nDo you want to download ".$fhversion." now (y/n)?\n";
	ReadMode 4;
	while (not defined ($inp = ReadKey(-1))) {
	}
	ReadMode 0;
	if($inp =~ /y/i) {
	    while (1 == 1) {
		print"Type absolute path to the folder in which FeynHiggs should be installed (e.g. /home/username, no worries a new FeynHiggs-folder will be created under this path!). Press ENTER if you would like to install FeynHiggs in the ggH folder.\n";
		$fhpath =  <STDIN>;
		chomp ($fhpath);
		if($fhpath eq "") {
		    $fhpath = $path;
		    mkdir($fhpath);
		    last;
		}
		if(-d $fhpath) {
		    last;
		} else
		{
		    print"\nError: Path does not exist.\n\n";
		}
	    }
	    chdir($fhpath) or die "$!";
	    print "\n**** Downloading FeynHiggs ****\n";
	    system("curl \"http://wwwth.mpp.mpg.de/members/heinemey/feynhiggs/newversion/$filename\" -o \"$filename\" 2>> $current_path/FH_curl.log");
	    if($? != 0) {
		print "Download error. Please check FH_curl.log for more information.\n";
		die;
	    }
	    print "Download complete.\n";
	    print "\n**** Untaring FeynHiggs tarball ****\n";
	    system("tar zxf $filename");
	    if($? != 0) {
		print "Error while untaring $filename.\n";
		die;
	    }
	    print "\n**** Compiling FeynHiggs ****\n";
	    chdir($fhversion) or die "$!"; 
	    system("./configure > $current_path/FH_conf.log");
	    if($? != 0) {
		print "Error while running configure. Please check FH_conf.log for more information.\n";
		die;
	    }
	    print "This may take some time.\n";	    
	    system("make > $current_path/FH_make.log 2>> $current_path/FH_make.log");
	    if($? != 0) {
		print "Error while compiling FeynHiggs. Please check FH_make.log for more information.\n";
	    } else {
		print "FeynHiggs has been successfully installed.\n\n";
	    }
	    system("make install >> $current_path/FH_make.log");
	    $fhpath = $fhpath."/".$fhversion."/build";
	} else {	    
	    while (1 == 1) {
		print"Type absolute path to lib folder of FeynHiggs (e.g. /home/username/FeynHiggs-X.X.X/x86_64-Darwin/lib). Press ENTER if you would like to link libFH.a by hand.\n";
		$fhpath =  <STDIN>;
		chomp ($fhpath);
		if($fhpath eq "") {
		    last;
		}
		if(-d $fhpath) {
		    last;
		} else
		{
		    print"\nError: Path does not exist.\n\n";
		}
	    }
	}
    } else {
	while (1 == 1) {
	    print"Downloading FeynHiggs is not possible at the moment. Type absolute path to lib folder of FeynHiggs (e.g. /home/username/FeynHiggs-X.X.X/x86_64-Darwin/lib). Press ENTER if you would like to link libFH.a by hand.\n";
	    $fhpath =  <STDIN>;
	    chomp ($fhpath);
	    if($fhpath eq "") {
		last;
	    }
	    if(-d $fhpath) {
		last;
	    } else
	    {
		print"\nError: Path does not exist.\n\n";
	    }
	}
    }
    chomp ($fhpath);
}
# SusHi:
if($#ARGV > 0) {
    $sushipath = $ARGV[1];
}
else{
    system("curl \"https://www.hepforge.org/archive/sushi/\" -o \"sushi\" 2> $current_path/SusHi_curl.log");
    if($? == 0) {
	$line = `grep 'SusHi' sushi | tail -1`;
	chomp ($line);
	system("rm sushi");
	$start=index($line,"SusHi-");
	$end=index($line,".tar.gz");
	$filename=substr($line,$start,$end+7-$start);
	$sushiversion=substr($line,$start,$end-$start);
	print"\nDo you want to download ".$sushiversion." now (y/n)?\n";
	ReadMode 4;
	while (not defined ($inp = ReadKey(-1))) {
	}
	ReadMode 0;
	if($inp =~ /y/i) {
	    chdir($current_path) or die "$!"; 
	    while (1 == 1) {
		print"Type absolute path to the folder in which SusHi should be installed (e.g. /home/username, no worries a new SusHi-folder will be created under this path!). Press ENTER if you would like to install SusHi in the ggH folder.\n";
		$sushipath =  <STDIN>;
		chomp ($sushipath);
		if($sushipath eq "") {
		    $sushipath = $path;
		    last;
		}
		if(-d $sushipath) {
		    last;
		} else
		{
		    print"\nError: Path does not exist.\n\n";
		}
	    }
	    chdir($sushipath) or die "$!";
	    print "\n**** Downloading SusHi ****\n";
	    system("curl \"https://www.hepforge.org/archive/sushi/$filename\" -o \"$filename\" 2>> $current_path/SusHi_curl.log");
	    if($? != 0) {
		print "Download error. Please check SusHi_curl.log for more information.\n";
		die;
	    }
	    print "Download complete.\n";
	    print "\n**** Untaring SusHi tarball ****\n";
	    system("tar zxf $filename");
	    if($? != 0) {
		print "Error while untaring $filename.\n";
		die;
	    }
	    chdir($sushiversion) or die "$!";
	    print("\nWould you like to add flag \'F77FLAGS=-mcmodel=large\' to compilation (y/n)?\n--> Solves relocation errors on some clusters, but leads to compilation errors on some Mac OS X systems.\n");
	    ReadMode 4;
	    while (not defined ($inp = ReadKey(-1))) {
	    }
	    ReadMode 0;
	    if($inp =~ /y/i) {
		$flag="F77FLAGS=-mcmodel=large";
	    }
	    else
	    {
		$flag=""
	    }
	    print "\n**** Compiling SusHi ****\n";
	    system("./configure $flag > $current_path/SusHi_conf.log");
	    if($? != 0) {
		print "Error while running configure. Please check SusHi_conf.log for more information.\n";
		die;
	    }
	    print "This may take some time.\n";
	    system("make predef=FH > $current_path/SusHi_make.log 2>> $current_path/SusHi_make.log");
	    if($? != 0) {
		print "Error while compiling SusHi. Please check SusHi_make.log for more information.\n";
	    } else {
		print "SusHi has been successfully installed.\n\n";
	    }
	    $sushipath = $sushipath."/".$sushiversion."/lib";
	} else {
	    while (1 == 1) {
		print"Type absolute path to lib folder of SusHi (e.g. /home/username/Sushi-X.X.X/lib). Press ENTER if you would like to link libsushiFH.a by hand.\n";
		$sushipath =  <STDIN>;
		chomp ($sushipath);
		if($sushipath eq "") {
		    last;
		}
		if(-d $sushipath) {
		    last;
		} else
		{
		    print"\nError: Path does not exist.\n\n";
		}
	    }
	}
    } else {
	while (1 == 1) {
	    print"Downloading SusHi is not possible at the moment. Please type absolute path to lib folder of SusHi (e.g. /home/username/Sushi-X.X.X/lib). Press ENTER if you would like to link libsushiFH.a by hand.\n";
	    $sushipath =  <STDIN>;
	    chomp ($sushipath);
	    if($sushipath eq "") {
		last;
	    }
	    if(-d $sushipath) {
		last;
	    } else
	    {
		print"\nError: Path does not exist.\n\n";
	    }
	}
    }
    chomp ($sushipath);
}

chdir($current_path) or die "$!"; 

$compare="\n";
chomp($compare);

if(-d $sushipath){
    unless($sushipath =~ /^\s*\//){
	$sushipath = $current_path."/".$sushipath;
    }
    if(-f $sushipath."/libsushiFH.a"){
	$sushilib = $sushipath."/libsushiFH.a";
	$libpath  = $path."/lib";
	$linkedsushi = $libpath."/libsushiFH.a";
	if(-f $linkedsushi){
	    system("rm $linkedsushi");
	}
	system("ln -ns $sushilib $libpath");
	print"INFO: libsushiFH.a linked from $sushipath\n";
    }
    else{
	print"WARNING: libsushiFH.a not found in path to lib folder of SusHi: $sushipath. You have to link libsushiFH.a by hand to the lib folder of your ggH output!\n"
    }
}
elsif($sushipath eq ""){
    print"INFO: You have to link libsushiFH.a by hand to the lib folder of your ggH output!\n";
}
else{
    print"WARNING: folder of SusHi $sushipath does not exist. You have to link libsushiFH.a by hand to the lib folder of your ggH output!\n"
}

if(-d $fhpath){
    unless($fhpath =~ /^\s*\//){
	$fhpath = $current_path."/".$fhpath;
    }
    if(-f $fhpath."/libFH.a"){
	$fhlib = $fhpath."/libFH.a";
	$libpath  = $path."/lib";
	$linkedFH = $libpath."/libFH.a";
	if(-f $linkedFH){
	    system("rm $linkedFH");
	}
	system("ln -ns $fhlib $libpath");
	print"INFO: libFH.a linked from $fhpath\n";
    }
    else{
	print"WARNING: libFH.a not found in path to lib folder of FeynHiggs: $fhpath. You have to link libFH.a by hand to the lib folder of your ggH output!\n";
    }
}
elsif($fhpath eq ""){
    print"INFO: You have to link libFH.a by hand to the lib folder of your ggH output!\n";
}
else{
    print"WARNING: folder of FeynHiggs $fhpath does not exist. You have to link libFH.a by hand to the lib folder of your ggH output!\n"
}


$script_path = abs_path($0);
$gghpath=substr($script_path, 0, -26);
system("cp -r $gghpath/* $path\n");
system("rm $path/set_up_ggH_MSSM_script.pl");
system("rm $path/GGH_README");
print"INFO: modified ggH files copied from $gghpath/* to $path\n";


