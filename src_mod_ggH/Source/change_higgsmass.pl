#!/usr/bin/perl -w
use Cwd 'abs_path';

if($#ARGV < 0) {die}
$mH=$ARGV[0];

$script_path = abs_path($0);
#print"$script_path\n";

$path=substr($script_path, 0, -26)."Cards/param_card.dat";
#print"$path\n";

changeparam($path,'mass',25,$mH);

#- {{{ sub changeparam:

sub changeparam {
#
#  changeparam($file,$block,$entry,$val)
#
#  In Block $block of an SLHA input file, change the line
#  $entry  <somevalue>
#  to
#  $entry  $val
#
#  Example:   changeparam('slha.in','MASS',25,120)
#  This will set the light Higgs mass to 120 GeV.
#
    my($file,$block,$entry,$val);
    my($failed);
    ($file,$block,$entry,$val) = @_;
    movefile("$file","$file.bak");
    open(FILEIN,"$file.bak") || die;
    open(FILEOUT,">$file") || die "Cannot open file $file.";
    $failed = 1;
  LOOP: while (<FILEIN>) {
      if ((/^Block +$block /i) || (/^Block +$block$/i)) {
	  $failed = 0;
	  print {FILEOUT} ($_);
	  while (<FILEIN>) {
	      if (/^[BD]/i) { 
		  print {FILEOUT} ($_);
		  next LOOP;
	      }
	      if (/^[ \t]+$entry[ \t]+.*\#/) {
		  s/^[ \t]+$entry[ \t]+.*\#/ $entry $val \#/;
	      }
	      print {FILEOUT} ($_);
	  }
      } 
	else {print {FILEOUT} ($_)}
    }
    if ($failed) {
	print("Error in changeparam() when looking for Block $block.\n");
	exit 1;
    }
    close(FILEIN);
    close(FILEOUT);
    unlink("$file.bak");
}

#- }}}
#- {{{ sub movefile:

sub movefile {
    my($fromfile,$tofile,$targetdir);
    ($fromfile,$tofile) = @_;
    if ($tofile =~ /\//) {
	($targetdir = $tofile) =~ s/\/[^\/]*$//;
	unless (-d $targetdir) { system("mkdir -p $targetdir") }
    }
    unless (-f $fromfile) {
	print("Error in movefile: Cannot find file $fromfile.\n");
	exit 1;
    }
    system("/bin/mv -f $fromfile $tofile");
}

#- }}}
