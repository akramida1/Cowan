#!perl
use strict;
# This program adds a given fixed value to all Eav parameters in the RCG input file ING11}

my ($S1,$S3,$i,$j,$Np,$Res,$Nline,$dE,$HF);

open(ING11,"<ing11.") or die('Error opening ING11 file.');
open(OUT,">tmp.") or die('Error creating output file.');

$S1 = shift;
$dE = $S1+0;
if ( ($S1 !~ /^-{0,1}\d*(\.\d*){0,1}$/) || !$dE) {
  print("\nUsage: \nadd_eav <shift_value>\nwhere <shift_value> is the quantity to add to all Eav parameters (in 1000 cm-1).\n\n");
  exit;
}

$Nline = 0;
while ( ($S1 = <ING11>) && defined($S1) ) {
  chomp $S1;
  $Nline++;
  if ((length($S1)<=71) || (substr($S1,70,1) !~ /H/i) || (substr($S1,9,1) eq '-') || (substr($S1,54,1) eq '/')) {
    print OUT "$S1\n" or die 'Error writing output file.';
    next;
  }
  print OUT substr($S1,0,20);
  $S3 = substr($S1,18,2);
  $S3 =~ s/^\s+//g;
  $Np = $S3+0;
  if ( ($S3 !~ /^-{0,1}\d+$/) || !$Np) {
    die("ING11 file format error in line $Nline: wrong number of parameters in columns 19-20.");
  }
  $S3 = substr($S1,20,9);
  $S3 =~ s/^\s+//g;
  $HF = $S3+0;
  if ( ($S3 !~ /^-{0,1}\d*(\.\d*){0,1}(E[+-]\d+){0,1}$/) || (!$HF && ($S3 ne '0')) ) {
    die("ING11 file format error in line $Nline: wrong EAV value.");
  }
  my $thousands = ($S3 =~ /\./);
  if ( $thousands ) {
    $HF += $dE;
  } else {
    $HF += $dE*10000;
  }

  my $s = sprintf("%d",$HF);
  if ( !$thousands && (($s > 999999999) || ($s < -99999999)) ) {
    $thousands = 1;
    $HF = $HF*0.0001;
  }
  if ( $thousands ) {
    $S3 = sprintf("%9.1f",$HF);
  } else {
    $S3 = sprintf("%9d",$HF);
  }

  print OUT $S3,substr($S1,29,length($S1)-29),"\n" or die 'Error writing output file.';
}
close(ING11);
close(OUT) || die "Error writing output file";

# Delete ING11.BAK if it exists ------------------------------------------
if ( -f 'ING11.BAK' ) {
  `del ING11.BAK`;
}
# Rename ING11 to ING11.BAK ----------------------------------------------
`copy ING11 ING11.BAK`;

# Rename temp file to $out
`copy tmp. ING11.`;
`del tmp.`;

print "Adding finished.\n";
