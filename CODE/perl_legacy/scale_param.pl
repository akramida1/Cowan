#!perl
# scale_param.pl

use strict;

my ($F_factor, $G_factor, $CI_factor, $group_factors) = (1,1,1,{});
my $param = 'CI=1';

# Read command-line parameters
my $in_file = shift;
my $out_file = shift;

my $num_groups = 0;
do {
  $param = shift;
  if ( $param =~ /^F *= *(\S+)$/i ) {
    $F_factor = $1;
  } elsif ( $param =~ /^G *= *(\S+)$/i ) {
    $G_factor = $1;
  } elsif ( $param =~ /^CI *= *(\S+)$/i ) {
    $CI_factor = $1;
  } elsif ( $param =~ /^GROUP *([0-9-]+) *= *(\S+)$/i ) {
    my $group_num = abs($1);
    $group_factors->{$group_num} = $2;
    $num_groups++;
  }
} while ( $param =~ /^([FG]|CI|GROUP *([0-9-]+)) *=/i );
my $conf_set = $param;
my $fix = shift;
my $fix_all = ($fix =~ /ALL/i);

# Check the command-line syntax
if ( !$in_file || !$out_file || (!$CI_factor && !$F_factor && !$G_factor && !$num_groups) || !$conf_set ) {
  print "\nUsage:\n";
  print "scale_param <inp_file> <out_file> F=<F_factor>|G=<G_factor>|CI=<CI_factor>|GROUP<number> <conf_set> [FIX]\n" .
    qq{where <inp_file> is a saved copy of RCEINP,
      <out_file> is a new copy of RCEINP to be created with scaled CI parameters,
      <XX_factor> is the factor by which to scale the corresponding group of parameters
        (any or all F=... G=... CI=... GROUP<number>=... may be given in any order),
      <conf_set> can be 1 for the first parity, 2 for the second parity, or 12 for both.
      Optional additional option FIX can be specified to fix all free Slater parameters 
      for the given conf. set(s).

};
  exit;
}

open(INP, $in_file) || die "Could not open $in_file";
open(OUT, ">$out_file") || die "Could not create $out_file";

# Read parameter values from the input file
my $param_section = 0;
while ( <INP> ) {
  chomp;
  if ( /PARAMETER/ ) {
    $param_section++;
  }
  if ( ($conf_set =~ /$param_section/) && /^(.{10} )( *[0-9-]{1,4})( +[0-9.-]{6,13})( +0\.0000)( +[^-]+-[^-]+)* *$/) {
    # The Slater parameters section
    my ($param, $flag, $value, $zeros, $configs) = ($1, $2, $3, $4, $5);
    my $group_num = abs($flag);
    $group_num = 0 unless $group_num < 100;

    if ( $fix_all ) {
      $flag =~ s/^ +//g;
      $flag = ($flag < 0) ? '-100' : ' 100';
    }
    if ( ($param =~ /^F/) && ($F_factor != 1)  ) {
      # Scale the F parameters by the given factor
      $value =~ s/^ +//g;
      $value = sprintf("%14.4f", $value * $F_factor);
      if ( $fix ) {
        $flag =~ s/^ +//g;
        $flag = ($flag < 0) ? '-100' : ' 100';
      }
    } elsif ( ($param =~ /^G/) && ($G_factor != 1)  ) {
      # Scale the G parameters by the given factor
      $value =~ s/^ +//g;
      $value = sprintf("%14.4f", $value * $G_factor);
      if ( $fix ) {
        $flag =~ s/^ +//g;
        $flag = ($flag < 0) ? '-100' : ' 100';
      }
    } elsif ( $configs && ($CI_factor != 1) ) {
      # Scale the CI parameters by the given factor
      $value =~ s/^ +//g;
      $value = sprintf("%14.4f", $value * $CI_factor);
      if ( $fix ) {
        $flag =~ s/^ +//g;
        $flag = ($flag < 0) ? '-100' : ' 100';
      }
    } elsif ( defined($group_factors->{$group_num}) && ($group_factors->{$group_num} != 1) ) {
      $value = sprintf("%14.4f", $value * $group_factors->{$group_num});
      if ( $fix ) {
        $flag =~ s/^ +//g;
        $flag = ($flag < 0) ? '-100' : ' 100';
      }
    }
    print OUT "$param$flag$value$zeros$configs\n";
  } else {
    print OUT "$_\n";
  }
}
close INP;
close OUT || die "Error writing $out_file";

print "Ok.\n";