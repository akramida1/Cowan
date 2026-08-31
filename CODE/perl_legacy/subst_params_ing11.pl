#!perl
# This program substitutes parameter values from one ING11 file into another one
# Substitution is done line by line. If the conf.name(s) on the line read
# from the input file matches those in the output file, the parameter line(s) for these
# configs are substituted in place of similar lines in the output file.
# The input ING11 file is assumed to have both parities, while the output ING11 is supposed
# to have only one parity set.
# The main purpose is to substitute fitted parameters to set up ING11 for
# autoionization calculations.

use strict;

my $inp = shift;
my $out = shift;
my $conf_set = shift;

if ( !$inp || !$out || ($conf_set !~ /^[12]$/) ) {
  print "\nUsage: \nsubst_params_ing11 <input_ING11> <output_ING11> <conf_set>\n";
  print " where <conf_set> can be 1 or 2, for conf. set #1 or #2 in  <input_ING11>.\n\n";
  exit;
}

my $out_tmp = "$out.tmp";
open (INP, "<$inp") || die "Could not open $inp for reading";
open (INP1, "<$out") || die "Could not open $out for reading";
open (OUT, ">$out_tmp") || die "Could not create $out_tmp";

my ($s,$s1);

# Read the first lines; copy the first line from <output_ING11> to temp file
$s = <INP>;
$s1 = <INP1>;
print OUT $s1;

# Skip to the first conf. parameter line in <input_ING11>
while ( defined($s = <INP>) && ($s =~ /^([SPDFGHIKLMNOQTUV]( \d|\d\d)  ){8}/i) ) {
  next;
}
$s = $s;
# Copy the shell specification lines from <output_ING11> to temp. file
while ( defined($s1 = <INP1>) && ($s1 =~ /^([SPDFGHIKLMNOQTUV]( \d|\d\d)  ){8}/i) ) {
  print OUT $s1;
  next;
}

my $par = 1;
if ( $conf_set == 2 ) {
  # Skip to the second conf. set in <input_ING11>
  # The end of the first set is determined by the end of the first CI section (minus in 10th position in <input_ING11>)
  # Skip to the beginning of the first CI section
  while ( defined($s = <INP>) && ($s !~ /^.{9}-.{8}( \d|\d\d)/) ) {
    next;
  }
  $s = $s;
  # Skip to the end of the first CI section
  while ( defined($s = <INP>) && ($s =~ /^(.{9}-.{8})( \d|\d\d)/) ) {
    my ($conf1, $n_par1) = ($1,$2);
    $n_par1 =~ s/^\s+|\s+$//g;
    if ( $n_par1 > 5 ) {
      # Skip additional parameter lines for this config
      my $n_add_p = $n_par1 - 5;
      my $n_lines = 1;
      {
        use integer;
        $n_lines = ($n_add_p - 1) / 7 + 1;
      }
      for ( my $i = 1; $i <= $n_lines; $i++ ) {
        $s = <INP>;
      }
    }

    next;
  }
}

# Read single-conf. parameters.
# Compare the line from <input_ING11> with the line in <output_ING11>.
# If the same config, write this conf section (maybe several lines) from <input_ING11> to temp file,
# otherwise, write this conf. section from <output_ING11> to temp file.

my ($conf1, $conf2, $n_par1, $n_par2);

do {
  if ( $s =~ /^(.{18})([ 0-9]{2})/ ) {
    ($conf1, $n_par1) = ($1,$2);
  }
  $n_par1 =~ s/^\s+|\s+$//g;
  while ( defined($s1) && ($s1 =~ /^(.{18})([ 0-9]{2})/) && (($conf2, $n_par2) = ($1,$2)) && ($conf2 ne $conf1) ) {
    print OUT $s1;
    $n_par2 =~ s/^\s+|\s+$//g;
    if ( $n_par2 > 5 ) {
      # Read and write additional parameter lines for this config from <output_ING11>
      my $n_add_p = $n_par2 - 5;
      my $n_lines = 1;
      {
        use integer;
        $n_lines = ($n_add_p - 1) / 7 + 1;
      }
      for ( my $i = 1; $i <= $n_lines; $i++ ) {
        $s1 = <INP1>;
        print OUT $s1;
      }
    }
    $s1 = <INP1>; # Read next line from <output_ING11>
  }
  if ( (!defined($s1)) || ($s1 =~ /^.{9}-.{8}( \d|\d\d)/) ) {
    die "Conf. $conf1 not found in $out";
  }
  $n_par2 =~ s/^\s+|\s+$//g;
  if ($n_par1 ne $n_par2) {
    die "Parameters number mismatch for conf. $conf1";
  }
  if ( $conf1 eq $conf2 ) {
    print OUT $s;
    if ( $n_par1 > 5 ) {
      # Read and write additional parameter lines for this config
      my $n_add_p = $n_par1 - 5;
      my $n_lines = 1;
      {
        use integer;
        $n_lines = ($n_add_p - 1) / 7 + 1;
      }
      for ( my $i = 1; $i <= $n_lines; $i++ ) {
        $s = <INP>;
        $s1 = <INP1>;
        print OUT $s;
      }
    }
  } else {
    die "Conf. $conf1 not found in $out";
  }
  # Read next lines
} while ( defined($s = <INP>) && defined($s1 = <INP1>) && ( $s !~ /^.{9}-.{8}( \d|\d\d)/ ) );

# Read/write CI sections
do {
  if ( $s =~ /^(.{9}-.{8})( \d|\d\d)/ ) {
    ($conf1, $n_par1) = ($1,$2);
  }
  $n_par1 =~ s/^\s+|\s+$//g;
  while ( defined($s1) && ($s1 =~ /^(.{18})( \d|\d\d)/) && (($conf2, $n_par2) = ($1,$2)) && ($conf2 ne $conf1) ) {
    print OUT $s1;
    $n_par2 =~ s/^\s+|\s+$//g;
    if ( $n_par2 > 5 ) {
      # Read and write additional parameter lines for this config from <output_ING11>
      my $n_add_p = $n_par2 - 5;
      my $n_lines = 1;
      {
        use integer;
        $n_lines = ($n_add_p - 1) / 7 + 1;
      }
      for ( my $i = 1; $i <= $n_lines; $i++ ) {
        $s1 = <INP1>;
        print OUT $s1;
      }
    }
    $s1 = <INP1>; # Read next line from <output_ING11>
  }
  if ( (!defined($s1)) || ($s1 !~ /^.{9}-.{8}( \d|\d\d)/) ) {
    die "Configs $conf1 not found in CI section of $out";
  }
  $n_par2 =~ s/^\s+|\s+$//g;
  if ( $conf1 eq $conf2 ) {
    if ($n_par1 ne $n_par2) {
      print "Warning: param. number mismatch for configs $conf1: $n_par1 in first file and $n_par2 in second. I skip this section.\n";
      print OUT $s1;
      if ( $n_par2 > 5 ) {
        # Read and write additional parameter lines for this config from <output_ING11>
        my $n_add_p = $n_par2 - 5;
        my $n_lines = 1;
        {
          use integer;
          $n_lines = ($n_add_p - 1) / 7 + 1;
        }
        for ( my $i = 1; $i <= $n_lines; $i++ ) {
          $s1 = <INP1>;
          print OUT $s1;
        }
      }
      if ( $n_par1 > 5 ) {
        # Skip additional parameter lines for this config in <inputING11>
        my $n_add_p = $n_par1 - 5;
        my $n_lines = 1;
        {
          use integer;
          $n_lines = ($n_add_p - 1) / 7 + 1;
        }
        for ( my $i = 1; $i <= $n_lines; $i++ ) {
          $s = <INP>;
        }
      }

    } else {
      print OUT $s;
      if ( $n_par1 > 5 ) {
        # Read and write additional parameter lines for this config
        my $n_add_p = $n_par1 - 5;
        my $n_lines = 1;
        {
          use integer;
          $n_lines = ($n_add_p - 1) / 7 + 1;
        }
        for ( my $i = 1; $i <= $n_lines; $i++ ) {
          $s = <INP>;
          $s1 = <INP1>;
          print OUT $s;
        }
      }
    }
  } else {
    die "Conf. $conf1 not found in $out";
  }
  # Read next lines
} while ( defined($s = <INP>) && defined($s1 = <INP1>) && ($s =~ /^.{9}-.{8}( \d|\d\d)/) && ($s1 =~ /^.{9}-.{8}( \d|\d\d)/) );

$s = $s;
# Read/write the rest of the <output_ING11> file
do {
  print OUT $s1;
} while ( defined($s1 = <INP1>) );

close INP;
close INP1;
close OUT || die "Error writing to $out_tmp";

# Rename temp file to $out
`del $out`;
rename $out_tmp, $out;

print "Ok.\n";
