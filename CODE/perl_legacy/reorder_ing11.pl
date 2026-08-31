#!perl
use strict;

my $out_file = shift;
my $new_ord_str = shift;
my @new_ord_arr = ();
if ( length($new_ord_str) > 1 ) {
  @new_ord_arr = split(//,$new_ord_str);
} else {
  $new_ord_arr[0] = $new_ord_str if $new_ord_str;
  for ( my $i = 1; $i <=7; $i++ ) {
    my $n = shift;
    last unless $n;
    push(@new_ord_arr, $n);
  }
}
my %new_ord_hash = {};
for ( my $i = 1; $i <=8; $i++ ) {
  my $n = $new_ord_arr[$i-1];
  last unless $n;
  $new_ord_hash{$n} = $i;
}
#new_order = 54231
#new_ord_arr = 53421
#new_ord_hash = (5=>1, 3=>2, 4=>3, 2=>4, 1=>5);
#new_order = 54231
my @new_order = ();
#for ( my $i = 1; $i <= 8; $i++ ) {
#  my $n = shift;
#  last unless $n;
#  push(@new_order, $n);
#}
foreach my $key (sort {$a<=>$b} keys %new_ord_hash) {
  push(@new_order, $new_ord_hash{$key}) unless !defined($new_ord_hash{$key});
}
if ( !$out_file ) {
  print "\nUsage: \nreorder_ing11 <out_file_name> <1st_shell_no><2nd_shell_no>[3rd_shell_no[...[8th_shell_no]]]]]]\n\n";
  exit;
}

open ING11, "<ING11" or die "Could not open input file ING11";
open OUT_FILE, ">$out_file" or die "Could not open output file " . $out_file;

my $s = '';
my $parity = 0;
my $spectrum_name = '';

# Initialize arrays
my @confs = ({},{});
my @EEt = ({},{});
my @EEk = ({},{});
my @params = ({},{});
my @CI = ([],[]);
my %confs36 = ();
my @shells = ({},{});
my @parities = ('','');
my $n_conf = '';
my %conf_nums = ();
my $last_shell_num = 1;
my $last_movable_shell = 1;

&read_in36();

# Read shell definition section of ING11
$s = &read_confs($s);

if ( $#new_order + 1 != $last_movable_shell)
{
  print "You must provide the new subshell numbers for the first $last_movable_shell occupied subshells.\n";
  exit;
}

my $i = 0;
foreach my $j (sort @new_order) {
  $i++;
  if ( $i != $j ) {
    print "Incorrect new subshell numbers. You must give all integers between 1 and $last_movable_shell.\n";
    exit;
  }
}

$s = &read_params($s);

&read_outg11();
&Reorder();

&write_shells();
&write_params();

do {
  chomp $s;
  print OUT_FILE "$s\n";
} while ( defined($s = <ING11>) ) ;

close ING11;
close OUT_FILE;

1;

############################################################################
sub read_in36		#3/16/2004 2:51PM A.Kramida
############################################################################
{
	open IN36, "<IN36" or die "Could not open input file IN36";
	my $s = '';
	my $parity = -1;

	$s = <IN36>; # Skip the first line

	# Read configurations
  while ( (defined($s = <IN36>)) && ($s =~ /^.{16}(.{12}) *\d* *(\d+)([^ 0-9])(\d*) *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1}/) ) {
		my $conf = $1;
		my ($n1,$s1,$w1,$n2,$s2,$w2,$n3,$s3,$w3,$n4,$s4,$w4,$n5,$s5,$w5,$n6,$s6,$w6,$n7,$s7,$w7,$n8,$s8,$w8,$n9,$s9,$w9) =
		  ($2,$3,$4, $6,$7,$8, $10,$11,$12, $14,$15,$16, $18,$19,$20, $22,$23,$24, $26,$27,$28, $30,$31,$32, $34,$35,$36);
  	$conf =~ s/^\s+|\s+$//g;

		$confs36{$conf} = [];
		for ( my $n_shell = 1; $n_shell <=9; $n_shell++ ) {
			my ($n,$sh,$w) = (0,'',0);
			eval("\$w = \$w$n_shell");
			eval("\$n = \$n$n_shell");
			eval("\$sh = \$s$n_shell");
			if ( ($w eq '') && ($n > 0) && $sh ) {
				$w = 1;
			}
			if ( $w > 0 ) {
				$s = lc($s);
				push(@{$confs36{$conf}},[$n,$sh,$w]);
			}
		}
	}

	my $delete_first = 1;
	foreach my $conf (keys %confs36) {
		if ( ($confs36{$conf}->[0]->[1] ne 'f') || ($confs36{$conf}->[0]->[2] != 14)  ) {
			$delete_first = 0;
		}
	}
	if ( $delete_first ) {
		foreach my $conf (keys %confs36) {
		  shift(@{$confs36{$conf}});
		}
	}
} # read_in36

############################################################################
sub read_confs    #9/30/2004 2:34PM A.Kramida
############################################################################
{
  my $s = '';
  my $par = 1;
  my $first_par = 0;
  my $nc = 0;  # Config. number
  while ( defined($s = <ING11>)) {
    chomp $s;
    if ($s =~ /^ {4}\d /) {
      # Read the first (options) line(s) of ING11 and write them to the output file
      print OUT_FILE "$s\n";
      next;
    }

    if ($s =~ /^(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)  (.{6})(.{12}) *([0-9.-]+) +([0-9.-]+)$/ ) {
      # Read configurations
      $nc++;

      # Store the shell filling numbers
      my ($s1,$w1,$s2,$w2,$s3,$w3,$s4,$w4,$s5,$w5,$s6,$w6,$s7,$w7,$s8,$w8,$sp_name,$conf_name,$Etot,$Ek) =
        ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20);

      $spectrum_name = $sp_name unless $spectrum_name;
      $conf_name =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces

      # Store the shell filling numbers
      my @conf_shells = ([lc($s1),$w1],[lc($s2),$w2],[lc($s3),$w3],[lc($s4),$w4],[lc($s5),$w5],[lc($s6),$w6],[lc($s7),$w7],[lc($s8),$w8]);

      # Calculate the parity of this config
      my $par_value = 1;
      for ( my $n_shell = 1; $n_shell <= 8; $n_shell++ ) {
        my ($L, $occup) = @{$conf_shells[$n_shell-1]};
        for ( my $i = 1; $i<=$occup; $i++ ) {
          $par_value *= &shell_parity($L);
        }
        # Determine the number of the last occupied subshell
        if ( $occup && ($last_shell_num < $n_shell) ) {
          $last_shell_num = $n_shell;
        }
      }

      if ( !$first_par ) {
        # Set parity value of the first configuration set
        $first_par = $par_value;
      }
      if ( $par_value != $first_par ) {
        # Switch to the next config set
        $first_par = $par_value;
        $par++;
        $nc = 1;
        if ( $par > 2 ) {
          die "Something wrong with ING11: more than 2 parities detected.";
        }
      }
      $shells[$par-1]->{$nc} = \@conf_shells;

      $EEt[$par-1]->{$nc} = $Etot;
      $EEk[$par-1]->{$nc} = $Ek;

      # Set the name of the config
      $confs[$par-1]->{$nc} = $conf_name;
      # Store the config parity and number for backwards reference
      $conf_nums{$conf_name} = [$par,$nc];
      if ( $par_value == -1 ) {
        $parities[$par-1] = 'o';
      } else {
        $parities[$par-1] = 'e';
      }

      # Find this conf in confs36 and replace shell names with those from IN36 (including the principal quantum #)
      my @shells36 = @{$confs36{$conf_name}};
      for ( my $i = 0; $i <= $#shells36; $i++ )
      {
        my ($n, $sh, $w) = @{$shells36[$i]};
        for ( my $j = 0; $j <= 7; $j++ ) {
          if ( ($sh eq $shells[$par-1]->{$nc}->[$j]->[0]) && ($w == $shells[$par-1]->{$nc}->[$j]->[1]) ) {
            $shells[$par-1]->{$nc}->[$j]->[0] = "$n$sh";
            last;
          }
        }
      }
    } else {
      last;
    }
    $par = $par;
  }

  # If the last occupied shell specified in ING11 has different principal quantum numbers
  # in different configs, set the last shell unmovable
  $last_movable_shell = $last_shell_num;
  my $prev_nL = '';
  for ( my $ip = 1; $ip <= 2; $ip++ ) {
    last unless $last_movable_shell == $last_shell_num;
    foreach my $nc (keys %{$shells[$ip-1]} ) {
      my ($nL, $occup) = @{$shells[$ip-1]->{$nc}->[$last_shell_num-1]};
      next unless $occup;
      if ( $prev_nL && ($prev_nL ne $nL) ) {
        $last_movable_shell--;
        last;
      }
      $prev_nL = $nL;
    }
  }
  return $s;
}	##read_confs

############################################################################
sub read_params   #9/30/2004 9:01PM A.Kramida
############################################################################
{
  my $s = shift;
  do {
    chomp $s;
	if ($s =~ /^2s2      -2p2/) {
	  $s = $s;
	}
    #if ( $s =~ /^(.{6})([^-]{12})(\d\d| \d) *([0-9.-]{1,9})(\d) *([0-9.-]{1,9})(\d) *([0-9.-]{1,9})(\d) *([0-9.-]{1,9})(\d) *([0-9.-]{1,9})(\d)([A-Z]{2}\d{8})\s*$/  ){
    if ( $s =~ /^(.{6})([^-]{12})(.{2})(.{9})(\d)(.{9})(\d)(.{9})(\d)(.{9})(\d)(.{9})(\d)([A-Z]{2}\d{8})\s*$/i  ){
      my ($spectr, $conf, $num_param, $p1, $t1, $p2, $t2, $p3, $t3, $p4, $t4, $p5, $t5, $scaling) =
        ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14);
      foreach ($conf,$num_param,$p1,$t1,$p2,$t2,$p3,$t3,$p4,$t4,$p5,$t5) {
        s/^\s+|\s+$//g;  # Strip leading and trailing spaces
      }

      my ($par,$nc) = ();
      if ( defined($conf_nums{$conf}->[0]) ) {
        ($par,$nc) = @{$conf_nums{$conf}}; # Get the number of this config
      }
      die "Config. name $conf not recognized in ING11" unless $nc;

      # Initialize the parameters hash for this config
      $params[$par-1]->{$nc} = {};

      $params[$par-1]->{$nc}->{'num_param'} = $num_param;
      $params[$par-1]->{$nc}->{'scaling'} = $scaling;

      my $n = ($num_param <=5) ? $num_param : 5;
      for ( my $j = 1; $j <= $n; $j++ ) {
        my ($p, $t) = ();
        my $expression = "(\$p, \$t) = (\$p$j,\$t$j)";
        eval($expression);
        $params[$par-1]->{$nc}->{"p$j"} = [$p,$t];
        if ( !defined($params[$par-1]->{$nc}->{$t}) ) {
          $params[$par-1]->{$nc}->{$t} = [];
        }
        push(@{$params[$par-1]->{$nc}->{$t}}, $p);
      }

      if ( $num_param > 5 ) {
        my $rest = $num_param - 5;
        # Determine how many additional lines of parameters there are for this config
        my $n_lines = 0;
        { use integer;
          $n_lines = ($rest - 1) / 7 + 1;
        }
        # Read additional parameter lines
        for ( my $i = 1; $i <= $n_lines; $i++ ) {
          $s = <ING11>;
          chomp $s;
          if ( $s =~ /^ +([0-9.e+-]{1,9})(\d)( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9e+.-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1} *$/i ) {
            my $num_to_read = ($rest >=7) ? 7 : $rest;

            my ($p1, $t1, $p2, $t2, $p3, $t3, $p4, $t4, $p5, $t5, $p6, $t6, $p7, $t7) =
               ($1, $2, $4, $5, $7, $8, $10, $11, $13, $14, $16, $17, $19, $20);

            my $par_num = $num_param - $rest;

            for ( my $j = 1; $j <= $num_to_read; $j++ ) {
              $par_num++;
              my ($p,$t) = ();
              my $expression = "(\$p,\$t) = (\$p$j,\$t$j)";
              eval($expression);

              $params[$par-1]->{$nc}->{"p$par_num"} = [$p,$t];
              if ( !defined($params[$par-1]->{$nc}->{$t}) ) {
                $params[$par-1]->{$nc}->{$t} = [];
              }
              push(@{$params[$par-1]->{$nc}->{$t}}, $p);
            }
            $rest -= $num_to_read;
          } else {
            die "Format error in ING11 in parameters for config $conf";
          }
        }
      }
    } elsif ( $s =~ /^(.{9})-(.{8})([ \d]{2})([ 0-9.e+-]{9})5([ 0-9.e+-]{9})5([ 0-9.e+-]{9})5([ 0-9.e+-]{9})5([ 0-9.e+-]{9})5([A-Z]{2}\d{8}) *$/i ) {
      # CI parameters section
      my ($conf1,$conf2,$num_param, $p1,$p2,$p3,$p4,$p5,$scaling) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
	  foreach ($conf1,$conf2,$num_param, $p1,$p2,$p3,$p4,$p5) {
        $_ =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces
      }

      my ($par1,$nc1,$par2,$nc2) = ();
      if ( defined($conf_nums{$conf1}->[0]) ) {
        ($par1,$nc1) = @{$conf_nums{$conf1}}; # Get the number of this config
      }
      die "Config. name $conf1 not recognized in CI section of ING11" unless $nc1;
      if ( defined($conf_nums{$conf2}->[0]) ) {
        ($par2,$nc2) = @{$conf_nums{$conf2}}; # Get the number of this config
      }
      die "Config. name $conf2 not recognized in CI section of ING11" unless $nc2;
      die "Config. parity mismatch for configs $conf1 and $conf2 in CI section of ING11" unless $par1 == $par2;

      if ( ($conf1 eq 'd10p2') && ($conf2 eq 'd9p2d' ) ) {
        $par1 = $par1;
      }
      $CI[$par1-1]->[$nc1-1] = [] unless defined($CI[$par1-1]->[$nc1-1]);
      $CI[$par1-1]->[$nc1-1]->[$nc2-1] = {} unless defined($CI[$par1-1]->[$nc1-1]->[$nc2-1]);
      $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'} = [] unless defined($CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'});
      $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'num_param'} = $num_param;

      my $n = ($num_param <= 5) ? $num_param : 5;
      for ( my $j = 1; $j <= $n; $j++ ) {
        my $p;
        my $expression = "\$p = \$p$j";
        eval($expression);

        $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$j-1] = $p;
      }
      if ( $num_param > 5 ) {
        # Since num. of params is a two-digit number,
        # There can be many additional lines of CI parameters
        my $rest = $num_param - 5;
        # Determine how many additional lines of parameters there are for this config
        my $n_lines = 0;
        { use integer;
          $n_lines = ($rest - 1) / 7 + 1;
        }
        # Read additional parameter lines
        for ( my $i = 1; $i <= $n_lines; $i++ ) {
          $s = <ING11>;
          chomp $s;
          if ( $s =~ /^ *([0-9.-]{5,9})5( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1} *$/i ) {
            my ($p1, $p2, $p3, $p4, $p5, $p6, $p7) =
              ($1, $3, $5, $7, $9, $11, $13);
            my $num_to_read = ($rest >=7) ? 7 : $rest;

            my $par_num = $num_param - $rest;

            for ( my $j = 1; $j <= $num_to_read; $j++ ) {
              $par_num++;
              my $p = 0;
              my $expression = "\$p = \$p$j";
              eval($expression);

              $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$j+4] = $p;
            }
            $rest -= $num_to_read;
          } else {
            die "Format error in ING11 in CI parameters for configs $conf1 and $conf2";
          }
        }
      }
    } else {
      # End of parameters section; stop here
      return $s;
    }
  } while ( defined($s = <ING11> ) );
  return $s;

} ##read_params

############################################################################
sub write_shells    #9/30/2004 12:57PM A.Kramida
############################################################################
{
  for ( my $par = 1; $par <= 2; $par++ ) {
    foreach my $nc ( sort {$a<=>$b} keys %{$shells[$par-1]} ) {
      my $s1 = '';
      for ( my $n_shell = 0; $n_shell <= 7; $n_shell++ ) {
        my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell]};
        $nL =~ s/^[0-9]+//g; # Strip principal quantum number
        $s1 .= uc($nL) . sprintf("%2d  ",$occup) if ($nL);
      }
      my $conf_name = $confs[$par-1]->{$nc};
      next unless ($nc && $s1);
      my $Ek = $EEk[$par-1]->{$nc};
      my $Et = $EEt[$par-1]->{$nc};
      $s1 .= "$spectrum_name" . sprintf("%-12s%13.3f%9.4f",$conf_name,$Et,$Ek);

      print OUT_FILE "$s1\n";
    }
  }
} ##write_shells

############################################################################
sub write_params    #9/30/2004 1:18PM A.Kramida
############################################################################
{
  for ( my $par = 1; $par <= 2; $par++ ) {
    my $scaling = '';
    # Write the Slater parameters for each config
    foreach my $nc ( sort {$a<=>$b} keys %{$params[$par-1]} ) {
      my $conf_name = $confs[$par-1]->{$nc};
      my $num_param = $params[$par-1]->{$nc}->{'num_param'};
      $scaling   = $params[$par-1]->{$nc}->{'scaling'};

      my $s1 = $spectrum_name . sprintf("%-12s%2d",$conf_name, $num_param);

      my $np = ($num_param <= 5) ? 5 : $num_param;
      for ( my $j = 1; $j <= $np; $j++ ) {
        my ($p, $t) = ('0','0');
        if ( $j <= $num_param ) {
          #if ( !defined($params[$par-1]->{$nc}->{"p$j"}) ) {
          #  $par = $par;
          #}
          ($p, $t) = @{$params[$par-1]->{$nc}->{"p$j"}};
        }

        if ( $p =~ /[e]/i ) {
          $p = sprintf("%10.2e",$p);
          $p =~ s/(e[+-])0/$1/;  # Remove the extra zero in the exponent: Perl has one more compared to Fortran
          if ( $p eq ' 0.00e+00' ) {
            $p = '        0';
          }
          $s1 .= "$p$t";
        } elsif ( $p =~ /[.]\d\d/ ) {
          $s1 .= sprintf("%9.2f%d",$p, $t);
        } elsif ( $p =~ /[.]\d/ ) {
          $s1 .= sprintf("%9.1f%d",$p, $t);
        } else {
          $s1 .= sprintf("%9d%d",$p, $t);
        }
        #if ( ($p >= -1e-3) && ($p <= 1e-3) && ($j > 1) ) {
        #  $p = sprintf("%10.2e",$p);
        #  $p =~ s/(e[+-])0/$1/;  # Remove the extra zero in the exponent: Perl has one more compared to Fortran
        #} elsif ( ($p >= -9999.99) && ($p <= 99999.99) ) {
        #  $p = sprintf("%9.0f", $p);
        #} else {
        #  if ( $j == 1 ) {
        #    $p = sprintf("%9.1f", $p);
        #  } else {
        #    $p = sprintf("%9.2f", $p);
        #  }
        #}
        #$s1 .= "$p$t";
        if ( $j == 5 ) {
          # Write the main parameter line for this config and reset $s1
          $s1 .= $scaling;
          print OUT_FILE "$s1\n";
          $s1 = '';
          next;
        }
        if ( (($j - 5) % 7 == 0) || (($j > 5) && ($j == $num_param)) ) {
          # Write the additional parameter line for this config and reset $s1
          print OUT_FILE "$s1\n";
          $s1 = '';
        }
      }
    }
    # Write the CI section for each parity
    my $num_c1 = $#{$CI[$par-1]} +  1;
    for (my $nc1=1; $nc1 <= $num_c1; $nc1++ ) {
      my $conf1 = $confs[$par-1]->{$nc1};
      my $num_c2 = $#{$CI[$par-1]->[$nc1-1]} +  1;
      for (my $nc2=1; $nc2 <= $num_c2; $nc2++ ) {
        my $conf2 = $confs[$par-1]->{$nc2};
        next unless defined($CI[$par-1]->[$nc1-1]->[$nc2-1]->{'num_param'});
        my $num_param = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'num_param'};
        my $np = ($num_param > 5) ? $num_param : 5;
        my $s1 = '';
        for ( my $i = 1; $i<=$np; $i++ ) {
          if ( $i == 1 ) {
            $s1 = sprintf("%-9s-%-8s%2d",$conf1,$conf2,$num_param);
          }
          my $p = ($i <= $num_param) ? $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$i-1] : 0;
          $s1 .= sprintf("%9.4f5",$p);
          if ( $i == 5 ) {
            # Write the main parameter line for this pair of configs and reset $s1
            $s1 .= $scaling;
            print OUT_FILE "$s1\n";
            $s1 = '';
            next;
          }
          if ( ($i == $np) || (($i > 5) && (($i-5) % 7 == 0)) ) {
            # Write the additional parameter line for this config and reset $s1
            print OUT_FILE "$s1\n";
            $s1 = '';
          }
        }
      }
    }
  }
} ##write_params

############################################################################
sub read_outg11   #9/30/2004 8:47PM A.Kramida
############################################################################
{
  open OUTG11, "<OUTG11" or die "Could not open input file OUTG11";
	my $s = '';
  my $spectrum = $spectrum_name;
  $spectrum =~ s/([+$().:^#])/\\$1/g; # Preceed special chars with escape symbol

  # Read parameters
  for ( my $par = 1; $par<=2; $par++ ) {
    while ( (defined($s = <OUTG11>)) && ($s !~ /^ *$spectrum([^-]{12}) *PARAMETER VALUES IN +([0-9.]+) /) ) {
      next;
    }
    last unless defined $s;
    do {
      $s =~ /^ *$spectrum([^-]{12}) *PARAMETER VALUES IN +([0-9.]+) /;
      my ($conf, $units) = ($1, $2);
      $conf =~ s/^\s+|\s+$//g;

      my ($par1,$nc) = ();
      if ( defined($conf_nums{$conf}->[0]) ) {
        ($par1,$nc) = @{$conf_nums{$conf}}; # Get the number of this config and its parity
      }
      die "Config. name $conf not recognized in OUTG11" unless $nc;
      die "Parity mismatch between OUTG11 and ING11 for config $conf" unless ($par1 == $par);

      my $num_param = $params[$par-1]->{$nc}->{'num_param'};

      $s = <OUTG11>;
      $s = <OUTG11>;

      if ( $s =~ /^.{39}(EAV) *([A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1} *$/  ) {
        my @pars = ($1,$2,$3,$4,$5);
        for ( my $i = 1; $i <= 5; $i ++ ) {
          my $p = $pars[$i-1];
          $p =~ s/^\s+|\s+$//g;
          # Store the parameter name
          $params[$par-1]->{$nc}->{"p$i"}->[2] = $p if $p;
          last if $i >= $num_param;
        }
      } else {
        die "Format error in the main parameter line of OUTG11 for parity $par, config $conf";
      }
      if ( $num_param > 5 ) {
        my $rest = $num_param - 5;
        my $n_lines = 0;
        { use integer;
          $n_lines = ($rest - 1) / 7 + 1;
        }
        for ( my $n_line = 1; $n_line <= $n_lines; $n_line++ ) {
          $s = <OUTG11>;
          if ( $s =~ /^( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\))( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1}( +[A-Z]{4,5} {0,1}\d{0,1}| +[FGT][0-9*]{0,1}\([^()]{2,3}\)){0,1} *$/  ) {
            my @pars = ($1,$2,$3,$4,$5,$6,$7);
            for ( my $i = 1; $i <= 7; $i ++ ) {
              my $p = $pars[$i-1];
              $p =~ s/^\s+|\s+$//g;
              # Store the parameter name
              my $j = 5 + ($n_line - 1) * 7 + $i;
              $params[$par-1]->{$nc}->{"p$j"}->[2] = $p if $p;
              last if $j >= $num_param;
            }
          } else {
            die "Format error in the additional parameter line $n_line of OUTG11 for parity $par, config $conf";
          }
        }

      }

      # Scroll to CI section or next parity
      while ( (defined($s = <OUTG11>)) && ($s !~ /^ *$spectrum(.{12}) *PARAMETER VALUES IN +([0-9.]+) /) ) {
        if ( $s =~ /ENERGY MATRIX|COUPLING|EIGEN|PURITY/ ) {
          # Switch to next parity
          last;
        }
        if ( $s =~ /^ (.{9})-(.{8}). *PARAMETER VALUES IN +([0-9.]+) / ) {
          # Read the CI parameter section
          my ($conf1, $conf2, $units) = ($1, $2, $3);
          $conf1 =~ s/^\s+|\s+$//g;
          $conf2 =~ s/^\s+|\s+$//g;
          if ( ($conf1 eq '3d65p') && ($conf2 eq 'd5s4p') ) {
            $s=$s;
          }
          my ($par1,$par2, $nc1, $nc2) = ();
          if ( defined($conf_nums{$conf1}->[0]) ) {
            ($par1,$nc1) = @{$conf_nums{$conf1}}; # Get the number of this config and its parity
          }
          if ( defined($conf_nums{$conf2}->[0]) ) {
            ($par2,$nc2) = @{$conf_nums{$conf2}}; # Get the number of this config and its parity
          }
          die "Config. name $conf1 not recognized in CI section of OUTG11" unless $nc1;
          die "Config. name $conf2 not recognized in CI section of OUTG11" unless $nc2;
          die "Parity mismatch between OUTG11 and ING11 for config $conf1 in CI section" unless ($par1 == $par);
          die "Parity mismatch between OUTG11 and ING11 for config $conf2 in CI section" unless ($par2 == $par);

          my $num_param = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'num_param'};

          $s = <OUTG11>;
          $s = <OUTG11>;
          $s = <OUTG11>;

          if ( $s =~ /^.{39}([^ ]{8})( +[^ ]{8}){0,1}( +[^ ]{8}){0,1}( +[^ ]{8}){0,1}( +[^ ]{8}){0,1} *$/  ) {
            my @pars = ($1,$2,$3,$4,$5);
            $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'param_names'} = [];
            for ( my $i = 1; $i <= 5; $i ++ ) {
              my $p = $pars[$i-1];
              $p =~ s/^\s+|\s+$//g;
              # Store the parameter name
              $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'param_names'}->[$i-1] = $p if $p;
              last if $i >= $num_param;
            }
          } else {
            die "Format error in the main parameter line of OUTG11 for parity $par, config $conf";
          }
          if ( $num_param > 5 ) {
            my $rest = $num_param - 5;
            my $n_lines = 0;
            { use integer;
              $n_lines = ($rest - 1) / 7 + 1;
            }
            for ( my $n_line = 1; $n_line <= $n_lines; $n_line++ ) {
              $s = <OUTG11>;
              if ( $s =~ /^( +[^ ]{8})( +[^ ]{8}){0,1}( +[^ ]{8}){0,1}( +[^ ]{8}){0,1}( +[^ ]{8}){0,1}( +[^ ]{8}){0,1}( +[^ ]{8}){0,1} *$/  ) {
                my @pars = ($1,$2,$3,$4,$5,$6,$7);
                for ( my $i = 1; $i <= 7; $i ++ ) {
                  my $p = $pars[$i-1];
                  $p =~ s/^\s+|\s+$//g;
                  # Store the parameter name
                  my $j = 5 + ($n_line - 1) * 7 + $i;
                  $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'param_names'}->[$j-1] = $p if $p;
                  last if $j >= $num_param;
                }
              } else {
                die "Format error in the additional parameter line $n_line of OUTG11 for parity $par, config $conf";
              }
            }

          }

        }
        next;
      }
      if ( !defined($s) || ($s =~ /ENERGY MATRIX|COUPLING|EIGEN|PURITY/) ) {
        # Switch to next parity
        next;
      }
      #elsif ($s =~ /^ *$spectrum(.{12}) *PARAMETER VALUES IN +([0-9.]+) /) {
      #  next;
      #}
    } while ( 1 ) ; #(defined($s = <OUTG11>)) && ($s =~ /^ *$spectrum(.{12}) *PARAMETER VALUES IN +([0-9.]+) /) );
    if ( $s =~ /ENERGY MATRIX|COUPLING|EIGEN|PURITY/ ) {
      # Switch to next parity
      next;
    }
  }
  close OUTG11;
} ##read_outg11

############################################################################
sub Reorder   #10/1/2004 8:00AM A.Kramida
############################################################################
{
  # Reorder parameter values and shell names according to the new shell sequence
  # given at the command line and stored in array @new_order

  # Reorder parameters
  for ( my $par = 1; $par <= 2; $par++ ) {
    # Reorder the Slater parameters for each config
    foreach my $nc ( sort {$a<=>$b} keys %{$params[$par-1]} ) {
      my $conf_name = $confs[$par-1]->{$nc};
      my $num_param = $params[$par-1]->{$nc}->{'num_param'};
      my %new_params = ();
      my %prev_shell = ('FG' => 0, 'ABGT' => 0);
      my $prev_type = '';
      for ( my $j = 1; $j <= $num_param; $j++ ) {
        my ($p, $t, $par_name) = @{$params[$par-1]->{$nc}->{"p$j"}};

        # Parameter types:
        # Eav      - 0
        # Fn(ll), ALPHA(l), BETA(l), GAMMA(l), Tn(l)  - 1
        # ZETA l   - 2
        # Fn(l'l") - 3
        # Gn(l'l") - 4
        #
        # Within each type, parameters are sorted in the following order:
        # type 1: l, F (rank), ALPHA, BETA, GAMMA, T(rank)
        # type 2: l
        # type 3: l', l", rank
        # type 4: l', l", rank
        #
        # For ALPHA, BETA, and GAMMA, l is not given in OUTG11 but it is implied;
        # For Tn(l), l is not given in OUTG11, but the shell name and occupation number are given.
        #
        my $key = $t;
        my $new_par_name = $par_name;
        if ( $t == 1 ) {
          if ( $par_name =~ /^([FG])(\d)\((\d){2}\)/ ) {
            my ($FG, $rank, $n_shell) = ($1,$2,$3);
            if ( ($n_shell > 1) && ($n_shell != $prev_shell{'FG'}) ) {
              $prev_shell{'FG'} = $n_shell - 1;
            }
            my $new_n_shell = $n_shell;
            if ( $n_shell <= $last_movable_shell ) {
              $new_n_shell = $new_order[$n_shell-1] - 1;
            }
            $key .= "$new_n_shell$FG$rank";
            $new_par_name = "$FG$rank\($new_n_shell$new_n_shell\)";
            if ( $prev_type =~ /A|B|G|T/ ) {
              $prev_type = 'FG';
              $prev_shell{'ABGT'}++;
            }

          } elsif ( $par_name =~ /^T(\d){0,1}\(([SPDFGH]) (\d)\)/ ) {
            my ($rank, $L, $occ) = ($1 + 0, $2, $3);
            $L = lc($L);
            if ( $prev_type =~ /^T(\d)\(([SPDFGH]) (\d)\)/ ) {
              my ($r1,$L1,$occ1) = ($1,$2,$3);
              if ( ($r1 > $rank) || (lc($L1) ne lc($L)) || ($occ1 != $occ) ) {
                $prev_shell{'ABGT'}++;
              }
            }
            $prev_type = $par_name;
            # Find this shell's number
            my $ns = 0;
            for ( my $n_shell = $prev_shell{'ABGT'}; $n_shell <= 7; $n_shell++ ) {
              my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell]};
              $nL =~ s/^\d+//; # strip the shell's number
              if ( ($nL eq $L) && ($occup == $occ) ) {
                $ns = $n_shell + 1;
                last;
              }
            }
            die "Shell designation is not recognized in T parameter for config. $conf_name in OUTG11" unless $ns;

            if ( ($ns > 1) && ($ns != $prev_shell{'ABGT'}) ) {
              $prev_shell{'ABGT'} = $ns - 1;
            }
            my $new_n_shell = $ns;
            if ( $ns <= $last_movable_shell ) {
              $new_n_shell = $new_order[$ns-1] - 1;
            }

            $key .= "$new_n_shell" . "T$rank";

          } else {
            # For ALPHA, BETA, GAMMA :
            # Find the first non-processed shell with equivalent electrons
            $par_name =~ /^(.)/;
            my $cur_type = $1; # The first letter of the parameter name
            if ( ($prev_type =~ /A|B|G|T/) && ($prev_type ge $cur_type) ) {
              $prev_shell{'ABGT'}++;
            }
            $prev_type = $cur_type;
            my $ns = 0;
            for ( my $n_shell = $prev_shell{'ABGT'}; $n_shell <= 7; $n_shell++ ) {
              my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell]};
              $nL =~ s/^\d+//; # strip the shell's number
              if ( ($nL eq 'p') && ($occup > 1) && ($occup < 5) ) {
                $ns = $n_shell + 1;
                last;
              }
              if ( ($nL eq 'd') && ($occup > 1) && ($occup < 9) ) {
                $ns = $n_shell + 1;
                last;
              }
              if ( ($nL eq 'f') && ($occup > 1) && ($occup < 14) ) {
                $ns = $n_shell + 1;
                last;
              }
              if ( ($nL =~ /fghiklmno/ ) && ($occup > 1) ) {
                $ns = $n_shell + 1;
                last;
              }
            }
            unless ($ns) {
              for ( my $n_shell = $prev_shell{'ABGT'}-1; $n_shell >= 0; $n_shell-- ) {
                my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell]};
                $nL =~ s/^\d+//; # strip the shell's number
                if ( ($nL eq 'p') && ($occup > 1) && ($occup < 5) ) {
                  $ns = $n_shell + 1;
                  last;
                }
                if ( ($nL eq 'd') && ($occup > 1) && ($occup < 9) ) {
                  $ns = $n_shell + 1;
                  last;
                }
                if ( ($nL eq 'f') && ($occup > 1) && ($occup < 14) ) {
                  $ns = $n_shell + 1;
                  last;
                }
                if ( ($nL =~ /fghiklmno/ ) && ($occup > 1) ) {
                  $ns = $n_shell + 1;
                  last;
                }
              }
            }
            unless ($ns) {
              die "Could not find shell with equivalent electrons for $par_name parameter\nfor config. $conf_name in OUTG11" ;
            }
            if ( ($ns > 1) && ($ns != $prev_shell{'ABGT'}) ) {
              $prev_shell{'ABGT'} = $ns - 1;
            }
            my $new_n_shell = $ns;
            if ( $ns <= $last_movable_shell ) {
              $new_n_shell = $new_order[$ns-1] - 1;
            }
            $key .= "$new_n_shell" . "H$par_name";
            # Parameter name does not change for ALPHA, BETA, GAMMA
          }
        } elsif ($t == 2) {
          # ZETA l
          my $ns = 0;
          if ( $par_name =~ /ZETA (\d)/ ) {
            $ns = $1;
          }
          unless ($ns) {
            die "Shell number not recognized for ZETA parameter of config. $conf_name in OUTG11";
          }

          my $new_n_shell = $ns;
          if ( $ns <= $last_movable_shell ) {
            $new_n_shell = $new_order[$ns-1] - 1;
          }
          $key .= "ZETA$new_n_shell";
          $new_par_name = "ZETA $new_n_shell";
        } elsif ( ($t == 3) || ($t == 4) ) {
          # Fn, Gn(l',l")
          if ( $par_name =~ /^([FG])(\d)\((\d)(\d)\)/ ) {
            my ($FG, $rank, $n1, $n2) = ($1, $2, $3, $4);
            my $new_n1 = $n1;
            if ( $n1 <= $last_movable_shell ) {
              $new_n1 = $new_order[$n1-1] - 1;
            }
            my $new_n2 = $n2;
            if ( $n2 <= $last_movable_shell ) {
              $new_n2 = $new_order[$n2-1] - 1;
            }
            if ( $new_n2 < $new_n1 ) {
              # Exchange $new_n2 and $new_n1 so that they go in increasing order
              my $n = $new_n1;
              $new_n1 = $new_n2;
              $new_n2 = $n;
            }
            $key .= "$new_n1$new_n2$rank";
            $new_par_name = "$FG$rank\($new_n1$new_n2\)";
          }
        } elsif ($t != 0) {
          die "Unrecognized parameter type for parameter $par_name of config $conf_name in ING11";
        }
        # Now we have the new sort key for each parameter.
        $new_params{$key} = [$p, $t, $new_par_name];
        # Substitute the old parameter values with new ones
        my $par_num = 0;
        foreach my $key (sort {$a cmp $b} keys %new_params) {
          $par_num++;
          $params[$par-1]->{$nc}->{"p$par_num"} = $new_params{$key};
        }
      }
    }

    # Reorder the CI section parameters for each pair of configs
    my $num_c1 = $#{$CI[$par-1]} + 1;
    for (my $nc1 = 1; $nc1 <= $num_c1; $nc1++ ) {
      #print "$nc1, $nc2\n";
      my $conf1 = $confs[$par-1]->{$nc1};
      my $num_c2 = $#{$CI[$par-1]->[$nc1-1]} + 1;
      for (my $nc2 = 1; $nc2 <= $num_c2; $nc2++ ) {
        my %new_params = ();
        my $conf2 = $confs[$par-1]->{$nc2};
        #print "$nc1, $nc2\n";
        next unless defined($CI[$par-1]->[$nc1-1]->[$nc2-1]->{'num_param'});
        my $num_param = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'num_param'};
        for ( my $i = 1; $i<=$num_param; $i++ ) {
          my $p = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$i-1];
          my $par_name = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'param_names'}->[$i-1];
          # Par. name format:
          # nc1 (A2), nc2 (A2), rank (I1), D or E (A1), shell numbers: L1, L2, L3, L4
          #
          # Parameter ordering in ING11:
          # 1) l1,l2,l3,l4
          # 2) First Direct, then Exchange
          # 3) rank
          if ( $par_name =~ /^(..)(\d)([DE])(\d)(\d)(\d)(\d)/ ) {
            my ($cc, $rank, $DE, $n1, $n2, $n3, $n4) = ($1, $2, $3, $4, $5, $6, $7);

            my ($class1,$transpose1) = &CI_class($par, $nc1, $nc2, $n1, $n2, $n3, $n4);

            my $new_n1 = $n1-1;
            if ( $n1 <= $last_movable_shell ) {
              $new_n1 = $new_order[$n1-1] - 1;
            }
            my $new_n2 = $n2-1;
            if ( $n2 <= $last_movable_shell ) {
              $new_n2 = $new_order[$n2-1] - 1;
            }
            my $new_n3 = $n3-1;
            if ( $n3 <= $last_movable_shell ) {
              $new_n3 = $new_order[$n3-1] - 1;
            }
            my $new_n4 = $n4-1;
            if ( $n4 <= $last_movable_shell ) {
              $new_n4 = $new_order[$n4-1] - 1;
            }

            my ($class2,$transpose2) = &CI_class($par, $nc1, $nc2, $new_n1+1, $new_n2+1, $new_n3+1, $new_n4+1);

            my ($nL1, $occup1) = @{$shells[$par-1]->{$nc1}->[$n1-1]};
            $nL1 =~ s/^\d+//g;
            my ($nL2, $occup2) = @{$shells[$par-1]->{$nc1}->[$n2-1]};
            $nL2 =~ s/^\d+//g;
            my ($nL3, $occup3) = @{$shells[$par-1]->{$nc2}->[$n3-1]};
            $nL3 =~ s/^\d+//g;
            my ($nL4, $occup4) = @{$shells[$par-1]->{$nc2}->[$n4-1]};
            $nL4 =~ s/^\d+//g;

            my $old_DE = $DE;
            if ( $new_n2 < $new_n1 ) {
              # Exchange $new_n2 and $new_n1 so that they go in increasing order
              my $n = $new_n1;
              $new_n1 = $new_n2;
              $new_n2 = $n;
              #$p = -$p;
              #if ( ($nL3 ne $nL4) &&
                #(($n1 eq $nL3) || ($nL1 eq $nL4) || ($nL2 eq $nL3) || ($nL2 eq $nL4)) ) {
              #if ( ($nL1 ne $nL2) && ($n3 != $n4) &&
              #  (($n1 == $n3) || ($n1 == $n4 ) || ($n2 == $n3) || ($n2 == $n4)) ){
              #    #$DE = ($DE eq 'E') ? 'D' : 'E';
              #}
              if ( ($class1 =~/^[6789]$|^10$/ ) ) {
                $DE = ($DE eq 'E') ? 'D' : 'E';
              }
            }
            if ( $new_n4 < $new_n3 ) {
              # Exchange $new_n4 and $new_n3 so that they go in increasing order
              my $n = $new_n3;
              $new_n3 = $new_n4;
              $new_n4 = $n;
              #$p = -$p;
              #if ( ($nL1 ne $nL2) &&
              #  (($nL1 eq $nL3) || ($nL1 eq $nL4) || ($nL2 eq $nL3) || ($nL2 eq $nL4)) ) {
              #if ( ($nL3 ne $nL4) && ($n1 != $n2) &&
              #  (($n1 == $n3) || ($n1 == $n4 ) || ($n2 == $n3) || ($n2 == $n4)) ){
                  #$DE = ($DE eq 'E') ? 'D' : 'E';
              #}
              if ( ($class1 =~/^[6789]$|^10$/) ) {
                $DE = ($DE eq 'E') ? 'D' : 'E';
              }
            }

            if ( $par_name =~ /^2F/ ) {
              $p = $p;
            }
            #if ( ($class1 != $class2) && (($class1 == 7) || ($class2 == 7)) && ($class2 != 0)) {
            #  $n1 = $n1;
            #  if ( $DE eq 'E' ) {
                #$p = -$p;
            #  }
            #}
            #if ( $class1 == $class2 ) {
            #  $DE = $old_DE;
            #}
            #if ( $transpose1 != $transpose2 ) {
              #$DE = ($DE eq 'E') ? 'D' : 'E';
            #}
            #if ( ($class1 != $class2) && (($class1 == 7) || ($class2 == 7)) ) {
              #$DE = ($DE eq 'E') ? 'D' : 'E';
            #}
            #print "$par\t$par_name\t$rank\t$old_DE\t$nL1$nL2$nL3$nL4\t$class1\t$class2\n";

            my $rank_par = 2 - $rank % 2;
            my $key = "$new_n1$new_n2$new_n3$new_n4$DE$rank";
            $new_params{$key} = $p;

          } else {
            die "Wrong format for CI parameter name $par_name in OUTG11, configs $conf1 and $conf2";
          }
        }

        # Store the reordered parameters in the same $CI hash
        my $par_num = 0;
        foreach my $key ( sort {$a cmp $b} keys %new_params ) {
          $par_num++;
          $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$par_num-1] = $new_params{$key};
        }
      }
    }
  }

  # Reorder shells
  for ( my $par = 1; $par <= 2; $par++ ) {
    my %new_shells = ();
    foreach my $nc ( sort {$a<=>$b} keys %{$shells[$par-1]} ) {
      $new_shells{$nc} = [];
      for ( my $n_shell = 0; $n_shell <= 7; $n_shell++ ) {
        my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell]};
        my $new_n_shell = $n_shell;
        if ( $n_shell+1 <= $last_movable_shell ) {
          $new_n_shell = $new_order[$n_shell] - 1;
        }

        $new_shells{$nc}->[$new_n_shell] = [$nL, $occup];
      }
    }
    # Replace old shells with new shells
    $shells[$par-1] = \%new_shells;
  }
} ##Reorder

############################################################################
sub shell_parity($)    #10/3/2004 9:32PM A.Kramida
############################################################################
{
  my $shell_name = shift;
  if ( $shell_name =~ /[pfhkmo]/ ) {
    return -1;
  } else {
    return 1;
  }
} ##shell_parity($)

############################################################################
sub CI_class($$$$$$$)    #10/3/2004 9:35PM A.Kramida
# parameters: parity number, N_conf1, N_conf2, ns11, ns12, ns21, ns22
############################################################################
{
  my ($par, $nc1, $nc2, @CI_shells) = @_;
  my ($n1, $n2, $n3, $n4) = @CI_shells;
  if ( !defined($shells[$par-1]->{$nc1}->[$n1-1]) ) {
    $n1 = $n1;
  }
  my ($nL1, $occup1) = @{$shells[$par-1]->{$nc1}->[$n1-1]};
  $nL1 =~ s/^\d+//g;
  my ($nL2, $occup2) = @{$shells[$par-1]->{$nc1}->[$n2-1]};
  $nL2 =~ s/^\d+//g;
  my ($nL3, $occup3) = @{$shells[$par-1]->{$nc2}->[$n3-1]};
  $nL3 =~ s/^\d+//g;
  my ($nL4, $occup4) = @{$shells[$par-1]->{$nc2}->[$n4-1]};
  $nL4 =~ s/^\d+//g;

  my ($i, $j) = ($n1 > $n3) ? (3,1) : (1,3);

  my $n_same_shell = 1;
  for ( my $k = 1; $k <=3; $k++ ) {
    for ( my $k1 = $k+1; $k1 <=4; $k1++ ) {
      if ( $CI_shells[$k-1] == $CI_shells[$k1-1] ) {
        $n_same_shell++;
      }
    }
  }
  my $class = 0;
  if ( $n_same_shell == 1 ) {
    if ($CI_shells[$i] < $CI_shells[$j-1]) {
      $class = 9;
    }
    if ($CI_shells[$j-1] < $CI_shells[$i]) {
      $class = 10;
    }
  } elsif ( $n_same_shell == 2 ) {
    if ( ($CI_shells[$i-1] == $CI_shells[$i]) && ($CI_shells[$i-1] < $CI_shells[$j-1]) ) {
      $class = 3;
    }
    if ( ($CI_shells[$j-1] == $CI_shells[$j]) && ($CI_shells[$j-1] < $CI_shells[$i]) ) {
      $class = 4;
    }
    if ( ($CI_shells[$j-1] == $CI_shells[$j]) && ($CI_shells[$j-1] > $CI_shells[$i]) ) {
      $class = 5;
    }
    if ( $CI_shells[$i-1] == $CI_shells[$j-1] ) {
      $class = 6;
    }
    if ( $CI_shells[$i] == $CI_shells[$j-1] ) {
      $class = 7;
    }
    if ( $CI_shells[$i] == $CI_shells[$j] ) {
      $class = ( $nL2 eq $nL4 ) ? 8 : 10;
    }
  } elsif ( $n_same_shell == 3 ) {
    $class = ($n1 == $n3) ? ($nL2 ne $nL4) ? 6 : 11 : 2;
  } elsif ( $n_same_shell == 4 ) {
    $class = 1;
  }

  my $transpose = 0;
  if ( ($class != 2) || ($n1 > $n3)  ) {
    if ( ($n3 != $n4) ||
         (($n1 != $n4) && ($n1 >= $n2) && (($n1 != $n3) || ($n2 > $n3)))
    ) {
      $transpose = 1;
    }
  }

  return ($class, $transpose);
} ##CI_class($$$$)
