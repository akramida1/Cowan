#!perl
use strict;
my $inp = shift;
my $out = shift;
my $param_str = shift;
if ( !$inp || !$out ) {
  print "\nUsage: \nsubst_exp_levs <input file> <output_file> [params=flags|noflags|Eav|F|G|Z|CI|none]\n\n";
  exit;
}
if ( !-f $inp ) {
  die "$inp file not found";
}
if ( !-f $out ) {
  die "$out file not found";
}
my ($noparams,$param_option, $lopt_levs) = ('','','');
$noparams = ($param_str eq 'NONE');    # To skip the parameter section
if ( $param_str =~ /^params *= *(flags|noflags|Eav|F|G|Z|CI|none)/i ) {
  $param_option = uc($1);
  $noparams = 1 if ($param_option eq 'NONE');
}
my @inp_levs = ({},{});
my @out_levs = ({},{});
my @inp_confs = ({},{});
my @out_confs = ({},{});
my @inp_params = ({},{});
my @out_params = ({},{});
my $max_unc = 1e6;

if ( $inp =~ /\.LEV$/i ) {
  $noparams = 1;
  $lopt_levs = 1;
}

&ReadRCE($inp,\@inp_levs,\@inp_params);
&ReadRCE($out,\@out_levs,\@out_params);
&Identify_RCE_levs();
&PrintOut();
&print_suspicious();
&print_duplicate();
&print_confs();

print "Done.\n";

############################################################################
sub ReadRCE($$$)   #9/28/2004 1:38PM A.Kramida
############################################################################
{
  my $file = shift;
  my $RCE_levs = shift;
  my $params = shift;
  open RCE, "<$file" or die "Could not open input file $file";

  print "Reading $file...\n";
  my $s = '';
  # Read levels of each parity
  for ( my $par = 0; $par<=1; $par++ ) {
    while ( (defined($s = <RCE>)) && ($s !~ /^ *([0-9-]+\.\d{0,6})([* ]+)(-{0,1}\d+\.\d{1,3}) +([0-9\/]+) +([0-9-]+) (.{6}) +(\S+) (\S+)( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) +(.{6}) +(\S+) (\S+)){0,1}/ ) ) {
      # Check if this file is actually an LOPT output file and not RCEINP/RCEOUT
      if ( $s =~ /^([^\t]*\t){9}[^\t]*$/ ) {
        $lopt_levs = 1;
        $noparams = 1;
        $max_unc = $param_str+0 if $param_str;
        close RCE;
        &Read_LOPT_Levs($file, $RCE_levs);
        return;
      }
      next;
    }
    my $Jprev = '';
    my $i = 0;
    do {
      $s =~ /^ *([0-9-]+\.\d{0,6})([* ]+)(-{0,1}\d+\.\d{1,3}) +(\d+)(\/(\d+))* +([0-9-]+) (.{6}) +(\S+) (\S+)( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}/;
      my ($Ee, $exp_c, $Ec, $J1,$J2, $p1, $conf1, $term1) = ($1,$2,$3,$4,$5,  $7,$8,"$9 $10");
      my ($A1,$A2,$A3,$A4);
      if ( $lopt_levs ) {
        $term1 = $10;
        $A1 = [$12,$13,$15];
        $A2 = [$17,$18,$20];
        $A3 = [$22,$23,$25];
        $A4 = [$27,$28,$30];
      } else {
        $A1 = [$12,$13,"$14 $15"];
        $A2 = [$17,$18,"$19 $20"];
        $A3 = [$22,$23,"$24 $25"];
        $A4 = [$27,$28,"$29 $30"];
      }
      $exp_c =~ s/\s+//g; # Trim spaces
      $conf1 =~ s/^\s+|\s+$//g;
      my $J = "$J1$J2";
      if ( $Jprev ne $J ) {
        $RCE_levs->[$par]->{$J} = {};
        $i = 0;
      }
      $Jprev = $J;
      $i++;
      $RCE_levs->[$par]->{$J}->{$i} = {};
      $RCE_levs->[$par]->{$J}->{$i}->{'Ee'} = $Ee;
      $RCE_levs->[$par]->{$J}->{$i}->{'exp_c'} = $exp_c;
      $RCE_levs->[$par]->{$J}->{$i}->{'Ec'} = $Ec;
      $RCE_levs->[$par]->{$J}->{$i}->{'vector'} = {"$conf1,$term1"=>$p1};

      for ( my $j = 1; $j<=4; $j++ ) {
        my $c = 0;
        my $expression = "\$c = \$A$j";
        eval($expression);
        my ($p, $conf, $term) = @{$c};

        # Store only components with amplitudes greater than 10
        if ( abs($c->[0]) >= 10  ) {
          $conf =~ s/^\s+|\s+$//g;
          $RCE_levs->[$par]->{$J}->{$i}->{'vector'}->{"$conf,$term"} = $p;
        } else {
          last;
        }
      }
    } while ( (defined($s = <RCE>)) && ($s !~ /PARAMETER/) );
    if ( !$noparams ) {
      my $conf = '';
      my $prev_conf = '';
      my $CI = 0;
      my $n_CI = 0;
      do {
        $s = <RCE>;
        if ((defined $s) && ($s !~ /^ {4}\d {4}\d {4}\d/) ) {
          if ( ($s =~ /^EAV (.{6}) ([ 0-9-]{4}) ([ .0-9-]{13}) (.+)$/) && !$CI  ) {
            my ($c, $flag, $value, $tail) = ($1,$2,$3,$4);
            $conf = $c;
            $params->[$par]->{$conf} = {};
            $params->[$par]->{$conf}->{'EAV'} = [$flag, $value, $tail];
          } elsif ( $s =~ /^.{63}-/ ) {
            $CI = 1;
            if ( $s =~ /^(.{10}) ([ 0-9-]{4}) ([ .0-9-]{13}) (.{27})(.{6}-.{2,6}) *$/ ) {
              my ($param, $flag, $value, $zeros, $confs) = ($1,$2,$3,$4, $5);
              if ( $confs ne $prev_conf ) {
                $n_CI = 0;
              }
              $prev_conf = $confs;
              $n_CI++;
              $params->[$par]->{$confs}->{$n_CI} = [$flag, $value, $zeros];
            } else {
              die "$file format error in CI section:\n$s";
            }
          } elsif ( ($s =~ /^(.{10}) ([ 0-9-]{4}) ([ .0-9-]{13}) (.+)$/) && !$CI  ) {
            my ($param, $flag, $value, $tail) = ($1,$2,$3,$4);
            $params->[$par]->{$conf}->{$param} = [$flag, $value, $tail];
          } else {
            die "$file format error on the following line:\n$s";
          }
        }
      } while ( (defined $s) && ($s !~ /^ {4}\d {4}\d {4}\d/) );
    } else {
      do {
        $s = <RCE>;
      } while ( (defined $s) && ($s !~ /^ {4}\d {4}\d {4}\d/) );
    }
  }
  close RCE;
} ##ReadRCE

############################################################################
sub Read_LOPT_Levs    #8/8/2006 8:47AM A.Kramida
############################################################################
{
  my $file = shift;
  my $levs_hash = shift;
  my $params = shift;
  open LEVS, "<$file" or die "Could not open input file $file";

  my $s = '';

  # Read levels of each parity
  while ( (defined($s = <LEVS>)) && ($s !~ /^([^\t]*\t){9}[^\t]*$/) ) {
    next;
  }
  my $J = '';
  my $i = 0;
  my $par = 0;
  my %parities = ();
  do {
    $s =~ /^([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/;
    my ($conf_term,$J,$Ee,$unc1,$unc2) = ($1,$2,$3,$4,$5);
    my ($conf1, $term1) = ($9,$10);
    chomp $term1;
    $Ee *= 0.001; # LOPT levels are in cm-1, RCE levels are in 10^3 cm-1

    if ( $unc2 <= $max_unc ) {
      my ($conf, $term) = ($conf_term,'');
      my $len = length($conf_term)-3;
      ($conf,$term) = (substr($conf_term,0,$len),substr($conf_term,$len,3));
      $conf =~ s/^\s+|\s+$//g;
      $term =~ s/^\s+|\s+$//g;
      my $exp_c = '';

      my $parity = ( $term =~ /[*]$/ ) ? 'o' : 'e';

      if ( $i == 0 ) {
        %parities = ($parity eq 'e') ? ('e' => 0, 'o' => 1) : ('e' => 1, 'o' => 0);
      }

      $par = $parities{$parity};
      if ( abs($Ee-3581.2) < 0.001 ) {
        $i = $i;
      }

      if ( !$J && ($J ne '0') ) {
        if ( $term1 =~ /^(\d+)([SPDFGHIKLMNOPQRSTUV])/ ) {
          my ($S,$L) = ($1,$2);
          $L = &get_L($L); # Convert to a number
          $S = ($S-1)*0.5;   # Convert from multiplicity to spin momentum
          for ( my $J1 = abs($L - $S); $J1 <= $L + $S; $J1 += 1.0 ) {
            $J = sprintf("%4.1f",$J1);
            $J =~ s/^\s+|\s+$//g;
            if ( $J =~ /^(\d+)\.5$/ ) {
              $J = ($1*2 + 1) . '/2';
            } else {
              $J = sprintf("%d",$J);
            }
            $i = &store_level($levs_hash, $par, $J, $Ee, $exp_c, $Ee, $conf1, $term1, 100);
          }
        }
      } else {
        $i = &store_level($levs_hash, $par, $J, $Ee, $exp_c, $Ee, $conf1, $term1, 100);
      }
    }

  } while ( defined($s = <LEVS>) );
  close LEVS;
} ##Read_LOPT_Levs

############################################################################
sub get_L($)    #8/8/2006 9:38AM A.Kramida
############################################################################
{
  # Convert a symbol for the L orbital momentum into a number
  my $L = shift;
  my $L_str = 'SPDFGHIKLMNOPQRSTUV';
  my $i = index($L_str,$L);
  return $i;
} ##get_L($)

############################################################################
sub store_level   #8/8/2006 9:31AM A.Kramida
############################################################################
{
  my ($levs_hash, $par, $J, $Ee, $exp_c, $Ec, $conf, $term, $amplitude) = @_;
  my $i = 0;
  if ( !defined($levs_hash->[$par]->{$J}) ) {
    $levs_hash->[$par]->{$J} = {};
    $i = 1;
  } else {
    my @Js = keys %{$levs_hash->[$par]->{$J}};
    $i = $#Js + 1;  # Number of already read levels with this J
    $i++;
  }
  $levs_hash->[$par]->{$J}->{$i} = {};
  $levs_hash->[$par]->{$J}->{$i}->{'Ee'} = $Ee;
  $levs_hash->[$par]->{$J}->{$i}->{'exp_c'} = $exp_c;
  $levs_hash->[$par]->{$J}->{$i}->{'Ec'} = $Ec;
  $levs_hash->[$par]->{$J}->{$i}->{'vector'} = {"$conf,$term"=>$amplitude};
  return $i;
} ##store_level

############################################################################
sub Identify_RCE_levs   #9/28/2004 4:15PM A.Kramida
############################################################################
{
  print "Identifying levels from $inp with levels from $out...\n";
  for ( my $par = 0; $par<=1; $par++ ) {
    foreach my $J ( keys %{$inp_levs[$par]} ) {
      foreach my $num_inp_lev (sort {$a<=>$b} keys %{$inp_levs[$par]->{$J}} ) {
        # Try to identify only the experimentally known levels in the input file
        next if $inp_levs[$par]->{$J}->{$num_inp_lev}->{'exp_c'};
        my $Ec = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'Ec'};
        my $Ee = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'Ee'};
        my $found = 0;

        # Go through the five RCE vector components and find corresponding RCG basis-state components
        my $minD = 1e17;
        my $best_match = 0;
        if ( ($J eq '1/2') && (abs($Ee -180.20192) < 0.0001) ) {
          $found = 0;
        }
        if ( abs($Ee-3581.2) < 0.001 ) {
          $found = 0;
        }
        foreach my $num_out_lev ( sort { my ($E1,$E2) = ($out_levs[$par]->{$J}->{$a}->{'Ec'}, $out_levs[$par]->{$J}->{$b}->{'Ec'});
                                        abs($E1-$Ec)<=>abs($E2-$Ec);
                                      } keys %{$out_levs[$par]->{$J}})
        {
          my $Ee_out = $out_levs[$par]->{$J}->{$num_out_lev}->{'Ee'};
          next unless $out_levs[$par]->{$J}->{$num_out_lev}->{'exp_c'} || (abs($Ee_out - $Ee) < 0.001) ; # Skip already mapped levels

          my ($D1,$D2) = (0,0);
          # Compare five amplitudes of the RCE vector with those of the RCG vector if these components are found there
          foreach my $conf_term_out (sort {my ($p1,$p2) = ($out_levs[$par]->{$J}->{$num_out_lev}->{'vector'}->{$a}, $out_levs[$par]->{$J}->{$num_out_lev}->{'vector'}->{$b});
                                        abs($p2)<=>abs($p2);
                                      } keys %{$out_levs[$par]->{$J}->{$num_out_lev}->{'vector'}} ) {
            my $p_out = $out_levs[$par]->{$J}->{$num_out_lev}->{'vector'}->{$conf_term_out};
            if (!defined($inp_levs[$par]->{$J}->{$num_inp_lev}->{'vector'}->{$conf_term_out})) {
              $D1 += $p_out*$p_out;
              $D2 += $p_out*$p_out;
              next;
            }
            my $p_inp = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'vector'}->{$conf_term_out};
            $D1 += ($p_inp-$p_out)*($p_inp-$p_out);
            $D2 += ($p_inp+$p_out)*($p_inp+$p_out);
            if ( (abs($p_inp) > 70) && (abs($p_out) > 70) ) {
              # Stop search if both vectors have more than 50% of the same conf and term
              $best_match = $num_out_lev;
              $found = 1;
              last;
            }
          }
          if ( !$found ) {
            # Add amplitudes of the inp vector that are not present in the out vector
            my $j = 0;
            foreach my $conf_term_inp ( keys %{$inp_levs[$par]->{$J}->{$num_inp_lev}->{'vector'}} ) {
              $j++;
              if (!defined($out_levs[$par]->{$J}->{$num_out_lev}->{'vector'}->{$conf_term_inp})) {
                my $p_inp = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'vector'}->{$conf_term_inp};
                $D1 += $p_inp*$p_inp;
                $D2 += $p_inp*$p_inp;
              }
            }
          } else {
            last;
          }

          if ( $D1 > $D2 ) {
            $D1 = $D2;
          }
          if ( $D1 < $minD ) {
            $minD = $D1;
            $best_match = $num_out_lev;
          }
        }
        if ( !$best_match ) {
          print "Warning! Could not identify $inp level \n      (parity = " .
            ($par+1) . ", J = $J, Ec = $Ec, Ee = $Ee) with a $out level.\n";
        } else {
          $out_levs[$par]->{$J}->{$best_match}->{'exp_c'} = '';
          $out_levs[$par]->{$J}->{$best_match}->{'Ee'} = $Ee;
          $inp_levs[$par]->{$J}->{$num_inp_lev}->{'map'} = $best_match;
        }
      }
    }
  }
} ##Identify_RCE_levs

############################################################################
sub PrintOut()   #9/28/2004 1:38PM A.Kramida
############################################################################
{
  my $file = $out;
  my $tmp_file = "tmp_$out";
  open RCE, "<$file" or die "Could not open input file $file";
  open TMP, ">$tmp_file" or die "Could not create temporary file $tmp_file";

  print "Writing $tmp_file...\n";
  my $s = '';
  # Read levels of each parity
  for ( my $par = 0; $par<=1; $par++ ) {
    while ( (defined($s = <RCE>)) && ($s !~ /^ *([0-9-]+\.\d{0,6})([* ]+)(\d+\.\d{1,3}) +([0-9\/]+) +([0-9-]+) (.{6}) +(\S+) (\S+)( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) +(.{6}) +(\S+) (\S+)){0,1}/ ) ) {
      print TMP $s;
      next;
    }
    my $Jprev = '';
    my $i = 0;
    do {
      $s =~ /^ *([0-9-]+\.\d{0,6})([* ]+)(\d+\.\d{1,3}) +(\d+)(\/(\d+))* (.*)$/;
      my ($Ee, $exp_c, $Ec, $J1,$J2, $vect) = ($1,$2,$3,$4,$5,  $7);
      my $J = "$J1$J2";
      if ( $Jprev ne $J ) {
        $i = 0;
      }
      $Jprev = $J;
      $i++;
      $exp_c = $out_levs[$par]->{$J}->{$i}->{'exp_c'};
      if ( $exp_c ) {
        print TMP $s;
      } else {
        $Ee = $out_levs[$par]->{$J}->{$i}->{'Ee'};

        $Ee = sprintf("%13.6f",$Ee);
        $Ec = sprintf("%10.3f",$Ec);
        $J = ($J =~ /\/2/) ? sprintf("%6s",$J) : sprintf("%4s",$J);
        print TMP "$Ee $Ec$J $vect\n";
      }
    } while ( (defined($s = <RCE>)) && ($s !~ /PARAMETER/) );

    if ( !$noparams ) {
      print TMP $s;
        # Substitute parameter values and flags from the input file into the output file
      my $conf = '';
      my $prev_conf = '';
      my $CI = 0;
      my $n_CI = 0;
      do {
        $s = <RCE>;
        if ( (defined $s) && ($s !~ /^ {4}\d {4}\d {4}\d/) ) {
          if ( ($s =~ /^EAV (.{6}) ([ 0-9-]{4}) ([ .0-9-]{13}) (.+)$/) && !$CI  ) {
            # Single-configuration EAV parameters
            my ($c, $flag, $value, $tail) = ($1,$2,$3,$4);
            $conf = $c;
            if ( defined $inp_params[$par]->{$conf}->{'EAV'} ) {
              my ($flag1, $value1, $tail1) = @{$inp_params[$par]->{$conf}->{'EAV'}};
              if ( ($param_option eq 'FLAGS') || $param_option && ($param_option ne 'EAV') ) {
                $value1 = $value; # leave the value unchanged
              }
              if ( ($param_option eq 'NOFLAGS') || $param_option && ($param_option ne 'FLAGS') ) {
                $flag1 = $flag;  # Leave the flag unchanged
              }
              print TMP "EAV $conf $flag1 $value1 $tail1\n";
            } else {
              print TMP $s;
            }
          } elsif ( $s =~ /^.{63}-/ ) {
            # CI parameters
            $CI = 1;
            if ( $s =~ /^(.{10}) ([ 0-9-]{4}) ([ .0-9-]{13}) (.{27})(.{6}-.{2,6}) *$/ ) {
              my ($param, $flag, $value, $zeros, $confs) = ($1,$2,$3,$4, $5);
              if ( $confs ne $prev_conf ) {
                $n_CI = 0;
              }
              $prev_conf = $confs;
              $n_CI++;
              if ( defined $inp_params[$par]->{$confs}->{$n_CI} ) {
                my ($flag1, $value1, $zeros1) = @{$inp_params[$par]->{$confs}->{$n_CI}};
                if ( ($param_option eq 'FLAGS' ) || $param_option && ($param_option ne 'CI') ) {
                  $value1 = $value; # leave the value unchanged
                }
                if ( ($param_option eq 'NOFLAGS') || $param_option && ($param_option ne 'FLAGS') ) {
                  $flag1 = $flag;  # Leave the flag unchanged
                }
                print TMP "$param $flag1 $value1 $zeros1$confs\n";
              } else {
                print TMP $s;
              }
            } else {
              die "$file format error in CI section:\n$s";
            }
          } elsif ( ($s =~ /^(.{10}) ([ 0-9-]{4}) ([ .0-9-]{13}) (.+)$/) && !$CI  ) {
            # Single-configuration parameters
            my ($param, $flag, $value, $tail) = ($1,$2,$3,$4);
            if ( defined $inp_params[$par]->{$conf}->{$param} ) {
              my ($flag1, $value1, $tail1) = @{$inp_params[$par]->{$conf}->{$param}};
              if ( ($param_option eq 'FLAGS') || $param_option && ($param !~ /$param_option/i) ) {
                $value1 = $value; # leave the value unchanged
              }
              if ( ($param_option eq 'NOFLAGS') || $param_option && ($param_option ne 'FLAGS') ) {
                $flag1 = $flag;  # Leave the flag unchanged
              }
              print TMP "$param $flag1 $value1 $tail1\n";
            } else {
              print TMP $s;
            }
          } else {
            die "$file format error on the following line:\n$s";
          }
        }
      } while ( (defined $s) && ($s !~ /^ {4}\d {4}\d {4}\d/) );
      print TMP $s;
    } else {
      print TMP $s;
      do {
        $s = <RCE>;
        print TMP $s;
      } while ( (defined $s) && ($s !~ /^ {4}\d {4}\d {4}\d/) );
    }
  }
  close RCE;
  close TMP;
  rename($tmp_file,$out) or die "Could not rename $tmp_file to $out";
} ##PrintOut

############################################################################
sub print_suspicious    #7/19/2006 11:19AM A.Kramida
############################################################################
{
  for ( my $par = 0; $par<=1; $par++ ) {
    my $num_susp = 0;
    foreach my $J ( keys %{$inp_levs[$par]} ) {
      foreach my $num_inp_lev (sort {$a<=>$b} keys %{$inp_levs[$par]->{$J}} ) {
        # Look only at the experimentally known levels in the input file
        next if $inp_levs[$par]->{$J}->{$num_inp_lev}->{'exp_c'};
        my $Ec = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'Ec'};
        my $Ee = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'Ee'};
        my $map = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'map'};
        my $inp_vector = $inp_levs[$par]->{$J}->{$num_inp_lev}->{'vector'};
        my $out_vector = $out_levs[$par]->{$J}->{$map}->{'vector'};
        my $Ec_out = $out_levs[$par]->{$J}->{$map}->{'Ec'};
        my $Ee_out = $out_levs[$par]->{$J}->{$map}->{'Ee'};
        my ($p_inp,$p_out, $t_inp,$t_out);

        my @v_inp = ();
        foreach my $conf_term (sort {my ($p1,$p2) = ($inp_vector->{$a},$inp_vector->{$b});
                                 abs($p2)<=>abs($p1);
                               } keys %{$inp_vector} )
        {
          $p_inp = $inp_vector->{$conf_term};
          $t_inp = $conf_term;
          push @v_inp, [$p_inp,$t_inp];
          #last;
        }
        my @v_out = ();
        foreach my $conf_term (sort {my ($p1,$p2) = ($out_vector->{$a},$out_vector->{$b});
                                 abs($p2)<=>abs($p1);
                               } keys %{$out_vector})
        {
          $p_out = $out_vector->{$conf_term};
          $t_out = $conf_term;
          push @v_out, [$p_out,$t_out];
          #last;
        }
        if ( (abs($v_inp[0]->[0]) > 70) && (abs($v_out[0]->[0]) > 70) ) {
          next;
        }
        if ( ($v_inp[0]->[1] eq  $v_out[0]->[1]) && ($v_inp[1]->[1] . '' eq  $v_out[1]->[1] . '') ) {
          next;
        }

        if ( !$num_susp ) {
          my $cs = $par+1;
          print "\nSuspicious levels, conf. set \#$cs:\n";
        } else {
          print "\n";
        }
        $num_susp++;
        print sprintf("J=%4s input : %13.6f%10.3f", $J, $Ee, $Ec), &vect_to_str(2,@v_inp), "\n";
        print sprintf("J=%4s output: %13.6f%10.3f", $J, $Ee_out, $Ec_out), &vect_to_str(2,@v_out),"\n";
      }
    }
  }
} ##print_suspicious

############################################################################
sub vect_to_str(@)    #7/20/2006 8:15AM A.Kramida
############################################################################
{
  my ($max_comp,@vect) = @_;
  my $n_comp = $#vect + 1;
  if ( $n_comp > $max_comp ) {
    $n_comp = $max_comp;
  }
  my $s = '';
  for ( my $i = 0; $i < $n_comp; $i++ ) {
    my ($p, $t) = @{$vect[$i]};
    my ($conf,$term) = split(/,/, $t);
    $conf = sprintf("%-7s",$conf);
    $term = sprintf("%-8s",$term);
    {use integer;
      $p = $p*$p/100;
    }
    $p = sprintf("%4d",$p);
    $s .= "$p% $conf$term";
  }
  return $s;
} ##vect_to_str($$$)

############################################################################
sub print_duplicate    #7/19/2006 11:19AM A.Kramida
############################################################################
{
  for ( my $par = 0; $par<=1; $par++ ) {
    my $num_dups = 0;
    foreach my $J ( keys %{$inp_levs[$par]} ) {
      my $prev_Ee = -1000000;
      my $prev_Ec = -1000000;
      my $prev_n = 0;
      my @prev_v = ();
      my $n_dup_E = 0;
      foreach my $num_out_lev (sort {$out_levs[$par]->{$J}->{$a}->{'Ee'}<=>$out_levs[$par]->{$J}->{$b}->{'Ee'}} keys %{$out_levs[$par]->{$J}} ) {
        # Look only at the experimentally known levels in the output file
        next if $out_levs[$par]->{$J}->{$num_out_lev}->{'exp_c'};
        my $Ec = $out_levs[$par]->{$J}->{$num_out_lev}->{'Ec'};
        my $Ee = $out_levs[$par]->{$J}->{$num_out_lev}->{'Ee'};
        my @v_out = ();
        if ( ($Ee == $prev_Ee) && ($Ec != $prev_Ec) ) {
          $n_dup_E++;
          my $out_vector = $out_levs[$par]->{$J}->{$num_out_lev}->{'vector'};
          foreach my $conf_term (sort {my ($p1,$p2) = ($out_vector->{$a},$out_vector->{$b});
                                  abs($p2)<=>abs($p1);
                                } keys %{$out_vector})
          {
            my $p = $out_vector->{$conf_term};
            my $t = $conf_term;
            push @v_out, [$p,$t];
          }
          if ( !$num_dups ) {
            my $cs = $par+1;
            print "\nDuplicate exp. energies, conf. set \#$cs:\n";
          } else {
            print "\n";
          }
          $num_dups++;
          if ( $n_dup_E == 1 ) {
            my $out_vector = $out_levs[$par]->{$J}->{$prev_n}->{'vector'};
            foreach my $conf_term (sort {my ($p1,$p2) = ($out_vector->{$a},$out_vector->{$b});
                                    abs($p2)<=>abs($p1);
                                  } keys %{$out_vector})
            {
              my $p = $out_vector->{$conf_term};
              my $t = $conf_term;
              push @prev_v, [$p,$t];
            }
            print sprintf("J=%4s%13.6f%10.3f", $J, $prev_Ee, $prev_Ec), &vect_to_str(2,@prev_v), "\n";
          }
          print sprintf("J=%4s%13.6f%10.3f", $J, $Ee, $Ec), &vect_to_str(2,@v_out),"\n";
        } else {
          $n_dup_E = 0;
          $prev_n = 0;
          @prev_v = ();
        }
        $prev_Ee = $Ee;
        $prev_Ec = $Ec;
        $prev_n = $num_out_lev;
      }
    }
  }
} ##print_duplicate

############################################################################
sub print_confs   #7/20/2006 11:24AM A.Kramida
############################################################################
{
  my @confs = ({},{});
  for ( my $par = 0; $par<=1; $par++ ) {
    foreach my $J ( keys %{$inp_levs[$par]} ) {
      foreach my $num_out_lev (sort {$out_levs[$par]->{$J}->{$a}->{'Ee'}<=>$out_levs[$par]->{$J}->{$b}->{'Ee'}} keys %{$out_levs[$par]->{$J}} ) {
        # Look only at the experimentally known levels in the output file
        next if $out_levs[$par]->{$J}->{$num_out_lev}->{'exp_c'};
        my $Ee = $out_levs[$par]->{$J}->{$num_out_lev}->{'Ee'};
        my $out_vector = $out_levs[$par]->{$J}->{$num_out_lev}->{'vector'};
        my ($p,$t);
        foreach my $conf_term (sort {my ($p1,$p2) = ($out_vector->{$a},$out_vector->{$b});
                                abs($p2)<=>abs($p1);
                              } keys %{$out_vector})
        {
          $p = $out_vector->{$conf_term};
          $t = $conf_term;
          last;
        }
        my ($conf,$term) = split(/,/,$t);
        if ( (!defined $confs[$par]->{$conf}) || ($Ee < $confs[$par]->{$conf})  ) {
          $confs[$par]->{$conf} = $Ee
        }
      }
    }
    # Print sorted list of experimentally known configs
    my $cs = $par+1;
    print "\nKnown configs, conf. set \#$cs:\n";
    foreach my $conf (sort {$confs[$par]->{$a}<=>$confs[$par]->{$b}} keys %{$confs[$par]}) {
      print "$conf\n";
    }
  }
} ##print_confs