#!perl
use strict;
# This program reads fitted parameters, eigenvectors, and experimental energies from RCEOUT,
# varies parameters randomly according to stdev given in OUTE, also randomly varies the E2
# transition integrals from ING11, substitutes these random trials into ING11, saves this
# new ING11 in 10 separate trial subdirectories, runs RCG in each of them, reads resulting
# eigenvectors and transition data from the corresponding OUTG11, identifies the eigenvectors
# with the original ones, calculates the dispersion and average shift of log(S) for each
# transition, and saves the line list with these data in the output file tp_M1E2.txt.
#
# Syntax: perl tp_M1E2.pl
#
# Sample ING11 with necessary RCG options is in C:\work\cowan\ti-like\fe5\Forbid
#
use vars qw{@parities @energies @map_RCG_RCE @RCE_lev @vectors @basis};
require 'conv_cowan.pl';
require 'vacair.pl';

  # Define global parameters
my $num_trials           = 0; #100;      #10000;  #Must be an even number
my $store_trials         = 100;      #10000;  # Number of trial results to store for each transition
my $R2_variance          = 0.15;   # Fractional standard deviation
my $vary_fitted_params   = 1;      # Set to 1 if fitted RCE parameters are to be varied
my $var_grouped_params   = 1;      # Set to 0 if parameter grouping is to be ignored
my $fixed_param_variance = 0.00;   # Set to 0 if fixed RCE parameters are not to be varied
my $mearged_list         = 1;      # set to 1 if a mearged M1+E2 list is desired;
                                   # otherwise, all M1 and E2 transition lists will be given separately
my $min_mixed_fraction   = 0.02;   # Transition will be omitted if its A-value is smaller than
                                   # this fraction of the total M1+E2 A-value
                                   # This setting has effect only if $meaged_list != 0.
my $min_BF               = 0.00010;  # Transitions with smaller radiative branching fraction will be omitted
my $out_tr_file          = "Pr10_tp_M1E2_mearged_100.txt";
my $print_distribution   = 0;      # Set to 0 if no printing of statistical distribution of A values is desired
my $print_dBF            = 1;      # Set to 0 to suppress printing stdev of branching fractions
my $print_cf_trials      = 0;      # Set to 0 to suppress printing trial data for CF
my $read_only            = 0;      # Set to even number N>0 if RCG calculations have already been made in N trial subdirectories
my $test_random          = 0;      # debugging randomization routines
my $sort_trials          = 1;      # Set to 0 if no sorting of trial data is desired.
                                   # WARNING: sorting trials destroys correspondence between A and cf in trials
                                   # as well as the branching fractions in trials.

$store_trials = $num_trials if $store_trials > $num_trials;
$read_only = $num_trials if $read_only > $num_trials;

# Read initial transition arrays from OUTG11
open OUT_TR_FILE, ">$out_tr_file" or die "Could not open output file $out_tr_file" unless $test_random;

my @transitions = ();
my @trans_keys = ();
my %sum_A_hash = ();
my %dev_hash_logA = ();
my %dev_hash_A = ();
my %dev_hash_cf = ();
my %dev_hash_param = ();
my %dev_hash_param_groups = ();
my %dev_hash_param_fixed = ();
my @map_shells = ({},{});
my @param_stats = ();
# Prepare the trial subdirectories
unless ( $test_random ) {
  for ( my $i = 1; $i <= $num_trials; $i++ ) {
    my $dir = "./trial$i";
    unless ( -d $dir ) {
      mkdir $dir or die "Unable to create directory $dir";
    }
    if ( $i > $read_only ) {
      `copy LEVELS1 $dir` if (-f 'LEVELS1');
    }
  }
}

# Read $Emax1,$Emax2,$ING_header, R2 values, and $ING_footer from ING11
my ($Emax1,$Emax2, $ING_header, $ING_footer, @R2) = &read_R2();
my @R2_0;
for ( my $i = 0; $i <= $#R2; $i++ ) {
  $R2_0[$i] = $R2[$i]->[1];
}

# Vary R2 and write the new ING11 files in trial subdirectories
my @R2_trial = ();
my $i1 = 0;
for ( my $i = 1; $i <= $num_trials/2; $i++ ) {
  &vary_R2();
  for ( my $k = 0; $k <= 1; $k++ ) {
    $i1++;
    if ( $i1 > $read_only ) {
      my $dir = "./trial$i1";
      chdir $dir unless $test_random;
      my $ix = $i1;
      for ( my $j = 0; $j <= $#R2; $j++ ) {
        $R2[$j]->[1] = $R2_trial[$j]->[$k];

        if ( ($R2_variance != 0) && $print_distribution && ($R2_0[$j] != 0) ) {
          my $k = ($R2[$j]->[1] / $R2_0[$j] - 1)/$R2_variance*10;
          push(@param_stats,$k/10);
          $k = &get_rounded_key($k);
          $dev_hash_param{$k}++;
        }

       # if ( $j != 2 ) {
       #   $R2[$j]->[1] = $R2_0[$j];
       # } else {
       #   $ix++ if $i1 >= ($num_trials/2 + 1);
       #   my $start_value = $R2_0[$j] - $R2_0[$j] * $R2_variance;
       #   $R2[$j]->[1] = $start_value + $R2_0[$j] * 2*$R2_variance/$num_trials * ($ix-1);
       # }
      }
      #print "$i1\t$ix\t" . $R2[0]->[1] . "\n";
      unless ( $test_random ) {
        &print_R2();
        chdir '../';
      }
    }
  }
}

open OUTG11, "<OUTG11" or die "Could not open input file OUTG11";

# Initialize global variables
&init_vars();
&read_RCG_options();
&read_in36();
# Read configurations from OUTG11 for the first parity only
&read_confs();
&read_ING11_params();
# Read the basis state definitions printed by CALCFC in OUTG11
my $s = &read_basis();
# Read Slater parameter labels and values of the first parity ------------
# Read eigenvalues and eigenvectors of the first parity
$s = &read_OUTG11_params(1,$s);
# Continue to read the OUTG11 file. Find and read the LS basis state labels printed by ENERGY
# If second parity present, read Slater parameter labels and values of the first parity,
# and read eigenvalues and eigenvectors of the second parity
$s = &read_basis_labels($s,1);

#&Reorder_Params();
&ReadOUTE();
my ($params,$params_CI) = &get_params();
my @params = @{$params->[0]->{'sq'}};
#my @CI = @{$params_CI->[0]};
&ReadRCE();
unless ($test_random) {
  &sort_vectors();
  &Identify_RCE_levs($Emax1,$Emax2);
  &ReadTransitions_M1E2(0,$s);
}
close (OUTG11);

# Run RCG and collect the new transition data in trial subdirectories
for ( my $i = 1; $i <= $num_trials; $i++ ) {
  unless ($test_random) {
    print "\nTrial $i:\n";
    my $dir = "./trial$i";
    chdir $dir;
  }
  &vary_params() if ( ($vary_fitted_params || $fixed_param_variance) && ($i > $read_only) );
  unless ( $test_random ) {
    `rcg` if $i > $read_only;
    open OUTG11, "<OUTG11" or die "Could not open input file OUTG11";
    $s = undef;
    if ( $vary_fitted_params || $fixed_param_variance ) {
      $s = &read_OUTG11_params(1,$s);
      $s = &read_basis_labels($s,0);
      &Identify_RCE_levs($Emax1,$Emax2);
    }
    &ReadTransitions_M1E2($i,$s);
    close (OUTG11);
    chdir '../';
  }
}

unless ( $test_random ) {
  # Calculate variations of S and wite the results
  #&print_trans($M1_trans,$E2_trans);
  my @types = ($mearged_list ? ('') : ('M1','E2'));
  for (my $i = 0; $i <= $#types; $i++) {
    &PrintTransitions($types[$i],$i);
  }
  close OUT_TR_FILE;
  print "Done.\n";
}

&print_distribution() if ($print_distribution && $store_trials);

############################################################################
sub read_R2() {  #11/19/2013 8:00AM
############################################################################
  my @R2 = ();
  my $ING_header = '';
  my $ING_footer = '';
  my $line_ING11 = 0;
  my $s;

  open(ING11,'<ing11') or die('Cannot open ING11 file.');
  my ($Emax1,$Emax2);
  for ( my $i = 1; $i <= 3; $i++ ) {
    $s = <ING11>;
    $ING_header .= $s;
    if ( $s =~ /^    0 {55}(.{10})(.{10})/ ) {
      ($Emax1,$Emax2) = ($1,$2);
      last;
    }
  }
  print "Reading E2 transition integrals...";
  while ( (defined($s = <ING11>)) && ($s !~ /^.{54}\/\/R2\/\//i) ) {
    $ING_header .= $s;
    next;
  }
  while ( $s =~ /^(.{38})(.{12})(.{4})\/\/R2\/\/(.{4})(.{6})HR(.{4})(.{4})/i ) {
    my ($s1,$r2,$s2,$s3,$r21,$n1,$n2) = ($1,$2,$3,$4,$5,$6,$7);
    $r2 =~ s/^\s+|\s+$//g;
    push(@R2,[$s1,$r2,$s2,$s3,$r21,$n1,$n2]);
    $s = <ING11>;
  }
  do {
    $ING_footer .= $s;
    $s = <ING11>;
  } while ( defined($s) );

  close(ING11);
  print "Done.\n";
  return ($Emax1,$Emax2,$ING_header, $ING_footer, @R2);
} ##read_R2()

############################################################################
sub vary_R2() {  #11/19/2013 10:05AM
############################################################################
  for ( my $i = 0; $i <= $#R2; $i++ ) {
    #my ($s1,$r2,$s2,$s3,$r21,$n1,$n2) = $R2[$i];
    my $r2 = $R2_0[$i];
    # Set $r2 to a normally distributed random number with an expectation value $r2
    # and a standard deviation $R2_variance*$r2
    $R2_trial[$i] = [];
    ($R2_trial[$i]->[0],$R2_trial[$i]->[1]) = &get_normal_random_trial($r2, abs($R2_variance*$r2));;
  }
} ##vary_R2()


sub getUniformPair() {
# This function returns two numbers (0 to 1)
# from a uniform distribution
  my $a = 1 - rand(1);  # Get a random value in the interval (0,1] (excluding 0 but including 1
  my $b = 1 - rand(1);  # Get a random value in the interval (0,1] (excluding 0 but including 1
  return ($a,$b);
}

sub getNormalPair($$) {
  # This function accepts two number (0 to 1] from a
  # uniform distribution and returns two numbers from
  # the standard normal distribution. (mean 0, variance 1)

  my ($a,$b) = @_;
  my $pi = 3.14159265359;

  # Box-Muller Transformation
  my $x = sqrt(-2 * log($a)) * cos(2*$pi*$b);
  my $y = sqrt(-2 * log($a)) * sin(2*$pi*$b);

  return ($x,$y);
}

sub scale($$$) {
# convert the standard normal distribution to our normal
# distribution with standard deviation sigma and expectation value mu

  my ($x0, $mu, $sigma) = @_;
  $x0 *= $sigma;
  $x0 += $mu;
  return $x0;
}

############################################################################
sub get_normal_random_trial($$) {   #11/19/2013 12:08PM
############################################################################
  my ($mean,$stdev) = @_;
  my ($a,$b);
  # obtain two random numbers
  while (!$a || !$b) { ($a,  $b) = getUniformPair(); }

  # transform to normal pair
  ($a,$b) = getNormalPair($a,$b);

  # transform to our distribution
  $a = scale($a,$mean,$stdev);
  $b = scale($b,$mean,$stdev);
  my $x = 1 - rand(1);

  return ($a , $b);
} ##get_normal_random_trial($$) {

############################################################################
sub print_R2() {    #11/19/2013 12:15PM
############################################################################
  open(OUT, ">ING11") or die "Unable to create ING11";
  print OUT $ING_header;
  for ( my $i = 0; $i <= $#R2; $i++ ) {
    my ($s1,$r2,$s2,$s3,$r21,$n1,$n2) = @{$R2[$i]};
    $r2 = sprintf("%12.5f",$r2);
    my $dL = length($r2) - 12;
    if ( $dL > 5 ) {
      $r2 = sprintf("%12.2e",$r2);
    } elsif ( $dL > 0 ) {
      my $n_dig = 5-$dL;
      $r2 = sprintf("%12.${n_dig}f",$r2);
    }
    print OUT $s1, $r2, $s2,'//R2//',$s3,$r21,'HR',$n1,$n2,"\n";
  }
  print OUT $ING_footer;
  close(OUT) or die "Unable to write to ING11";
} ##print_R2() {

############################################################################
sub S_from_gA($$$) {    #11/19/2013 3:01PM
############################################################################
  my ($gA,$WL,$type) = @_;
  my $S = 0;
  if (($type eq '') || ($type eq 'E1')) {
    $S = $gA*$WL*$WL*$WL / 2.0261269E+18;
  } elsif ($type eq 'M1') {
    $S = $gA*$WL*$WL*$WL / 269735e8;
  } elsif ($type eq 'E2') {
    $S = $gA*$WL*$WL*$WL*$WL*$WL / 1.11995E+18;
  } elsif ($type eq 'M2') {
    $S = $gA*$WL*$WL*$WL*$WL*$WL / 14909714e6;
  }
  return $S;
} ##S_from_gA($$$) {

############################################################################
sub gA_from_S($$$) {    #11/19/2013 3:59PM
############################################################################
  my ($S,$WL,$type) = @_;
  my $gA = 0;
  if (($type eq '') || ($type eq 'E1')) {
    $gA = $S*2.0261269E+18/($WL*$WL*$WL);
  } elsif ($type eq 'M1') {
    $gA = $S*269735e8/($WL*$WL*$WL);
  } elsif ($type eq 'E2') {
    $gA = $S*1.11995E+18/($WL*$WL*$WL*$WL*$WL);
  } elsif ($type eq 'M2') {
    $gA = $S*14909714e6/($WL*$WL*$WL*$WL*$WL);
  }
  return $gA;
} ##gA_from_S($$$) {

############################################################################
sub ReadTransitions_M1E2($$) {  #12/10/2013 11:33AM A.Kramida
############################################################################
  my ($num_trial,$s) = @_;
  print "Reading transitions...";

  #@map_shells = ({},{});
  my $sp_data = undef;
  while ( !defined($sp_data = &spectrum_header($s)) && defined($s = <OUTG11>) && !defined($sp_data = &spectrum_header($s)) ) {
    next;
  }
  while ( 1 ) {
    last unless defined $s;
    if ($s =~ /^ *\d{0,5}$/) {
      $s = <OUTG11>;
      last unless defined $s;
    }
    last unless defined($sp_data = &spectrum_header($s));
    my ($trans_type,$spectrum_type) = @{$sp_data};
    last unless $trans_type;
    while ( (defined($s = <OUTG11>)) && ($s !~ /0           E      J   CONF              EP      JP  CONFP           DELTA E   LAMBDA/) ) {
      next;
    }
    $s = <OUTG11>; # skip the blank line

    #while ( (defined($s = <OUTG11>)) && ($s !~ /^. +\* \* \* +([0-9.-]+) +([0-9.]+) (.{3}) (.{8})  (.{8}) (.{8}) +\* \* \*/) ) {
    #  next;
    #}
    last unless defined $s;
    while ( defined($s = <OUTG11>) && ($s =~ /^.{5}(.{11})(.{5})(.{3})(.{9})(.{13})(.{5})(.{3})(.{9})(.{13})(.{12}).{25} (.{10}) (.{9})/) ) {
      my ($E,$J,$nc,$term,$EP,$JP,$ncP,$termP,$dE,$lambda,$gA,$cf) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
      my ($n1,$n2);
      if ($s =~ /^.{5}(.{11})(.{5})(.{3})(.{9})(.{13})(.{5})(.{3})(.{9})(.{13})(.{12}).{25} (.{10}) (.{9}) (.{5}) (.{5})/) {
        ($n1,$n2) = ($13,$14);
      }
      foreach ($E,$J,$nc,$term,$EP,$JP,$ncP,$termP,$dE,$lambda,$gA,$cf,$n1,$n2) {
        $_ =~ s/^\s+|\s+$//g;
      }
      $n1 += 0;
      $n2 += 0;
      $lambda = 1e5/$dE if ( ((length($dE) > length($lambda)) && ($dE > 0)) || ($lambda =~ /[*]/) );

      #if ( ($J eq '4.0') && ($E eq '49.4372') && ($JP eq '5.0') && ($EP eq '69.3490') ) {
      #  $J = $J;
      #}

      my ($n_lev1, $E1, $lead_c_no, $lead_term, $n_lev2, $E21, $lead_c_no2, $lead_term2);
      if (!$n1) {
        # Identify the initial level
        ($n_lev1, $E1, $lead_c_no, $lead_term) = &FindLev(1, $J, $E, $nc, $term);
        if ( !$lead_c_no ) {
          die "Initial level not identified: J= $J, E= $E, level $nc $term";
        }
        if ( ($lead_c_no != $nc) || ($lead_term ne $term)  ) {
          #die "First level parameters mismatch: J = $J, conf. set $first_par; in transitions section, E = $E, conf_no = $nc, term = $term\n" .
        #  "In eigenvectors section: E = $E1, conf_no = $lead_c_no, term = $lead_term";
        }
        # Identify the final level
        ($n_lev2, $E21, $lead_c_no2, $lead_term2) = &FindLev(1, $JP, $EP, $ncP, $termP);
        if ( !$lead_c_no2 ) {
          die "Final level not identified: J= $JP, E= $EP, level $ncP $termP";
        }
        if ( ($lead_c_no2 != $ncP) || ($lead_term2 ne $termP)  ) {
          #die "\nFinal level parameters mismatch: J = $J2, conf. set $second_par; in transitions section, E = $E2, conf_no = $nc2, term = $term2;\n" .
          #  "In eigenvectors section: E = $E21, conf_no = $lead_c_no2, term = $lead_term2;";
        }
      } else {
        $n_lev1 = $n1;
        $n_lev2 = $n2;
        $E1 = $energies[0]->{$J}->{$n_lev1}->[0];
        $E21 = $energies[0]->{$JP}->{$n_lev2}->[0];
        my $E1_0 = $E1;
        if (abs($E1 - $E) > 0.03) {
          $n_lev1 = $n2;
          $n_lev2 = $n1;
          $E1 = $energies[0]->{$J}->{$n_lev1}->[0];
          $E21 = $energies[0]->{$JP}->{$n_lev2}->[0];
          if (abs($E1 - $E) > 0.03) {
            print "Error: too large energy differnce between identified lower RCG levels, J=$J: expected $E, found $E1 or $E1_0";
            exit(1);
          }
          if (abs($E21 - $EP) > 0.03) {
            print "Error: too large energy differnce between identified upper RCG levels, J=$JP: expected $EP, found $E21";
            exit(1);
          }
        }
        $lead_c_no = $nc;
        $lead_term = $term;
      }


      my ($par1, $par2, $E_init, $E_fin, $J_init, $J_fin, $n_init, $n_fin) = (1, 1, $EP, $E, $JP, $J, $n_lev2, $n_lev1);
      # Substitute experimental energies if known, and correct the gA value
      my $num_RCE_E1 = $map_RCG_RCE[$par1-1]->{$J_init}->{$n_init};
      my $RCE_data = $RCE_lev[$par1-1]->{$J_init}->[$num_RCE_E1-1];
      my $Ec1 = $RCE_data->{'Ec'};
      my $Ee1 = $RCE_data->{'Ee'};
      my $ec1 = $RCE_data->{'exp_c'};  # Star means "no experimental value"

      my $num_RCE_E2 = $map_RCG_RCE[$par2-1]->{$J_fin}->{$n_fin};
      $RCE_data = $RCE_lev[$par2-1]->{$J_fin}->[$num_RCE_E2-1];
      my $Ec2 = $RCE_data->{'Ec'};
      my $Ee2 = $RCE_data->{'Ee'};
      my $ec2 = $RCE_data->{'exp_c'};  # Star means "no experimental value"

      if ( ($n_init == $n_fin) && ($J_init == $J_fin) ) {
        die "Same level $n_init, J=$J_init, E=$E_init, type=$trans_type appears as lower and upper, wl=$lambda, dE=$dE";
      }

      if ( ((($ec1 eq '*') || ($ec2 eq '*')) && ($Ec2 > $Ec1)) || (($ec1 ne '*') && ($ec2 ne '*') && ($Ee2 > $Ee1)) ) {
        # Swap the levels, so that E_init is always the greatest
        my @tmp = ($n_init, $par1, $J_init, $E_init,$Ec1,$Ee1,$ec1);
        ($n_init, $par1, $J_init, $E_init, $Ec1, $Ee1, $ec1) = ($n_fin, $par2, $J_fin, $E_fin, $Ec2, $Ee2, $ec2);
        ($n_fin, $par2, $J_fin, $E_fin, $Ec2, $Ee2, $ec2) = @tmp;
      }
      $E_init = $Ec1;
      $E_fin = $Ec2;
      my $num_RCE_E_init = $map_RCG_RCE[$par1-1]->{$J_init}->{$n_init};
      my $num_RCE_E_fin = $map_RCG_RCE[$par1-1]->{$J_fin}->{$n_fin};

      $map_shells[$par1-1]->{$J_init} = {} unless defined $map_shells[$par1-1]->{$J_init};
      $map_shells[$par2-1]->{$J_fin} = {} unless defined $map_shells[$par2-1]->{$J_fin};
	  if (!defined($map_shells[$par1-1]->{$J_init}->{$num_RCE_E_init})) {
        my $shells = &get_leading_LS_term($par1,$J_init,$n_init); # Take the leading term as the level designation
        $map_shells[$par1-1]->{$J_init}->{$num_RCE_E_init} = $shells;
	  }
	  if (!defined($map_shells[$par2-1]->{$J_fin}->{$num_RCE_E_fin})) {
        my $shells2 = &get_leading_LS_term($par2,$J_fin,$n_fin); # Take the leading term as the level designation
        $map_shells[$par2-1]->{$J_fin}->{$num_RCE_E_fin} = $shells2;
      }
      my $dEc = abs($Ec1 - $Ec2);
      my $dEe = ( (($ec1 eq '*') || ($ec2 eq '*')) ? $dEc : abs($Ee1 - $Ee2));
      my $lambdaE = 1e5/$dEe;

      # Store the transition data
      $transitions[$par1-1]->{$J_init} = {} unless defined $transitions[$par1-1]->{$J_init};
      $trans_keys[$par1-1]->{$J_init} = {} unless defined $trans_keys[$par1-1]->{$J_init};
      $transitions[$par1-1]->{$J_init}->{$num_RCE_E_init} = [] unless defined $transitions[$par1-1]->{$J_init}->{$num_RCE_E_init};
      $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init} = {} unless defined $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init};

      my $S = &S_from_gA($gA,$lambda,$trans_type);
      my $gA_corr = ( (($ec1 eq '*') || ($ec2 eq '*')) ? $gA : &gA_from_S($S,$lambdaE,$trans_type) );
      my $A = $gA_corr / (2*$J_init + 1);
      my $logA = log($A);
      my $logA2 = $logA*$logA;

      if ( !$num_trial ) {
        push(@{$transitions[$par1-1]->{$J_init}->{$num_RCE_E_init}},
          [$par2, $J_fin, $n_fin, $dEc, 1e5/$dEc, $gA_corr, $S, $cf, $trans_type,
            [1,$logA,$logA2,$cf,$cf*$cf],[],[]]);
        my $tr_index = $#{$transitions[$par1-1]->{$J_init}->{$num_RCE_E_init}};
        my $key = "$par2$trans_type$J_fin$num_RCE_E_fin";
        $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init}->{$key} = $tr_index;

        $sum_A_hash{$J_init} = {[]} unless defined($sum_A_hash{$J_init});
        # Accumulate the total radiative decay rate for the upper level
        $sum_A_hash{$J_init}->{$num_RCE_E_init}->[0] += $gA_corr/(2*$J_init +1);
      } else {
        # Locate the transition $n_init -> $n_fin and accumulate statistics for it
        my $key = "$par2$trans_type$J_fin$num_RCE_E_fin";
        my $tr_index = undef;
        $tr_index = $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init}->{$key} if defined $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init};
        $sum_A_hash{$J_init}->{$num_RCE_E_init}->[$num_trial] += $gA_corr/(2*$J_init +1);
        if ( defined $tr_index ) {
          my $trans = $transitions[$par1-1]->{$J_init}->{$num_RCE_E_init}->[$tr_index];
          #if ( ref($trans->[9]) ne 'ARRAY' ) {
          #  $key = $key;
          #}
          my ($n, $sum_logA,$sum_logA2,$sum_cf,$sum_cf2) = @{$trans->[9]};
          $n++;
          $sum_logA += $logA;
          $sum_logA2 += $logA2;
          $sum_cf += $cf;
          $sum_cf2 += $cf*$cf;
          $trans->[9] = [$n, $sum_logA,$sum_logA2,$sum_cf,$sum_cf2];
          if ($num_trial <= $store_trials) {
            push(@{$trans->[10]},$A);
            push(@{$trans->[11]},$cf);
          }
        } else {
          # Add this transitin to the list and initialize its statistics
          $transitions[$par1-1]->{$J_init}->{$num_RCE_E_init} = [] unless defined $transitions[$par1-1]->{$J_init}->{$num_RCE_E_init};
          $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init} = {} unless defined $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init};
          push(@{$transitions[$par1-1]->{$J_init}->{$num_RCE_E_init}}, [$par2, $J_fin, $num_RCE_E_fin, $dE, $lambda, '', '', '', $trans_type,[1,$logA,$logA2,$cf,$cf*$cf],[$A],[$cf]]);
          my $tr_index = $#{$transitions[$par1-1]->{$J_init}->{$num_RCE_E_init}};
          my $key = "$par2$trans_type$J_fin$num_RCE_E_fin";
          $trans_keys[$par1-1]->{$J_init}->{$num_RCE_E_init}->{$key} = $tr_index;
        }
      }
    }
  }
  print "Done.\n";
} ##ReadTransitions_M1E2()

############################################################################
sub spectrum_header($) {    #12/10/2013 1:42PM
############################################################################
  my $s = shift;
  my $spectrum_data = undef;
  if ( ($s =~ /^. +([^ ]+) +([^ ]+) +SPECTRUM +[()]ENERGIES IN UNITS OF/) ) {
    my $spectrum_type = "$1 $2";
    my $trans_type='';
    if ( $spectrum_type =~ /MAG DIP/i ) {
      $trans_type = 'M1';
    } elsif ( $spectrum_type =~ /ELEC QUD/i ) {
      $trans_type = 'E2';
    }
    $spectrum_data = [$trans_type,$spectrum_type];
  }
  return $spectrum_data;
} ##spectrum_header($)

############################################################################
sub vary_params() {   #12/10/2013 4:16PM
############################################################################
  my $param_count = $#params;
  unless ( $test_random ) {
    # Create a dummy RCEOUT file in the current trial directory
    # and call s11.bat to transfer the varied parameters from it to ING11
    open(RCEOUT,">RCEOUT") or die "Could not create RCEOUT";
    print RCEOUT " PARAMETER FLAG      VALUE     MAX.VALUE      DENOM\n";
  }
  my %groups = ();
  for ( my $i = 0; $i <= $param_count; $i++ ) {
    my $k = -100;
    my ($par_name, $flag, $value, $sdx) = @{$params[$i]};
    my $a = $value;
    if ( $vary_fitted_params && ($sdx != 0) ) {
      if ( $var_grouped_params && (abs($flag) != 100) && ($value != 0) ) {
        if ( !defined($groups{$flag}) ) {
          # Vary the group ratio and store the coefficient in the $groups hash
          $groups{$flag} = &get_normal_random_trial(1,abs($sdx/$value));;
          if ( $print_distribution ) {
            $k = ($groups{$flag} - 1)/abs($sdx/$value)*10;
            push(@param_stats,$k/10);
            $k = &get_rounded_key($k);
            $dev_hash_param{$k}++;
            $dev_hash_param_groups{$flag} = {} unless defined $dev_hash_param_groups{$flag};
            $dev_hash_param_groups{$flag}->{$k}++;
          }
        }
        my $ratio = $groups{$flag};
        $a *= $ratio;
      } else {
        ($a) = &get_normal_random_trial($value,abs($sdx));
        if ( $print_distribution ) {
          $k = ($a - $value)/$sdx*10;
          push(@param_stats,$k/10);
          $k = &get_rounded_key($k);
          $dev_hash_param{$k}++;
        }
      }
    } elsif ( ($value != 0) && abs($flag == 100) && ($fixed_param_variance != 0) ) {
      ($a) =  &get_normal_random_trial($value,abs($value * $fixed_param_variance));
      if ( $print_distribution && ($value != 0) ) {
        $k = ($a/$value - 1)/$fixed_param_variance*10;
        push(@param_stats,$k/10);
        $k = &get_rounded_key($k);
        $dev_hash_param{$k}++;
        $dev_hash_param_fixed{$k}++;
      }
    }
    #if ( $k == 0 ) {
    #  $k = $k;
    #}
    # use only one of the returned random values
    print RCEOUT sprintf("%-10s%5d%14.6f%14.6f                    - \n", $par_name, $flag, $a, 0) unless $test_random;
  }
  unless ( $test_random ) {
    close RCEOUT or die "Error writing RCEOUT";
    `s11`;
    `del *.bak`;
  }
} ##vary_params()

############################################################################
sub PrintTransitions($;$) {    #12/10/2013 4:56PM
############################################################################
  my ($type,$num_type) = @_; # $type = M1 or E2, or empty if a total list is needed
                             # $num_type == 0 signals to create a new output file;
                             # $num_type != 0 signals to append an existing output file
  # Create sorted hash of lines
  my %lines = ();
  print "Sorting $type lines...\n";
  my $par = 1;
  foreach my $J (sort {$a<=>$b} keys %{$RCE_lev[$par-1]}) {
    for ( my $num_RCE_E1 = 1; $num_RCE_E1 <= $#{$RCE_lev[$par-1]->{$J}}+1; $num_RCE_E1++ ) {
      #my $num_RCE_E1 = $map_RCG_RCE[$par-1]->{$J}->{$n_lev1};
      next unless defined $transitions[$par-1]->{$J}->{$num_RCE_E1};
      #my ($E) = @{$energies[$par-1]->{$J}->{$n_lev1}};
      my $AR_tot = $sum_A_hash{$J}->{$num_RCE_E1}->[0];
      my @trans = @{$transitions[$par-1]->{$J}->{$num_RCE_E1}};
      my $num_trans = $#trans;

      my $RCE_data = $RCE_lev[$par-1]->{$J}->[$num_RCE_E1-1];
      my $Ee1 = $RCE_data->{'Ee'};
      my $Ec1 = $RCE_data->{'Ec'};
      my $ec1 = $RCE_data->{'exp_c'};  # Star means "no experimental value"

      for ( my $i = 0; $i <= $num_trans; $i++) {
        my ($second_par, $J2, $num_RCE_E2, $dE, $lambda, $gA, $S, $cf, $trans_type, $stats_accum, $trial_A_data, $trial_cf_data) = @{$trans[$i]};
        my $E2_frac = (($trans_type eq 'M1') ? 0 : 1);
        if ( ($J eq '5.0') && ($num_RCE_E1 == 10) && ($J2 eq '4.0') && ($num_RCE_E2 == 12) ) {
          $J = $J;
        }
        next if ( $type && ($type ne $trans_type) );
        my ($n_trials,$sum_logA,$sum_logA2,$sum_cf,$sum_cf2) = @{$stats_accum};
        # Mean ln(A)
        my $A_mean  = $sum_logA/$n_trials;
        # Stdev of ln(A)
        my $stdev_A = $sum_logA2/$n_trials - $A_mean*$A_mean;
        $stdev_A = (($stdev_A > 0 ) ? sqrt($stdev_A) : 0);
        # Mean A
        $A_mean = exp($A_mean);
        # Relative stdev of S in %
        $stdev_A = (exp($stdev_A) - 1) * 100;
        # Mean cf
        my $cf_mean = $sum_cf/$n_trials;
        # stdev of cf
        my $stdev_cf = $sum_cf2/$n_trials - $cf_mean*$cf_mean;
        $stdev_cf = (($stdev_cf > 0 ) ? sqrt($stdev_cf) : 0);

        #my $num_RCE_E2 = $map_RCG_RCE[$second_par-1]->{$J2}->{$n_lev2};
        $RCE_data = $RCE_lev[$second_par-1]->{$J2}->[$num_RCE_E2-1];
        my $Ee2 = $RCE_data->{'Ee'};
        my $Ec2 = $RCE_data->{'Ec'};
        my $ec2 = $RCE_data->{'exp_c'};  # Star means "no experimental value"

        my $dEe = (($ec1 || $ec2) ? abs($Ec1-$Ec2) : abs($Ee1 - $Ee2));
        my $lambdaE = 1e5/$dEe;

        # Mean S and gA from A_mean
        my $gA_mean = $A_mean * (2*$J + 1);
        my $S_mean = &S_from_gA($gA_mean,$lambdaE,$trans_type);

        my $key = sprintf("%12.3f", $lambdaE);
        my $E2_frac = (($trans_type eq 'M1') ? 0 : 1);
        my $stdev_frac = 0;
        my $continue = 1;
        while ( defined($lines{$key}) && $continue ) {
          my ($px, $Jx, $n_lev1x, $p2x, $J2x, $n_lev2x, $gAx, $Sx, $cfx,$gA_meanx,$S_meanx,$stdev_Ax,$cf_meanx,$stdev_cfx, $t_x, $E2_frac_x, $stdev_frac_x, $AR_tot_x, $trial_A_data_x, $trial_cf_data_x) = @{$lines{$key}};
          if ( ("$J $num_RCE_E1" eq "$Jx $n_lev1x") && ("$J2 $num_RCE_E2" eq "$J2x $n_lev2x") ) {
            if ( ($t_x eq $trans_type) || ($t_x eq 'M1+E2') ) {
              die "Duplicate transitions: J=$J, n_lev1=$num_RCE_E1, J2=$J2, n_lev2=$num_RCE_E2, lambda=$lambda, type=$trans_type";
            } elsif ($type eq '') {
              # Merge the two?
              # 'x' is always 'M1'
              if ( $t_x ne 'M1' ) {
                my @tmp = ($par, $J, $num_RCE_E1, $second_par, $J2, $num_RCE_E2, $gA, $S, $cf,$gA_mean,$S_mean,$stdev_A,$cf_mean,$stdev_cf, $trans_type,$E2_frac,$stdev_frac, $AR_tot,$trial_A_data,$trial_cf_data);
                ($par, $J, $num_RCE_E1, $second_par, $J2, $num_RCE_E2, $gA, $S, $cf,$gA_mean,$S_mean,$stdev_A,$cf_mean,$stdev_cf, $trans_type,$E2_frac,$stdev_frac, $AR_tot,$trial_A_data,$trial_cf_data) =
                  ($px, $Jx, $n_lev1x, $p2x, $J2x, $n_lev2x, $gAx, $Sx, $cfx,$gA_meanx,$S_meanx,$stdev_Ax,$cf_meanx,$stdev_cfx, $t_x, $E2_frac_x, $stdev_frac_x, $AR_tot_x, $trial_A_data_x,$trial_cf_data_x);
                ($px, $Jx, $n_lev1x, $p2x, $J2x, $n_lev2x, $gAx, $Sx, $cfx,$gA_meanx,$S_meanx,$stdev_Ax,$cf_meanx,$stdev_cfx, $t_x, $E2_frac_x, $stdev_frac_x, $AR_tot_x, $trial_A_data_x, $trial_cf_data_x) =
                  @tmp;
              }
              my ($gA1,$gAx1) = (($gA ne '') ? ($gA,$gAx) : ($gA_mean,$gA_meanx));
              my $gAtot = $gA1 + $gAx1;
              $E2_frac = $gA1/$gAtot;
              my $stdev_logA = log($stdev_A/100+1);
              my $stdev_logAx = log($stdev_Ax/100+1);
              my $r = $gAx1/$gA1;
              $stdev_logA = sqrt($stdev_logA*$stdev_logA + $r*$r*$stdev_logAx*$stdev_logAx)/(1 + $r);
              $stdev_A = (exp($stdev_logA) - 1)*100;
              $gA_mean += $gA_meanx;
              $cf = (abs($cf)*$gA1+abs($cfx)*$gAx1)/$gAtot if $cf ne '';
              $stdev_cf = $stdev_cf*$gA1/$gAtot;
              $stdev_cfx = $stdev_cfx*$gAx1/$gAtot;
              $stdev_cf = sqrt($stdev_cf*$stdev_cf + $stdev_cfx*$stdev_cfx);
              $cf_mean = (abs($cf_mean)*$gA1+abs($cf_meanx)*$gAx1)/$gAtot;
              $gA = $gAtot if $gA ne '';
              my $n_trials = $#{$trial_A_data};
              my $n_trials_x = $#{$trial_A_data_x};
              #$n_trials_x = 0 if $n_trials <= 0;
              my $n_t = $n_trials;
              $n_t = $n_trials_x if ( $n_trials_x > $n_trials );
              my @A_data = ();
              my @cf_data = ();
              my $n1 = 0;
              my $sum_frac = 0;
              my $sum_frac2 = 0;
              for ( my $i = 0; $i <= $n_t; $i++ ) {
                my ($A_this, $cf_this, $Ax, $cfx) = (0,0,0,0);
                if ( $i <= $n_trials ) {
                  $A_this = $trial_A_data->[$i];
                  $cf_this = abs($trial_cf_data->[$i]);
                }
                if ($i <= $n_trials_x) {
                  $Ax = $trial_A_data_x->[$i];
                  $cfx = abs($trial_cf_data_x->[$i]);
                  if ($i <= $n_trials) {
                    my $E2_frac_trial = $A_this/($Ax + $A_this);
                    $n1++;
                    $sum_frac += $E2_frac_trial;
                    $sum_frac2 += $E2_frac_trial*$E2_frac_trial;
                  }
                }
                my $A_tot = $A_this + $Ax;
                if ( $A_tot == 0) {
                  $A_tot = $A_tot;
                }
                $cf_mean = ($cf_this*$A_this + $cfx*$Ax)/$A_tot;
                #if ( $cf_mean > 1 ) {
                #  $cf_mean = $cf_mean;
                #}
                push(@A_data, $A_tot);
                push(@cf_data,$cf_mean);
              }
              my $mean_frac = ($n1 <= 0 ? $E2_frac : $sum_frac/$n1);
              $stdev_frac = ($n1 <= 0 ? 0 : $sum_frac2 / $n1 - $mean_frac*$mean_frac);
              $stdev_frac = 0 if $stdev_frac < 0;
              my $d_frac = $E2_frac - $mean_frac;
              $stdev_frac = sqrt($stdev_frac + $d_frac*$d_frac);
              if ( $E2_frac + 2*$stdev_frac < $min_mixed_fraction ) {
                # Ignore the current transition and leave only the one already stored in lines hash ('x' = 'M1')
                my $E2_frac_x = $E2_frac;
                ($par, $J, $num_RCE_E1, $second_par, $J2, $num_RCE_E2, $gA, $S, $cf,$gA_mean,$S_mean,$stdev_A,$cf_mean,$stdev_cf, $trans_type,$E2_frac,$AR_tot,$trial_A_data,$trial_cf_data) =
                  ($px, $Jx, $n_lev1x, $p2x, $J2x, $n_lev2x, $gAx, $Sx, $cfx,$gA_meanx,$S_meanx,$stdev_Ax,$cf_meanx,$stdev_cfx, $t_x, $E2_frac_x, $AR_tot_x, $trial_A_data_x, $trial_cf_data_x);
                #next;
              } elsif ( 1 - $E2_frac + 2*$stdev_frac < $min_mixed_fraction ) {
                # Leave the key as is; This will
                # replace the stored transition in line hash with the current new transition
                #$key = $key;
              } else {
                # Merge the two types in one transition
                $S = '';
                $S_mean = '';
                $trans_type = 'M1+E2';
                $trial_A_data = \@A_data;
                $trial_cf_data = \@cf_data;
              }
              $continue = 0;
            }
          } else {
            $continue = 1;
            $key .= '1';
          }
        }
        $lines{$key} = [$par, $J, $num_RCE_E1, $second_par, $J2, $num_RCE_E2, $gA, $S, $cf,$gA_mean,$S_mean,$stdev_A,$cf_mean,$stdev_cf, $trans_type,$E2_frac,$stdev_frac, $AR_tot,$trial_A_data,$trial_cf_data];
      }
    }
  }

  # Print lines
  print "Printing $type lines...\n";
  # Print header
  my $delim = "\t";
  unless ( $num_type ) {
    # Print file header
    print OUT_TR_FILE join($delim,'type','conf1','t1','J1','conf2','t2','J2','lande1','lande2','E1','E2','wl_c(A)','E1_exp','E2_exp','wl_exp(A)','A(s-1)', 'S','cf','sum A','BF');
    if ( $mearged_list ) {
      print OUT_TR_FILE join($delim, '','E2_frac');
      if ( $num_trials ) {
        print OUT_TR_FILE join($delim, '','d_E2frac');
      }
    }
    if ( $num_trials ) {
      print OUT_TR_FILE join($delim, '','d_BF','A_mean','S_mean','std_A%','cf_mean','std_cf');
    }
    if ( $store_trials ) {
      my @a = ();
      for ( my $i = 1; $i < $store_trials; $i++ ) {
        push(@a, $i+1);
      }
      print OUT_TR_FILE join($delim,'','Trial A data',@a);
      print OUT_TR_FILE join($delim,'','Trial cf data',@a) if $print_cf_trials;
    }
    print OUT_TR_FILE "\n";
  }

  # Print lines data
  foreach my $wl (sort {$a<=>$b} keys %lines) {
    my ($par, $J, $num_RCE_E1, $second_par, $J2, $num_RCE_E2, $gA, $S, $cf,$gA_mean,$S_mean,$stdev_A,$cf_mean,$stdev_cf, $trans_type,$E2_frac,$stdev_frac, $AR_tot,$trial_A_data,$trial_cf_data)
      = @{$lines{$wl}};
    my $RCE_data = $RCE_lev[$par-1]->{$J}->[$num_RCE_E1-1];
    my $ec1 = $RCE_data->{'exp_c'};  # Star means "no experimental value"
    my $Ee1 = ($ec1) ? '' : $RCE_data->{'Ee'};
    my $Ec1 = $RCE_data->{'Ec'};
    my $lande1 = $RCE_data->{'lande'};

    my $AR_tot = $sum_A_hash{$J}->{$num_RCE_E1}->[0];
    if ( $AR_tot == 0 ) {
      my $knz = 0;
      for (my $i = 0; $i < $num_trials; $i++) {
        if ( $i <= $#{$trial_A_data} ) {
          my $AR_tot_trial = $sum_A_hash{$J}->{$num_RCE_E1}->[$i];
          if ( $AR_tot_trial > 0 ) {
            $AR_tot += $AR_tot_trial;
            $knz++;
          }
        }
      }
      $AR_tot /= $knz if ($knz > 0);
    }
    my $BF = (($AR_tot != 0) ? ($gA > 0 ? $gA : $gA_mean)/(2*$J+1)/$AR_tot : 1);
    next if $BF < $min_BF;
    my $shells = $map_shells[$par-1]->{$J}->{$num_RCE_E1};

    # Substitute Eexp from RCE if available
    $RCE_data = $RCE_lev[$second_par-1]->{$J2}->[$num_RCE_E2-1];
    my $ec2 = $RCE_data->{'exp_c'};  # Star in $ec2 means "no experimental value"
    my $Ee2 = ($ec2) ? '' : $RCE_data->{'Ee'};
    my $Ec2 = $RCE_data->{'Ec'};
    my $lande2 = $RCE_data->{'lande'};

    my $shells2 = $map_shells[$par-1]->{$J2}->{$num_RCE_E2}; # Take the leading term as the level designation

    my $dE = abs($Ec1 - $Ec2);
    my ($wl_out,$wl_out_e) = (1e5/$dE, '');

    if ( !$ec1 && !$ec2) {
      $dE = abs($Ee1 - $Ee2);
    }
    if ( $dE == 0 ) {
      next;
    }

    if ( ($dE > 5) && ($dE < 50) ) {
      $wl_out = sprintf("%12.4f",&Lair($wl_out));
      if ( !$ec1 && !$ec2) {
        $wl_out_e = sprintf("%12.4f",&Lair(1e5/$dE));
      }
    } else {
      $wl_out = sprintf("%12.4f",$wl_out);
      if ( !$ec1 && !$ec2) {
        $wl_out_e = sprintf("%12.4f",1e5/$dE);
      }
    }
    $wl_out =~ s/^\s+//g;
    $wl_out_e =~ s/^\s+//g;
    my $AR_out = sprintf("%12.5e",$AR_tot);
    my $A_out = (($gA ne '') ? sprintf("%12.5e",$gA/(2*$J+1)) : '');
    my $A_mean = $gA_mean/(2*$J+1);
    my $A_mean_out = sprintf("%12.5e",$A_mean);
    my ($S_out,$S_mean_out) = ('','');
    if ( $trans_type ne 'M1+M2' ) {
      $S_out = sprintf("%11.4e",$S) if ($S ne '');
      $S_mean_out = sprintf("%11.4e",$S_mean) if ($S_mean ne '');
    }
	if ($shells2 !~ /$delim/) {
	  $shells2 .= "$delim";
	}
	if ($shells !~ /$delim/) {
	  $shells .= "$delim";
	}
    print OUT_TR_FILE join($delim,$trans_type,$shells2,$J2,$shells,$J,$lande2,$lande1,$Ec2,$Ec1,$wl_out,$Ee2,$Ee1,$wl_out_e,$A_out,$S_out, $cf,$AR_out,$BF);
    if ( $mearged_list ) {
      print OUT_TR_FILE join($delim, '',$E2_frac);
    }
    my $dBF = '';
    if ( $num_trials ) {
      if ( ($gA ne '') || ($gA_mean > 0) ) {
        my $stdev_logA = log($stdev_A/100 + 1);
        my $d_logA = log($gA_mean/($gA || $gA_mean));
        $stdev_logA = sqrt($stdev_logA*$stdev_logA + $d_logA*$d_logA);
        $stdev_A = (exp($stdev_logA) - 1)*100;
        my $cf_init = ($gA ne '') ? $cf : $cf_mean;
        my $d_cf = $cf_init - $cf_mean;
        $stdev_cf = sqrt($stdev_cf*$stdev_cf + $d_cf*$d_cf);
        my $A_init = ($gA || $gA_mean) / (2*$J + 1);
        my @BF_trials = ();
        if ( $print_dBF ) {
          for (my $i = 0; $i < $store_trials; $i++) {
            my $A = '';
            if ( $i <= $#{$trial_A_data} ) {
              $A = $trial_A_data->[$i];
              my $Atot = $sum_A_hash{$J}->{$num_RCE_E1}->[$i+1];
              if ($Atot != 0) {
                push(@BF_trials,$A/$Atot);
              }
            }
          }
          $dBF = &stdev($BF,@BF_trials);
        }
        if ( $print_distribution && $store_trials && ($gA ne '')) {
          my $Aw = 0;
          for (my $i = 0; $i < $store_trials; $i++) {
            my ($A,$cf1) = ('','');
            if ( $i <= $#{$trial_A_data} ) {
              $A = $trial_A_data->[$i];
              $cf1 = $trial_cf_data->[$i];
              my $dev_A = log($A/$A_init)/$stdev_logA; # Relative deviation in terms of standard deviation
              my $k = $dev_A * 10;
              $k = &get_rounded_key($k);
              $dev_hash_logA{$k}++;

              if ( $stdev_cf > 0 ) {
                $k = ($cf1 - $cf_init)/$stdev_cf * 10;
                $k = &get_rounded_key($k);
                $dev_hash_cf{$k}++;
              }

              my $dA = $A - $A_init;
              $Aw += $dA*$dA;
            }
          }
          $Aw = sqrt($Aw/($#{$trial_A_data} + 1));
          for (my $i = 0; $i < $store_trials; $i++) {
            if ( $i <= $#{$trial_A_data} ) {
              my $A = $trial_A_data->[$i];
              my $dev_A = ($A - $A_init)/$Aw; # Relative deviation in terms of standard deviation
              my $k = $dev_A * 10;
              $k = &get_rounded_key($k);
              $dev_hash_A{$k}++;
            }
          }
        }
      }
      my $stdev_A_out = sprintf("%10.2f",$stdev_A);
      my $cf_mean_out = sprintf("%10.7f",$cf_mean);
      my $stdev_cf_out = sprintf("%10.7f",$stdev_cf);
      if ( $mearged_list ) {
        print OUT_TR_FILE join($delim, '',$stdev_frac);
      }
      print OUT_TR_FILE join($delim, '',$dBF,$A_mean_out,$S_mean_out,$stdev_A_out,$cf_mean_out,$stdev_cf_out);
    }
    if ( $store_trials ) {
      my @dummy = ();
      for (my $i = $#{$trial_A_data} + 1; $i < $store_trials; $i++) {
        push(@dummy, '');
      }
      foreach (@{$trial_A_data}) {
        $_ = sprintf("%12.5e", $_);
      }
      foreach (@{$trial_cf_data}) {
        $_ = sprintf("%10.7e", $_);
      }
      if ( $sort_trials ) {
        print OUT_TR_FILE join($delim, '',(sort {$a<=>$b} @{$trial_A_data}),@dummy);
        print OUT_TR_FILE join($delim, '',(sort {$a<=>$b} @{$trial_cf_data}), @dummy) if $print_cf_trials;
      } else {
        print OUT_TR_FILE join($delim, '',(@{$trial_A_data}),@dummy);
        print OUT_TR_FILE join($delim, '',(@{$trial_cf_data}), @dummy) if $print_cf_trials;
      }
    }
    print OUT_TR_FILE "\n";
  }
} ##PrintTransitions()


############################################################################
sub print_distribution() {    #12/18/2013 11:04AM
############################################################################
  print "\nDistribution of logarithmic deviations:\n";
  print "d/stdev\tN\n";
  foreach my $k (sort {$a<=>$b} keys %dev_hash_logA) {
    print join("\t", $k*0.1, $dev_hash_logA{$k}), "\n";
  }

  print "\nDistribution of straight deviations:\n";
  print "d/stdev\tN\n";
  foreach my $k (sort {$a<=>$b} keys %dev_hash_A) {
    print join("\t", $k*0.1, $dev_hash_A{$k}), "\n";
  }

  print "\nDistribution of CF deviations:\n";
  print "d/stdev\tN\n";
  foreach my $k (sort {$a<=>$b} keys %dev_hash_cf) {
    print join("\t", $k*0.1, $dev_hash_cf{$k}), "\n";
  }

  print "\nDistribution of input parameter deviations:\n";
  my $kurt = &Kurtosis(0,@param_stats);
  print "Kurtosis = $kurt\n";
  print "d/stdev\tN\n";
  foreach my $k (sort {$a<=>$b} keys %dev_hash_param) {
    print join("\t", $k*0.1, $dev_hash_param{$k}), "\n";
  }

  return;

  if ( $var_grouped_params ) {
    print "\nDistribution of input parameter deviations by group:\n";
    foreach my $group (sort {abs($a)<=>abs($b)} keys %dev_hash_param_groups) {
      print "\nGroup $group:\n";
      print "d/stdev\tN\n";
      my %hash = %{$dev_hash_param_groups{$group}};
      foreach my $k (sort {$a<=>$b} keys %hash) {
        print join("\t", $k*0.1, $hash{$k}), "\n";
      }
    }
  }

  print "\nDistribution of fixed-input-parameter deviations:\n";
  print "d/stdev\tN\n";
  foreach my $k (sort {$a<=>$b} keys %dev_hash_param_fixed) {
    print join("\t", $k*0.1, $dev_hash_param_fixed{$k}), "\n";
  }

} ##print_distribution()

############################################################################
sub get_rounded_key($) {    #12/18/2013 5:03PM
############################################################################
  my $k = shift;
  my $sign = (($k < 0) ? -1 : 1);
  $k = abs($k) + 0.5;
  { use integer;
    $k += 0;
  }
  return $k*$sign;
} ##get_rounded_key($)

############################################################################
sub Kurtosis(@) {  #12/19/2013 4:02PM
############################################################################
  # As defined in Excel, see http://office.microsoft.com/en-us/excel-help/kurt-HP005209150.aspx
  my ($mean, @data) = @_;
  my $n = $#data + 1;
  return '' if ($n<4);
  # Compute variance === standard deviation squared
  my $std2 = 0;
  foreach (@data) {
    $std2 += $_*$_;
  }
  $std2 /= $n;
  $std2 -= $mean*$mean;
  $std2 = 0 if $std2 < 0;
  return '' if ($std2 == 0);
  # Compute Kurtoise
  my $K = 0;
  foreach (@data) {
    my $d = ($_ - $mean);
    $d = $d*$d/$std2;
    $K += $d*$d;
  }
  my $n_23 = ($n-2)*($n-3);
  my $n_1 = $n-1;
  $K = ($K * $n*($n+1)/$n_1 - 3*$n_1*$n_1)/$n_23;
  return $K;
} ##Kurtosis(@)

############################################################################
sub stdev(@) {    #12/19/2013 4:26PM
############################################################################
  my ($mean, @data) = @_;
  my $n = $#data + 1;
  return '' if ($n < 1);
  # Compute standard deviation
  my $std2 = 0;
  my $av = 0;
  foreach (@data) {
    $av += $_;
    $std2 += $_*$_;
  }
  $av /= $n;
  $std2 /= $n;
  $std2 -= $av*$av;
  $std2 = 0 if $std2 < 0;
  my $d_av = $mean - $av;
  return sqrt($std2 + $d_av*$d_av);
} ##stdev(@)

############################################################################
sub scale_A($) {    #12/19/2013 4:31PM
############################################################################
  my $A = shift;
  return exp(ln($A)/3); # Power of 1/3
} ##scale_A($)

############################################################################
sub scale_A_back($) {    #12/19/2013 4:31PM
############################################################################
  my $A = shift;
  return $A*$A*$A;
} ##scale_A_back($)

############################################################################
sub average(@) {    #12/20/2013 12:12PM
############################################################################
  my @arr = @_;
  my $n = 0;
  my $sum = 0;
  foreach my $a (@arr) {
    next unless defined $a;
    $n++;
    $sum += $a;
  }
  return ($n ? $sum/$n : undef);
} ##average(@)

############################################################################
sub skew(@) {   #12/20/2013 12:54PM
############################################################################
# As defined in Excel, see http://office.microsoft.com/en-us/windows-sharepoint-services-help/skew-function-HA001161067.aspx?CTT=1
  my ($mean, @data) = @_;
  my $n = $#data + 1;
  return '' if ($n<3);
  # Compute variance === standard deviation squared
  my $std2 = 0;
  foreach (@data) {
    $std2 += $_*$_;
  }
  $std2 /= $n;
  $std2 -= $mean*$mean;
  $std2 = 0 if $std2 < 0;
  return '' if ($std2 == 0);
  # Compute skewness
  my $std = sqrt($std2);
  my $K = 0;
  foreach (@data) {
    my $d = ($_ - $mean);
    $d = $d*$d*$d/($std2*$std);
    $K += $d;
  }
  $K *= $n/(($n-1)*$n-2);
  return $K;
} ##skew(@)
