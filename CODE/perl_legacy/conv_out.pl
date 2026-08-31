#!perl
use strict;
use vars qw{@parities @energies @map_RCG_RCE @RCE_lev @vectors @basis %L_moment %conf_nums @confs};
require 'conv_cowan.pl';
require 'vacair.pl';

# Determine if files for IDEN must be created.
my $IDEN = 0;
my $truncate = 0;
my $num_spliced = 0;
my $linear_intensity_scale = 0;
my $convert2log = 0;
do {
    $num_spliced = 0;
    for (my $i = 0; $i <= $#ARGV; $i++) {
      if ( $ARGV[$i] =~ /truncate/ ) {
        $truncate = 1;
        splice(@ARGV,$i,1);
        $num_spliced++
      } elsif ( $ARGV[$i] =~ /IDEN/ ) {
        $IDEN = 1;
        splice(@ARGV,$i,1);
        $num_spliced++
      } elsif ( $ARGV[$i] =~ /linear/ ) {
        $linear_intensity_scale = 1;
        splice(@ARGV,$i,1);
        $num_spliced++
      } elsif ( $ARGV[$i] =~ /convert2log/ ) {
        $convert2log = 1;
        splice(@ARGV,$i,1);
        $num_spliced++
      }
    }
} while $num_spliced;
if ($linear_intensity_scale && $convert2log) {
  print("Options 'linear' and 'convert2log' are incompatible. Only one of them can be specified.");
  exit;
}

# Read input parameters from remaining command-line arguments
my $upper_lev_bound=1e20;
my $out_lev_file = shift;
my $out_tr_file = shift;
if ($out_tr_file =~ /^u([0-9.e+-]+)$/) {
  $upper_lev_bound = $1+0;
  $out_tr_file = shift;
}
my $no_RCE = uc(shift);
my $sp_type = uc(shift);
#my $DR_channel = shift;
my $t_flight = 0.114; # Beam-foil spectroscopy: ions time of flight across the entrance slit of the spectrograph (ns)
if ( $sp_type =~ /^BF:([0-9.]+)$/ ) {
  $t_flight = $1;
}
my $T = shift;
my $scale_int = shift;
my $ai_set1_file = shift;
my $ai_set2_file = shift;

my $ev_cm = 8.06554393734921; # CODATA >= 2018 -- Conversion factor from eV to cm-1

#my $T = 19284.49; #Effective temperature for Boltzmann distribution in units of 1000 cm-1
#my $T = 5000.0;  # W42+
        # Ne VIII: 1671.60
        # Ne IX:   1928.447


if ( !$out_lev_file || !$out_tr_file || ($no_RCE !~ /RCE/) || ($sp_type !~ /LTE|LTEABS|BF:[0-9.]+|DR|AUGER|CM/) ||
  ($T !~ /^[0-9.Ee+-]+$/) || ($T == 0) || ($scale_int !~ /^[0-9.Ee+-]+$/) || ($scale_int < 1) || ($ai_set1_file && ! $ai_set2_file)
) {
  print "\nUsage:\n";
  print "conv_out.bat <out_lev_file_name> <out_trans_file_name> [no]RCE <LTE|LTEABS|BF:<time_of_flight_ns>|DR|AUGER|CM> " .
    "<eff_temperature (10^3 cm-1)> <max. scaled. intensity> [<AI_even_file_name> <AI_odd_file_name> [IDEN [truncate [linear] [convert2log]]]]\n";
  print "For Beam-Foil spectra, time of flight in nanoseconds is needed.\n" .
    "This is the time needed for ions in the beam to travel across the viewing aperture of the spectrometer.\n" .
    "CM is 'Cascade Matrix' option.\n" .
    "In all modes, effective temperature parameter is required, although it is used only in the LTE and LTEABS options.\n" .
    "RCE           - fetch experimental energies from RCEOUT and calculate Rits wavelengths with them ('noRCE' ignores RCE output)\n" .
    "LTE           - gives relative emission intensities modeled by Boltzmann populations of upper levels of transitions.\n" .
    "LTEABS        - gives relative reduced absorption intensity modeled by Boltzmann populations of lower levels of transitions.\n" .
    "IDEN          - if included, creates input files for the IDEN2 code, ENLEV.DAT and TRANS.DAT (provided that the files dlv.dat and numset.dat are present).\n" .
    "truncate      - if included, truncates the list of levels for IDEN2 by deleting the levels not connected to any others by calculated transitions.\n" .
    "linear        - if included, denotes that the calculated intensities for IDEN2 in TRANS.DAT must be on a linear scale (default is log scale).\n" .
    "convert2log   - if included, denotes that a new DLV.DAT must be created with intensities converted to a log scale, I_log = 10*ln(I).\n\n";
  exit;
}

$no_RCE = ($no_RCE =~ /^NO/) ? 1 : 0;

open OUTG11, "<OUTG11" or die "Could not open input file OUTG11";
open OUT_LEV_FILE, ">$out_lev_file" or die "Could not open output file " . $out_lev_file;

# variables for IDEN
my @lines = ();
my @numset = ();
my %IDEN_lev_map = ();
my %TRANSDAT_map = ();
my $Athresh = 0.1;
my $I_obs_max = 0;

if ( $IDEN ) {
  &ReadLines();
  &ReadNumset();
  &ConvertNumset();
}

# Global parameters
my $num_printed_components = 5;          # Number of printed eigenvector components
my $min_printed_percentage = 1.5;        # Omit 3rd, 4th, 5th, etc. eigenvector components with percentage less than that
my $min_printed_second_percentage = 1.5; # Omit second eigenvector components with percentage less than that
my $max_scaled_intensity = $scale_int;
my $s = '';
my $parity = 0;
my ($J, $J_prev) = (-1,1000);
my $start = 0;

my $delim = "\t"; # Delimitor to use in the output files
$delim = ',' if ($sp_type =~ /^CM$/i); # Use comma as delimitor in the 'CM' mode

# Initialize arrays
my @lev_nums = ({},{});
my @transitions = ({},{});
my %lev_num_hash = ();
my @AI_rates = ({},{});
my @AI_cores = ([],[]);
my $n_conf = '';
my %lines = ();

# # Initialize variables
my %E_map = ();
my %E_reverse_map = ();
my %E_map_ENLEV = ();
my %E_reverse_map_ENLEV = ();
my $reverse_order = 0;
my $num_levs = 0;

&init_vars();

# Read RCG options
&read_RCG_options();
&read_in36();

# Start processing ...

# Read configurations from OUTG11 for the first parity and second parity, if present
&read_confs();

&read_ING11_params();

# Read the basis state definitions printed by CALCFC in OUTG11
$s = &read_basis();


# Read Slater parameter labels and values of the first parity ------------
# Read eigenvalues and eigenvectors of the first parity
$s = &read_OUTG11_params(1,$s);

# Continue to read the OUTG11 file. Find and read the LS basis state labels printed by ENERGY
# If second parity present, read Slater parameter labels and values of the first parity,
# and read eigenvalues and eigenvectors of the second parity
&read_basis_labels($s,1);

&sort_vectors();

if ( !$no_RCE ) {
  &ReadRCE();
  &Identify_RCE_levs();
}

if ( $ai_set1_file && $ai_set2_file ) {
  &ReadAI($ai_set1_file, 1);
  &ReadAI($ai_set2_file, 2);
}

my $spectrum = &ReadTransitions();

# Print the energy levels in LS and JJ coupling
print OUT_LEV_FILE "Cpl${delim}par${delim}Lev_\#${delim}" .
  "Ee${delim}e_c${delim}Ec${delim}J"
  . (!$no_RCE ? "${delim}Lande_g" : '')
  . ($ai_set1_file && $ai_set2_file ? "${delim}Eai${delim}AA_tot" : '')
  . "${delim}AR_tot${delim}n${delim}l${delim}occ${delim}S${delim}L${delim}" .
  "\%1${delim}conf1${delim}term1${delim}\%2${delim}conf2${delim}term2\n";

foreach my $cpl  ('LS','JJ' ) {

  last if ( ($cpl ne 'LS') && ($sp_type =~ /^CM$/i) );

  print "Printing $cpl levels...";

  # Create a map of energy values for energy sorting
  for ( $parity = 1; $parity <=2; $parity++ ) {
    my $i = 0;
    foreach my $J ( keys %{$energies[$parity-1]} ) {
      foreach my $num_e ( keys %{$energies[$parity-1]->{$J}}) {
        $i++;
        my ($E,$Eai,$AA_tot, $AR_tot, $DR) = @{$energies[$parity-1]->{$J}->{$num_e}};
        my $key = "$parity${delim}$E${delim}$i";
        $E_map{$key} = [$J,$num_e,$Eai,$AA_tot,$AR_tot];
        $E_reverse_map{$parity} = {} unless defined $E_reverse_map{$parity};
        $E_reverse_map{$parity}->{$J} = {} unless defined $E_reverse_map{$parity}->{$J};
        $E_reverse_map{$parity}->{$J}->{$num_e} = [$key];
      }
    }
  }

  my $lev_num = 0;
  foreach my $key (sort {&sort_by_par_E($a,$b)} keys %E_map) {
    my ($parity,$E,$k) = split(/${delim}/, $key);
    next if ($E+0 > $upper_lev_bound);
    my $par_code = $parities[$parity-1];
    $lev_num++;
    my ($J,$num_e,$Eai,$AA_tot,$AR_tot) = @{$E_map{$key}};
    push(@{$E_map{$key}},$lev_num);
    push(@{$E_reverse_map{$parity}->{$J}->{$num_e}},$lev_num);
    $lev_nums[$parity-1] = {} unless defined $lev_nums[$parity-1];
    $lev_nums[$parity-1]->{$J} = {} unless defined $lev_nums[$parity-1]->{$J};
    $lev_nums[$parity-1]->{$J}->{$num_e} = $lev_num;
    $lev_num_hash{$lev_num} = [$parity,$J,$num_e];

    my ($num_RCE_E,$exp_c,$Ee, $Ec, $lande_g) = (0, '*',$E,$E,'');

    # Substitute Eexp and Ec from RCE if available
    if ( !$no_RCE ) {
      #$num_RCE_E = $map_RCG_RCE[$parity-1]->{$J}->{$num_e};
      #$exp_c = $RCE_lev[$parity-1]->{$J}->{$num_RCE_E}->{'exp_c'};
      #$Ec = $RCE_lev[$parity-1]->{$J}->{$num_RCE_E}->{'Ec'};
      #$Ee = ($exp_c) ? $Ec : $RCE_lev[$parity-1]->{$J}->{$num_RCE_E}->{'Ee'};
      #$lande_g = $RCE_lev[$parity-1]->{$J}->{$num_RCE_E}->{'lande'};
      $num_RCE_E = $map_RCG_RCE[$parity-1]->{$J}->{$num_e};
      my $hash = $RCE_lev[$parity-1]->{$J}->[$num_RCE_E-1];
      $exp_c = $hash->{'exp_c'};
      $Ec = $hash->{'Ec'};
      $Ee = ($exp_c) ? $hash->{'Ec'} : $hash->{'Ee'};
      $lande_g = $hash->{'lande'};
    }

    print OUT_LEV_FILE "$cpl${delim}$par_code${delim}$lev_num${delim}" .
#      sprintf("%12.4f${delim}%s${delim}%12.3f${delim}%4.1f${delim}%12.4f${delim}%9.2e${delim}%9.2e",$Ee,$exp_c,$E,$J,$Eai,$AA_tot,$AR_tot);
      sprintf("%12.6f${delim}%s${delim}%12.4f${delim}%4.1f",$Ee,$exp_c,$Ec,$J) .
      (!$no_RCE ? "$delim$lande_g" : '') .
      ($ai_set1_file && $ai_set2_file ? sprintf("${delim}%12.4f${delim}%9.2e",$Eai,$AA_tot) : '') .
      sprintf("${delim}%9.2e", $AR_tot);
    my $i = 0;
    my $need_CR = 1;
    foreach my $num_bas (sort {abs($vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$b})<=>
                              abs($vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$a})}
                          keys %{$vectors[$parity-1]->{$cpl}->{$J}->{$num_e}}) {
      $i++;
      if ( $i > $num_printed_components ) {
        print OUT_LEV_FILE "\n";
        $need_CR = 0;
        last;
      }
      my $A = $vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$num_bas};
      $A *= $A*100;
      if ( ($A < $min_printed_percentage) && (($i > 2) || ($A < $min_printed_second_percentage)) ) {
        print OUT_LEV_FILE "\n";
        $need_CR = 0;
        last;
      }
      $A += 0.5; # For proper rounding
      $A = sprintf("%2d",$A);
      $A =~ s/^ +//g;
      my $shells= $basis[$parity-1]->{$cpl}->{$J}->{$num_bas}->{'sh'};
      # Correct the final parity
      if ($par_code eq 'e') {
        $shells =~ s/[*]$//;
      } else {
        $shells .= '*' unless $shells =~ /[*]$/;
        # Bug fix for odd repeating terms of dn shells
        $shells =~ s/^(.+)\t([^*]+)[*]([^*]+)[*]$/$1\t$2$3*/;
      }

      if ( $i==1 ) {
        my $outer_shell_nl_occ = $basis[$parity-1]->{$cpl}->{$J}->{$num_bas}->{'lastsh'};
        my ($n,$l,$occup) = (ref($outer_shell_nl_occ) eq 'ARRAY') ? @{$outer_shell_nl_occ} : ('','','');

        my ($L,$S) = ('','');
        if ( $cpl eq 'LS' ) {

          my $SL_final = '';
          if ($shells =~ /${delim}([^${delim}]+)$/) {
            $SL_final = $1;
            $SL_final =~ s/[*]//g;
            $SL_final =~ s/\d+$//g; # Strip additional quantum numbers that may be present for d and f shells
            if ( $SL_final =~ /^(\d+)([^0-9]+)$/ ) {
              ($S,$L) = ($1,$2);
              $S = sprintf("%4.1f",($S-1)/2);
              $L = $L_moment{lc($L)}; # Convert the letter code to the integer quantum number
            }
          }
        }
        $l = $L_moment{lc($l)}; # Convert the shell's letter code to the integer quantum number
        print OUT_LEV_FILE "${delim}$n${delim}$l${delim}$occup${delim}$S${delim}$L";
      }

      print OUT_LEV_FILE "${delim}$A${delim}$shells";
    }
    if ( $need_CR ) {
      print OUT_LEV_FILE "\n";
    }
  }

  $num_levs = $lev_num;
  print "Done.\n";
}
close OUT_LEV_FILE;

&PrintTransitions($spectrum, $T, $max_scaled_intensity);

close (OUTG11);

if ( $IDEN ) {
  &PrintENLEV();
  &WriteTRANS_DAT();
}

1;

############################################################################
sub sort_by_par_E($$) {   #02/06/2017 4:20PM
############################################################################
  my ($a,$b) = @_;
  my ($p1,$E1,$k1,$p2,$E2,$k2) = (split(/${delim}/,$a),split(/${delim}/,$b));
  my ($J1,$J2,$num_e);
  if (!$no_RCE) {
    ($J1,$num_e) = @{$E_map{$a}};
    my $num_RCE_E = $map_RCG_RCE[$p1-1]->{$J1}->{$num_e};
    my $hash = $RCE_lev[$p1-1]->{$J1}->[$num_RCE_E-1];
    $E1 = $hash->{'Ec'};
    ($J2,$num_e) = @{$E_map{$b}};
    $num_RCE_E = $map_RCG_RCE[$p2-1]->{$J2}->{$num_e};
    $hash = $RCE_lev[$p2-1]->{$J2}->[$num_RCE_E-1];
    $E2 = $hash->{'Ec'};
  }
  if ( $reverse_order ) {
    if ($E2 != $E1) {
      return ($E2 <=> $E1);
    } elsif ($p1 != $p2) {
      return ($p2 <=> $p1);
    } else {
      return ($J2 <=> $J1);
    }
  } else {
    if ($sp_type !~ /^CM$/i) {
      $E1 += $p1*1e6;
      $E2 += $p2*1e6;
    }
    #return (($sp_type !~/^CM$/i)?$p1*1e6:0)+$E1 <=> (($sp_type !~/^CM$/i)?$p2*1e6:0)+$E2;
    if ($E2 != $E1) {
      return ($E1 <=> $E2);
    } elsif ($p1 != $p2) {
      return ($p1 <=> $p2);
    } else {
      return ($J1 <=> $J2);
    }
  }
} ##sort_by_par_E($$)

############################################################################
sub ReadAI($$)    #1/3/2005 10:47PM A.Kramida
############################################################################
{
  my ($ai_filename, $par) = @_;
  open AI, "<$ai_filename" or die "Could not open input file $ai_filename";
  my $s = <AI>; # Read the header line
  chomp $s;
  my @AI_channels = split(/${delim}/,$s);
  $AI_cores[$par-1] = \@AI_channels;
  my $num_channels = $#AI_channels - 3;
  foreach my $J1 (sort {$a<=>$b} keys %{$energies[$par-1]}) {
    foreach my $num_state (sort {$a<=>$b} keys %{$energies[$par-1]->{$J1}}) {
      $s = <AI>;
      if ( !defined($s) ) {
        die "Unexpected end of file $ai_filename";
      }
      chomp $s;
      my @line = split(/${delim}/, $s);
      my ($J, $E, $AA_tot) = ($line[0], $line[1], $line[$#line]);
      my $DR = 0;
      if ( $s=~ /:/ ) {
        $E = $E;
      }
      for ( my $i = 0; $i <= $num_channels; $i++ ) {
        my $AI_rate = $line[$i+2];
        #if ($AI_channels[$i+2] =~ /$DR_channel/) {
        if ( $i == 0 ) {
          # The first AI channel must be to the ground state of the ionized atom.
          # The DR rate is proportional to AI rate to this channel.
          $DR += $AI_rate;
        }
        if ( $AI_rate > 0 ) {
          $AI_rates[$par-1]->{$J} = {} unless defined($AI_rates[$par-1]->{$J});
          $AI_rates[$par-1]->{$J}->{$num_state} = {} unless defined($AI_rates[$par-1]->{$J}->{$num_state});
          $AI_rates[$par-1]->{$J}->{$num_state}->{$AI_channels[$i+2]} = $AI_rate;
        }
      }
      if ( $J ne $J1 ) {
        die "J value mismatch between OUTG11 and $ai_filename on the following line:\n$s";
      }
      push(@{$energies[$par-1]->{$J}->{$num_state}}, $E, $AA_tot, 0, $DR);
    }
  }
} ##ReadAI($$)

############################################################################
sub ReadTransitions()   #1/5/2005 12:02PM A.Kramida
############################################################################
{
  my $s;
  my $spectrum = '';
  print "Reading transitions...";
  close(OUTG11);
  open OUTG11, "<OUTG11" or die "Could not open input file OUTG11";
  my %trans_hash = ();
  my $line_no = 0;

  while ( 1 ) {
    while ( (defined($s = <OUTG11>)) && ($s !~ /^. +([^ ]+) +([^ ]+) +SPECTRUM +[()]ENERGIES IN UNITS OF/) ) {
      $line_no++;
      next;
    }
    $line_no++;
    last unless defined $s;
    if ($line_no > 1090866) {
      $s = $s;
    }
    $s =~ /^. +([^ ]+) +([^ ]+) +SPECTRUM +[()]ENERGIES IN UNITS OF/;
    my $spectrum_type = "$1 $2";
    my $same_parity = ($spectrum_type =~ /MAG DIP|ELEC QUD/) ? 1 : 0;
    my $fmt_old = '^. +\* \* \* +([0-9.-]+) +([0-9.]+) (.{3}) (.{8})  (.{8}) (.{8}) +\* \* \*';
    my $fmt_new = '^. +\* \* \* +([0-9.-]+) +([0-9.]+) (.{3}) (.{8}) (.{6}) (.{5}) (.{8}) +\* \* \*';
    my $wl_order_pattern = '^[0-9]           E      J   CONF              EP      JP  CONFP           DELTA E   LAMBDA';
    while ( (defined($s = <OUTG11>)) && (($s !~ /$fmt_old/) && ($s !~ /$fmt_new/)) ) {
      if ($s =~ /$wl_order_pattern/) {
        last;
      }
      $line_no++;
      if ($s =~ /^. +([^ ]+) +([^ ]+) +SPECTRUM +[()]ENERGIES IN UNITS OF/) {
        $line_no = $line_no;
      }
      next;
    }
    $line_no++;
    last unless defined $s;
    next if ($s =~ /$wl_order_pattern/);
# 112     44.164  0.5100 (3D) 4D       37.123   2693.757  -2.972 9.808E+05 -0.05    54
#1218     56.407  2.0  4 (2S) 1D        0.0462166579.439  -3.864 1.941E-01 -0.03     8
#1245     61.616  1.0 12 (2D) 1P        0.007***********  -5.441 1.107E-04  0.22     8
    my $fmt_fin_old = '^ *(\d+|\*+) +([0-9.-]+) +([0-9]+\.[0-9]) *(\d+) (.{8}) +([0-9]+\.[0-9]{3}) *([0-9.*]+) +([0-9.-]+) ([0-9.E+-]+) +([0-9.-]+)';
    my $fmt_fin_new = '^ *(\d+|\*+) +([0-9.-]+) +([0-9]+\.[0-9]) *(\d+) (.{8}) +([0-9]+\.[0-9]{3}) *([0-9.*]+) +([0-9.-]+) ([0-9.E+-]+) +([0-9.-]+) + ([0-9]+)';
    my $fmt_init = $fmt_old;
    my $fmt_fin = $fmt_fin_old;
    my $new_format = 0;
    if ( $s =~ /$fmt_new/ ) {
      $fmt_init = $fmt_new;
      $fmt_fin = $fmt_fin_new;
      $new_format = 1;
    }
    while ($s =~ /$fmt_init/) {
      my ($E, $J, $nc, $term, $n1, $spectrum_name, $conf) = ($new_format ?
        ($1,$2,$3,$4,$5,$6,$7) : ($1,$2,$3,$4,$7,$5,$6));
      $E =~ s/^\s+|\s+$//g;
      next if ($E+0 > $upper_lev_bound);
      $nc =~ s/^\s+|\s+$//g;
      $term =~ s/^\s+|\s+$//g;
      $conf =~ s/^\s+|\s+$//g;
      $spectrum_name =~ s/^\s+|\s+$//g;
      $n1 =~ s/^\s+|\s+$//g;
      $spectrum = $spectrum_name unless $spectrum;

      die "Unknown conf. name $conf in OUTG11 $spectrum_type section" unless defined($conf_nums{$conf});
      my ($first_par,$nc1) = @{$conf_nums{$conf}};
      die "Conf. number mismatch in OUTG11 $spectrum_type section: read $nc, should be $nc1" unless $nc == $nc1;

      my $second_par = ($same_parity) ? $first_par : 3 - $first_par;

      # Identify the initial level
      my ($n_lev1, $E1, $lead_c_no, $lead_term) = ($new_format ? &FindLev1($n1,$first_par, $J, $E, $nc, $term, $line_no, 1) :
        &FindLev($first_par, $J, $E, $nc, $term));

      if ( !$lead_c_no ) {
        die "First level $n1 not identified in OUTG11:\nJ= $J, E=$E, level $conf $term";
      }
      if ( ($lead_c_no != $nc) || ($lead_term ne $term)  ) {
        #die "First level parameters mismatch: J = $J, conf. set $first_par; in transitions section, E = $E, conf_no = $nc, term = $term\n" .
        #  "In eigenvectors section: E = $E1, conf_no = $lead_c_no, term = $lead_term";
      }

      $s = <OUTG11>; # Skip one blank line
      $line_no++;
      chomp $s;
      # Read transitions from this level
      while ( (defined($s = <OUTG11>)) && ($s !~ /^\s*$/) && ($s !~ /SUMFI/) && ($s !~ /^ +NO. OF LINES/) ) {
        $line_no++;
        if ($s !~ /$fmt_fin/i) {
          next;
        }
        $s =~ /$fmt_fin/i;
        my ($E2, $J2, $nc2, $term2, $dE, $lambda, $log_gf, $gA, $cf,$n2) = ($2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
        $n1 =~ s/^\s+|\s+$//g;
        $E2 =~ s/^\s+|\s+$//g;
        next if ($E2+0 > $upper_lev_bound);
        $nc2 =~ s/^\s+|\s+$//g;
        $term2 =~ s/^\s+|\s+$//g;
        $n2 =~ s/^\s+|\s+$//g;

        if (($line_no == 1090940) || ($s =~/    1      0.097  4.0  1 (2F) 3F       36.689   2725.617  -9.102 7.099E-01  0.00     1/)) {
          $s = $s;
        }
        
        # Identify the final level
        my ($n_lev2, $E21, $lead_c_no2, $lead_term2) = ($new_format ? &FindLev1($n2, $second_par, $J2, $E2, $nc2, $term2, $line_no, 2) :
          &FindLev($second_par, $J2, $E2, $nc2, $term2));

        if ( !$lead_c_no2 ) {
          die "Final level $n2 not identified in OUTG11:\nJ= $J2, E=$E2, level $nc2 $term2";
        }
        if ( ($lead_c_no2 != $nc2) || ($lead_term2 ne $term2)  ) {
          #die "\nFinal level parameters mismatch: J = $J2, conf. set $second_par; in transitions section, E = $E2, conf_no = $nc2, term = $term2;\n" .
          #  "In eigenvectors section: E = $E21, conf_no = $lead_c_no2, term = $lead_term2;";
        }

        my ($par1, $par2, $E_init, $E_fin, $J_init, $J_fin, $n_init, $n_fin) = ($first_par, $second_par, $E, $E2, $J, $J2, $n_lev1, $n_lev2);
        if ( $E < $E2 ) {
          # Swap the levels
          my ($n_tmp, $par_tmp, $J_tmp, $E_tmp) = ($n_init, $par1, $J_init, $E_init);
          ($n_init, $par1, $J_init, $E_init) = ($n_fin, $par2, $J_fin, $E_fin);
          ($n_fin, $par2, $J_fin, $E_fin) = ($n_tmp, $par_tmp, $J_tmp, $E_tmp);
        }
        # Store the transition data
        $transitions[$par1-1]->{$J_init} = {} unless defined $transitions[$par1-1]->{$J_init};
        $transitions[$par1-1]->{$J_init}->{$n_init} = [] unless defined $transitions[$par1-1]->{$J_init}->{$n_init};
        if ($lambda =~ /\*/) {
          if ($E_init < $E_fin) {
            my ($n_tmp, $par_tmp, $J_tmp, $E_tmp) = ($n_init, $par1, $J_init, $E_init);
              ($n_init, $par1, $J_init, $E_init) = ($n_fin, $par2, $J_fin, $E_fin);
              ($n_fin, $par2, $J_fin, $E_fin) = ($n_tmp, $par_tmp, $J_tmp, $E_tmp);
          }
          my $dE = $E_init - $E_fin;
          next if ($dE == 0);
          $lambda = 1e5/$dE;
        }
        my $I = $gA/$lambda;
        if ( $sp_type =~ /LTEABS/ ) {
          my $S = &SfromgA($gA,$lambda,'');
          my $gf = &gffromS($S,$lambda,'');
          $I = $gf*$lambda*$lambda; # Reduced absorption intensity; see https://www.nist.gov/pml/atomic-spectroscopy-spectral-lines#node172
          $I *= exp(-$E_fin/$T); # Boltzmann-population absorption intensity
        } elsif ( $sp_type =~ /LTE/ ) {
          $I *= exp(-$E_init/$T); # Boltzmann-population emission intensity
        }
        my $tr_key = "$par1 $J_init $n_init $par2 $J_fin $n_fin";
        if (defined($trans_hash{$tr_key})) {
          my ($gA_prev,$ind) = @{$trans_hash{$tr_key}};
          if ($gA >= $gA_prev) {
            # Replace the old gA with the new one.
            # Ignore the smaller gA of the M1 and E2 transitions between the same levels. 
            $trans_hash{$tr_key} = [$gA, $ind];
            $transitions[$par1-1]->{$J_init}->{$n_init}->[$ind] = [$par2, $J_fin, $n_fin, $dE, $lambda, $gA, $cf, $I];
            $energies[$par1-1]->{$J_init}->{$n_init}->[3] += ($gA - $gA_prev)/(2*$J_init +1);
          }
        } else {
          push(@{$transitions[$par1-1]->{$J_init}->{$n_init}}, [$par2, $J_fin, $n_fin, $dE, $lambda, $gA, $cf, $I]);
          my $ind = scalar @{$transitions[$par1-1]->{$J_init}->{$n_init}} - 1;
          $trans_hash{$tr_key} = [$gA, $ind];
        }
          
        if ( !defined($energies[$par1-1]->{$J_init}->{$n_init}->[3]) ) {
          $energies[$par1-1]->{$J_init}->{$n_init}->[3] = 0;
        }
        # Accumulate the total radiative decay rate for the upper level
        $energies[$par1-1]->{$J_init}->{$n_init}->[3] += $gA/(2*$J_init +1);

      }
      $line_no++;
      $s = <OUTG11>; # Skip the SUMFI line
      $line_no++;
      $s = <OUTG11>; # Read next initial level line
      $line_no++;
      $s = <OUTG11>; # Skip one blank line
      $line_no++;
    }
    next;  # Proceed to next initial level
  }
  # Scale the intensities
  # First step: find the max. intensity
  my $Imax = 0;
  if ( $sp_type !~ /AUGER/ ) {
    for (my $par = 1; $par <= 2; $par++) {
      foreach my $J (sort {$a<=>$b} keys %{$energies[$par-1]}) {
        foreach my $n_lev1 (sort {$a<=>$b} keys %{$energies[$par-1]->{$J}}) {
          next unless defined $transitions[$par-1]->{$J}->{$n_lev1};
          my ($E,$E1,$AA_tot, $AR_tot, $DR) = @{$energies[$par-1]->{$J}->{$n_lev1}};
          my $lft_ns = 1e9/($AA_tot + $AR_tot); # Lifetime of the upper level in nanoseconds
          my $rad_branch = $AR_tot/($AA_tot + $AR_tot);

          my @trans = @{$transitions[$par-1]->{$J}->{$n_lev1}};
          my $num_trans = $#trans;
          for ( my $i = 0; $i <= $num_trans; $i++) {
            my ($second_par, $J2, $n_lev2, $dE, $lambda, $gA, $cf, $I) = @{$trans[$i]};
            # Substitute Ec from RCE output and recalculate lambda if RCE output is available
            if ( !$no_RCE ) {
              my $num_RCE_E = $map_RCG_RCE[$par-1]->{$J}->{$n_lev1};
              my $hash = $RCE_lev[$par-1]->{$J}->[$num_RCE_E-1];
              my $Ec1 = $hash->{'Ec'};

              $num_RCE_E = $map_RCG_RCE[$second_par-1]->{$J2}->{$n_lev2};
              $hash = $RCE_lev[$second_par-1]->{$J2}->[$num_RCE_E-1];
              my $Ec2 = $hash->{'Ec'};
              $lambda = ($Ec2 == $Ec1) ? 9.9999e7 : 1e5/abs($Ec2-$Ec1);
              $dE = abs($Ec2-$Ec1);
              #if (($dE > 5.0) && ($dE < 50.0)) {
              #  $lambda = &Lair(1e8/$dE);
              #}
              $trans[$i]->[4] = $lambda;
            }

            if ( $sp_type =~ /BF/ ) {
              # Suppose that the decrease in population vs. Boltzmann is proportional
              # to the branching fraction of radiative vs autoionization decay rates
              # my $rad_branch = $AR_tot/($AA_tot + $AR_tot);
              $I *= $lft_ns * (1-exp(-$t_flight/$lft_ns)); # * $rad_branch;

            } elsif ($sp_type =~ /DR/) {
              $I *= $lft_ns * $DR;
            } elsif ($sp_type =~ /AUGER/) {
            } else {
              $I *= $rad_branch;
            }
            $transitions[$par-1]->{$J}->{$n_lev1}->[$i]->[7] = $I;
            if ( $I > $Imax ) {
              $Imax = $I;
            }
          }
        }
      }
    }
    # Scale the intensities
    # Second step: divide all intensities by the max. intensity and multiply by max. scaled intensity
    for (my $par = 1; $par <= 2; $par++) {
      foreach my $J (sort {$a<=>$b} keys %{$energies[$par-1]}) {
        foreach my $n_lev1 (sort {$a<=>$b} keys %{$energies[$par-1]->{$J}}) {
          next unless defined $transitions[$par-1]->{$J}->{$n_lev1};
          my ($E,$E1,$AA_tot, $AR_tot) = @{$energies[$par-1]->{$J}->{$n_lev1}};
          my @trans = @{$transitions[$par-1]->{$J}->{$n_lev1}};
          my $num_trans = $#trans;
          for ( my $i = 0; $i <= $num_trans; $i++) {
            my ($second_par, $J2, $n_lev2, $dE, $lambda, $gA, $cf, $I) = @{$trans[$i]};
            # Suppose that the decrease in population vs. Boltzmann is proportional
            # to the branching fraction of total radiative decay vs autoionization decay rates
            $I *= $max_scaled_intensity/$Imax;
            $transitions[$par-1]->{$J}->{$n_lev1}->[$i]->[7] = $I;
          }
        }
      }
    }
  } else {
    @transitions = ([],[]);
    for (my $par = 1; $par <= 2; $par++) {
      my @AI_channels = @{$AI_cores[$par-1]};
      my $num_cores = $#AI_channels-3;
      $transitions[$par-1] = {};
      foreach my $J (sort {$a<=>$b} keys %{$energies[$par-1]}) {
        $transitions[$par-1]->{$J} = {};
        foreach my $n_lev1 (sort {$a<=>$b} keys %{$energies[$par-1]->{$J}}) {
          $transitions[$par-1]->{$J}->{$n_lev1} = {};
          my ($E,$E1,$AA_tot, $AR_tot) = @{$energies[$par-1]->{$J}->{$n_lev1}};
          next unless ( $AA_tot > 0 );
          for ( my $num_channel = 0; $num_channel <= $num_cores; $num_channel++ ) {
            next unless defined($AI_rates[$par-1]->{$J}->{$n_lev1}->{$AI_channels[$num_channel+2]});
            my $AI_rate = $AI_rates[$par-1]->{$J}->{$n_lev1}->{$AI_channels[$num_channel+2]};
            next unless ($AI_rate > 0);
            my $I = (2*$J + 1) * $AI_rate/($AA_tot + $AR_tot);
            $transitions[$par-1]->{$J}->{$n_lev1}->{$num_channel+2} = [$I,$AI_rate];
            if ( $I > $Imax ) {
              $Imax = $I;
            }
          }
        }
      }
    }
    # Scale intensities
    for (my $par = 1; $par <= 2; $par++) {
      my @AI_channels = @{$AI_cores[$par-1]};
      my $num_cores = $#AI_channels-3;
      foreach my $J (sort {$a<=>$b} keys %{$energies[$par-1]}) {
        foreach my $n_lev1 (sort {$a<=>$b} keys %{$energies[$par-1]->{$J}}) {
          next unless defined $transitions[$par-1]->{$J}->{$n_lev1};
          my ($E,$E1,$AA_tot, $AR_tot) = @{$energies[$par-1]->{$J}->{$n_lev1}};
          next unless ($AA_tot > 0);
          for ( my $num_channel = 0; $num_channel <= $num_cores; $num_channel++ ) {
            next unless defined($transitions[$par-1]->{$J}->{$n_lev1}->{$num_channel+2});
            #my ($AI_rate, $E_kin) = @{$AI_rates[$par-1]->{$J}->{$n_lev1}}->{$AI_channels[$num_channel+2]};
            my ($I, $AI_rate) = @{$transitions[$par-1]->{$J}->{$n_lev1}->{$num_channel+2}};
            $I *= $max_scaled_intensity/$Imax;
            $transitions[$par-1]->{$J}->{$n_lev1}->{$num_channel+2}->[0] = $I;
          }
        }
      }
    }
  }
  print "Done.\n";
  return $spectrum;
} ##ReadTransitions()

############################################################################
sub PrintTransitions($$$)   #1/6/2005 8:44AM A.Kramida
############################################################################
{
  my ($spectrum, $T, $max_scaled_intensity) = @_;

  # Create sorted hash of lines
  print "Sorting lines..." unless ($sp_type =~ /^CM$/i);
  if ( $sp_type =~ /AUGER/i ) {
    for (my $par = 1; $par <= 2; $par++) {
      my @AI_channels = @{$AI_cores[$par-1]};
      my $num_cores = $#AI_channels-3;
      foreach my $J (sort {$a<=>$b} keys %{$energies[$par-1]}) {
        foreach my $n_lev1 (sort {$a<=>$b} keys %{$energies[$par-1]->{$J}}) {
          next unless defined $transitions[$par-1]->{$J}->{$n_lev1};
          my ($E,$E1,$AA_tot, $AR_tot, $DR) = @{$energies[$par-1]->{$J}->{$n_lev1}};
          #my $AI_rate = $AI_rates[$par-1]->{$J}->{$n_lev1};
          next unless ($AA_tot > 0);
          for ( my $num_channel = 0; $num_channel <= $num_cores; $num_channel++ ) {
            next unless defined($transitions[$par-1]->{$J}->{$n_lev1}->{$num_channel+2});
            my ($I,$AI_rate) = @{$transitions[$par-1]->{$J}->{$n_lev1}->{$num_channel+2}};
            next unless $I >= 0.5;
            my $core = $AI_channels[$num_channel+2];
            my ($core_name, $E2) = split(/:/,$core);
            my $E_kin = $E - $E2;
            my $key = sprintf("%10.3f", $E_kin *1000 / $ev_cm); # Translate to eV
            while ( defined($lines{$key}) ) {
              $key .= '1';
            }
            $lines{$key} = [$par, $J, $n_lev1, 1e9/($AA_tot+$AR_tot), $core_name, $E2, $I, $AI_rate];
          }
        }
      }
    }
  } elsif ( $sp_type !~ /^CM$/i ) {
    for (my $par = 1; $par <= 2; $par++) {
      foreach my $J (sort {$a<=>$b} keys %{$energies[$par-1]}) {
        foreach my $n_lev1 (sort {$a<=>$b} keys %{$energies[$par-1]->{$J}}) {
          next unless defined $transitions[$par-1]->{$J}->{$n_lev1};
          my ($E,$E1,$AA_tot, $AR_tot, $DR) = @{$energies[$par-1]->{$J}->{$n_lev1}};
          my @trans = @{$transitions[$par-1]->{$J}->{$n_lev1}};
          my $num_trans = $#trans;
          for ( my $i = 0; $i <= $num_trans; $i++) {
            my ($second_par, $J2, $n_lev2, $dE, $lambda, $gA, $cf, $I) = @{$trans[$i]};
            if ( $I >= 0.5 ) {
              my $key = sprintf("%13.4f", $lambda);
              while ( defined($lines{$key}) ) {
                $key .= '1';
              }
              $lines{$key} = [$par, $J, $n_lev1, 1e9/($AA_tot+$AR_tot), $second_par, $J2, $n_lev2, $gA, $cf, $I,$AA_tot,$AR_tot];
            }
          }
        }
      }
    }
  } else {

    print "Printing Cascades Matrix...";
    open OUT_TR_FILE, ">$out_tr_file" or die "Could not open output file " . $out_tr_file;

    my $num_str = '';
    for ( my $i = 1; $i <= $num_levs; $i++ ) {
      $num_str .= ($num_str ? $delim : '') . $i;
    }
    print OUT_TR_FILE "CM${delim}$num_str\n";

    foreach my $key ( sort {my ($J1,$num_e1,$Eai1,$AA_tot1,$AR_tot1,$lev_n1) = @{$E_map{$a}};
                            my ($J2,$num_e2,$Eai2,$AA_tot2,$AR_tot2,$lev_n2) = @{$E_map{$b}};
                            $lev_n1<=>$lev_n2}
                     keys %E_map ) {
      my ($par,$E,$k) = split(/${delim}/, $key);
      my ($J,$num_e,$Eai,$AA_tot,$AR_tot,$lev_n1) = @{$E_map{$key}};

      $num_str = '';
      my @A = ();
      for ( my $i = 1; $i <= $num_levs; $i++ ) {
        push(@A,0);
        $num_str .= ($num_str ? $delim : '') . $i;
      }

      if (defined $transitions[$par-1]->{$J}->{$num_e}) {
        my @trans = @{$transitions[$par-1]->{$J}->{$num_e}};
        my $num_trans = $#trans;
        for ( my $i = 0; $i <= $num_trans; $i++) {
          my ($second_par, $J2, $n_lev2, $dE, $lambda, $gA, $cf, $I) = @{$trans[$i]};
          my ($key2,$lev_num2) = @{$E_reverse_map{$second_par}->{$J2}->{$n_lev2}};
          #my ($par2,$E2,$k2) = split(/${delim}/, $key2);
          $A[$lev_num2-1] = sprintf("%10.3e",$gA/(2*$J+1));
        }
      }
      print OUT_TR_FILE "$lev_n1${delim}",join($delim,@A),"\n";
    }
  }
  print "Done.\n";

  if ( $sp_type =~ /^CM$/i ) {
    return;
  }

  # Print lines
  print "Printing lines...";
  open OUT_TR_FILE, ">$out_tr_file" or die "Could not open output file " . $out_tr_file;

       # Print transition file headers
  if ( $sp_type =~ /LTE/ ) {
    my $LTE_type = ($sp_type =~ /LTEABS/) ? 'LTE absorption' : 'LTE';
    print OUT_TR_FILE "Spectrum${delim}$spectrum${delim}${delim}$LTE_type eff. temperature $T kK${delim}${delim}${delim}${delim}Max. intens.${delim}${delim}$max_scaled_intensity\n";
    if ( $no_RCE ) {
      if ( $ai_set1_file && $ai_set2_file ) {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','tot AI','sum A2i','A21','BRrad') . "\n";
      } else {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','sum A2i','A21','BRrad') . "\n";
      }
    } else {
      if ( $ai_set1_file && $ai_set2_file ) {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','E1_exp','E2_exp','wl_exp(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','tot AI','sum A2i','A21','BRrad','lande1','lande2') . "\n";
      } else {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','E1_exp','E2_exp','wl_exp(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','sum A2i','A21','BRrad','lande1','lande2') . "\n";
      }
    }
  } elsif ( $sp_type =~ /AUGER/ ) {
    print OUT_TR_FILE "Spectrum${delim}$spectrum${delim}${delim}AUGER electron spectrum${delim}${delim}${delim}${delim}Max. intens.${delim}${delim}$max_scaled_intensity\n";
    if ( $no_RCE ) {
      print OUT_TR_FILE join($delim,"up lev${delim}",'J_up','AI_channel','E_up','E_lo','dE(eV)','id_up','Icalc','AI_rate(s-1)', 'lft(ns)') . "\n";
    } else {
      print OUT_TR_FILE join($delim,"up lev${delim}",'J_up','AI_channel','E_up','E_lo','dE(eV)','id_up','Eup_exp','dE_exp(eV)','Icalc', 'AI_rate(s-1)', 'lft(ns)') . "\n";
    }
  } else {
    print OUT_TR_FILE "Spectrum${delim}$spectrum${delim}${delim}$sp_type${delim}${delim}${delim}${delim}Max. intens.${delim}${delim}$max_scaled_intensity\n";
    if ( $no_RCE ) {
      if ( $ai_set1_file && $ai_set2_file ) {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','tot AI','sum A2i','A21','BRrad') . "\n";
      } else {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','sum A2i','A21','BRrad') . "\n";
      }
    } else {
      if ( $ai_set1_file && $ai_set2_file ) {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','E1_exp','E2_exp','wl_exp(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','tot AI','sum A2i','A21','BRrad') . "\n";
      } else {
        print OUT_TR_FILE join($delim,"level1${delim}",'J1',"level2${delim}",'J2','E1','E2','lid1','lid2','wl(A)','E1_exp','E2_exp','wl_exp(A)','Icalc','g2A21', 'cf','lft1(ns)','lft2(ns)','sum A2i','A21','BRrad') . "\n";
      }
    }
  }

  # Print lines
  foreach my $wl (sort keys %lines) {
    if ( $sp_type !~ /AUGER/ ) {
      my ($par, $J, $n_lev1, $lft_ns, $second_par, $J2, $n_lev2, $gA, $cf, $I,$AA_tot, $AR_tot) = @{$lines{$wl}};
      my ($E,$E1,$AA_tot, $AR_tot,$DR) = @{$energies[$par-1]->{$J}->{$n_lev1}};
      my $lid1 = $lev_nums[$par-1]->{$J}->{$n_lev1};
      my $cpl = "LS";
      my $shells = &get_leading_term($par,$cpl,$J,$n_lev1); # Take the leading term as the level designation
      my ($E2,$E21,$AA_tot2, $AR_tot2,$DR2) = @{$energies[$second_par-1]->{$J2}->{$n_lev2}};
      my $lid2 = $lev_nums[$second_par-1]->{$J2}->{$n_lev2};
      my ($Ee1, $Ee2, $Ec1, $Ec2, $ec1, $ec2, $lande1, $lande2) = ('','',$E,$E2,' ',' ','','');

      # Substitute Eexp from RCE if available
      if ( !$no_RCE ) {
        my $num_RCE_E = $map_RCG_RCE[$par-1]->{$J}->{$n_lev1};
        my $hash = $RCE_lev[$par-1]->{$J}->[$num_RCE_E-1];
        $ec1 = $hash->{'exp_c'};  # Star in $ec1 means "no experimental value"
        $Ee1 = ($ec1) ? '' : $hash->{'Ee'};
        $Ec1 = $hash->{'Ec'};
        $lande1 = $hash->{'lande'};

        $num_RCE_E = $map_RCG_RCE[$second_par-1]->{$J2}->{$n_lev2};
        $hash = $RCE_lev[$second_par-1]->{$J2}->[$num_RCE_E-1];
        $ec2 = $hash->{'exp_c'};  # Star in $ec2 means "no experimental value"
        $Ee2 = ($ec2) ? '' : $hash->{'Ee'};
        $Ec2 = $hash->{'Ec'};
        $lande2 = $hash->{'lande'};
      }

      my $lft2_ns = 1e27;
      if ( $AA_tot2 + $AR_tot2 > 1e-19 ) {
        $lft2_ns = 1e9/($AA_tot2 + $AR_tot2); # Lifetime of the lower level in nanoseconds
      }

      my $shells2 = &get_leading_term($second_par,$cpl,$J2,$n_lev2); # Take the leading term as the level designation

      #my $dE = abs($E - $E2);
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
      my $I_out = sprintf("%10.0f",$I);
      my $AA_out = sprintf("%10.2e",$AA_tot);
      my $AR_out = sprintf("%10.2e",$AR_tot);
      my $BR_out = sprintf("%5.3f",$gA/(2*$J+1)/($AA_tot+$AR_tot));
      my $A_out = sprintf("%10.2e",$gA/(2*$J+1));
      my $lft_out = sprintf("%10.2e",$lft_ns);
      my $lft2_out = sprintf("%10.2e",$lft2_ns);

      if ( $no_RCE ) {
        if ( $ai_set1_file && $ai_set2_file ) {
          print OUT_TR_FILE join($delim,$shells2,$J2,$shells,$J,$Ec2,$Ec1,$lid2,$lid1,$wl_out,$I_out,$gA, $cf,$lft2_out,$lft_out,$AA_out,$AR_out,$A_out,$BR_out) . "\n";
        } else {
          print OUT_TR_FILE join($delim,$shells2,$J2,$shells,$J,$Ec2,$Ec1,$lid2,$lid1,$wl_out,$I_out,$gA, $cf,$lft2_out,$lft_out,$AR_out,$A_out,$BR_out) . "\n";
        }
      } else {
        if ( $ai_set1_file && $ai_set2_file ) {
          print OUT_TR_FILE join($delim,$shells2,$J2,$shells,$J,$Ec2,$Ec1,$lid2,$lid1,$wl_out,$Ee2,$Ee1,$wl_out_e,$I_out,$gA, $cf,$lft2_out,$lft_out,$AA_out,$AR_out,$A_out,$BR_out)
          . ( $sp_type =~ /LTE/ ? "$delim$lande2$delim$lande1" : '')
          . "\n";
        } else {
          print OUT_TR_FILE join($delim,$shells2,$J2,$shells,$J,$Ec2,$Ec1,$lid2,$lid1,$wl_out,$Ee2,$Ee1,$wl_out_e,$I_out,$gA, $cf,$lft2_out,$lft_out,$AR_out,$A_out,$BR_out)
          . ( $sp_type =~ /LTE/ ? "$delim$lande2$delim$lande1" : '')
          . "\n";
        }
      }

    } else {

        # Auger transitions
      my ($par, $J, $n_lev1, $lft_ns, $AI_channel, $E2, $I, $AI_rate) = @{$lines{$wl}};
      my ($E,$E1,$AA_tot, $AR_tot,$lid1) = @{$energies[$par-1]->{$J}->{$n_lev1}};
      my $par_code = $parities[$par-1];
      my $cpl = "LS";
      my ($Ee1, $ec1) = ('',' ');
      #my ($Ec1, $Ec2) = ($E,$E1);

      # Substitute Eexp from RCE if available
      if ( !$no_RCE ) {
        my $num_RCE_E = $map_RCG_RCE[$par-1]->{$J}->{$n_lev1};
        #if (ref($RCE_lev[$par-1]->{$J}) ne 'HASH') {
        #  $J=$J;
        #}
        #$ec1 = $RCE_lev[$par-1]->{$J}->{$num_RCE_E}->{'exp_c'};  # Star in $ec1 means "no experimental value"
        #$Ee1 = ($ec1) ? '' : $RCE_lev[$par-1]->{$J}->{$num_RCE_E}->{'Ee'};
        my $hash = $RCE_lev[$par-1]->{$J}->[$num_RCE_E-1];
        $ec1 = $hash->{'exp_c'};  # Star in $ec1 means "no experimental value"
        $Ee1 = ($ec1) ? '' : $hash->{'Ee'};
        #$Ec1 = $hash->{'Ec'};
      }

      my $shells = &get_leading_term($par,$cpl,$J,$n_lev1); # Take the leading term as the level designation

      my $dE = abs($E - $E2);
      my ($wl_out,$wl_out_e) = ($dE*1000.0/$ev_cm, '');
      my $wl_out = sprintf("%10.3f",$wl_out);
      if ( !$ec1 ) {
        $wl_out_e = sprintf("%12.4f",abs($Ee1-$E2)*1000/$ev_cm);
      }
      $wl_out =~ s/^\s+//g;
      $wl_out_e =~ s/^\s+//g;
      my $I_out = sprintf("%10.0f${delim}%9.2e",$I,$AI_rate);
      my $lft_out = sprintf("%10.2e",$lft_ns);
      if ( $no_RCE ) {
        print OUT_TR_FILE join($delim,$shells,$J,$AI_channel,$E,$E2,$wl_out,$lid1,$I_out,$lft_out) . "\n";
      } else {
        print OUT_TR_FILE join($delim,$shells,$J,$AI_channel,$E,$E2,$wl_out,$lid1,$Ee1,$wl_out_e,$I_out,$lft_out) . "\n";
      }
    }
  }
  close OUT_TR_FILE;
  print "Done.\n";
} ##PrintTransitions

############################################################################
sub get_leading_term    #8/10/2006 1:34PM A.Kramida
############################################################################
{
  my ($par,$cpl,$J,$n_lev) = @_;
  # Take the leading term as the level designation
  my $par_code = $parities[$par-1];
  my $shells = '';
  if ($n_lev == 7) {
    $n_lev = $n_lev;
  }
  foreach my $num_bas (sort {abs($vectors[$par-1]->{$cpl}->{$J}->{$n_lev}->{$b})<=>
                            abs($vectors[$par-1]->{$cpl}->{$J}->{$n_lev}->{$a})}
                        keys %{$vectors[$par-1]->{$cpl}->{$J}->{$n_lev}}) {
    $shells = $basis[$par-1]->{$cpl}->{$J}->{$num_bas}->{'sh'};
    # Correct the final parity
    if ($par_code eq 'e') {
      $shells =~ s/[*]$//;
    } else {
      $shells .= '*' unless $shells =~ /[*]$/;
      $shells =~ s/\*(\d+\*)$/$1/;
    }
    last;
  }
  return $shells;
} ##get_leading_term

# ############################################################################
# sub fix_trailing_index($$$$)    #3/30/2006 11:48AM A.Kramida
# ############################################################################
# {
#   my ($parity, $J, $basis_state_num, $term1) = @_;
#     #my $fill_str = $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'shells'};
#     #my $final_term_label = '';
#     #if ( $fill_str =~ /\t([^*]+)[*]*$/ ) {
#     #  $final_term_label = $1;
#     #}
#     #my ($bare_term, $index) = ('','');
#     #if ( $term1 =~ /(..)([^A-Z])*$/ ) {
#     #  ($bare_term, $index) = ($1,$2);
#     #}
#     #if ( $final_term_label ne $bare_term ) {
#     #  die "Final term label mismatch for basis state number $basis_state_num, J = $J: $final_term_label from genealogy, $bare_term in energy matrix.";
#     #}
#     #if ( $index =~ /\d/ ) {
#       # Add the missing Nielson & Koster index for the dn and fn configs
#     #  $fill_str .= $index;
#     #  $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'shells'} = $fill_str;
#     #}
#     my %complete_df = ('d' => 10, 'f' => 14);
#     if ( $term1 =~ /([0-9])$/ ) {
#       my $seniority = $1;
#       my $sh = $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'sh'};
#       if ( $sh =~ /4f5/ ) {
#         $sh = $sh;
#       }
#       foreach my $df ('d','f') {
#         if ( $sh !~ /$df(\d+)/ ) {
#           next;
#         }
#         my $occup = 0;
#         my ($df_shell, $df_replace) = ('','');
#         my $complete_shell = $complete_df{$df};
#         for ( my $i = 2; $i<= $complete_shell; $i++ ) {
#           if ( ($sh =~ /$df$i\.\((..)\)/) && ($i > $occup) ) {
#             # Incomplete d or f shell followed by another incomplete shell
#             $occup = $i;
#             $df_shell = "$df$i\.\\($1\\)";
#             $df_replace = "$df$i.($1$seniority)";
#           } elsif ( ($sh =~ /$df$i\.<(..)>/) && ($i > $occup) ) {
#             # Incomplete d or f shell followed by another incomplete shell
#             $occup = $i;
#             $df_shell = "$df$i\.\<$1\>";
#             $df_replace = "$df$i.($1$seniority)";
#           } elsif ( ($sh =~ /$df$i( |.<)(..)(>{0,1})\t/) && ($i > $occup) ) {
#             # Incomplete d or f shell which is the last open shell
#             $occup = $i;
#             $df_shell = "$df$i$1$2$3";
#             $df_replace = "$df$i$1$2$seniority$3";
#           } elsif ( ($sh =~ /$df$i\.(\d+s2|\d+p6|\d+d10|\d+f14)(\S*)( |.<)(..)(>{0,1})\t/) && ($i > $occup) ) {
#             # Incomplete d or f shell followed by a complete shell
#             $occup = $i;
#             $df_shell = "$df$i.$1$2$3$4$5";
#             my ($last_sh, $fin_term,$end) = ("$1$2","$3$4",$5);
#             $df_shell =~ s/\(/\\(/g;
#             $df_shell =~ s/\)/\\)/g;
#             $df_shell =~ s/\./\\./g;
#             $df_replace = "$df$i.$last_sh$fin_term$seniority$end";
#           } elsif ( ($sh =~ /$df$i\t([^\t]+)$/) && ($i > $occup) ) {
#             # Incomplete d or f shell which is the last open shell, without intermediate term
#             $occup = $i;
#             $df_shell = "$df$i\t$1";
#             my $fin_term = $1;
#             if ( $fin_term =~ /^(.+)([*]){0,1}$/ ) {
#               $fin_term = $1;
#               my $p = $2;
#               $df_replace = "$df$i\t$fin_term$seniority$p";
#             }
#           } elsif ( ($sh =~ /$df$i\.(\d+s2|\d+p6|\d+d10|\d+f14)(\S*)\t([^\t*]+)([*]{0,1})$/) && ($i > $occup) ) {
#             # Incomplete d or f shell followed by a complete shell, without intermediate term
#             $occup = $i;
#             $df_shell = "$df$i.$1$2\t$3$4";
#             my ($last_sh, $fin_term,$end) = ("$1$2",$3,$4);
#             $df_shell =~ s/\(/\\(/g;
#             $df_shell =~ s/\)/\\)/g;
#             $df_shell =~ s/\./\\./g;
#             $df_shell =~ s/\*/\\*/g;
#             $df_replace = "$df$i.$last_sh\t$fin_term$seniority$end";
#           }
#         }
#         if ( $occup > 0 ) {
#           $sh =~ s/$df_shell/$df_replace/;
#           $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'sh'} = $sh;
#         }
#       }
#     }
# } ##fix_trailing_index($$$$)

############################################################################
sub Create_IDEN_lev_map() {   #01/16/2024 8:13AM
############################################################################
  my $cpl = 'LS';

  # Re-create maps of energy values for energy sorting
  #%E_map = ();
  #%E_reverse_map = ();
  #for (my $parity = 1; $parity <= 2; $parity++) {
  #  my $i = 0;
  #  foreach my $J ( keys %{$energies[$parity-1]} ) {
  #    foreach my $num_e ( keys %{$energies[$parity-1]->{$J}}) {
  #      $i++;
  #      my ($E,$Eai,$AA_tot, $AR_tot, $DR) = @{$energies[$parity-1]->{$J}->{$num_e}};
  #      my $key = "$parity${delim}$E${delim}$i";
  #      $E_map{$key} = [$J,$num_e,$Eai,$AA_tot,$AR_tot];
  #      $E_reverse_map{$parity} = {} unless defined $E_reverse_map{$parity};
  #      $E_reverse_map{$parity}->{$J} = {} unless defined $E_reverse_map{$parity}->{$J};
  #      $E_reverse_map{$parity}->{$J}->{$num_e} = [$key];
  #    }
  #  }
  #}
  $reverse_order = 1; # To enforce reversed energy sorting regardless of parity
  my $i = 0;
  foreach my $key (sort {&sort_by_par_E($a,$b)} keys %E_map) {
    my ($parity,$E,$k) = split(/${delim}/, $key);
    #my $par_code = $parities[$parity-1];
    $i++;
    #my $E = $energies[$parity-1]->{$J}->{$num_e};
    my ($J,$num_e,$Eai,$AA_tot,$AR_tot) = @{$E_map{$key}};
    $E_map_ENLEV{$key} = $i;
    $E_reverse_map_ENLEV{$parity} = {} unless defined $E_reverse_map_ENLEV{$parity};
    $E_reverse_map_ENLEV{$parity}->{$J} = {} unless defined $E_reverse_map_ENLEV{$parity}->{$J};
    $E_reverse_map_ENLEV{$parity}->{$J}->{$num_e} = [$key,$i];
  }
  if (($i > 9999) && !$truncate) {
    print("\nNumber of levels $i is greater than maximum number of levels handled by IDEN2.\nAttempting to truncate the list...\n");
    $truncate = 1;
  }
  #open(ENLEV, ">ENLEV.DAT") or die ("Could not create ENLEV.DAT");
  foreach my $key (sort {$E_map_ENLEV{$a} <=> $E_map_ENLEV{$b}} keys %E_map_ENLEV) {
    my $lev_num_ENLEV = $E_map_ENLEV{$key};
    my ($parity,$E,$k) = split(/${delim}/, $key);
    my ($J,$num_e,$Eai,$AA_tot,$AR_tot) = @{$E_map{$key}};
    my $par_code = $parities[$parity-1];
    my ($J,$num_e,$Eai,$AA_tot,$AR_tot) = @{$E_map{$key}};

    my ($num_RCE_E,$exp_c,$Ee, $Ec) = (0, '*',$E,$E);

    # Substitute Eexp and Ec from RCE if available
    if ( !$no_RCE ) {
      $num_RCE_E = $map_RCG_RCE[$parity-1]->{$J}->{$num_e};
      my $hash = $RCE_lev[$parity-1]->{$J}->[$num_RCE_E-1];
      $exp_c = $hash->{'exp_c'};
      $Ec = $hash->{'Ec'};
      $Ee = ($exp_c) ? $hash->{'Ec'} : $hash->{'Ee'};
    }
    my $IDEN_label;
    foreach my $num_bas (sort {abs($vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$b})<=>
                              abs($vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$a})}
                          keys %{$vectors[$parity-1]->{$cpl}->{$J}->{$num_e}}) {
      my $term   = $basis[$parity-1]->{$cpl}->{$J}->{$num_bas}->{'label'};
      $term =~ s/[ ()>\]\}]//g;
      my $c_num  = $basis[$parity-1]->{$cpl}->{$J}->{$num_bas}->{'cn'};
      my $conf   = $confs[$parity-1]->{$c_num};
      my $p = ($par_code eq 'e' ? '_' : '~');
      $IDEN_label = sprintf("%-5s%1s%-5s", $conf, $p, $term);
      last;
    }

    $IDEN_lev_map{$lev_num_ENLEV} = [$parity,$E,$num_e,$J,$IDEN_label,$Ec,$Ee,$exp_c];

    # print ENLEV
    #   sprintf("%4d%12.3f%10.3f%12.3f%2s%12.3f%5.1f /%11s/\n", $lev_num_ENLEV,
    #     $Ec*1000, ($exp_c ? 5000 : 10), $Ee*1000, ($exp_c ? ' ' : '*'),
    #     ($Ee-$Ec)*1000, $J, $IDEN_label);
  }

  # print "Done.\n";
  # close ENLEV;
} ##Create_IDEN_lev_map()

############################################################################
sub PrintENLEV($) {   #01/16/2024 8:05AM
############################################################################
  my $force_print = shift;
  if (!$force_print) {
    print "Printing ENLEV.DAT...";
    &Create_IDEN_lev_map();
  }
  if (!$truncate || $force_print) {
    open(ENLEV, ">ENLEV.DAT") or die ("Could not create ENLEV.DAT");
    foreach my $lev_num_ENLEV (sort {$a<=>$b} keys %IDEN_lev_map) {
      my ($parity,$E,$num_e,$J,$IDEN_label,$Ec,$Ee,$exp_c) = @{$IDEN_lev_map{$lev_num_ENLEV}};

      print ENLEV
        sprintf("%4d%12.3f%10.3f%12.3f%2s%12.3f%5.1f /%11s/\n", $lev_num_ENLEV,
          $Ec*1000, ($exp_c ? 5000 : 10), $Ee*1000, ($exp_c ? ' ' : '*'),
          ($Ee-$Ec)*1000, $J, $IDEN_label);
    }
    print "Done printing ENLEV.DAT.\n";
    close ENLEV;
  }
} ##PrintENLEV()


############################################################################
sub ReadLines() {  #02/27/2013 9:56AM
############################################################################
  open(DLV, "<dlv.dat") or die "Could not open DLV.DAT for reading";
  my $i = 0;
  my $s;
  while ( defined($s = <DLV>) ) {
    $i++;
#     1011.9795  /          /       0.0015
#  25    486071.610      205.7310  /          /       0.007      1
    chomp($s);
    if ( $s =~ /^(.{5})(.{14})(.{34})(.{1,7})(.*)$/ ) {
      my ($intens,$wn,$mid,$wl_unc,$tail) = ($1,$2,$3,$4,$5);
      $intens =~ s/^\s+|\s+$//g;
      $wn =~ s/^\s+|\s+$//g;
      $wl_unc =~ s/^\s+|\s+$//g;
      $intens += 0;
      $wn += 0;
      $wl_unc += 0;
      $wl_unc = abs($wl_unc);
      $I_obs_max = abs($intens) if $I_obs_max < abs($intens);  # Ignore the minus sign in intensities
      if ( $wn <= 0 ) {
        die("Error in wavenumber in DLV.DAT, line $i");
      }
      push(@lines,[abs($intens),$wn,$wl_unc,[],$mid,$tail]);
    } else {
      die("Format error in DLV.DAT, line $i");
    }
  }
  close DLV;
  if ($convert2log) {
    &WriteDLV();
  }
} ##ReadLines()

############################################################################
sub WriteDLV() {  #01/16/2024 1:25PM
############################################################################
  my $Im_log_sc = 10*log($I_obs_max);
  my $dI = 999 - $Im_log_sc;
  die "Maximum intensity in $I_obs_max in DLV.DAT is too large for conversion to log scale." if $dI < 0;
  `copy dlv.dat dlv.dat.bak`;   # Back up the original dlv.dat
  open(DLV, ">dlv.dat") or die "Could not open DLV.DAT for writing";
  for (my $lin_num = 1; $lin_num <= $#lines + 1; $lin_num++) {
    my ($intens, $wn,$wl_unc,$trans,$mid,$tail) = @{$lines[$lin_num-1]};
    my $I = 10*log($intens);
    $lines[$lin_num-1]->[0] = $I;
    my $wl_unc_str = ' ';
    if ($wl_unc > 0) {
      $wl_unc_str = sprintf("%7.4f", $wl_unc);
    } else {
      $wl_unc_str = ' ';
    }
    print DLV sprintf("%5d%14.3f%34s%7s%s\n", $I, $wn,$mid,$wl_unc_str,$tail) or die "Error writing new DLV.DAT.";
  }
  close(DLV) or die "Could not create a new DLV.DAT.";
} ##WriteDLV()

############################################################################
sub ReadNumset() {   #02/27/2013 10:46AM
############################################################################
  open(NUMSET, "<numset.dat") or die "Could not open NUMSET.DAT for reading";
  my $i = 0;
  my $s;
  while ( defined($s = <NUMSET>) ) {
    $i++;
    next if ( $s =~ /^[\$]/);
    #if ( $s =~ /^.(.{5})(.{10})(.{7})/ ) {
      #my ($group_start, $uncert) = ($1,$3);
      #$group_start =~ s/^\s+|\s+$//g;
      #$uncert =~ s/^\s+|\s+$//g;
      #$group_start += 0;
      #$group_start = 1 unless $group_start;
      #$uncert += 0;
      #if ( ($group_start <= 0) || ($group_start > $#lines + 1) ) {
      #  die "Group start number out of range in NUMSET.DAT, line $i";
      #}
    if ( $s =~ /^([ a-z]{5})(.{10})(.{7})/ ) {
      my ($group_lbl, $wl_start, $uncert) = ($1,$2,$3);
      $group_lbl =~ s/^\s+|\s+$//;
      $wl_start =~ s/^\s+|\s+$//;
      $uncert =~ s/^\s+|\s+$//;
      $wl_start += 0;
      $uncert += 0;
      if ( $wl_start <=0 ) {
        die "Error in starting wavelength value in NUMSET.DAT, line $i (must be > 0)";
      }
      if ( $uncert <=0 ) {
        die "Error in group uncertainty value in NUMSET.DAT, line $i (must be > 0)";
      }
      #push(@numset,[$group_start, $uncert]);
      push(@numset,[$group_lbl, $wl_start, $uncert]);
    } else {
      die "Format error in NUMSET.DAT, line $i";
    }
  }
  # $i = $i;
  close(NUMSET);
} ##ReadNumset()

############################################################################
sub AssignLines() {  #02/27/2013 2:12PM
############################################################################
  foreach my $wl (sort keys %lines) {
    my ($par, $J, $n_lev1, $lft_ns, $second_par, $J2, $n_lev2, $gA, $cf, $I,$AA_tot, $AR_tot) = @{$lines{$wl}};
    my ($E,$E1,$AA_tot, $AR_tot) = @{$energies[$par-1]->{$J}->{$n_lev1}};
    my $cpl = "LS";
    #my $shells = &get_leading_term($par,$cpl,$J,$n_lev1); # Take the leading term as the level designation
    my ($E2,$E21,$AA_tot2, $AR_tot2) = @{$energies[$second_par-1]->{$J2}->{$n_lev2}};
    my ($Ee1, $Ee2, $Ec1, $Ec2, $ec1, $ec2, $lande1, $lande2) = ($E,$E2,$E,$E2,' ',' ','','');

    # Substitute Eexp from RCE if available
    if ( !$no_RCE ) {
      my $num_RCE_E = $map_RCG_RCE[$par-1]->{$J}->{$n_lev1};
      my $hash = $RCE_lev[$par-1]->{$J}->[$num_RCE_E-1];
      $ec1 = $hash->{'exp_c'};
      $Ec1 = $hash->{'Ec'};
      $Ee1 = ($ec1) ? $hash->{'Ec'} : $hash->{'Ee'};
      $lande1 = $hash->{'lande'};

      $num_RCE_E = $map_RCG_RCE[$second_par-1]->{$J2}->{$n_lev2};
      $hash = $RCE_lev[$second_par-1]->{$J2}->[$num_RCE_E-1];
      $ec2 = $hash->{'exp_c'};
      $Ec2 = $hash->{'Ec'};
      $Ee2 = ($ec2) ? $hash->{'Ec'} : $hash->{'Ee'};
      $lande2 = $hash->{'lande'};
    }

    next unless ( !$ec1 && !$ec2);

    my $dE = abs($E - $E2);
    my $dEexp = abs($Ee1 - $Ee2);
    next if ( ($dE == 0) || ($dEexp == 0) );

    my ($wl_c,$wl_e);
    if ( ($dEexp > 5) && ($dEexp < 50) ) {
      $wl_c = &Lair(1e5/$dE);
      $wl_e = &Lair(1e5/$dEexp);
    } else {
      $wl_c = 1e5/$dE;
      $wl_e = 1e5/$dEexp;
    }
    $dEexp = $dEexp*1000;
    my ($ns_unc,$ns_start,$ns_start_wn) = (0,0,0);
    if ( abs($dEexp-300224.1)<0.1  ) {
      $ns_unc = $ns_unc;
    }
    foreach my $NS (@numset) {
      my ($start_num,$unc) = @{$NS};
      my ($start_wn, $start_wn_unc) = &get_wn_unc($start_num,$unc);
      if ( $dEexp <= $start_wn + $start_wn_unc) {
        $ns_unc = $unc;
        $ns_start = $start_num;
        $ns_start_wn = $start_wn;
        last;
      }
    }
    my $lin_num = ($ns_start ? &find_line($dEexp,$ns_start,$ns_unc) : 0);
    if ( $lin_num ) {
      my @assigned_trans = @{$lines[$lin_num-1]->[3]};
      push(@assigned_trans,$wl);
      # Sort assigned transitions in the order of decreasing calculated intensity
      @assigned_trans = sort { $lines{$b}->[9] <=> $lines{$a}->[9] } @assigned_trans;
      my @at1 = ();
      my $I0 = 0;
      # Trancate the assigned transitions array by deleting
      # transitions with calculated intensity less than $Athresh times the max intensity
      foreach my $at (@assigned_trans) {
        my $I = $lines{$at}->[9];
        $I0 = $I unless $I0;
        last if $I/$I0 < $Athresh;
        push(@at1, $at);
      }
      $lines[$lin_num-1]->[3] = \@at1;
    }
  }

  for (my $lin_num = 1; $lin_num <= $#lines + 1; $lin_num++) {
    my @assigned_trans = @{$lines[$lin_num-1]->[3]};
    foreach my $wl (@assigned_trans) {
      push(@{$lines{$wl}}, $lin_num);
    }
  }
} ##AssignLines()

############################################################################
sub get_wn_unc($$) {   #02/27/2013 3:51PM
############################################################################
  my ($line_num,$default_unc) = @_;
  my $line_wn = $lines[$line_num-1]->[1];
  my $line_wl_unc = $lines[$line_num-1]->[2];
  my $wl_unc = $default_unc;
  $wl_unc = $line_wl_unc if $line_wl_unc;
  if ( !$wl_unc ) {
    die 'No default uncertainty in NUMSET.DAT and empty uncertainty in DLV.DAT';
  }
  my $wn_unc = $wl_unc/(1e8/$line_wn)*$line_wn;
  return ($line_wn,$wn_unc);
} ##get_wn_unc($$)

sub ConvertNumset() {
  for (my $ns = 0; $ns <= $#numset; $ns++) {
    my ($range_lbl, $start_wl, $default_unc) = @{$numset[$ns]};
    my $start_wn = 1e8/((($start_wl > 2000) && ($start_wl < 20000)) ? Lvac($start_wl) : $start_wl);
    my $start_num = 0;
    for (my $lin_num = 1; $lin_num <= $#lines + 1; $lin_num++) {
      my $line_wn = $lines[$lin_num-1]->[1];
      my $line_wl_unc = $lines[$lin_num-1]->[2];
      my $dE0 = abs($line_wn-$start_wn);
      if ( $line_wn > $start_wn) {
        my $line_wl_vac = 1e8/$line_wn;
        my $start_wn_unc = $line_wn*(($line_wl_unc > $default_unc) ? $line_wl_unc : $default_unc)/$line_wl_vac;
        if ( !$start_wn_unc ) {
          die 'No default uncertainty in NUMSET.DAT (range $range_lbl) and empty uncertainty in DLV.DAT (line no. $lin_num)';
        }
        if ( $dE0 <= $start_wn_unc) {
          $start_num = $lin_num;
          last;
        } 
      } else {
        $start_num = $lin_num;
        last;
      }
    }
    $start_num = $#lines + 1 unless $start_num;
    my $line_wl_unc = $lines[$start_num-1]->[2];
    if (!$line_wl_unc) {
      $lines[$start_num-1]->[2] = $default_unc;
    } else {
      $default_unc = $line_wl_unc;
    }
    $start_wn = $lines[$start_num-1]->[1];
    $numset[$ns] = [$start_num, $default_unc];
  }
}

############################################################################
sub find_line($$$) {  #02/27/2013 2:56PM
############################################################################
  my ($wn, $start_num, $default_unc) = @_;
  $start_num = 1 unless $start_num;
  my ($start_wn, $start_wn_unc) = &get_wn_unc($start_num,$default_unc);
  my $dE0 = abs($wn-$start_wn);
  my $end_num = $#lines + 1;
  if ( $wn > $start_wn) {
    if ( $dE0 <= $start_wn_unc) {
      return $start_num;
    } else {
      return 0;
    }
  } else {
    my $dn = $end_num - $start_num;
    while ( $dn > 0 ) {
      { use integer;
        $dn = ($end_num - $start_num)/2;
      }
      last unless $dn;
      my $n1 = $start_num + $dn;
      my ($wn1, $wn_unc1) = &get_wn_unc($n1,$default_unc);
      if ( $wn1 < $wn ) {
        $end_num = $n1;
      } else {
        $start_num = $n1;
        my ($start_wn, $start_wn_unc) = &get_wn_unc($start_num,$default_unc);
        $dE0 = abs($wn-$start_wn);
        if ( $dE0 <= $start_wn_unc) {
          return $start_num;
        }
      }
    }
    ($start_wn, $start_wn_unc) = &get_wn_unc($start_num,$default_unc);
    $dE0 = abs($wn - $start_wn);
    my $n = $start_num;
    my $wn_unc = $start_wn_unc;
    my ($wn1, $wn_unc1) = &get_wn_unc($end_num,$default_unc);
    my $dE1 = abs($wn - $wn1);
    if ( ($dE1 <= $wn_unc1) && (($dE1 < $dE0) || ($dE0 > $start_wn_unc)) ) {
      $n = $end_num;
      $dE0 = $dE1;
      $wn_unc = $wn_unc1;
    }
    $n = 0 if ( $dE0 > $wn_unc);
    return $n;
  }
} ##find_line($$$)

############################################################################
sub Create_TRANSDAT_map() {  #01/16/2024  9:22AM
############################################################################
  # open(TRANSDAT, ">TRANS.DAT") or die "Could not create TRANS.DAT";
  # print('Writing trans. probabilities to TRANS.DAT ... ');

  &AssignLines();
  
  my $Im_log_sc = 10*log($scale_int);  # Perl's log is natural logarithm
  my $dI = 0;
  if (!$convert2log) {
    $dI = $I_obs_max - $Im_log_sc;
  } else {
    $dI = 10*log($I_obs_max) - $Im_log_sc;
  }
  
  if ((($I_obs_max > 99999) && $linear_intensity_scale) || ($dI > 99999 && !$linear_intensity_scale)) {
    print("\nMaximum intensity in dlv.dat, $I_obs_max, is too large for IDEN2. Cannot create TRANS.DAT.");
    exit(1);
  }
  
  my ($upp_lev_num, $upp_lev_E) = (0,0);
  foreach my $wl (sort {
      my ($par_a, $J_a, $n_lev1_a, $lft_ns_a, $second_par_a, $J2_a, $n_lev2_a, $gA_a, $cf_a, $I_a,$AA_tot_a, $AR_tot_a) = @{$lines{$a}};
      my $IDEN_lev_num1_a = $E_reverse_map_ENLEV{$par_a}->{$J_a}->{$n_lev1_a}->[1];
      my $IDEN_lev_num2_a = $E_reverse_map_ENLEV{$second_par_a}->{$J2_a}->{$n_lev2_a}->[1];
      my ($par_b, $J_b, $n_lev1_b, $lft_ns_b, $second_par_b, $J2_b, $n_lev2_b, $gA_b, $cf_b, $I_b,$AA_tot_b, $AR_tot_b) = @{$lines{$b}};
      my $IDEN_lev_num1_b = $E_reverse_map_ENLEV{$par_b}->{$J_b}->{$n_lev1_b}->[1];
      my $IDEN_lev_num2_b = $E_reverse_map_ENLEV{$second_par_b}->{$J2_b}->{$n_lev2_b}->[1];
      $IDEN_lev_num1_a*1000 + $IDEN_lev_num2_a <=> $IDEN_lev_num1_b*1000 + $IDEN_lev_num2_b;
    } keys %lines) {
    my ($par, $J, $n_lev1, $lft_ns, $second_par, $J2, $n_lev2, $gA, $cf, $I,$AA_tot, $AR_tot, $lin_num) = @{$lines{$wl}};
    my ($E,$E1,$AA_tot, $AR_tot) = @{$energies[$par-1]->{$J}->{$n_lev1}};
    my $IDEN_lev_num1 = $E_reverse_map_ENLEV{$par}->{$J}->{$n_lev1}->[1];
    my $IDEN_lev_num2 = $E_reverse_map_ENLEV{$second_par}->{$J2}->{$n_lev2}->[1];
    my $cpl = "LS";
    #my $shells = &get_leading_term($par,$cpl,$J,$n_lev1); # Take the leading term as the level designation
    my ($E2,$E21,$AA_tot2, $AR_tot2) = @{$energies[$second_par-1]->{$J2}->{$n_lev2}};
    my ($Ee1, $Ee2, $Ec1, $Ec2, $ec1, $ec2) = ($E,$E2,$E,$E2,' ',' ');

    # Substitute Eexp from RCE if available
    if ( !$no_RCE ) {
      my $num_RCE_E = $map_RCG_RCE[$par-1]->{$J}->{$n_lev1};
      my $hash = $RCE_lev[$par-1]->{$J}->[$num_RCE_E-1];
      $ec1 = $hash->{'exp_c'};
      $Ec1 = $hash->{'Ec'};
      $Ee1 = ($ec1) ? $hash->{'Ec'} : $hash->{'Ee'};

      $num_RCE_E = $map_RCG_RCE[$second_par-1]->{$J2}->{$n_lev2};
      my $hash = $RCE_lev[$second_par-1]->{$J2}->[$num_RCE_E-1];
      $ec2 = $hash->{'exp_c'};
      $Ec2 = $hash->{'Ec'};
      $Ee2 = ($ec2) ? $hash->{'Ec'} : $hash->{'Ee'};
    }

    if ( $IDEN_lev_num1 != $upp_lev_num) {
      for ( my $i = $upp_lev_num + 1; $i <= $IDEN_lev_num1; $i++) {
        my ($parity,$E,$num_e,$J,$IDEN_label) = @{$IDEN_lev_map{$i}};
        my ($num_RCE_E,$exp_c,$Ee, $Ec) = (0, '*',$E,$E);
        # Substitute Eexp and Ec from RCE if available
        if ( !$no_RCE ) {
          $num_RCE_E = $map_RCG_RCE[$parity-1]->{$J}->{$num_e};
          my $hash = $RCE_lev[$parity-1]->{$J}->[$num_RCE_E-1];
          $exp_c = $hash->{'exp_c'};
          $Ec = $hash->{'Ec'};
          $Ee = ($exp_c) ? $hash->{'Ec'} : $hash->{'Ee'};
        }

        $upp_lev_E = $Ee;
        #writeln(f3,'$',i:4,'    J=',(PL^.g-1)/2:4:1,'  ',LName(PL),abs(PL^.E)*1000:13:3,' ',Estr[PL^.Fixed],' ');
#$   1    J= 3.5  f2   _3H2F   1247270.996
        # print TRANSDAT sprintf("\$%4d    J=%4.1f  %11s%13.3f%2s \n", $i, $J, $IDEN_label, abs($Ee)*1000, ($exp_c ? ' ' : '*'));
      }
      $upp_lev_num = $IDEN_lev_num1;
    }
    next if $I < 1;
    if (!$linear_intensity_scale) {
      $I = 10*log($I) + $dI;  # Perl's log is natural logarithm
    } else {
      $I *= $I_obs_max / $scale_int;
    }
    $lin_num += 0;
    my $wnCalc = abs($Ee1-$Ee2)*1000;
    my ($intens,$wn,$wl_unc) = $lin_num ? @{$lines[$lin_num-1]} : (0, 0, 0);
    my $dWn = 0;
    if ( $lin_num ) {
      $dWn = $wn - $wnCalc;
    }
    #writeln(f3,'+',ToLevel^.N:4,gA:5:0,abs(ToLevel^.E)*1000:12:3,
    #  ' ',Estr[ToLevel^.Fixed],WnCalc:13:3,LInt:6,WnObs:12:3,dWn:11:3,LNo:6);
#+ 103 2391  153863.007 *  1093407.990     0       0.000      0.000     0
    # print TRANSDAT sprintf("+%4d%5d%12.3f%2s%13.3f%6d%12.3f%11.3f%6d\n",
    #   $IDEN_lev_num2, $I, $Ee2*1000, ($ec2 ? ' ' : '*'), $wnCalc, $intens, $wn, $dWn, $lin_num);
    my $key = ("$IDEN_lev_num1.$IDEN_lev_num2") + 0;
    $TRANSDAT_map{$key} = [$IDEN_lev_num1, $IDEN_lev_num2, $I, $Ee2, $ec2, $wnCalc, $intens, $wn, $dWn, $lin_num];
  }
  # close TRANSDAT or die "Error closing TRANS.DAT file, possibly due to insufficient space on disk.";
  #print "Done.\n";
} ##Create_TRANSDAT_map()

sub TruncateIDEN() {
  my %IDEN_lev_map_new = ();
  my %TRANSDAT_map_new = ();
  my %lev_numbers_hash = ();
  foreach my $key (sort {$a<=>$b} keys %TRANSDAT_map) {
    my ($IDEN_lev_num1, $IDEN_lev_num2, $I, $Ee2, $ec2, $wnCalc, $intens, $wn, $dWn, $lin_num) = @{$TRANSDAT_map{$key}};
    $lev_numbers_hash{$IDEN_lev_num1} = 1;
    $lev_numbers_hash{$IDEN_lev_num2} = 1;
  }
  my $new_lev_no = 0;
  foreach my $old_no (sort {$a<=>$b} keys %lev_numbers_hash) {
    $new_lev_no++;
    $lev_numbers_hash{$old_no} = $new_lev_no;
    $IDEN_lev_map_new{$new_lev_no} = $IDEN_lev_map{$old_no};
  }
  if ($new_lev_no > 9999) {
    print("After truncation, the number of levels $new_lev_no is still greater than 9999. Aborting creation of files for IDEN.");
    exit(1);
  }
  foreach my $key (sort {$a<=>$b} keys %TRANSDAT_map) {
    my ($IDEN_lev_num1, $IDEN_lev_num2, $I, $Ee2, $ec2, $wnCalc, $intens, $wn, $dWn, $lin_num) = @{$TRANSDAT_map{$key}};
    my $IDEN_lev_num1_new = $lev_numbers_hash{$IDEN_lev_num1};
    my $IDEN_lev_num2_new = $lev_numbers_hash{$IDEN_lev_num2};
    my $key_new = "$IDEN_lev_num1_new.$IDEN_lev_num2_new" + 0;
    $TRANSDAT_map_new{$key} = [$IDEN_lev_num1_new, $IDEN_lev_num2_new, $I, $Ee2, $ec2, $wnCalc, $intens, $wn, $dWn, $lin_num];
  }
  %IDEN_lev_map = %IDEN_lev_map_new;
  %TRANSDAT_map = %TRANSDAT_map_new;
} ##TruncateIDEN()

############################################################################
sub WriteTRANS_DAT() {  #01/16/2024 8:25AM
############################################################################
  open(TRANSDAT, ">TRANS.DAT") or die "Could not create TRANS.DAT";
  print('Writing trans. probabilities to TRANS.DAT ... ');
  &Create_TRANSDAT_map();
  if ($truncate) {
    &TruncateIDEN();
    &PrintENLEV(1);
  }
  my ($upp_lev_num, $upp_lev_E) = (0,0);
  foreach my $key (sort {$a<=>$b} keys %TRANSDAT_map) {
    my ($IDEN_lev_num1, $IDEN_lev_num2, $I, $Ee2, $ec2, $wnCalc, $intens, $wn, $dWn, $lin_num) = @{$TRANSDAT_map{$key}};
    if ( $IDEN_lev_num1 != $upp_lev_num) {
      for ( my $i = $upp_lev_num + 1; $i <= $IDEN_lev_num1; $i++) {
        my ($parity,$E,$num_e,$J,$IDEN_label,$Ec,$Ee,$exp_c) = @{$IDEN_lev_map{$i}};
        #writeln(f3,'$',i:4,'    J=',(PL^.g-1)/2:4:1,'  ',LName(PL),abs(PL^.E)*1000:13:3,' ',Estr[PL^.Fixed],' ');
#$   1    J= 3.5  f2   _3H2F   1247270.996
        if (!$truncate || ($i == $IDEN_lev_num1)) {
          print TRANSDAT sprintf("\$%4d    J=%4.1f  %11s%13.3f%2s \n", $i, $J, $IDEN_label, abs($Ee)*1000, ($exp_c ? ' ' : '*'));
        }
      }
      $upp_lev_num = $IDEN_lev_num1;
    }
    next if $I < 1;
    #writeln(f3,'+',ToLevel^.N:4,gA:5:0,abs(ToLevel^.E)*1000:12:3,
    #  ' ',Estr[ToLevel^.Fixed],WnCalc:13:3,LInt:6,WnObs:12:3,dWn:11:3,LNo:6);
#+ 103 2391  153863.007 *  1093407.990     0       0.000      0.000     0
    print TRANSDAT sprintf("+%4d%5d%12.3f%2s%13.3f%6d%12.3f%11.3f%6d\n",
      $IDEN_lev_num2, $I, $Ee2*1000, ($ec2 ? ' ' : '*'), $wnCalc, $intens, $wn, $dWn, $lin_num);
  }
  close TRANSDAT or die "Error closing TRANS.DAT file, possibly due to insufficient space on disk.";
  print "Done.\n";
} ##WriteTRANSDAT()


############################################################################
sub complete_shells($) {    #04/17/2013 12:43PM
############################################################################
  my $sh = shift;
  my $n = 2*(2*$L_moment{$sh} + 1);
  return $n;
} ##complete_shells($) {

############################################################################
sub gffromS($$$) {    #07/14/2017 9:45AM
############################################################################
  my ($S,$WL,$type) = @_;
  my $gf = 0;
  if (($type eq '') || ($type eq 'E1')) {
    $gf = $S * 303.75568885954 / $WL;
  } elsif ($type eq 'M1') {
    $gf = $S * 0.0040438504 / $WL;
  } elsif ($type eq 'E2') {
    $gf = $S * 167.84224 / ($WL*$WL*$WL);
  } elsif ($type eq 'M2') {
    $gf = $S * 0.002235255 / ($WL*$WL*$WL);
  } elsif ($type eq 'E3') {
    $gf = $S * 47.140897 / ($WL*$WL*$WL*$WL*$WL);
  } elsif ($type eq 'M3') {
    $gf = $S * 0.000627579 / ($WL*$WL*$WL*$WL*$WL);
  }
  return $gf;
} ##gffromS($$$)

sub SfromgA($$$) {    #07/14/2017 9:45AM
############################################################################
  my ($gA,$WL,$type) = @_;
  my $S = 0;
  if (($type eq '') || ($type eq 'E1')) {
    $S = $gA * $WL*$WL*$WL/2.0261269E+18;
  } elsif ($type eq 'M1') {
    $S = $gA * $WL*$WL*$WL/269735e8;
  } elsif ($type eq 'E2') {
    $S = $gA * $WL*$WL*$WL*$WL*$WL/1.11995E+18;
  } elsif ($type eq 'M2') {
    $S = $gA * $WL*$WL*$WL*$WL*$WL/14909714e6;
  } elsif ($type eq 'E3') {
    $S = $gA * $WL*$WL*$WL*$WL*$WL*$WL*$WL/3.14441E+17;
  } elsif ($type eq 'M3') {
    $S = $gA * $WL*$WL*$WL*$WL*$WL*$WL*$WL/4.1861e12;
  }
  return $S;
} ##SfromgA($$$) {

1;
