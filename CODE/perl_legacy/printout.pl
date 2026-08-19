#!perl
use strict;
use vars qw{@parities @energies @map_RCG_RCE @RCE_lev @vectors @basis};
require 'conv_cowan.pl';

# Global parameters
my $num_printed_components = 5;          # Number of printed eigenvector components
my $min_printed_percentage = 1.5;        # Omit 3rd, 4th, 5th, etc. eigenvector components with percentage less than that
my $min_printed_second_percentage = 1.5; # Omit second eigenvector components with percentage less than that

# Read command-line parameters
my $out_file = shift;
my $no_RCE = shift;
if ( !$out_file ) {
  print "\nUsage: \nprintout <out_levels_file_name> [noRCE] [param_out_file_name]\n\n";
  exit;
}
my $param_file = '';
if ( $no_RCE !~ /(no)*RCE$/i ) {
  $param_file = $no_RCE;
  $no_RCE = shift;
} else {
  $param_file = shift;
}

$no_RCE = '' unless $no_RCE =~ /noRCE$/i;
$param_file = 'params.txt' unless $param_file; # By default, assume name 'params.txt' for the parameters output file

# Open input and output files
open OUTG11, "<OUTG11" or die "Could not open input file OUTG11";
open OUT_FILE, ">$out_file" or die "Could not open output file " . $out_file;
open PARAM_FILE, ">$param_file" or die "Could not open parameters output file " . $param_file;

# Initialize global variables
&init_vars();

# Read RCG options
&read_RCG_options();

&read_in36();

# Start processing ...
#my ($s, $J_prev);
my $s;
#my $start = 0;

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

close (OUTG11);

&Reorder_Params();
&ReadOUTE();
&ReadHF();
&write_params();

if ( !$no_RCE ) {
  $no_RCE = &ReadRCE();
  unless ($no_RCE) {
    &sort_vectors();
    &Identify_RCE_levs();
  }
}

# Print levels file header
print OUT_FILE "cpl\tpar\tEe\te_c\tEc\tEe-Ec\tJ\tLande_g\t\%1\tconf1\tterm1\t\%2\tconf2\tterm2\t\%3\tconf3\tterm3\t\%4\tconf4\tterm4\t\%5\tconf5\tterm5\n";
# Print the energy levels in LS and JJ coupling

foreach my $cpl  ('LS','JJ' ) {
  print "Printing $cpl levels...";
  # Create a map of energy values for energy sorting
  for ( my $parity = 1; $parity <=2; $parity++ ) {
    my %E_map = ();
    my $i = 0;
    foreach my $J ( keys %{$energies[$parity-1]} ) {
      foreach my $num_e ( keys %{$energies[$parity-1]->{$J}}) {
        $i++;
        my $E = $energies[$parity-1]->{$J}->{$num_e}->[0];
        my $key = "$E\t$i";
        $E_map{$key} = [$J,$num_e];
      }
    }
    my $par_code = $parities[$parity-1];

    foreach my $key (sort {my ($E1,$k1,$E2,$k2) = (split(/\t/,$a),split(/\t/,$b)); $E1<=>$E2}
                     keys %E_map) {
      my ($E,$k) = split(/\t/, $key);
      my ($J,$num_e) = @{$E_map{$key}};

      my ($num_RCE_E,$exp_c,$Ee,$Ec,$lande_g) = (0, '*',$E,$E,'');

      # Substitute Eexp from RCE if available
      if ( !$no_RCE ) {
        $num_RCE_E = $map_RCG_RCE[$parity-1]->{$J}->{$num_e};
        my $hash = $RCE_lev[$parity-1]->{$J}->[$num_RCE_E-1];
        $exp_c = $hash->{'exp_c'};
        $Ec = $hash->{'Ec'};
        $Ee = ($exp_c) ? $hash->{'Ec'}
                        : $hash->{'Ee'};
        $lande_g = $hash->{'lande'};
      }

      my $dif_E = $exp_c ? '' : sprintf("%12.4f",$Ee-$Ec);
      $dif_E =~ s/^\s+|\s+$//g;

      print OUT_FILE "$cpl\t$par_code\t" .
        sprintf("%12.6f\t%s\t%12.4f\t%s\t%4.1f\t%6.3f",$Ee,$exp_c,$Ec,$dif_E,$J,$lande_g);
      my $i = 0;
      my $need_CR = 1;
      if (($parity == 1) && ($Ec eq '14.0000')) {
        $i = $i;
      }
      foreach my $num_bas (sort {abs($vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$b})<=>
                                abs($vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$a})}
                            keys %{$vectors[$parity-1]->{$cpl}->{$J}->{$num_e}}) {
        $i++;
        if ( $i > $num_printed_components ) {
          print OUT_FILE "\n";
          $need_CR = 0;
          last;
        }
        my $A = $vectors[$parity-1]->{$cpl}->{$J}->{$num_e}->{$num_bas};
        $A *= $A*100;
        if ( ($A < $min_printed_percentage) && (($i > 2) || ($A < $min_printed_second_percentage)) ) {
          print OUT_FILE "\n";
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
		  $shells =~ s/[*]+/*/g;
        }

        print OUT_FILE "\t$A\t$shells";
      }
      if ( $need_CR ) {
        print OUT_FILE "\n";
      }
    }
  }
  print "Done.\n";
}
close OUT_FILE;

1;
