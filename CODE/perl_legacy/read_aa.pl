#!perl
use strict;

# Read parameters from the command line
my $cores = shift; # Comma-delimited list of cores with their energies relative to the ground state of ionized atom in units of KK (10^3 cm-1).
                   # Sample (target ion = He-like): 1s2:0,1s2s:120.1,1s2p:125.7
                   # Sample (target ion = N-like Ne IV): p3.+4S:0,p3.+2D:41.25,p3.+2P:62.44
                   # No spaces are allowed in this list.

my $IP = shift;    # IP in 10^3 cm-1

if ( !$cores && !$IP ) {
  print "\nUsage: \nread_aa <cores_list> <IP[10^3 cm-1]>\n";
  print "where <cores_list> is a comma-delimited list of core names with energies relative to the ground state\n";
  print "of the next ion in units of 10^3 cm-1. Core names may be Unix-type patterns, to sum over several terms. Examples:\n";
  print "For AI to a He-like ion: <cores_list> = 1s2:0,1s2s:120.1,1s2p:125.7\n";
  print "For AI to N-like Ne IV:  <cores_list> = p3.+4S:0,p3.+2D:41.25,p3.+2P:62.44\n";
  print "No spaces are allowed in this list of cores.\n\n";
  exit;
}

my %cores = ();
my $Ground_State_Core = '';
my @cores = split(/,/, $cores);
my $num_cores = $#cores;
for ( my $i = 0; $i <= $num_cores; $i++ ) {
  # In the parameters, the energies of the core states must be given
  # relative to the ground state of the ionized atom.
  # For example, if we are calculating AI of Ne VII, the core energies
  # must be given relative to the ground state of Ne VIII.
  my $core = $cores[$i];
  my ($core_name, $core_energy) = split(/:/, $core);
  $cores{$core_name} = $core_energy;
  if ( $core_energy == 0 ) {
    $Ground_State_Core = $core;
  }
}

foreach my $core (keys %cores) {
  # Recalculate the core energies relative to the ground state of the
  # atom for which we are calculating the AI rates.
  $cores{$core} += $IP;
}

$IP = -100000 unless defined $IP;

my $s="";
my $J_str = '';
my $J_prev = -1;
my %confs = ();

open OUTG11, "<OUTG11" or die "Cannot open OUTG11";

my %AA = ();
my %eigenvalues = ();
my %EK = ();

my $exit = 0;
my $confs_read = 1;
do {
  if ( !$confs_read ) {
    while ( (defined($s = <OUTG11>)) && ($s !~ /^.{98}(.{9}) +([0-9.-]+) +([0-9.]+)/) ) {
      next;
    }
    # Read conf. names and total energies
    my $prev_spectrum = '';
    my $E_min1 = 100000000;
    while ( (defined($s = <OUTG11>)) && ($s =~ /^.{92}(.{6})(.{9}) +([0-9.-]+) +([0-9.]+)/) ) {
      my ($spectrum, $conf, $Et, $Ek) = ($1,$2,$3,$4);
      if ( $prev_spectrum && ($spectrum ne $prev_spectrum) ) {
        $E_min1 = $Et;
      }
      $prev_spectrum = $spectrum;
      $conf =~ s/^\s|\s+//g;
      $confs{$conf} = [$Et, $Ek];
      next;
    }
    # Find the lowest config
    my $E_min = 100000000;
    my $lowest_conf = '';
    foreach my $conf (sort {$confs{$a}->[0]<=>$confs{$b}->[0]} keys %confs) {
      $E_min = $confs{conf}->[0];
      $lowest_conf = $conf;
      last;
    }
    # Find ionization energies
    foreach my $conf (keys %confs) {
      my $Et = $confs{$conf}->[0];
      if ( !$confs{$conf}->[1] ) {
        $confs{$conf}->[1] = $Et - $E_min - $IP;
      } else {
        $confs{$conf}->[1] = $E_min1 - $Et;
      }
    }
    $confs_read = 1;
  }

  while ( (defined($s = <OUTG11>)) && ($s !~ /EIGENVALUES      [()]J=([^()]+)[()]/) ) {
    next;
  }
  if ( $s =~ /EIGENVALUES      [()]J=([^()]+)[()]/ ) {
    $J_str = $1;
    $J_str =~ s/^\s+|\s+$//g;
    $exit = 1 if ($J_str <= $J_prev);
    $eigenvalues{$J_str} = [] unless $exit;
    $EK{$J_str} = {} unless $exit;
  }
  $J_prev = $J_str;
  if ( !$exit ) {
    if ( defined($s) ) {
      chomp $s;
      #print "$s\n";
    }

    do {
      # Read eigenvalues
      while ( (defined($s = <OUTG11>)) && ($s !~ /^.{28}( *[0-9-]+[.][0-9]+){1,11}/) ) {
        last if ($s =~ /AVEEIG|CONFIG/);
        next;
      }
      if ( defined($s) && ($s =~/^.{28}(.{9})(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}/) ){
        my @arr = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
        chomp $s;
        $s =~ s/^\s+|\s+$//g;
        #my @arr = split(/\s+/, $s);
        #push @{$eigenvalues{$J_str}}, @arr;
        for ( my $i = 1; $i<=11;$i++ ) {
          last unless defined $arr[$i-1];
          $arr[$i-1] =~ s/^\s+|\s+$//g;
          push @{$eigenvalues{$J_str}}, $arr[$i-1];
        }
        #print "$s\n";
      }
    } while ( defined($s) && ($s !~ /CONFIG/) );
    # Read kinetic energies
    while ( (defined($s = <OUTG11>)) && ($s !~ /KE FOR CORE/) ) {
      last if ($s =~ /EIGENVECTORS/);
      next;
    }
    while ( $s =~ /KE FOR CORE.{8}(.{17})(.{9})(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}/ ) {
      my ($core,@arr) = ($1,0,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
      $core =~ s/ +/ /g;
      $core =~ s/^\s+|\s+$//g;
      $EK{$J_str}->{$core} = [];
      for ( my $i = 1; $i<=11;$i++ ) {
        last unless defined $arr[$i];
        $arr[$i] =~ s/^\s+|\s+$//g;
        push @{$EK{$J_str}->{$core}}, $arr[$i];
      }
      while (defined($s = <OUTG11>) && ($s =~/^.{28}(.{9})(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}/) && ($s !~ /KE FOR CORE|EIGENVECTORS/) ){
        my @arr = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
        for ( my $i = 1; $i<=11;$i++ ) {
          last unless defined $arr[$i-1];
          $arr[$i-1] =~ s/^\s+|\s+$//g;
          push @{$EK{$J_str}->{$core}}, $arr[$i-1];
        }
      }
    }

    do {
      while ( (defined($s = <OUTG11>)) && ($s !~ /AA[()]UNITS/) ) {
        last if ($s =~ /AVEEIG/);
        next;
      }
      if ( defined($s) && ($s !~ /AVEEIG/) ) {
        chomp $s;
        #print "$s\n";
        my @arr = ();
        my $prev_target_state = '';
        while ( (defined($s = <OUTG11>)) && ($s =~ /[*]{2}AA[*]{2}(.{22})(.+)$/) ) {
          my $target_state = $1;
          $s = $2;
          $target_state =~ s/^\s+|\s+$//g;
          $target_state =~ s/\s+/ /g;
          $s =~ s/^\s+|\s+$//g;
          if ( $target_state ) {
            my $must_add = 0; # Some target states may have the same LS label in OUTG11
            if ( @arr and defined($arr[0]) and $prev_target_state) {
              if ( !defined($AA{$prev_target_state}->{$J_str}) ) {
                $AA{$prev_target_state} = {} unless defined $AA{$prev_target_state};
                $AA{$prev_target_state}->{$J_str} = [];
                $must_add = 0;
              } else {
                $must_add = 1;
              }
              if ( !$must_add ) {
                push @{$AA{$prev_target_state}->{$J_str}}, @arr;
              } else {
                # A target state with the same label was already filled in.
                # Add new AA values to the previuos ones.
                my $num_val = $#arr;
                for ( my $j = 0; $j <= $num_val; $j++ ) {
                  $AA{$prev_target_state}->{$J_str}->[$j] += $arr[$j];
                }
              }
            }
            @arr = ();
            $prev_target_state = $target_state;
          }
          my @arr1 = split(/\s+/, $s);
          push @arr, @arr1;

          #chomp $s;
          #print "$s\n";
        }
        if ( @arr and defined($arr[0]) and $prev_target_state) {
          $AA{$prev_target_state} = {} unless defined $AA{$prev_target_state};
          $AA{$prev_target_state}->{$J_str} = [] unless defined $AA{$prev_target_state}->{$J_str};
          push @{$AA{$prev_target_state}->{$J_str}}, @arr;
        }
        $prev_target_state = '';
      }
    } while ( defined($s) && ($s !~ /AVEEIG/) );
  }

} while ( defined($s) && ($s !~ /SPECTRUM/ && !$exit) );

my $s_out = "J\tEnergy";
foreach my $core (sort {$cores{$a}<=>$cores{$b}} keys %cores) {
  # Header line: print core names and energies, divided by colon
  $s_out .= "\t" . $core . ':' . $cores{$core};
}
$s_out .= "\tAA_tot";
print "$s_out\n";

# Find correction for the kinetic energies
#my $EK_min = 10000000;
#foreach my $J_str (sort {$a<=>$b} keys %EK) {
#  foreach my $core (keys %{$EK{$J_str}}) {
#    next unless ($core =~ /$Ground_State_Core/);
#    my @EK = @{$EK{$J_str}->{$core}};
#    my $num_eigenvalues = $#EK;
#    for (my $i = 0; $i<= $num_eigenvalues; $i++) {
#      if ( $EK[$i] < $EK_min ) {
#        $EK_min = $EK[$i];
#      }
#    }
#  }
#}
#my $dEkin = ($IP > 0) ? $IP + $EK_min : 0; # Note that $EK_min is negative

foreach my $J_str (sort {$a<=>$b} keys %eigenvalues) {
  my @eigenvalues = @{$eigenvalues{$J_str}};
  my $num_eigenvalues = $#eigenvalues;
  for (my $i = 0; $i<= $num_eigenvalues; $i++) {
    my $E = $eigenvalues[$i];
    next if $E <-5000;
    $s_out = "$J_str\t" . $E;
    my $AA_tot = 0;
    #if ( $E eq '7210.57' ) {
    #  $E = $E;
    #}
    foreach my $core (sort {$cores{$a}<=>$cores{$b}} keys %cores) {
      my $AI_rate = 0;
      my $IP_core = $cores{$core};
      foreach my $target (sort keys %AA) {
        next unless ($target =~ /$core/);
        my $AA = defined($AA{$target}->{$J_str}->[$i]) ?
           $AA{$target}->{$J_str}->[$i] :
           '';
        my @EK = defined($EK{$J_str}->{$target}) ? @{$EK{$J_str}->{$target}} : ();
        $AA = 0 if ($AA && ($E < $IP_core));
        #$AA = 0 if ( $AA && (!defined($EK[$i]) || ($EK[$i] < 0)) );
        $AI_rate += $AA;
      }
      $s_out .= sprintf("\t%9.2e",$AI_rate);
      $AA_tot += $AI_rate;
    }
    $s_out .= sprintf("\t%9.2e\n",$AA_tot);
    print $s_out;
  }
  next;
}

1;
