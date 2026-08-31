#!perl
use strict;
# This program substitutes fitted parameters from RCEOUT into the RCG input file ING11}
# Syntax: perl s11.pl [2] ------------------------------------------------
# If the optional parameter 2 is omitted, fitted parameters are substituted only for the first-parity
# configuration set, otherwise the substitution is done for both parities.

my $param = shift;

open(ING11,'<ing11') or die('Cannot open ING11 file.');
open(RCEOUT,'<RCEOUT') or die('Cannot open RCEOUT file.');
open(OUT,'>tmp') or die('Cannot create temporary output file.');

my $line_ING11 = 0;
my $line_RCEOUT = 0;
my $s1 = '';
my $start;
my $Eav_shift = 0;
do {
  $s1 = &readING11();
  if ( !$start && ($s1 =~ /^([ 0]{5}).{15}([ .e+\d-]){1,10}/i) && ($1 == 0) && $2) {
    $Eav_shift = $2;
  }
  if ($s1 =~ /^.{70}H/i)  {
    $start = 1;
  } else {
    print OUT $s1,"\n";
  }
} while ( !eof(ING11) && !$start );
        # Read fitted parameters for configurations of 1st parity
$s1 = &readFittedPars($s1,$Eav_shift);
if ($param =~ /2/) {
        # Read fitted parameters for configurations of 2nd parity
  $s1 = &readFittedPars($s1,$Eav_shift);
}

do {
  print OUT $s1,"\n";
  if (!eof(ING11)) {
    $s1 = &readING11();
  } else {
    $s1 = '';
  }
} while ( !eof(ING11) && $s1 );
print OUT $s1,"\n" if $s1;
close(ING11);
close(RCEOUT);
close(OUT) or die('Error closing temporary file.');
if ( -f 'ing11.bak' ) {
  `del ing11.bak`;
}
`ren ING11 ing11.bak`;
`ren tmp ing11`;
print "Substitution done.\n";

sub readING11() {
  my $s = <ING11> or die('Error reading ING11 file.');
  chomp $s;
  #if ( $s =~ /d56p     -d5p3/ ) {
  #  $s = $s;
  #}
  $line_ING11++;
  return $s;
}

sub readRCEOUT() {
  my $s = <RCEOUT> or die('Error reading RCEOUT file.');
  chomp $s;
  $line_RCEOUT++;
  return $s;
}

sub readFittedPars($$) {
  my $s1 = shift;
  my $Eav_shift = shift;

  my $s2 = '';
  do {
    $s2 = readRCEOUT();
  } while ( defined($s2) && ($s2 !~ /^ PARAMETER FLAG/i) );
  #$s2 = readRCEOUT();
  #die('Error in files ING11 or PARVALS.') if (eof(ING11) || eof(PARVALS));
  my ($num_conf,$CI_count,$stop,$skip_first);
  do {
    while ( $s1 && (($s1 !~ /^.{70}H/i) || ($s1 =~ /\/R2\//)) ) {
      print OUT $s1,"\n" ;
      $s1 = readING11();
    }
    my $CI = ( $s1 =~ /^.{9}-.{8}/ );
    # Configuration or CI block started
    my @params = ();
    my ($spectrum,$cur_conf,$cur_conf_trimmed,$num_params,$tail);
    my ($cur_conf1,$cur_conf2);
    if ( $s1 =~ /^(.{6})(.{12})(.{2})(.{9})(.)(.{9})(.)(.{9})(.)(.{9})(.)(.{9})(.)(.+)$/ ) {
      my ($sp,$conf,$np1,$p1,$t1,$p2,$t2,$p3,$t3,$p4,$t4,$p5,$t5,$tt) = ($1,$2,$3, $4,$5, $6,$7, $8,$9, $10,$11, $12,$13, $14);
      foreach ($np1,$p1,$t1,$p2,$t2,$p3,$t3,$p4,$t4,$p5,$t5) {
        s/^\s+|\s+$//g;
      }
      $num_params = $np1;
      if ( $CI ) {
        $s1 =~ /^(.{9})-(.{8})/;
        ($cur_conf1,$cur_conf2) = ($1,$2);
        $tail = $tt;
      } else {
        ($spectrum,$cur_conf,$tail) = ($sp,$conf,$tt) ;
        $cur_conf_trimmed = $cur_conf;
        $cur_conf_trimmed =~ s/^\s+|\s+$//g;
      }
      #@params = ([$p1,$t1],[$p2,$t2],[$p3,$t3],[$p4,$t4],[$p5,$t5]);
      $np1 = $num_params;
      $np1 = 5 if $num_params > 5;
      for ( my $i = 1; $i <= $np1; $i++ ) {
        my $cmd = 'push(@params,[$p' . $i .',$t' . $i .'])';
        eval($cmd);
      }
      if ( $num_params > 5 ) {
        my $need_params = $num_params - 5;
        while ( $need_params > 0 ) {
          $s1 = readING11();
          if ( $s1 =~ /^(.{10})(.{10}){0,1}(.{10}){0,1}(.{10}){0,1}(.{10}){0,1}(.{10}){0,1}(.{10}){0,1}/ ) {
            my @params1 = ($1,$2,$3,$4,$5,$6,$7);
            foreach (@params1) {
              s/^\s+|\s+$//g;
            }
            my $np = $need_params;
            $np = 7 if $np > 7;
            for ( my $i = 1; $i <= $np; $i++ ) {
              $params1[$i-1] =~ /^(.+)(.)$/;
              my ($p,$t) = ($1, $2);
              die ("Parameter format error in ING11, line $line_ING11.") unless (defined($p) && defined($t));
              push(@params,[$p,$t]);
            }
            $need_params -= $np;
          } else {
            die ("Parameter format error in ING11, line $line_ING11.");
          }
        }
      #} else {
      #  $s1 = readING11();
      }
    } else {
      die ("Parameter format error in ING11, line $line_ING11.");
    }

    my $np = $#params + 1;
    my $np1 = $np;
    $np1 = 5 if $np1 < 5;
    my @params_RCE = ();
    my $i1 = 0;
    for ( my $i = 1; $i <= $np1; $i++ ) {
      $s2 = readRCEOUT() unless ( (($i==1) && $skip_first) || ($i > $np) );
      $i1++;
      my ($par_name,$par_value,$conf1,$conf2);
      if ( $i <= $np ) {
        if ( $CI ) {
          if ( $s2 =~ /^(.{11})(.{4})(.{14})(.{28})(.{6})-(.{0,6})/ ) {
            my ($pn,$flag,$pv,$c1,$c2) = ($1,$2,$3,$5,$6);
            foreach ($pn,$flag,$pv,$c1,$c2) {
              s/^\s+|\s+$//g;
            }
            ($par_name,$par_value,$conf1,$conf2) = ($pn,$pv,$c1,$c2);
          } else {
            die("RCEOUT CI parameter format mismatch for configurations $cur_conf1,$cur_conf2, RCEOUT line $line_RCEOUT.");
          }
        } else {
          if ( $s2 =~ /^(.{11})(.{4})(.{14})/ ) {
            my ($pn,$flag,$pv) = ($1,$2,$3);
            $pv =~ s/^\s+|\s+$//g;
            ($par_name,$par_value) = ($pn,$pv);
          } else {
            die("RCEOUT parameter format mismatch for configuration $cur_conf_trimmed, RCEOUT line $line_RCEOUT.");
          }
        }
      } else {
        $par_value = '0';
      }
      if ( $i == 1 ) {
        if ( $CI ) {
          my ($c1, $cc1, $c2, $cc2) = ($conf1, $cur_conf1, $conf2, $cur_conf2);
          foreach ($c1, $cc1, $c2, $cc2) {
            s/^\s+|\s+$//g;   # Trim leading and trailing spaces
          }
          if ( ($cc1 !~ /^$c1/) || ($cc2 !~ /^$c2/) ) {
            die("Configuration names mismatch in CI block, ING11 line $line_ING11: $cc1-$cc2, RCEOUT line $line_RCEOUT: $c1-$2.");
          }
          print OUT "$cur_conf1-$cur_conf2", sprintf("%2d", $np);
        } else {
          my $conf = '';
          if ($par_name =~ /^EAV (.{7})$/i) {
            $conf = $1;
          }
          my ($c, $cc) = ($conf, $cur_conf);
          foreach ($c, $cc) {
            s/^\s+|\s+$//g;   # Trim leading and trailing spaces
          }
          if ( !$c || ($cc !~ /^$c/) ) {
            die("Configuration name mismatch, ING11 line $line_ING11: $cc, RCEOUT line $line_RCEOUT: $c.");
          }
          $par_value -= $Eav_shift;
          print OUT $spectrum, $cur_conf, sprintf("%2d", $np);
        }
#          if (VP.GE.-1D-3 .and. VP.LE.1D-3 .and. I.GT.1) then
#            write(CVPAR(I),27) VP
#   27       FORMAT (1PE9.2)
#          else
#c           If parameter is too large, use floating-point format
#c           to avoid overflow
#            if (VP .GE. -9999.99D0 .AND. VP .LE. 99999.99D0) then
#              IVP=10000.D0*VP
#              write(CVPAR(I),2701) IVP
# 2701         FORMAT (I9)
#            else
#              if (i.eq.1) then
#                write(CVPAR(I),2702) VP
#              else
#                write(CVPAR(I),2703) VP
#              endif
# 2702         FORMAT (F9.1)
# 2703         FORMAT (F9.2)
#            endif
#          endif
      }

      my ($pv,$par_type);
      if ( $i <= $np ) {
        $par_type = $params[$i-1]->[1];
        if ( $CI ) {
          $pv = sprintf("%9.4f", $par_value);
        } else {
          if ( ($par_value >= -1e-3) && ($par_value <= 1e-3) && ($i > 1) ) {
            $pv = sprintf("%10.2e",$par_value);
            $pv =~ s/(e[+-])0/$1/;
          } elsif ( ($par_value >= -9999.99) && ($par_value <= 99999.99) ) {
            $pv = sprintf("%9.0f", $par_value*10000.0);
          } else {
            if ( $i == 1 ) {
              $pv = sprintf("%9.1f", $par_value);
            } else {
              $pv = sprintf("%9.2f", $par_value);
            }
          }
        }
      } else {
        if ( $CI ) {
          $par_type = 5;
          $pv = '   0.0000';
        } else {
          $par_type = 0;
          $pv = '        0';
        }
      }
      print OUT $pv,$par_type;
      if ( $i == 5 ) {
        $i1 = 0;
        print OUT "$tail\n";
      } elsif ( ($i1 == 7) || ($i == $np1) ) {
        $i1 = 0;
        print OUT "\n";
      }
    }
    $i1 = 0;
    $s1 = readING11();
    if ( !eof(RCEOUT) ) {
      $s2 = readRCEOUT();
      $skip_first = 1;
    } else {
      $stop=1;
    }
    if ( $s2 =~ /^    2    \d    \d    \d[ .\d-]{10}    \d    \d[ .\d-]{10}/ ) {
      $stop = 1;
    }
    #if ($s1 =~ /\/R2\//) {
    #  $s1 = $s1;
    #}
  } while ( ($s1 =~ /^.{70}H/i) && ($s1 !~ /\/R2\//) && !eof(ING11) && $s2 && !$stop);
  return $s1;
}

