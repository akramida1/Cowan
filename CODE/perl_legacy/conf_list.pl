#!perl
use strict;

my %confs36 = ();
my @conf_keys = ({},{});
my %L_moment = ('s' => 0, 'p' => 1, 'd' => 2, 'f' => 3, 'g' => 4, 'h' => 5,
                'i' => 6, 'k' => 7, 'l' => 8, 'm' => 9, 'n' =>10, 'o' =>11);
&read_in36();

foreach my $conf (keys %confs36) {
  my $key = '';
  my $key1 = '';
  my $par = 1;
  my $last_shell = &get_last_shell($conf);
  my $n_max = 0;
  for ( my $n_shell = 0; $n_shell <= $last_shell; $n_shell++ ) {
    #next unless defined($confs36{$conf}->[$n_shell]);
    my ($n,$sh,$w) = @{$confs36{$conf}->[$n_shell]};
    #next unless $w;
    $key .= ($n * 10000) . ' ' . $L_moment{$sh} .' ' . (99 - $w);
    if ( ($sh =~ /[pfhkmo]/) && ($w % 2 > 0) ) {
      $par *= -1;
    }
    $key1 .= (1000 + $n);
    $n_max = $n if ($n > $n_max);
  }
  #$key = $key1 . ' ' . $key;
  $key = (1000 + $n_max) . ' ' . $key;
  $par = 1 - ($par + 1) / 2;
  $conf_keys[$par]->{$key} = $conf;
  #print "$conf\t$key\n";
}

# print sorted conf list for each parity
for ( my $par = 0; $par <= 1; $par++ ) {
  my $s = '';
  foreach my $key (sort {$a cmp $b} keys %{$conf_keys[$par]}) {
    my $conf = $conf_keys[$par]->{$key};
    $s .= ', ' if ($s);
    my $w_prev = 0;
    for ( my $n_shell = 0; $n_shell <=7; $n_shell++ ) {
      next unless defined($confs36{$conf}->[$n_shell]);
      my ($n,$sh,$w) = @{$confs36{$conf}->[$n_shell]};
      if ( $w ) {
        $s .= '.' unless (!$w_prev);
        $w_prev = $w;
        if ( $w == 1 ) {
          $w = '';
        }
        $s .= $n . $sh . $w;
      }
    }
  }
  print "$s\n";
}

############################################################################
sub read_in36		#3/16/2004 2:51PM A.Kramida
############################################################################
{
	open IN36, "<IN36" or die "Could not open input file IN36";
	my $s = '';
	my $parity = -1;

	$s = <IN36>; # Skip the first line

	# Read configurations
  while ( (defined($s = <IN36>)) && ($s =~ /^.{16}(.{9}).{4} *(\d+)([^ 0-9])(\d*) *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1}/) ) {
		my $conf = $1;
		my ($n1,$s1,$w1,$n2,$s2,$w2,$n3,$s3,$w3,$n4,$s4,$w4,$n5,$s5,$w5,$n6,$s6,$w6,$n7,$s7,$w7,$n8,$s8,$w8,$n9,$s9,$w9) =
		  ($2,$3,$4, $6,$7,$8, $10,$11,$12, $14,$15,$16, $18,$19,$20, $22,$23,$24, $26,$27,$28, $30,$31,$32, $34,$35,$36);
  	$conf =~ s/^\s+|\s+$//g;

    if ( defined($confs36{$conf}) ) {
      die "Duplicate configuration name $conf in IN36";
    }
    $confs36{$conf} = [];
		for ( my $n_shell = 1; $n_shell <=9; $n_shell++ ) {
			my ($n,$sh,$w) = (0,'',0);
			eval("\$w = \$w$n_shell");
			eval("\$n = \$n$n_shell");
			eval("\$sh = \$s$n_shell");
			if ( ($w eq '') && ($n > 0) && $sh ) {
				$w = 1;
			}
      #if ( $w > 0 ) {
        #$s = lc($s);
        $sh = lc($sh);
				push(@{$confs36{$conf}},[$n,$sh,$w]);
      #}
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

	close IN36;
	return 1;
}	##read_in36

############################################################################
sub get_last_shell($)		#3/16/2004 2:04PM A.Kramida
############################################################################
{
  my $conf = shift;
  my $last_shell = 0;
  for ( my $ns = 7; $ns>=0; $ns-- ) {
    next unless defined($confs36{$conf}->[$ns]);
    my ($n, $sh, $occup) = @{$confs36{$conf}->[$ns]};
    if ( $occup > 0 ) {
      $last_shell = $ns;
      last;
		}
	}
	return $last_shell;
}	##get_last_shell