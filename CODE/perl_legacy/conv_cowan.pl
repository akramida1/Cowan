# include file for printout.pl and conv_out.pl
#
use strict;
use vars qw{@parities @energies @map_RCG_RCE @RCE_lev @vectors @basis %L_moment %conf_nums @confs};

my (%confs36,@shells,@params,@sdx,@param_flags,@param_group_counts,@Emax,
  @param_scaling,@max_cyc_no,@CI,@vectors_sorted, @vectors_index, @basis_hash,
  %shell_order_in36,%shell_order_in36_back,%shell_ord_seq,@last_shells,
  %no_genealogy_shells,@LS_templates,@JJ_templates,$lowest_complete_shell,
  @term_labels_RCG, @term_labels_map, @num_states_in_blocks, $num_blocks, %open_shells);
%open_shells = ();

############################################################################
sub fill_LS_shell(@) {  #9/27/2004 12:26PM A.Kramida
############################################################################
  my ($parity,$state_num,$J,$n_conf, @p) = @_;
  foreach (@p) {
    $_ =~ s/^\s+|\s+$//; # strip leading and trailing spaces
  }
  $n_conf += 0;

  $basis[$parity-1]->{'LS'}->{$J}->{$state_num} = {};
  my $L_format = $p[11]; # Will be non-empty if > 4 shells
  # Current and accumulated term symbols
  my @fill = ();
  for ( my $i = 0; $i <= 7; $i++ ) {
    $fill[$i] = [$p[2*$i],$p[2*$i+1]];
  }
  #if ( ($J eq '0.0') and ($state_num == 8) ) {
  #  $n_conf = $n_conf;
  #}

  my ($last_shell, $last_shell_reordered, $last_shell_summation_order) = @{$last_shells[$parity-1]->[$n_conf-1]};  # 0,1,...

  my $fill_str = '';
  my $acc_parity = 0;
  my $first_filled = -1;
  my $reordered = 0;
  my $final_term = '';
  my ($LS_curr,$LS_accum, $LS_curr_reordered, $LS_accum_reordered) = ('1S','1S','1S','1S');
  my $LS_last = '1S';

  # Determine the accumulated parity after adding each shell in the order of summation given in OUTG11
  my @accum_par = ();
  my $acc_par = 0;
  for ( my $ns = 0; $ns <= $last_shell; $ns++ ) {
    my ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    $acc_par++ if $sh_par_char;
    $accum_par[$ns] = ( $acc_par % 2 > 0 ) ? '*' : '';
  }

  my @shells_summation_order = ();
  #my @shells_display_order = ();
  my $shell_ind = 0;
  #my %shell_ord_sequence = ();
  my $prev_accum_par_char = '';
  my $accum_par_char = '';
  my $final_term = '';
  my $desig = '';
  my $desig_no_shell_num = '';
  my $desig_has_genealogy = 0;
  my $nL_prev = 0;
  my $last_shell_in_conf = 0;
  for ( my $ns = 0; $ns <= $last_shell; $ns++ ) {
    # Get the new sequential order of shells from the hash
    my $n_shell_display_order = $ns;
    if ( $shell_ord_seq{$ns+1} ) {
      $n_shell_display_order = $shell_ord_seq{$ns+1} - 1;
    #  $shell_ord_seq{$ns} = $n_shell_display_order;
    #} else {
    #  $n_shell_display_order = -1;
    }
    $reordered = 1 unless ($ns == $n_shell_display_order);

    ($LS_curr,$LS_accum) = @{$fill[$ns]};
    my ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$ns]};

    my $L = $sh;
    $L =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
    my $nL = $sh + 0; #$sh*10000 + $L_moment{$L};

    if ( !$occup || (($sh eq $lowest_complete_shell) && !$reordered) ) {
      $nL_prev = $nL;
      next;
    }

    $accum_par_char = (($prev_accum_par_char && !$sh_par_char) || (!$prev_accum_par_char && $sh_par_char)) ? '*' : '';

    my $last_shell_in_this_conf = ($reordered ? $last_shell_reordered : $last_shell_summation_order);
    my $complete_shell = ($occup == &complete_shells($L));
    if ( $complete_shell && ($ns != $last_shell_in_this_conf)) {
      $LS_curr = '';
      $LS_accum = '' unless ( $ns == $last_shell_summation_order );
      $last_shell_in_conf = $ns unless $last_shell_in_conf;
    } else {
      $LS_curr .= $sh_par_char;
      $LS_accum .= $accum_par_char;
      $last_shell_in_conf = $ns;
    }
    if ( $ns == $last_shell_summation_order ) {
      $final_term = $LS_accum;
      $LS_accum = '';
    }
    #$shell_ord_seq{$n_shell_display_order} = $shell_ind;
    $prev_accum_par_char = $accum_par_char;
    $occup = '' if $occup == 1;
    $desig .= ($desig ? '.' : '');
    my $added_desig = '';
    $desig .= $sh . $occup;
    $added_desig .= $sh . $occup;
    if ( !$complete_shell || ($nL > $nL_prev) ) {
      $desig_no_shell_num .= ($desig_no_shell_num ? '.' : '');
      $desig_no_shell_num .= $L . $occup;
    }

    if ( ($shell_ind || ($ns < $last_shell) ) &&
      (&next_shell_requires_genealogy($parity,$n_conf,$ns,$last_shell_summation_order,$desig_no_shell_num)
        && !defined($no_genealogy_shells{$desig_no_shell_num})
        || ($reordered && (($ns == $last_shell) || ($ns == $last_shell_summation_order)) && $desig_has_genealogy)
        || (($ns == $last_shell) && !defined($no_genealogy_shells{$L . $occup}))
      )
    ) {
      $desig .= ".<" . $LS_curr . ">" if $LS_curr;
      $added_desig .= ".<" . $LS_curr . ">" if $LS_curr;
      $desig_has_genealogy = 1;
    }
    if ( $shell_ind && &next_shell_requires_accum_LS($parity,$n_conf,$ns,$last_shell_summation_order) ) {
      $desig .= ".(" . $LS_accum . ")" if $LS_accum;
      $added_desig .= ".(" . $LS_accum . ")" if $LS_accum;
      $desig_has_genealogy = 1;
    }
    my $shell_hash = {'LS_curr' => $LS_curr, 'LS_accum' => $LS_accum, 'name' => $sh . $occup,
      'code_name' => $L . $occup,'occup'=>$occup, 'sh'=>$sh, 'L' => $L, 'desig' => $desig, 'added_desig' => $added_desig};
    $shells_summation_order[$ns] = $shell_hash;
    $shell_ind++;
    $nL_prev = $nL;
  }

  #my $num_shells = $#shells_summation_order + 1;
  if ( $reordered ) {
    $last_shell_in_conf = 0;
    for ( my $ns = 0; $ns <= 7; $ns++ ) {
      next unless ($shell_ord_seq{$ns+1});
      my $shell_ind = $shell_ord_seq{$ns+1} - 1;

      $last_shell_in_conf = $shell_ind;

      next unless defined $shells_summation_order[$shell_ind];
      my $shell_hash = $shells_summation_order[$shell_ind];
      my $added_desig = $shell_hash->{'added_desig'};
      $fill_str .= ($fill_str ? '.' : '') . $added_desig if $added_desig;
    }
  } else {
    $fill_str = $desig;
  }

  if ( !$final_term ) {
    $final_term = '1S';
  }

  $fill_str .= "\t$final_term";

  #if ( $fill_str =~ /s0/ ) {
  #  $fill_str = $fill_str;
  #}
  my ($sh,$occ_last,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$last_shell_in_conf]};
  my $L_last = $sh;
  $L_last =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
  my $n_last = $sh + 0;

  $basis[$parity-1]->{'LS'}->{$J}->{$state_num}->{'cn'} = $n_conf;
  $basis[$parity-1]->{'LS'}->{$J}->{$state_num}->{'sh'} = $fill_str;
  $basis[$parity-1]->{'LS'}->{$J}->{$state_num}->{'lastsh'} = [$n_last,$L_last,$occ_last];

  return $L_format;
} ##fill_LS_shell(@)

############################################################################
sub fill_JJ_shell(@) {  #3/16/2004 1:54PM A.Kramida
############################################################################
  my ($parity,$state_num,$n_conf,$J,@p) = @_;
  my ($last_shell, $last_shell_reordered, $last_shell_summation_order) = @{$last_shells[$parity-1]->[$n_conf-1]};  # 0,1,...
  my @fill = ();
  for ( my $i = 0; $i <= 7; $i++ ) {
    my $j = $i*3;
    $fill[$i] = [$p[$j],$p[$j+1],$p[$j+2]];
  }

  my $fill_str = '';
  my $acc_parity = 0;
  my $acc_par_char = '';
  my $first_filled = -1;
  my $reordered = 0;
  my $J_last = '0';
  my $J_prev_to_last = '0';
  my $J_last_shell = '0';
  my ($nL_last, $nL_prev) = (0,0);
  my $last_shell_in_conf = 0;
  my $dot = '';
  for (my $ns = 0; $ns <= 7; $ns++ ) {
    my ($n_shell, $sh, $occup, $sh_par_char) = ();
    if ($ns <= $last_shell) {
      # Get the new sequential order of shells from the hash
      $n_shell = $ns;
      if ( $shell_ord_seq{$ns+1} ) {
        $n_shell = $shell_ord_seq{$ns+1} - 1;
      }
      $reordered = 1 if ($ns != $n_shell);

      ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$n_shell]};
      if ( $sh_par_char ) {
        $acc_parity++;
      }
      if ( $acc_parity % 2 > 0 ) {
        $acc_par_char = '*';
      }
    }
    my $L = $sh;
    $L =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code

    my ($sh1,$occup1) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    if ( $occup1 && ($occup1 != &complete_shells($sh1)) ) {
      $nL_prev = $nL_last;
      $nL_last = $sh1;
      $nL_last =~ s/[a-z]//g;
      $sh1 =~ s/^\d+//g;
      $nL_last = $nL_last*10000 + $L_moment{$sh1};

      my $LS = $fill[$ns]->[0];
      my $J1 = &get_J($fill[$ns]->[1]);
      my $J2 = &get_J($fill[$ns]->[2]);
      $J_prev_to_last = $J_last;
      $J_last = $J2;  # Not re-ordered !  These are tracked in order to determine the final term.
      $J_last_shell = $J1;
    }
    if ( $ns > $last_shell ) {
      next;
    }
    if ( $occup > 0 ) {
      # Get the J value for the shell and the accumulated J value
      my $LS_reordered = $fill[$n_shell]->[0];
      my $J1_reordered = &get_J($fill[$n_shell]->[1]);
      my $J2_reordered = &get_J($fill[$n_shell]->[2]);

      if ( $occup == &complete_shells($L) ) {
        # Move the previously given parent term and total J outside the last
        # complete shell
        if ( ($fill_str =~ /^(.+)\.\(([^()]+)\)\.\(([^()]+)\)$/) && !$reordered ) {
          $fill_str = "$1.$sh$occup.($2).($3)";
          $dot = '.';
        } elsif ($ns != $last_shell) {
          if ($sh ne $lowest_complete_shell ) {
            # Skip the 1s2 in 1s2.2s, but include the 2s2 in 2s2.3s. Give no intermediate term.
            $fill_str .= "$dot$sh$occup";
            $dot = '.';
          }
        } elsif ( $sh ne $lowest_complete_shell ) {
          $fill_str .= "$dot$sh$occup.($LS_reordered$sh_par_char<$J1_reordered>).($J2_reordered)";
          $dot = '.';
        } else {
          $fill_str .= "$sh$occup";
          $dot = '.';
        }
        $last_shell_in_conf = $ns unless $last_shell_in_conf;
      } else {
        $first_filled = $n_shell unless $first_filled >= 0;

        $fill_str .= "$dot$sh";
        $dot = '.';

        if ( $occup > 1 ) {
          $fill_str .= $occup;
        }

        if ( ($LS_reordered ne '') && ($reordered || (($sh ne 's') || ($fill_str =~ /\.\(.+\)/))) ) {
          if ( $sh ne 's' ) {
            $fill_str .= ".($LS_reordered$sh_par_char<$J1_reordered>)";
          }
          if ( (($n_shell != $first_filled) && ($ns != $last_shell)) || ($reordered && $n_shell != $last_shell_reordered)) {
            $fill_str .= ".($J2_reordered)";
          }
        }
        $last_shell_in_conf = $ns;
      }
    }
    if ( ($ns == $last_shell) && !$reordered) {
      # If, due to complete last shells, the previous parent term has
      # moved up tp the very end, detach it and convert into the final term
      if ($fill_str =~ /^([^()]+)\.\(([^()<>]+)<([^<>]+)>\)\.\(([^()]*)\)$/ ) {
        $fill_str = "$1.($2<$3>)";
      }
      if ( $fill_str =~ /^([^()]+)\.\(([^()<>]+)<([^<>]+)>\)$/ ) {
        $fill_str = "$1.($2<$3>)";
      }
      # If there are only two intermediate shell J's and no intermediate accumulating J's,
      # make the final JJ term out of the two intermediate J's
      if ( $fill_str =~ /^([^<> ]*)<([^<> ]*)>([^<> ]*)<([^<> ]*)>\)$/ ) {
        # Detach the last shell's LSJ if it is a singly occupied shell or a shell with one hole
        $fill_str =~ s/([spdfghiklmno]|p5|d9|f13)\.\([^.()]+\)$/$1/;
      }
      # Take the last intermediate accumulating J,
      # make the final JJ term out of it and the last shell's J
      my $fs1 = $fill_str;
      my @Js = ();
      my $nJ = 0;
      my $J11 = '';
      my $term = '';
      while ( $fs1 =~ /^(.+)\.\(([0-9\/]+)\)\./ ) {
        $nJ++;
        $J11 = $2;
        $term = $J11 unless $term; # The regex finds the last intermediate J first
        $Js[$nJ - 1] = $J11;
        # Replace round parentheses with angular brackets
        $fs1 =~ s/\.\($J11\)\./.<$J11>./;
      }
      if ( $fill_str =~ /^(.+)<([^.<>()]+)>\)$/ ) {
        # Detach the last shell's LSJ if it is a singly occupied shell or a shell with one hole
        $fs1 =~ s/([spdfghiklmno]|p5|d9|f13)\.\([^.()]+\)$/$1/;
        $fill_str = $fs1;
      }

    }
  }
  # Remove the unnecessary intermediate states
  $fill_str =~ s/s\.\(2S<1\/2>\)\./s./g;

  # Append the final term to the state name
  my $final_term = ($nL_last > $nL_prev) ? "($J_prev_to_last,$J_last_shell)$acc_par_char"
                                         : "($J_last_shell,$J_prev_to_last)$acc_par_char";
  $fill_str .= "\t$final_term";

  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num}->{'cn'} = $n_conf;
  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num}->{'sh'} = $fill_str;

  my ($sh,$occ_last,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$last_shell_in_conf]};
  my $L_last = $sh;
  $L_last =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
  my $n_last = $sh + 0;
  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num}->{'lastsh'} = [$n_last,$L_last,$occ_last];

}  ##fill_JJ_shell

sub only_complete_shells($$) {
  my ($parity,$n_conf) = @_;
  my $only_complete_shells = 1;
  for ( my $ns = 7; $ns>=0; $ns-- ) {
    my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    $L =~ s/^\d//g;
    if ( ($occup > 0) && ($occup != &complete_shells($L)) ) {
      $only_complete_shells = 0;
      last;
    }
  }
  return $only_complete_shells;
}

sub fill_open_shells() {
  my %all_shells = ();
  for (my $parity = 1; $parity <= 2; $parity++) {
    for ( my $ns = 7; $ns>=0; $ns-- ) {
      foreach my $n_conf (keys %{$shells[$parity-1]}) {
        my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
        $all_shells{$L} = $occup;
      }
    }
  }
  for (my $parity = 1; $parity <= 2; $parity++) {
    for ( my $ns = 7; $ns>=0; $ns-- ) {
      foreach my $n_conf (keys %{$shells[$parity-1]}) {
        my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
        if (($occup != $all_shells{$L}) || ($occup != &complete_shells($L))) {
          $open_shells{$ns} = 1;
        }
      }
    }
  }
}

############################################################################
sub get_last_shell_summation_order($$) { #3/16/2004 2:04PM A.Kramida
############################################################################
  my ($parity,$n_conf) = @_;
  my $last_shell = 0;
  #if (($parity == 1) && ($n_conf == 4)) {
  #  $n_conf = $n_conf;
  #}
  my $only_complete_shells = &only_complete_shells($parity,$n_conf);
  for ( my $ns = 7; $ns>=0; $ns-- ) {
    my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    $L =~ s/^\d//g;
    if ( ($occup > 0) && (($occup != &complete_shells($L)) || $only_complete_shells) ) {
      $last_shell = $ns;
      last;
    }
  }
  return $last_shell;
} ##get_last_shell

############################################################################
sub get_last_shell($$) { #3/16/2004 2:04PM A.Kramida
############################################################################
  my ($parity,$n_conf) = @_;
  my $last_shell = 0;
  #if (($parity == 2) && ($n_conf == 1)) {
  #if (($parity == 1) && ($n_conf == 4)) {
  #  $n_conf = $n_conf;
  #}
  my $only_complete_shells = &only_complete_shells($parity,$n_conf);
  for ( my $ns = 7; $ns>=0; $ns-- ) {
    my $n_shell = $ns;
    if ( $shell_ord_seq{$ns+1} ) {
      $n_shell = $shell_ord_seq{$ns+1} - 1;
    }
    my $occup = $shells[$parity-1]->{$n_conf}->[$n_shell]->[1];
    if ( $occup > 0 ) {
      $last_shell = $ns;
      last;
    }
  }
  return $last_shell;
}  ##get_last_shell

############################################################################
sub get_last_shell_reordered($$) {  #3/16/2004 2:04PM A.Kramida
############################################################################
  my ($parity,$n_conf) = @_;
  my $last_shell = 0;
  my $only_complete_shells = &only_complete_shells($parity,$n_conf);
  for ( my $ns = 7; $ns>=0; $ns-- ) {
    my $n_shell = $ns;
    if ( $shell_ord_seq{$ns+1} ) {
      $n_shell = $shell_ord_seq{$ns+1} - 1;
    }
    my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    $L =~ s/^\d//g;
    if ( ($occup > 0) && (($occup != &complete_shells($L)) || $only_complete_shells) ) {
      $last_shell = $n_shell;
      last;
    }
  }
  return $last_shell;
} ##get_last_shell_reordered

############################################################################
sub next_shell_requires_accum_LS($$$$) {  #03/09/2010 11:20AM
############################################################################
  my ($parity,$n_conf,$ns,$last_shell) = @_;
  my $requires = 0;
  for ( my $ns1 = $ns+1; $ns1 <= $last_shell; $ns1++ ) {
    my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns1]};
    next unless $occup;
    $L =~ s/^\d//g;
    next if ( $occup >= &complete_shells($L) );
    $requires = 1;
    last;
  }
  return $requires;
} ##next_shell_requires_accum_LS($$$$)

############################################################################
sub next_shell_requires_genealogy($$$$$) {  #03/09/2010 11:20AM
############################################################################
  my ($parity,$n_conf,$ns,$last_shell,$curr_desig_no_shell_num) = @_;
  my $requires = 0;
  for ( my $ns1 = $ns+1; $ns1 <= $last_shell; $ns1++ ) {
    my ($L,$occup) = @{$shells[$parity-1]->{$n_conf}->[$ns1]};
    $L =~ s/^\d//g;
    next unless $occup;
    next if ( $occup >= &complete_shells($L) );
    my $next_desig = $curr_desig_no_shell_num . ($curr_desig_no_shell_num ? '.' : '') . $L . $occup;
    $requires = !defined($no_genealogy_shells{$next_desig});
    last;
  }
  return $requires;
} ##next_shell_requires_accum_LS($$$$)

############################################################################
sub read_confs() {   #3/16/2004 2:34PM A.Kramida
############################################################################
  my $s;
  while ( (defined($s = <OUTG11>)) && ($s !~ /0  K   J           CONFIGURATION/) ) {
    next;
  }
  $s = <OUTG11>; # Skip the blank line

  for ( my $k = 1; $k <= 2; $k++ ) {
    while ( (defined($s = <OUTG11>)) && ($s =~ /^ +(\d+) +(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\d+) +([1-]+).{25}(.{12})/) ) {

      # Store the shell filling numbers
      my ($par,$nc,$s1,$w1,$s2,$w2,$s3,$w3,$s4,$w4,$s5,$w5,$s6,$w6,$s7,$w7,$s8,$w8,$n_elec,$par_code,$conf_name) =
      ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21);

      $conf_name =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces

      if ( ($conf_name eq 'dp8p') ) {
        $nc = $nc;
      }
      # Set the name of the config
      $confs[$par-1]->{$nc} = $conf_name;
      # Store the config number for backwards reference
      $conf_nums{$conf_name} = [$par,$nc];
      $parities[$par-1] = ( $par_code eq '-1' ) ? 'o' : 'e';

      # Store the shell filling numbers
      $shells[$par-1]->{$nc} = [[lc($s1),$w1],[lc($s2),$w2],[lc($s3),$w3],[lc($s4),$w4],[lc($s5),$w5],[lc($s6),$w6],[lc($s7),$w7],[lc($s8),$w8]];

      #my $key_prev = '';
      # Add parity codes for each shell
      foreach ( @{$shells[$par-1]->{$nc}} ) {
        my ($sh,$occup) = @{$_};
        my $parity_code = &get_shell_parity($sh,$occup);
        my $parity_char = &get_parity_char($parity_code);
        push(@{$_},$parity_char);
      }

        # Find this conf in confs36 and replace shell names with those from IN36 (including the principal quantum #)
      if ( !defined{$confs36{$conf_name}} || (ref($confs36{$conf_name}) ne 'ARRAY')) {
        die "Configuration $conf_name not found in IN36";
      }
      my @shells36 = @{$confs36{$conf_name}};
      #my %subshells36 = ();
      #for ( my $i = 0; $i <= $#shells36; $i++ ) {
      #  my ($n, $sh, $w) = @{$shells36[$i]};
      #  $subshells36{"$n$sh"} = $w;
      #}
      my %subshells = ();
      #my $key_prev;
      #my $last_shell = &get_last_shell($par,$nc);
      for ( my $i = 0; $i <= $#shells36; $i++ ) {
        my ($n, $sh, $w) = @{$shells36[$i]};
        for ( my $j = 0; $j <= 7; $j++ ) {
          if ( ($sh eq $shells[$par-1]->{$nc}->[$j]->[0]) && ($w == $shells[$par-1]->{$nc}->[$j]->[1])
               && ($shells[$par-1]->{$nc}->[$j]->[0] !~ /^\d/) && ! $subshells{"$n$sh"}
          ) {

            my $key = $n . $L_moment{$sh};
            #next if (($j > $i) && &has_same_shell_greater_n($i, $n, $sh, $w, @shells36));
            #next if $key lt $key_prev;

            $shells[$par-1]->{$nc}->[$j]->[0] = "$n$sh";
            $subshells{"$n$sh"} = $j+1;

            if ( defined($shell_order_in36_back{$j+1}) ) {
              my $w1 = $shell_order_in36_back{$j+1};
              next if ($w1 != 0 ); # This shell was already mapped using another config. where it is non-empty
            }

            $shell_order_in36{$key} = $j+1;
            $shell_order_in36_back{$j+1} = $w;
            #$key_prev = $key;

            last;
          }
        }
      }
    }
  }

  # Renumber all shells in the order of increasing n and L
  my @shells = sort {$a<=>$b} values %shell_order_in36;
  my $i = -1;
  foreach my $key (sort {$a<=>$b} keys %shell_order_in36) {
    $i++;
    $shell_ord_seq{$shells[$i]} = $shell_order_in36{$key};
  }

  $lowest_complete_shell = &get_lowest_complete_shell();
}  ##read_confs

############################################################################
sub has_same_shell_greater_n(@) {   #12/05/2013 5:34PM
############################################################################
  my ($i, $n, $sh, $w, @shells36) = @_;
  for ( my $j = $i+1; $j <= $#shells36; $j++ ) {
    my ($n1, $sh1, $w1) = @{$shells36[$j]};
    return 1 if "$sh1$w1" eq "$sh$w";
  }
  return 0;
} ##has_same_shell_greater_n(@)

############################################################################
sub read_in36() {    #3/16/2004 2:51PM A.Kramida
############################################################################
  open IN36, "<IN36" or die "Could not open input file IN36";
  my $s = '';
  my $parity = -1;

  $s = <IN36>; # Skip the first line

   # Read configurations
  while ( (defined($s = <IN36>)) && ($s =~ /^.{16}(.{9}).{4} *\d* *(\d+)([^ 0-9])(\d*) *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1} *((\d+)([^ 0-9])(\d*)){0,1}/) ) {
    my $conf = $1;
    my ($n1,$s1,$w1,$n2,$s2,$w2,$n3,$s3,$w3,$n4,$s4,$w4,$n5,$s5,$w5,$n6,$s6,$w6,$n7,$s7,$w7,$n8,$s8,$w8,$n9,$s9,$w9) =
       ($2,$3,$4, $6,$7,$8, $10,$11,$12, $14,$15,$16, $18,$19,$20, $22,$23,$24, $26,$27,$28, $30,$31,$32, $34,$35,$36);
    $conf =~ s/^\s+|\s+$//g;

    #if ( $conf eq 'd24s2' ) {
    #  $conf = $conf;
    #}
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
        $sh = lc($sh);
        push(@{$confs36{$conf}},[$n,$sh,$w]);
      }
    }
  }

  my $delete_first = 1;
  my @confs36 = keys %confs36;
  while ( $delete_first && $#confs36) {
    foreach my $conf (@confs36) {
      #if ( ($confs36{$conf}->[0]->[1] ne 'f') || ($confs36{$conf}->[0]->[2] != 14) ) {
      my ($n,$sh,$occup) = @{$confs36{$conf}->[0]};
      my $complete_shell = ($occup == &complete_shells($sh));
      if ( !$complete_shell ) {
        $delete_first = 0;
      }
    }

    if ( $delete_first ) {
      foreach my $conf (keys %confs36) {
        shift(@{$confs36{$conf}});
      }
    }
    @confs36 = keys %confs36;
  }

  close IN36;
  return 1;
} ##read_in36

############################################################################
sub get_conf_name($) {   #6/13/2006 11:10AM A.Kramida
############################################################################
  my $conf = shift;
  my @shells = @{$confs36{$conf}};
  my $conf_name = '';
  my $i = 0;
  my $last_filled_shell = '';
  foreach my $shell (@shells) {
    my ($n,$sh,$w) = @{$shell};
    $i++;
    if ( $w > 0 ) {
      my $w1 = ($w > 1 ? $w : '');
      $last_filled_shell = "$n$sh$w1"
    }
    if ( $w && (($w != &complete_shells($sh)) || ($i > 1) ) ) {
      $w = '' unless $w > 1;
      $conf_name .= ($conf_name ? '.' : '') . "$n$sh$w";
    }
  }
  $conf_name = $last_filled_shell unless $conf_name;
  return $conf_name;
} ##get_conf_name($)

############################################################################
sub get_J($) {   #9/27/2004 3:05PM A.Kramida
############################################################################
  my $J1 = shift;
  if ( $J1 =~/[.]5/ ) {
    $J1 = $J1*2 . '/2';
  } else {
    $J1 =~ s/[.]0$//;
  }
  return $J1;
} ##get_J($)

############################################################################
sub read_energies($$) {  #9/27/2004 8:31PM A.Kramida
############################################################################
  my ($J, $par) = @_;
  my $num_state = 0;
  print "Reading energies: parity $par, J= $J...\n";
  $energies[$par-1]->{$J} = {};
  my $s;
  while ( defined($s = <OUTG11>) && ($s !~ /CONFIG. NO./) ) {
    if ( ($s =~ /^(\s*)$/) && !defined($s = <OUTG11>)  ) {
      # Empty string (or with spaces) and could not read another one
      die "Unexpected end of OUTG11 in EIGENVALUES section, J=$J, parity $par";
    }
    if ( $s =~ /CONFIG. NO./ ) {
      last;
    }
    if ( $s =~ /^\s{28}(.{9})(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1} *$/) {
      my @Earr = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
      my $E = '';
      my $num_E = 0;
      for ( my $i = 0; $i<=10; $i++ ) {
        $E = $Earr[$i];
        if ( defined($E) && ($E =~ / *(-{0,1}\d+\.\d+) */) ) {
          $Earr[$i] = $1;
          $num_E = $i;
        } else {
          last;
        }
      }
      if ( defined($E) && ($E !~ / *(-{0,1}\d+\.\d+) */) ) {
        next;
      }
      for ( my $i = 0; $i <= $num_E; $i++ ) {
        $num_state++;
        $energies[$par-1]->{$J}->{$num_state} = [$Earr[$i]];
      }
    }
  }
} ##read_energies($$)

############################################################################
sub read_vectors($$$) {  #9/27/2004 8:57PM A.Kramida
############################################################################
  my ($J, $cpl1, $par) = @_;
  my $num_state = 0;
  my $num_section = 0;
  my $cpl = 'unknown';

  my $last_bs_num = 0;
  my $s;

  if ( $cpl1 eq '' ) {
    $term_labels_RCG[$par-1]->{$J} = {};
    $term_labels_map[$par-1]->{$J} = {};

    # Read RCG configuration indexes for term labels
    $num_blocks = 0;
    while ( defined($s = <OUTG11>) && ($s !~ /EIGENVECTORS/) ) {
      if ( $s =~ /^(\s*)$/ ) {
        next;
      }
#                               10       18        6       11       12       12       12       12       13        4        7
      if ( $s =~ /^\s{25}(.{9})(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1} *$/) {
        $num_blocks++;
        $num_states_in_blocks[$num_blocks-1] = [];
        my @arr = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
        my $conf_index = '';
        my $num_i = 0;
        for ( my $i = 0; $i<=10; $i++ ) {
          $conf_index = $arr[$i];
          if ( defined($conf_index) && ($conf_index =~ / *(\d+) */) ) {
            $arr[$i] = $1;
            $num_i = $i;
          } else {
            last;
          }
        }
        if ( defined($conf_index) && ($conf_index !~ / *(\d+) */) ) {
          next;
        }
        for ( my $i = 0; $i <= $num_i; $i++ ) {
          $num_state++;
          $num_states_in_blocks[$num_blocks-1]->[$i] = $num_state;
          my $conf_num = $arr[$i];
          $conf_num = s/^\s+|$\s+//g;
          $term_labels_RCG[$par-1]->{$J}->{$num_state} = [$arr[$i]];
        }
      }
    }
    #if ( $J eq '1.5' ) {
    #  $J = $J;
    #}
  } else {
    $cpl = $cpl1;
  }

  if ( $s =~ /EIGENVECTORS   \( *(\S+) COUPLING/ ) {
    $cpl = $1;
  }

  print "Reading vectors for $cpl coupling...\n";

  # Initialize vectors hash for this J and coupling
  $vectors[$par-1]->{$cpl}->{$J} = {};

  while ( defined($s = <OUTG11>) && ($s !~ /PURITY=/) && ($s !~ /TIME=/) ) {
    if ( ($s =~ /^(\s*)$/) && !defined($s = <OUTG11>)  ) {
      # Empty string (or with spaces) and could not read another one
      die "Unexpected end of OUTG11 in EIGENVECTORS section, J=$J, parity $par";
    }
    if ( ($s =~ /PURITY=/) || ($s =~ /TIME=/) ) {
      last;
    }

    # Read blocks of eigenvector data
    # $row is the block number
    $num_state = 0;
    my $read_term_labels = ($cpl1 ? 0 : 1);

    my @state_nums = keys %{$energies[$par-1]->{$J}};
    my $num_states = $#state_nums + 1;
    for ( my $row = 1; $row <= $num_blocks; $row++ ) {
      my $rx = '^\s{16}(.{5}).{9}(.{8})( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1}( .{8}){0,1} *$';
      if ( $s !~ /rx/ ) {
        while ( defined($s = <OUTG11>) && ($s =~ /^\s*$/) ) {
          next;
        }
      }
      if ( $read_term_labels ) {
        # $s has the RCG configuration labels; read the next string with term labels
        my $s1 = <OUTG11>;
        # $s1 has the RCG term labels
  #                   1         3p5      p44p     p44p     p44p     p44p     p44p     p44p     p45p     p44f     sp5d     p45p
  #                             (1S) 2P  (3P) 4P  (3P) 4D  (3P) 2P  (3P) 2S  (1D) 2P  (1S) 2P  (3P) 4P  (3P) 4D  (3P) 4P  (3P) 4D
        if ( $s =~ /$rx/) {
          my @arr_conf = ($2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
          my $row_num = $1;
          $row_num =~ s/^\s+|\s+//g;
          if ( $row_num != $row ) {
            die("Row number mismatch while reading eigenvectors for J=$J of parity $par: expected $row, read $row_num.");
          }
          $s1 =~ /^\s{30}(.{9})(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1}(.{9}){0,1} *$/;
          my @arr_term = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
          foreach (@arr_term) {
            $_ =~s /^\s+|\s+//g;
          }
          my $conf_label = '';
          my $term_label = '';
          my $num_i = 0;
          for ( my $i = 0; $i<=10; $i++ ) {
            $conf_label = $arr_conf[$i];
            $term_label = $arr_term[$i];
            if ( defined($conf_label) && ($conf_label =~ /^ *(.+) *$/) ) {
              $arr_conf[$i] = $1;
              $term_label =~ s/^\s+|\s+$//g;
              if ( $term_label eq '' ) {
                die "Empty LS term label in OUTG11, parity $par, J=$J, block $row";
              }
              $arr_term[$i] = $term_label;
              $num_i = $i;
            } else {
              last;
            }
          }
          if ( defined($conf_label) && ($conf_label !~ /^ *(.+) *$/) ) {
            next;
          }
          #if ( ($row == 2) ) {
          #  $row = $row;
          #}
          for ( my $i = 0; $i <= $num_i; $i++ ) {
            $num_state = $num_states_in_blocks[$row-1]->[$i];
            push(@{$term_labels_RCG[$par-1]->{$J}->{$num_state}}, $arr_term[$i]);
          }
        } else {
          die("OUTG11 format error in eigenvector block $row for J=$J of parity $par");
        }

        #$s = <OUTG11>;
        #$read_term_labels = 0;
      }
# 1:d10s2 e     (1S) 1S   1   0.99602 -0.07021 -0.05286  0.00197 -0.00070  0.00071  0.00421 -0.00425  0.01279  0.00454  0.00035
# 1:f4          (5D) 5D   1   0.71052  0.45211 -0.52783  0.04458 -0.10010  0.01233 -0.00003 -0.00014 -0.00226  0.00048  0.00112
# 2:f35p        (2G) 3F1 30  -0.06121 -0.13791 -0.02937  0.46859 -0.28947 -0.06323 -0.03097  0.00823  0.11806  0.01837  0.01562
# 2:d9f5        (2P) 3P1119  -0.00004 -0.00012 -0.00014  0.00014  0.00012  0.00012 -0.00011  0.00006 -0.01127 -0.01782  0.00743
# 1:f4          (5D 0) 0   1  0.71052  0.45211 -0.52783  0.04458 -0.10010  0.01233 -0.00003 -0.00014 -0.00226  0.00048  0.00112
# 1:f4          (3D 6) 6  11  0.00349  0.00070  0.00114  0.00010 -0.00004 -0.00014  0.00017 -0.00028  0.00009  0.00000 -0.00004
      $rx = '^.{24}\s*(\d+|[*]+) +(-*\d+\.\d+) *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)* *(-*\d+\.\d+)*$';
      while ( $s !~ /$rx/ ) {
        $s = <OUTG11>;
      }
      my $states_read = 0;
      while ( $s =~ /$rx/ ) {
        $states_read++;
        my ($basis_state,@A_arr) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
        #foreach (@A_arr) {
        #  $_ =~s /^\s+|\s+//g;
        #}
        if ( $basis_state == 1 ) {
          #$num_section++;
          $last_bs_num = 0;
        } elsif ( $basis_state eq '***' ) {
          $basis_state = $last_bs_num + 1;
        }
        $last_bs_num = $basis_state;

        #$num_state = ($num_section - 1) * 11;
        my $states_in_row = $#{$num_states_in_blocks[$row-1]};
        #for ( my $i = 0; $i<=10; $i++ ) {
          #$num_state++;    # Number of the energy state
        for ( my $i = 0; $i <= $states_in_row; $i++ ) {
          $num_state = $num_states_in_blocks[$row-1]->[$i];

          my $A = $A_arr[$i];
          # If it is the last, incomplete section of the eigenvectors printout,
          # some of the $1,$2,$3,...,$10 variables will be empty,
          # and the last non-empty energy will be in $11
          #if ( !defined($A) || ($A eq '') )  {
          #  $num_state--;
          #  next;
          #}
          $A = $A + 0;
          if ( $basis_state == 1 ) {
            # On the first line of the eigenvectors section, init the hashes for each energy state
            $vectors[$par-1]->{$cpl}->{$J}->{$num_state} = {};
          }
          # Store only components that have percentage greater than 0.0001%
          if ( abs($A) >= 0.001 ) {
            $vectors[$par-1]->{$cpl}->{$J}->{$num_state}->{$basis_state} = $A;
          }
        }
        $s = <OUTG11>;
        last unless defined($s);
      }
      if ( $states_read != $num_states ) {
        die "Number of states mismatch in eigenvector block $row of OUTG11, parity $par, $cpl coupling, J=$J: expected $num_states, read $states_read";
      }
    }

    if ( $cpl1 eq '' ) {
      foreach my $ns ( keys %{$term_labels_RCG[$par-1]->{$J}} ) {
        my $key = join("\t",@{$term_labels_RCG[$par-1]->{$J}->{$ns}});
        #if ( length($key) < 9 ) {
        #  $key = $key;
        #}
        $term_labels_map[$par-1]->{$J}->{$key} = $ns;
      }
    }
    while ( defined($s = <OUTG11>) ) {
      if ( ($s =~ /PURITY=/) || ($s =~ /TIME=/) ) {
        last;
      }
    }
    last;
  }
  return $cpl;
} ##read_vectors($$$)

############################################################################
sub ReadRCE {  #9/28/2004 1:38PM A.Kramida
############################################################################
  if ( !(-f "RCEOUT") ) {
    return 1;
  }
  open RCEOUT, "<RCEOUT" or die "Could not open input file RCEOUT";

  print "Reading RCEOUT...\n";
  my $s;
  # Build basis key hash
  my @basis_key = ({},{});
  for ( my $par = 0; $par<=1; $par++ ) {
    foreach my $J (keys %{$basis[$par]->{'LS'}}) {
      $basis_key[$par]->{$J} = {} unless defined $basis_key[$par]->{$J};
      foreach my $state_num ( keys %{$basis[$par]->{'LS'}->{$J}} ) {
        my $c_no  = $basis[$par]->{'LS'}->{$J}->{$state_num}->{'cn'};
        my $c_name= $confs[$par]->{$c_no};
        # Truncate conf. name to max 5 symbols with no trailing spaces
        # because RCE does so
        if ( $c_name =~ /^(.{6})/ ) {
          $c_name = $1;
          $c_name =~ s/\s+$//g;
        }
        my $label = $basis[$par]->{'LS'}->{$J}->{$state_num}->{'label'};
        $basis_key[$par]->{$J}->{$c_name} = {} unless defined $basis_key[$par]->{$J}->{$c_name};
        my $existing = $basis_key[$par]->{$J}->{$c_name}->{$label};
        if ( $existing ) {
          $par++;
          die "Duplicate label for LS basis state, par=$par; J=$J; conf=$c_name; term=$label; state num=$state_num and $existing";
        }
        $basis_key[$par]->{$J}->{$c_name}->{$label} = $state_num;
      }
    }
  }

  # Read levels of each parity
  for ( my $par = 0; $par<=1; $par++ ) {
    my $regex = '^ *([0-9-]+\.\d{0,6})([* ]+)(-{0,1}\d+\.\d{1,3}) +(\d+)(\/2)*([ 0-9-]{4}) (.{6}) +(\S{4})( *\S+)( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}( *([0-9-]+) (.{6}) +(\S{4})( *\S+)){0,1}';
    #while ( (defined($s = <RCEOUT>)) && ($s !~ /^ *([0-9-]+\.\d{0,6})([* ]+)(-{0,1}\d+\.\d{1,3}) +([0-9\/]+) +([0-9-]+) (.{6}) +(\S+) (\S+)( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) +(.{6}) +(\S+) (\S+)){0,1}/ ) ) {
    while ( (defined($s = <RCEOUT>)) && ($s !~ /$regex/ ) ) {
      next;
    }
    my $Jprev = '';
    my $i = 0;
    do {
      #$s =~ /^ *([0-9-]+\.\d{0,6})([* ]+)(-{0,1}\d+\.\d{1,3}) +(\d+)(\/(\d+))* +([0-9-]+) (.{6}) +(\S+) (\S+)( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}( *([0-9-]+) (.{6}) +(\S+) (\S+)){0,1}/;
      $s =~ /$regex/;
      my ($Ee, $exp_c, $Ec, $J1,$J2, $p1, $conf1, $term1) = ($1,$2,$3,$4,$5,$6,$7,"$8$9");
			$p1 =~ s/^\s+|\s+$//g;
      my $A2 = [$11,$12,"$13$14"];
      my $A3 = [$16,$17,"$18$19"];
      my $A4 = [$21,$22,"$23$24"];
      my $A5 = [$26,$27,"$28$29"];
      my $A6 = [$31,$32,"$33$34"];
      my $A7 = [$36,$37,"$38$39"];
      my $A8 = [$41,$42,"$43$44"];
      my $A9 = [$46,$47,"$48$49"];
      my $A10 = [$51,$52,"$53$54"];
      $exp_c =~ s/\s+//g; # Trim spaces
      $conf1 =~ s/^\s+|\s+$//g;
      my $J = $J1;
      # Make J look the same as in RCG
      if ( $J2 ) {  # Half-integer J
        $J = $J1/2;
        $J = sprintf("%4.1f",$J);
        $J =~ s/ //g;
      } else {
        $J .= '.0';
      }
      if ( $Jprev ne $J ) {
        $RCE_lev[$par]->{$J} = [];
        $i = 0;
      }
      $Jprev = $J;
      $i++;
      $RCE_lev[$par]->{$J}->[$i-1] = {};
      $RCE_lev[$par]->{$J}->[$i-1]->{'Ee'} = $Ee;
      $RCE_lev[$par]->{$J}->[$i-1]->{'exp_c'} = $exp_c;
      $RCE_lev[$par]->{$J}->[$i-1]->{'Ec'} = $Ec;
      $RCE_lev[$par]->{$J}->[$i-1]->{'v'} = [[$p1,$conf1,$term1]];

      my @A_arr = ($A2,$A3,$A4,$A5,$A6,$A7,$A8,$A9,$A10);
      for ( my $j = 0; $j<=8; $j++ ) {
        #my $c = 0;
        #my $expression = "\$c = \$A$j";
        #eval($expression);
        my $c = $A_arr[$j];

        # Store only components with amplitudes greater than 10
        if ( abs($c->[0]) >= 10  ) {
          $c->[1] =~ s/^\s+|\s+$//g;
          push(@{$RCE_lev[$par]->{$J}->[$i-1]->{'v'}}, $c);
        } else {
          last;
        }
      }
      my $RCE_vector = $RCE_lev[$par]->{$J}->[$i-1]->{'v'};
      for ( my $j = 0; $j <= 9; $j++ ) {
        last unless defined($RCE_vector->[$j]);
        my ($p, $conf, $term) = @{$RCE_vector->[$j]};
        last unless $conf;
        #my $n_state = 0;
        #foreach my $state_num ( keys %{$basis[$par]->{'LS'}->{$J}} ) {
        #  my $c_no  = $basis[$par]->{'LS'}->{$J}->{$state_num}->{'cn'};
        #  my $c_name= $confs[$par]->{$c_no};
        #  # Truncate conf. name to max 5 symbols with no trailing spaces
        #  # because RCE does so
        #  if ( $c_name =~ /^(.{6})/ ) {
        #    $c_name = $1;
        #    $c_name =~ s/\s+$//g;
        #  }
        #  my $label = $basis[$par]->{'LS'}->{$J}->{$state_num}->{'label'};
        #  if ( ($c_name eq $conf) && ($label eq $term) ) {
        #    $n_state = $state_num;
        #    last;
        #  }
        #}
        my $n_state = $basis_key[$par]->{$J}->{$conf}->{$term} if defined $basis_key[$par]->{$J}->{$conf};
        if ( !$n_state ) {
          $par = $par+1;
          die "Could not identify the RCE basis state (parity = $par, J = $J, $conf $term) with RCG basis states."
        }
        $RCE_lev[$par]->{$J}->[$i-1]->{'v'}->[$j]->[3] = $n_state;
      }
    } while ( (defined($s = <RCEOUT>)) && ($s !~ /PARAMETER/) );
  }
  close RCEOUT;

  if (-f 'LEVELS1') {
    open LEVELS1, "<LEVELS1" or die "Could not open input file LEVELS1";
    print "Reading LEVELS1...\n";
    # Read levels of first parity only and get the Lande g value for each level
    my $num_parities = 1;
    $num_parities++ if defined($RCE_lev[1]);
    for ( my $par = 0; $par <= $num_parities-1; $par++ ) {
      my $j = 0;
#247010.577000 247051.4959-40.918902  1.267 2p2     J= 0.5      70% 2p2    (3P) 2P       30% 2p2    (3P) 4P        1% 2p2    (1S) 2S     
#      while ( (defined($s = <LEVELS1>)) && ($s !~ /^ *([0-9-]+\.\d{0,6})([* ]+)([0-9-]+\.\d{1,4}) +[.0-9-]+ +([.0-9-]+) .{6} +J= *([.0-9]+) +\d+[%] .{6} +[(]\d+[SPDFGHIKLMNOQRT][)]/ ) ) {
	    my $regex = '^([ 0-9-]{6}\.[\d ]{6})([ *])([ 0-9-]{6}\.\d{4})[ .0-9-*]{10} +([.0-9-]+) .{6} +J= *([.0-9]+) +\d+[%] .{6} +[(]\d+[SPDFGHIKLMNOQRT][)]';
      while ( (defined($s = <LEVELS1>)) && ($s !~ /$regex/ ) ) {
        next;
      }
      $j++;
      my $Jprev = '';
      my $i = 0;
      my $num_levs = 0;
      next unless defined $s;
      do {
#        $s =~ /^ *([0-9-]+\.\d{0,6})([* ]+)([0-9-]+\.\d{1,4}) +[.0-9-]+ +([.0-9-]+) .{6} +J= *([.0-9]+) +/;
        $s =~ /$regex/;
        my ($Ee, $exp_c, $Ec, $lande_g, $J) = ($1,$2,$3,$4,$5);
        foreach ($Ee, $exp_c, $Ec, $lande_g) {
          $_ =~ s/^\s+|\s+$//g;
        }
        if ( !defined($RCE_lev[$par]->{$J}) ) {
          die "J value $J in LEVELS1 is not defined in RCEOUT.";
        }
        $exp_c =~ s/\s+//g;                                     # Trim spaces
        if ( $Jprev ne $J ) {
          $i = 0;
          my @xx = @{$RCE_lev[$par]->{$J}};
          $num_levs = $#xx + 1;
        }
        $Jprev = $J;
        $i++;
        if ( $i > $num_levs ) {
          die "Number of levels with J=$J in LEVELS1 is greater than in RCEOUT.";
        }
        if ( $RCE_lev[$par]->{$J}->[$i-1]->{'Ee'} ne $Ee) {
          $par++;
          die "Energy of level $i in parity $par with J=$J in LEVELS1, $Ee, differs from value in RCEOUT, ".$RCE_lev[$par]->{$J}->[$i-1]->{'Ee'};
        }
        if ( $RCE_lev[$par]->{$J}->[$i-1]->{'exp_c'} ne $exp_c) {die "Exp_Th flag of level $i with J=$J in LEVELS1, '$exp_c', differs from value in RCEOUT, ".$RCE_lev[$par]->{$J}->[$i-1]->{'exp_c'};}
        if ( abs($RCE_lev[$par]->{$J}->[$i-1]->{'Ec'} - $Ec) > 0.0006) {die "Calc. energy of level $i with J=$J in LEVELS1, $Ec, differs from value in RCEOUT, ".$RCE_lev[$par]->{$J}->[$i-1]->{'Ec'};}
        $RCE_lev[$par]->{$J}->[$i-1]->{'Ec'} = $Ec;            # Ec has one more digit in LEVELS1 compared to RCEOUT
        $RCE_lev[$par]->{$J}->[$i-1]->{'lande'} = $lande_g;
#      } while ( (defined($s = <LEVELS1>)) && ($s =~ /^ *([0-9-]+\.\d{0,6})([* ]+)([0-9-]+\.\d{1,4}) +[.0-9-]+ +([.0-9-]+) .{6} +J= *([.0-9]+) +/ ) && $j++);
      } while ( (defined($s = <LEVELS1>)) && ($s =~ /^([ 0-9-]{6}\.[\d ]{6})([* ])([ 0-9-]{6}\.\d{4})[ *.0-9-]{10} +([.0-9-]+) .{6} +J= *([.0-9]+) +/ ) && $j++);
      unless ( $par ) {
        # Skip second coupling scheme
        my $j1 = 0;
        do {
          last unless defined($s = <LEVELS1>);
          $j1++
#        } while ( ($j1 < $j) && ($s =~ /^ *([0-9-]+\.\d{0,6})([* ]+)([0-9-]+\.\d{1,4}) +[.0-9-]+ +([.0-9-]+) .{6} +J= *([.0-9]+) +/) );
        } while ( ($j1 < $j) && ($s =~ /^([ 0-9-]{6}\.[\d ]{6})([* ])([ 0-9-]{6}\.\d{4})[ *.0-9-]{10} +([.0-9-]+) .{6} +J= *([.0-9]+) +/) );
      }
    }
    close LEVELS1;
  }
  return 0;
} ##ReadRCE

############################################################################
sub sort_vectors() {    #12/06/2013 4:04PM
############################################################################
  print "Sorting RCG vectors...\n";
  for ( my $par = 0; $par<=1; $par++ ) {
    $vectors_sorted[$par] = {};
    $vectors_index[$par] = {};
    #$vectors_J_index[$par] = {};
    $basis_hash[$par] = {};
    foreach my $J ( keys %{$energies[$par]} ) {
      #if ( $J eq '2.0' ) {
      #  $J = $J;
      #}
      $vectors_sorted[$par]->{'LS'}->{$J} = [];
      $vectors_index[$par]->{'LS'}->{$J} = [];
      #$vectors_J_index[$par]->{'LS'}->{$J} = [];
      $basis_hash[$par]->{$J} = [];
      foreach my $num_RCG_e ( keys %{$energies[$par]->{$J}}) {
        #if (!$num_RCG_e) {
        #  my $x = 0;
        #}
        my $sorted = [];
        my $j = 0;
        # Sort each vector in the order of decreading absolute value of amplitude,
        # scale amplitudes by a factor of 100 and retain only up to 10 leading components
        # in the sorted array
        foreach my $num_state ( sort {abs($vectors[$par]->{'LS'}->{$J}->{$num_RCG_e}->{$b})<=>abs($vectors[$par]->{'LS'}->{$J}->{$num_RCG_e}->{$a})} keys %{$vectors[$par]->{'LS'}->{$J}->{$num_RCG_e}} ) {
          my $A = $vectors[$par]->{'LS'}->{$J}->{$num_RCG_e}->{$num_state} * 100;
          push(@{$sorted}, [$A,$num_state]);
          $vectors_index[$par]->{'LS'}->{$J}->[$num_state-1] = [] unless defined $vectors_index[$par]->{'LS'}->{$J}->[$num_state-1];
          push(@{$vectors_index[$par]->{'LS'}->{$J}->[$num_state-1]}, $num_RCG_e) if $j <= 59;
          #$vectors_J_index[$par]->{'LS'}->{$J}->[$num_RCG_e-1] = [] unless defined $vectors_J_index[$par]->{'LS'}->{$J}->[$num_RCG_e-1];
          #push(@{$vectors_J_index[$par]->{'LS'}->{$J}->[$num_RCG_e-1]}, $num_state) if $j <= 9;
          my $n_conf = $basis[$par]->{'LS'}->{$J}->{$num_state}->{'cn'};
          my $label = $basis[$par]->{'LS'}->{$J}->{$num_state}->{'label'};
          $basis_hash[$par]->{$J}->[$n_conf-1] = {} unless defined $basis_hash[$par]->{$J}->[$n_conf-1];
          $basis_hash[$par]->{$J}->[$n_conf-1]->{$label} = [] unless defined $basis_hash[$par]->{$J}->[$n_conf-1]->{$label};
          push(@{$basis_hash[$par]->{$J}->[$n_conf-1]->{$label}}, [$A,$num_RCG_e,$num_state]);
          $j++;
          last if ( ($j > 9) || ($A == 0) );
        }
        $vectors_sorted[$par]->{'LS'}->{$J}->[$num_RCG_e-1] = $sorted;
      }
      # Sort basis hash so that for each conf. num and term label point to an array sorted
      # in decreasing order by absolute value of eigenvector component
      for (my $nc = 0; $nc <= $#{$basis_hash[$par]->{$J}}; $nc++) {
        foreach my $label (keys %{$basis_hash[$par]->{$J}->[$nc]}) {
          my @sorted = sort {abs($b->[0])<=>abs($a->[0])} @{$basis_hash[$par]->{$J}->[$nc]->{$label}};
          $basis_hash[$par]->{$J}->[$nc]->{$label} = \@sorted;
        }
      }
    }
  }
} ##sort_vectors()

############################################################################
sub Identify_RCE_levs() {  #9/28/2004 4:15PM A.Kramida
############################################################################
  my $Em1 = shift;
  my $Em2 = shift;
  @Emax = ($Em1+0, $Em2+0);
  print "Identifying RCE levels with RCG levels...\n";
  for ( my $par = 0; $par<=1; $par++ ) {
    foreach my $J ( keys %{$energies[$par]} ) {
      $map_RCG_RCE[$par]->{$J} = {};
      
      # Determine the maximum calculated energy for this $J, $par that has an experimental energy
      my $Ece_max = undef;
      for (my $num_RCE_E = 1; $num_RCE_E <= $#{$RCE_lev[$par]->{$J}}+1; $num_RCE_E++ ) {
        if ($RCE_lev[$par]->{$J}->[$num_RCE_E-1]->{'exp_c'} ne '*') {
          if (!defined($Ece_max) || ($RCE_lev[$par]->{$J}->[$num_RCE_E-1]->{'Ec'} > $Ece_max)) {
            $Ece_max = $RCE_lev[$par]->{$J}->[$num_RCE_E-1]->{'Ec'};
          }
        }
      }
      $Ece_max += 0.005 if defined($Ece_max);

      for (my $num_RCE_E = 1; $num_RCE_E <= $#{$RCE_lev[$par]->{$J}}+1; $num_RCE_E++ ) {
        if ( ($par == 1) && ($J eq '3.5' ) && (($num_RCE_E == 1766) || ($num_RCE_E == 1768)) )  {
          $par = $par;
        }
        my $Ec = $RCE_lev[$par]->{$J}->[$num_RCE_E-1]->{'Ec'};
        #if ( ($par == 1) && ($J eq '3.5' ) && (abs($Ec - (-2123.3607)) <0.01) ) {
        #  $par = $par;
        #}

        my $best_match_RCG = 0;
        if (defined($Ece_max) && ($Ec < $Ece_max)) {
          my @nums_RCG_e = ();
          my $RCE_vector = $RCE_lev[$par]->{$J}->[$num_RCE_E-1]->{'v'};
          my $hash = {};
          #if ( ($par == 0) && ($J eq '2.0') ) {
          #  $Ec = $Ec;
          #}
          for ( my $j = 0; $j <= 9; $j++ ) {
            last unless defined($RCE_vector->[$j]) && defined($RCE_vector->[$j]->[3]);
            my ($p, $conf, $term, $num_state) = @{$RCE_vector->[$j]};
            last unless $conf;

            if ($num_state) {
              # Add all RCG eigenstates that have this eigenvector component
              #if ( !defined($vectors_index[$par]->{'LS'}->{$J}->[$num_state-1]) ) {
              #  $num_state = $num_state;
              #}
              my @RCG_nums = @{$vectors_index[$par]->{'LS'}->{$J}->[$num_state-1]};
              foreach my $num_RCG_e (@RCG_nums) {
                if ( !defined($hash->{$num_RCG_e}) ) {
                  $hash->{$num_RCG_e} = 1;
                  push(@nums_RCG_e,$num_RCG_e);
                }
              }
            }
          }
          # Go through the 10 RCE vector components and find corresponding RCG basis-state components
          my $minD = 1e17;
          
          foreach my $num_RCG_e (sort { my ($E1,$E2) = ($energies[$par]->{$J}->{$a}->[0],
                                  $energies[$par]->{$J}->{$b}->[0]);
                                  abs($E1-$Ec)<=>abs($E2-$Ec);
                                } @nums_RCG_e) {
            next if defined($map_RCG_RCE[$par]->{$J}->{$num_RCG_e}); # Skip already mapped RCG states
            my ($D1,$D2) = (0,0);
            # Compare ten amplitudes of the RCE vector with those of the RCG vector if these components are found there
            my $RCE_vector = $RCE_lev[$par]->{$J}->[$num_RCE_E-1]->{'v'};
            for ( my $j = 0; $j <= 9; $j++ ) {
              last unless defined($RCE_vector->[$j]) && defined($RCE_vector->[$j]->[3]);
              my ($p, $conf, $term, $num_state) = @{$RCE_vector->[$j]};
              last unless $conf;

              my $A = 0;
              if ($num_state) {
                # If the RCG basis state was identified, use its amplitude
                $A = $vectors[$par]->{'LS'}->{$J}->{$num_RCG_e}->{$num_state} * 100;
              }
              $D1 += ($A-$p)*($A-$p);
              $D2 += ($A+$p)*($A+$p);
            }
            # Add amplitudes of the RCG vector that are not present in the RCE one
            my $j = 0;
            my $max_RCG_comp = $#{$vectors_sorted[$par]->{'LS'}->{$J}->[$num_RCG_e-1]};
            $max_RCG_comp = 9 if $max_RCG_comp > 9;
            for ( my $m = 0; $m <= $max_RCG_comp; $m++ ) {
              my ($A,$num_state) = @{$vectors_sorted[$par]->{'LS'}->{$J}->[$num_RCG_e-1]->[$m]};
              $j++;
              my $found_in_RCE = -1;
              for ( my $k = 0; $k <= 9; $k++ ) {
                last unless defined($RCE_vector->[$k]) && defined($RCE_vector->[$k]->[3]);
                my ($p, $conf, $term, $ns) = @{$RCE_vector->[$k]};
                last unless $conf;
                next unless ($ns == $num_state);
                $found_in_RCE = $k;
                last;
              }
              next unless ($found_in_RCE < 0);

              # If the RCE basis state was not identified, add the RCG amplitude to the "distance"
              #my $A = $vectors[$par]->{'LS'}->{$J}->{$num_RCG_e}->{$num_state} * 100;

              $D1 += $A*$A;
              $D2 += $A*$A;
              last if $j >= 4;
            }

            if ( $D1 > $D2 ) {
              $D1 = $D2;
            }
            if ( $D1 < $minD ) {
              $minD = $D1;
              $best_match_RCG = $num_RCG_e;
            }
          }
          if ( !$best_match_RCG && (($Emax[$par] == 0) || ($Ec <= $Emax[$par])) ) {
            $par++;
            die "Could not identify RCE state (parity = $par, J = $J, Ec = $Ec) with an RCG eigenstate.";
          }
        } else {
          # Assume that the calculated energies in RCEOUT are sorted in exactly the same order as in OUTG11
          $best_match_RCG = $num_RCE_E;
          if (defined($map_RCG_RCE[$par]->{$J}->{$best_match_RCG})) {
            die "RCE state # $num_RCE_E (parity = $par, J = $J, Ec = $Ec) has already been mapped to RCG eigenstate # $map_RCG_RCE[$par]->{$J}->{$best_match_RCG}.";
          }
        }
        $map_RCG_RCE[$par]->{$J}->{$best_match_RCG} = $num_RCE_E;
      }
    }
  }
} ##Identify_RCE_levs

############################################################################
sub read_ING11_params() {  #4/25/2005 1:58PM A.Kramida
############################################################################
  open ING11, "<ING11" or return;
  my $s = '';
  my $Eav_shift = 0;
  my $param_count = 0;
  my $parity = -1;
  my $started = 0;
  while ( defined($s = <ING11> )) {
    if ( !$started && ($s =~ /^([ 0]{5}).{15}([ .e+\d-]){10}/i) && ($1 == 0) && $2) {
      # Rescale card
      $Eav_shift = $2;
      $Eav_shift =~ s/^\s+|\s+$//g;
    }
    if (!$started && ($s !~ /^(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)  (.{6})(.{12}) *([0-9.-]+) +([0-9.-]+)$/ )) {
      # Optional input cards and first main line with program options
      next;
    }
    $started = 1;
    if ($s =~ /^(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)  (.{6})(.{12}) *([0-9.-]+) +([0-9.-]+)$|\/R[12]\// ) {
      # Shell definitions or E1 or E2 matrix elements
      next;
    } else {
      chomp $s;
      #if ( $s =~ /^(.{6})([^-]{12})(\d\d| \d) *([0-9.e+-]{1,9})(\d) *([0-9.e+-]{1,9})(\d) *([0-9.e+-]{1,9})(\d) *([0-9.e+-]{1,9})(\d) *([0-9.e+-]{1,9})(\d)([A-Z]{2}\d{8})\s*$/i  ){
      if ( $s =~ /^(.{6})([^-]{12})(\d\d| \d)([ 0-9.e+-]{9})(\d)([ 0-9.e+-]{9})(\d)([ 0-9.e+-]{9})(\d)([ 0-9.e+-]{9})(\d)([ 0-9.e+-]{9})(\d)([A-Z]{2}\d{8})\s*$/i  ){
        my ($spectr, $conf, $num_param, $p1, $t1, $p2, $t2, $p3, $t3, $p4, $t4, $p5, $t5, $scaling) =
          ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14);
        $p1 =~ s/^\s+|\s$//g;
        $p2 =~ s/^\s+|\s$//g;
        $p3 =~ s/^\s+|\s$//g;
        $p4 =~ s/^\s+|\s$//g;
        $p5 =~ s/^\s+|\s$//g;

        $conf =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces

        my ($par,$nc) = ();
        if ( defined($conf_nums{$conf}->[0]) ) {
          ($par,$nc) = @{$conf_nums{$conf}}; # Get the number of this config
        }
        die "Config. name $conf not recognized in ING11" unless $nc;

        if ( $par != $parity ) {
          $param_count = 0;
          $parity = $par;
        }

        # Initialize the parameters hash for this config
        $params[$par-1]->{$nc} = {};

        $params[$par-1]->{$nc}->{'np'} = $num_param;
        $params[$par-1]->{$nc}->{'scaling'} = $scaling;

        my $n = ($num_param <=5) ? $num_param : 5;
        for ( my $j = 1; $j <= $n; $j++ ) {
          my ($p, $t) = ();
          $param_count++;
          my $expression = "(\$p, \$t) = (\$p$j,\$t$j)";
          eval($expression);
          if ( $Eav_shift && ($j == 1) ) {
            if ( $p =~ /[.e]/ ) {
              $p += $Eav_shift;
            } else {
              $p += sprintf("%d",$Eav_shift*10000);
            }
          }
          $params[$par-1]->{$nc}->{"p$j"} = [$p,$t,'',$param_count];
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
            if ( $s =~ /^ +([0-9.e+-]{1,9})(\d)( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1}( +([0-9.e+-]{1,9})(\d)){0,1} *$/i ) {
              my $num_to_read = ($rest >=7) ? 7 : $rest;

              my ($p1, $t1, $p2, $t2, $p3, $t3, $p4, $t4, $p5, $t5, $p6, $t6, $p7, $t7) =
                ($1, $2, $4, $5, $7, $8, $10, $11, $13, $14, $16, $17, $19, $20);

              my $par_num = $num_param - $rest;

              for ( my $j = 1; $j <= $num_to_read; $j++ ) {
                $par_num++;
                $param_count++;
                my ($p,$t) = ();
                my $expression = "(\$p,\$t) = (\$p$j,\$t$j)";
                eval($expression);

                $params[$par-1]->{$nc}->{"p$par_num"} = [$p,$t,'',$param_count];
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
      } elsif ( $s =~ /^(.{9})-(.{8})([ 0-9]{2}) *([0-9.e+-]{5,9})5 *([0-9.e+-]{5,9})5 *([0-9.e+-]{5,9})5 *([0-9.e+-]{5,9})5 *([0-9.e+-]{5,9})5([A-Z]{2}\d{8}) *$/i ) {
#2p4      -2p6      11496.28845   0.00005   0.00005   0.00005   0.00005HR95999595
        # CI parameters section
        my ($conf1,$conf2,$num_param, $p1,$p2,$p3,$p4,$p5,$scaling) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
        $num_param =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces
        $conf1 =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces
        $conf2 =~ s/^\s+|\s+$//g;  # Strip leading and trailing spaces

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

        $CI[$par1-1]->[$nc1-1] = [] unless defined($CI[$par1-1]->[$nc1-1]);
        $CI[$par1-1]->[$nc1-1]->[$nc2-1] = {} unless defined($CI[$par1-1]->[$nc1-1]->[$nc2-1]);
        $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'} = [] unless defined($CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'});
        $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'np'} = $num_param;

        my $n = ($num_param <= 5) ? $num_param : 5;
        for ( my $j = 1; $j <= $n; $j++ ) {
          $param_count++;
          my $p;
          my $expression = "\$p = \$p$j";
          eval($expression);

          $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$j-1] = [$p,$param_count];
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
            if ( $s =~ /^ *([0-9.e+-]{5,9})5( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1}( *([0-9.e+-]{5,9})5){0,1} *$/i ) {
              my ($p1, $p2, $p3, $p4, $p5, $p6, $p7) =
                ($1, $3, $5, $7, $9, $11, $13);
              my $num_to_read = ($rest >=7) ? 7 : $rest;

              my $par_num = $num_param - $rest;

              for ( my $j = 1; $j <= $num_to_read; $j++ ) {
                $param_count++;
                $par_num++;
                my $p = 0;
                my $expression = "\$p = \$p$j";
                eval($expression);

                $CI[$par1-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$j+4] = [$p,$param_count];
              }
              $rest -= $num_to_read;
            } else {
              die "Format error in ING11 in CI parameters for configs $conf1 and $conf2";
            }
          }
        }
      } else {
        # End of parameters section; stop here
        close ING11;
        return;
      }
    }
  }
  close ING11;
  return;
} ##read_ING11_params()

############################################################################
sub read_OUTG11_params($$) {   #4/25/2005 2:33PM A.Kramida
############################################################################
  # Read Slater parameter labels and values
  my $par = shift; # parity number (1 or 2)
  my $s = shift; # Last string read from OUTG11
  if ( !$s ) {
    # To read the first parity, skip to the parameter section in OUTG11;
    # Otherwise, we are already at the beginning of the second-parity parameter section in OUTG11

    while ( (defined($s = <OUTG11>)) && ($s !~ /^.{7}([^-]{12}) *PARAMETER VALUES IN +([0-9.]+) /) ) {
      next;
    }
  }
  my $spectrum = '';
  do {
    $s =~ /^ (.{6})([^-]{12}) *PARAMETER VALUES IN +([0-9.]+) /;
    $spectrum = $1 unless $spectrum;

    #PARAMETER VALUES IN  1000.0 CM-1 (HR TIMES 0.85)
    if ( $s =~ /PARAMETER VALUES IN[^(]+\(\S+ TIMES ([.0-9]+)\)/i ) {
      $param_scaling[4] = $1; # Set the CI scaling factor
    }

    my ($conf, $units) = ($2, $3);

    $spectrum =~ s/^\s+|\s+$//g;
    $conf =~ s/^\s+|\s+$//g;

    my ($par1,$nc) = ();
    if ( defined($conf_nums{$conf}->[0]) ) {
      ($par1,$nc) = @{$conf_nums{$conf}}; # Get the number of this config and its parity
    }
    die "Config. name $conf not recognized in OUTG11" unless $nc;
    die "Parity mismatch between OUTG11 and ING11 for config $conf" unless ($par1 == $par);

    my $num_param = $params[$par-1]->{$nc}->{'np'};

    $s = <OUTG11>;
    $s = <OUTG11>;

    chomp $s;
    $s .= '               ';
    if ( $s =~ /^.{39}(EAV)   (.{15}){0,1}(.{15}){0,1}(.{15}){0,1}(.{15}){0,1} *$/  ) {
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
        chomp $s;
        $s .= '               ';
        if ( $s =~ /^.{7}(.{15})(.{15}){0,1}(.{15}){0,1}(.{15}){0,1}(.{15}){0,1}(.{15}){0,1}(.{15}){0,1} *$/  ) {
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
    while ( (defined($s = <OUTG11>)) && ($s !~ /^ +$spectrum(.{12}) *PARAMETER VALUES IN +([0-9.]+) /i) ) {
      #PARAMETER VALUES IN  1000.0 CM-1 (HR TIMES 0.85)
      if ( $s =~ /PARAMETER VALUES IN[^(]+\(\S+ TIMES ([.0-9]+)\)/i ) {
        $param_scaling[4] = $1; # Set the CI scaling factor
      }

      if ( $s =~ /ENERGY MATRIX|COUPLING|EIGEN|PURITY/ ) {
        # Switch to next parity
        last;
      }

      if ( $s =~ /^ (.{9})-(.{8}). *PARAMETER VALUES IN +([0-9.]+) / ) {
        # Read the CI parameter section
        my ($conf1, $conf2, $units) = ($1, $2, $3);
        $conf1 =~ s/^\s+|\s+$//g;
        $conf2 =~ s/^\s+|\s+$//g;
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

        my $num_param = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'np'};

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
          die "Format error in the main parameter line of CI section of OUTG11 for parity $par, configs $conf1, $conf2";
        }
        if ( $num_param > 5 ) {
          my $rest = $num_param - 5;
          my $n_lines = 0;
          { use integer;
            $n_lines = ($rest - 1) / 7 + 1;
          }
          for ( my $n_line = 1; $n_line <= $n_lines; $n_line++ ) {
            $s = <OUTG11>;
            #if ( $s =~ /^( +[A-Za-z0-9*]{8})( +[A-Za-z0-9*]{8}){0,1}( +[A-Za-z0-9*]{8}){0,1}( +[A-Za-z0-9*]{8}){0,1}( +[A-Za-z0-9*]{8}){0,1}( +[A-Za-z0-9*]{8}){0,1}( +[A-Za-z0-9*]{8}){0,1} *$/  ) {
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
      # End of parity
      return $s;
    }
  } while ( 1 ) ;
  if ( $s =~ /ENERGY MATRIX|COUPLING|EIGEN|PURITY/ ) {
    # End of parity
    return $s;
  }
} ##read_OUTG11_params()

############################################################################
sub Reorder_Params() {  #4/25/2005 8:00AM A.Kramida
############################################################################
  # Reorder parameter values and rename them according to their type and shell names

  # Reorder parameters
  for ( my $par = 1; $par <= 2; $par++ ) {
    # Reorder the Slater parameters for each config
    foreach my $nc ( sort {$a<=>$b} keys %{$params[$par-1]} ) {
      my $conf_name = $confs[$par-1]->{$nc};
      if ($conf_name eq 'dp8p') {
          $nc=$nc;
      }

      my $num_param = $params[$par-1]->{$nc}->{'np'};
      my %new_params = ();
      my %prev_shell = ('FG' => 0, 'ABG' => 0, 'ABGT' => {});
      my $prev_type = '';
      for ( my $j = 1; $j <= $num_param; $j++ ) {
        my ($p, $t, $par_name, $par_ind) = @{$params[$par-1]->{$nc}->{"p$j"}};
        if ($j == 3) {
          $j=$j;
        }

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
        if ( $t == 0 ) {
          $new_par_name = "Eav\($conf_name\)";
        } elsif ( $t == 1 ) {
          if ( $par_name =~ /^([FG])(\d)\((\d){2}\)/ ) {
            my ($FG, $rank, $n_shell) = ($1,$2,$3);
            if ( ($n_shell > 1) && ($n_shell != $prev_shell{'FG'} + 1) ) {
              $prev_shell{'FG'} = $n_shell - 1;
            }
            my $shell_name = $shells[$par-1]->{$nc}->[$n_shell-1]->[0];
            $key .= "$n_shell$FG$rank";
            $new_par_name = "$FG$rank\($shell_name,$shell_name\)";
            if ( $prev_type =~ /A|B|G|T/ ) {
              $prev_type = 'FG';
              my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell-1]};
              $nL =~ s/^\d+//; # strip the shell's number
              $prev_shell{'ABGT'}->{"$nL$occup"} = $n_shell-1;
              $prev_shell{'ABG'} = $n_shell-1;
            }

          } elsif ( $par_name =~ /^T(\d){0,1}\(([SPDFGH]) (\d)\)/ ) {
            my ($rank, $L, $occ) = ($1 + 0, $2, $3);
            $L = lc($L);
            if ( $prev_type =~ /^T(\d{0,1})\(([SPDFGH]) (\d)\)/ ) {
              my ($r1,$L1,$occ1) = ($1+0,$2,$3);
              if ( ($r1 > $rank) || (uc($L1) ne uc($L)) || ($occ1 != $occ) ) {
                $prev_shell{'ABGT'}->{"$L1$occ1"} = $prev_shell{'ABGT'}->{"$L1$occ1"} ? $prev_shell{'ABGT'}->{"$L1$occ1"} + 1 : 0;
                $prev_shell{'ABG'} = $prev_shell{'ABGT'}->{"$L1$occ1"};
              }
            }
            $prev_type = $par_name;
            # Find this shell's number
            my $ns = 0;
            for ( my $n_shell = $prev_shell{'ABGT'}->{"$L$occ"} + 0; $n_shell <= 7; $n_shell++ ) {
              my ($nL, $occup) = @{$shells[$par-1]->{$nc}->[$n_shell]};
              $nL =~ s/^\d+//; # strip the shell's number
              if ( ($nL eq $L) && ($occup == $occ) ) {
                $ns = $n_shell + 1;
                last;
              }
            }
            die "Shell designation is not recognized in T parameter for config. $conf_name of parity $par in OUTG11" unless $ns;

            if ( ($ns > 1) && ($ns != $prev_shell{'ABGT'}->{"$L$occ"} + 1) ) {
              $prev_shell{'ABGT'}->{"$L$occ"} = $ns - 1;
              $prev_shell{'ABG'} = $ns - 1;
            }
            my $n_shell = $ns;

            my $shell_name = $shells[$par-1]->{$nc}->[$n_shell-1]->[0];
            $key .= "$n_shell" . "T$rank";
            $new_par_name = ($rank) ? "T$rank\($shell_name\)" : "T\($shell_name\)";

          } else {
            # For ALPHA, BETA, GAMMA :
            # Find the first non-processed shell with equivalent electrons
            $par_name =~ /^(.)/;
            my $cur_type = $1; # The first letter of the parameter name
            if ( ($prev_type =~ /A|B|G|T/) && ($prev_type ge $cur_type) ) {
              #$prev_shell{'ABG'}++;
            }
            $prev_type = $cur_type;
            my $ns = 0;
            for ( my $n_shell = $prev_shell{'ABG'}; $n_shell <= 7; $n_shell++ ) {
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
              if ( ($nL =~ /ghiklmno/ ) && ($occup > 1) ) {
                $ns = $n_shell + 1;
                last;
              }
            }
            if ( !$ns ) {
              die "Could not find shell with equivalent electrons for $par_name parameter\nfor config. $conf_name of parity $par in OUTG11" unless $ns;
            }
            if ( ($ns > 1) && ($ns != $prev_shell{'ABG'} + 1) ) {
              $prev_shell{'ABG'} = $ns - 1;
            }
            my $n_shell = $ns;
            $key .= "$n_shell" . "H$par_name";
            my $shell_name = $shells[$par-1]->{$nc}->[$n_shell-1]->[0];
            my %ABG = ('A'=>'ALPHA', 'B' => 'BETA', 'G' => 'GAMMA');
            $new_par_name = $ABG{$cur_type} . "\($shell_name\)";
          }
        } elsif ($t == 2) {
          # ZETA l
          my $ns = 0;
          if ( $par_name =~ /ZETA (\d)/ ) {
            $ns = $1;
          }
          #if ( !$ns ) {
          #  $ns = $ns;
          #}
          unless ($ns) {
            die "Shell number not recognized for ZETA parameter of config. $conf_name of parity $par in OUTG11";
          }

          my $n_shell = $ns;
          $key .= "ZETA$n_shell";
          my $shell_name = $shells[$par-1]->{$nc}->[$n_shell-1]->[0];
          $new_par_name = "ZETA\($shell_name\)";
        } elsif ( ($t == 3) || ($t == 4) ) {
          # Fn, Gn(l',l")
          if ( $par_name =~ /^([FG])(\d)\((\d)(\d)\)/ ) {
            my ($FG, $rank, $n1, $n2) = ($1, $2, $3, $4);
            my $new_n1 = $n1;
            my $new_n2 = $n2;
            if ( $new_n2 < $new_n1 ) {
              # Exchange $new_n2 and $new_n1 so that they go in increasing order
              my $n = $new_n1;
              $new_n1 = $new_n2;
              $new_n2 = $n;
            }
            $key .= "$new_n1$new_n2$rank";
            $new_par_name = "$FG$rank\($new_n1$new_n2\)";
            my $shell_name1 = $shells[$par-1]->{$nc}->[$new_n1-1]->[0];
            my $shell_name2 = $shells[$par-1]->{$nc}->[$new_n2-1]->[0];
            $new_par_name = "$FG$rank\($shell_name1,$shell_name2\)";
          }
        } elsif ($t != 0) {
          die "Unrecognized parameter type for parameter $par_name of config $conf_name of parity $par in ING11";
        }
        # Now we have the new sort key for each parameter.
        $new_params{$key} = [$p, $t, $new_par_name, $par_ind];
      }
      # Substitute the old parameter values with new ones
      my $par_num = 0;
      foreach my $key (sort {$a cmp $b} keys %new_params) {
        $par_num++;
        $params[$par-1]->{$nc}->{"p$par_num"} = $new_params{$key};
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
        next unless defined($CI[$par-1]->[$nc1-1]->[$nc2-1]->{'np'});
        my $num_param = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'np'};
        for ( my $i = 1; $i<=$num_param; $i++ ) {
          my ($p,$p_num) = @{$CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$i-1]};
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

            #my ($class1,$transpose1) = &CI_class($par, $nc1, $nc2, $n1, $n2, $n3, $n4);

            my $new_n1 = $n1-1;
            my $new_n2 = $n2-1;
            my $new_n3 = $n3-1;
            my $new_n4 = $n4-1;

            #my ($class2,$transpose2) = &CI_class($par, $nc1, $nc2, $new_n1+1, $new_n2+1, $new_n3+1, $new_n4+1);

            my ($nL1, $occup1) = @{$shells[$par-1]->{$nc1}->[$n1-1]};
            my ($nL2, $occup2) = @{$shells[$par-1]->{$nc1}->[$n2-1]};
            my ($nL3, $occup3) = @{$shells[$par-1]->{$nc2}->[$n3-1]};
            my ($nL4, $occup4) = @{$shells[$par-1]->{$nc2}->[$n4-1]};
            my ($L1,$L2,$L3,$L4) = ($nL1,$nL2,$nL3,$nL4);
            my ($new_nL1,$new_nL2,$new_nL3,$new_nL4) = ($nL1,$nL2,$nL3,$nL4);
            $L1 =~ s/^\d+//g;
            $L2 =~ s/^\d+//g;
            $L3 =~ s/^\d+//g;
            $L4 =~ s/^\d+//g;

            my $old_DE = $DE;
            if ( $new_n2 < $new_n1 ) {
              # Exchange $new_n2 and $new_n1 so that they go in increasing order
              my $n = $new_n1;
              $new_n1 = $new_n2;
              $new_n2 = $n;
              $new_nL1 = $nL2;
              $new_nL2 = $nL1;
              #if ( ($class1 =~/^[6789]$|^10$/ ) ) {
              #  $DE = ($DE eq 'E') ? 'D' : 'E';
              #}
            }
            if ( $new_n4 < $new_n3 ) {
              # Exchange $new_n4 and $new_n3 so that they go in increasing order
              my $n = $new_n3;
              $new_n3 = $new_n4;
              $new_n4 = $n;
              $new_nL3 = $nL4;
              $new_nL4 = $nL3;
              #if ( ($class1 =~/^[6789]$|^10$/) ) {
              #  $DE = ($DE eq 'E') ? 'D' : 'E';
              #}
            }

            my $rank_par = 2 - $rank % 2;
            my $key = "$new_n1$new_n2$new_n3$new_n4$DE$rank";
            my $new_par_name = "R$DE$rank\($new_nL1$new_nL2,$new_nL3$new_nL4\)";
            $new_params{$key} = [$p,$new_par_name];

          } else {
            die "Wrong format for CI parameter name $par_name in OUTG11, configs $conf1 and $conf2";
          }
        }

        # Store the reordered parameters in the same $CI hash
        my $par_num = 0;
        foreach my $key ( sort {$a cmp $b} keys %new_params ) {
          $par_num++;
          my ($p,$p_num) = @{$CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$par_num-1]};
          $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$par_num-1] = [$new_params{$key},$p_num];
        }
      }
    }
  }

} ##Reorder_Params


############################################################################
sub write_params {   #4/25/2005 1:18PM A.Kramida
############################################################################
  print PARAM_FILE "parity\tConfiguration\t\tparameter\tLSF\tUncert.\tGroup\tHF\tLSF/HF\n";

  my %scale = ('F'=>0, 'Z'=>1, 'f'=>2, 'G'=>3, 'R'=>4);
  for ( my $par = 1; $par <= 2; $par++ ) {
    my $par_code = $parities[$par-1];
    my $scaling = '';
    my $group_cnt = $param_group_counts[$par-1]; # Hash of counts of params in each linked group
    my $max_cyc = $max_cyc_no[$par-1]; # Number of iterations in RCE, 0 if not fitted

    # Write the Slater parameters for each config
    foreach my $nc ( sort {$a<=>$b} keys %{$params[$par-1]} ) {
      next unless ($nc+0);
      my $conf_name = $confs[$par-1]->{$nc};
      my $np = $params[$par-1]->{$nc}->{'np'};
      for ( my $j = 1; $j <= $np; $j++ ) {
        my ($p, $t,$par_name,$par_num) = ('0','0','',0);
        ($p, $t, $par_name, $par_num) = @{$params[$par-1]->{$nc}->{"p$j"}};
        if ( $p !~ /[.]/ ) {
          $p /= 10000;
        }

        if ( ($p==0) && defined($params[$par-1]->{'sq'}->[$par_num-1]) && $params[$par-1]->{'sq'}->[$par_num-1]->[1] == -100 ) {
          # Skip "illegal" or effective parameters (alpfa, beta, ...) if they were not used in the fitting
          #next;
        }
        my ($group, $sdx, $hf, $hf_ratio, $fixed) = ('','','','','');
        if ( $max_cyc && defined($params[$par-1]->{'sq'}->[$par_num-1]) ) {
          # Fitting was done
          my ($par_name, $flag, $value, $sdx1) = @{$params[$par-1]->{'sq'}->[$par_num-1]};
          # $sdx = abs($sdx1) if (($flag+0) && (abs($flag) < 99));
          $sdx = abs($sdx1);
					if (abs($flag) >= 99) {
						$fixed = 'fixed';
					}
          if ( $group_cnt->{$flag} > 1 ) {
            $group = abs($flag);
          }
        }
        if ( defined($params[$par-1]->{'HF'}->[$par_num-1]) ) {
          my ($par_name, $flag, $value) = @{$params[$par-1]->{'HF'}->[$par_num-1]};
          if ( $par_name =~ /^([FGZ])/) {
            my $par_type = $1;
            if ( $par_type eq 'F' ) {
              if ( ($par_name =~ /\((\d)(\d)\)/) && ($1 != $2) ) {
                $par_type = 'f';
              }
            }
            my $scale_factor = $param_scaling[$scale{$par_type}];
            $value /= $scale_factor;
          }
          $hf = $value;
          if ( $hf+0 ) {
            $hf_ratio = $p/$hf;
          }
        }
        foreach ($p, $hf, $sdx) {
          next unless (($_+0) || (/[.]/)); # Scale by 1000 and round to one decimal place if non-zero or has a decimal point
          $_ = sprintf("%18.3f",$_*1000);
          s/^\s+|\s+$//g;
        }
        $hf_ratio = ($hf+0 ? $p/$hf : '') if $hf_ratio;  # relcalculate from rounded values
        $hf_ratio = sprintf("%15.8f",$hf_ratio) if $hf_ratio;
        $hf_ratio =~ s/^\s+|\s+$//g;
        # $sdx = 'fixed' if (!$sdx && $max_cyc);
        $sdx = $fixed if ($fixed && $max_cyc);

        #my $cn = ($j == 1) ? &get_conf_name($conf_name) : '';
        my $cn = &get_conf_name($conf_name);
        #if ( ($par_name =~ /^EAV/i) && ($cn eq '1s.2s.6d') ) {
        #  $p = $p;
        #}
        $par_name = 'Eav' if $par_name =~ /^EAV/i;
        $par_name = lc($par_name) if $par_name =~ /ZETA/;
        print PARAM_FILE "$par_code\t$cn\t\t$par_name\t$p\t$sdx\t$group\t$hf\t$hf_ratio\n";
      }
    }
    # Write the CI section for each parity
    my $num_c1 = $#{$CI[$par-1]} +  1;
    for (my $nc1=1; $nc1 <= $num_c1; $nc1++ ) {
      my $conf1 = $confs[$par-1]->{$nc1};
      my $num_c2 = $#{$CI[$par-1]->[$nc1-1]} +  1;
      for (my $nc2=1; $nc2 <= $num_c2; $nc2++ ) {
        my $conf2 = $confs[$par-1]->{$nc2};
        next unless defined($CI[$par-1]->[$nc1-1]->[$nc2-1]->{'np'});
        my $np = $CI[$par-1]->[$nc1-1]->[$nc2-1]->{'np'};
        for ( my $i = 1; $i<=$np; $i++ ) {
          my ($p_hash,$par_num) = @{$CI[$par-1]->[$nc1-1]->[$nc2-1]->{'params'}->[$i-1]};
          my ($p,$par_name) = @{$p_hash};

          my ($group, $sdx, $hf, $hf_ratio) = ('','','','');
          if ( $max_cyc && defined($params[$par-1]->{'sq'}->[$par_num-1]) ) {
            # Fitting was done
            my ($par_name, $flag, $value, $sdx1) = @{$params[$par-1]->{'sq'}->[$par_num-1]};
            $sdx = abs($sdx1) if (($flag+0) && (abs($flag) < 99));
            if ( $group_cnt->{$flag} > 1 ) {
              $group = abs($flag);
            }
          }
          if ( defined($params[$par-1]->{'HF'}->[$par_num-1]) ) {
            my ($par_name, $flag, $value) = @{$params[$par-1]->{'HF'}->[$par_num-1]};
            my $scale_factor = $param_scaling[4];
            $value /= $scale_factor;

            $hf = $value;
            if ( $hf+0 ) {
              $hf_ratio = $p/$hf;
            }
          }
          foreach ($p, $hf, $sdx) {
            next unless (($_+0) || (/[.]/)); # Do it if non-zero or has a decimal point
            $_ = sprintf("%15.1f",$_*1000);
            s/^\s+|\s+$//g;
          }
          $hf_ratio = sprintf("%15.8f",$hf_ratio) if $hf_ratio;
          $hf_ratio =~ s/^\s+|\s+$//g;
          $sdx = 'fixed' if (!$sdx && $max_cyc);

          #my $cn1 = ($i == 1) ? &get_conf_name($conf1) : '';
          #my $cn2 = ($i == 1) ? &get_conf_name($conf2) : '';
          my $cn1 = &get_conf_name($conf1);
          my $cn2 = &get_conf_name($conf2);
          print PARAM_FILE "$par_code\t$cn1\t$cn2\t$par_name\t$p\t$sdx\t$group\t$hf\t$hf_ratio\n";
        }
      }
    }
  }
  close PARAM_FILE;
} ##write_params

############################################################################
sub fix_trailing_index($$$$) {   #3/30/2006 11:48AM A.Kramida
############################################################################
  my ($parity, $J, $basis_state_num, $term1) = @_;
  my %complete_df = ('d' => 10, 'f' => 14);
  if ( $term1 =~ /([0-9])$/ ) {
    my $seniority = $1;
    my $sh = $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'sh'};
    #if ( $sh =~ /4f3.\(2H/ ) {
    #  $sh = $sh;
    #}
    foreach my $df ('d','f') {
      if ( $sh !~ /$df(\d+)/ ) {
        next;
      }
      my $occup = 0;
      my ($df_shell, $df_replace) = ('','');
      my $complete_shell = $complete_df{$df};
      for ( my $i = 3; $i<= $complete_shell-3; $i++ ) {
        #if ( ($sh =~ /$df$i\./) && ($sh =~ /d7/) ) {
        #  $sh = $sh;
        #}
        #if ( ($sh =~ /$df$i\.\((..)\)/) && ($i > $occup) ) {
        if ( ($sh =~ /$df$i\.\((\d[A-Z][*]*)\)/) && ($i > $occup) ) {
          # Incomplete d or f shell followed by another incomplete shell
          $occup = $i;
          $df_shell = "$df$i\.\\($1\\)";
          $df_replace = "$df$i.($1$seniority)";
          $df_shell =~ s/[*]/\\*/g;
        } elsif ( ($sh =~ /$df$i\.<(\d[A-Z][*]*)>/) && ($i > $occup) ) {
          # Incomplete d or f shell followed by another incomplete shell
          $occup = $i;
          $df_shell = "$df$i\.\<$1\>";
          $df_replace = "$df$i.($1$seniority)";
          $df_shell =~ s/[*]/\\*/g;
        } elsif ( ($sh =~ /$df$i( |.<)(\d[A-Z][*]*)(>{0,1})\t/) && ($i > $occup) ) {
          # Incomplete d or f shell which is the last open shell
          $occup = $i;
          $df_shell = "$df$i$1$2$3";
          $df_replace = "$df$i$1$2$seniority$3";
          $df_shell =~ s/[*]/\\*/g;
        } elsif ( ($sh =~ /$df$i\.(\d+s2|\d+p6|\d+d10|\d+f14)(\S*)( |.<)(\d[A-Z][*]*)(>{0,1})\t/) && ($i > $occup) ) {
          # Incomplete d or f shell followed by a complete shell
          $occup = $i;
          $df_shell = "$df$i.$1$2$3$4$5";
          my ($last_sh, $fin_term,$end) = ("$1$2","$3$4",$5);
          $df_shell =~ s/\(/\\(/g;
          $df_shell =~ s/\)/\\)/g;
          $df_shell =~ s/\./\\./g;
          $df_shell =~ s/[*]/\\*/g;
          $df_replace = "$df$i.$last_sh$fin_term$seniority$end";
        } elsif ( ($sh =~ /$df$i\t([^\t]+)$/) && ($i > $occup) ) {
          # Incomplete d or f shell which is the last open shell, without intermediate term
          $occup = $i;
          $df_shell = "$df$i\t$1";
          my $fin_term = $1;
          if ( $fin_term =~ /^(.+)([*]){0,1}$/ ) {
            $fin_term = $1;
            my $p = $2;
            $df_replace = "$df$i\t$fin_term$seniority$p";
          }
        } elsif ( ($sh =~ /$df$i\.(\d+s2|\d+p6|\d+d10|\d+f14)(\S*)\t([^\t*]+)([*]{0,1})$/) && ($i > $occup) ) {
          # Incomplete d or f shell followed by a complete shell, without intermediate term
          $occup = $i;
          $df_shell = "$df$i.$1$2\t$3$4";
          my ($last_sh, $fin_term,$end) = ("$1$2",$3,$4);
          $df_shell =~ s/\(/\\(/g;
          $df_shell =~ s/\)/\\)/g;
          $df_shell =~ s/\./\\./g;
          $df_shell =~ s/\*/\\*/g;
          $df_replace = "$df$i.$last_sh\t$fin_term$seniority$end";
        }
      }
      if ( $occup > 0 ) {
        $df_replace =~ s/(\*)(\d)/$2$1/;   # Interchange, e.g., '*1' to make it '1*'
        $sh =~ s/$df_shell/$df_replace/;
        $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'sh'} = $sh;
      }
    }
  }
} ##fix_trailing_index($$$$)

############################################################################
sub ReadOUTE() {   #6/12/2006 5:06PM A.Kramida
############################################################################
  my $file = 'OUTE';
  if ( !-f $file ) {
    return;
  }
  open(OUTE, "<$file") || return;

  print "Reading OUTE...\n";
  my $parity = 0;
  my $cyc_no = -1;
  my $param_count = 0;
  my $s;
  while ( (defined($s = <OUTE>)) && ($s !~ /CYC NO=/) ) {
    next;
  }
  if ( !defined($s) && $parity) {
    $max_cyc_no[$parity-1] = $cyc_no;
  }

  while ( $s =~ /CYC NO= *(\d+) / ) {
    my $cur_cyc = $1;
    if ( $cur_cyc <= $cyc_no ) {
      if ( $parity ) {
        $max_cyc_no[$parity-1] = $cyc_no;
      }
      $parity++;
    } else {
      @param_group_counts[$parity-1] = {}; # Reset all group counts to zeros for each cycle
    }
    $parity++ unless $parity;
    $cyc_no = $cur_cyc;
    $param_count = 0;  # Reset the parameter count to zero in the beginning of each cycle

    # Skip three lines
    <OUTE>;
    <OUTE>;
    <OUTE>;
    # Read parameter names, flags, values, and sdx (standard deviation)
    while ( (defined($s = <OUTE>)) && ($s !~ /Iteration .+finished/i) ) {
      if ( $s =~ /(.{15})(.{4})(.{13}).{39}(.{13})/ ) {
        my ($par_name, $flag, $value, $sdx) = ($1, $2, $3, $4);
        $par_name =~ s/^\s+|\s+$//g;
        $flag =~ s/^\s+|\s+$//g;
        $value =~ s/^\s+|\s+$//g;
        $sdx =~ s/^\s+|\s+$//g;
        $param_count++;
        $params[$parity-1]->{'sq'} = [] unless defined $params[$parity-1]->{'sq'};
        $params[$parity-1]->{'sq'}->[$param_count-1] = [$par_name, $flag, $value, $sdx];
        if ( ($flag +0) && (abs($flag) < 99) ) {
          #if ( !defined(@param_group_counts[$parity-1]) ) {
          #  $parity = $parity;
          #}
          @param_group_counts[$parity-1]->{$flag} = 0 unless defined(@param_group_counts[$parity-1]->{$flag});
          @param_group_counts[$parity-1]->{$flag}++;
        }
      }

      next;
    }
    while ( (defined($s = <OUTE>)) && ($s !~ /CYC NO=/) ) {
      next;
    }
    if ( !defined($s) && $parity) {
      $max_cyc_no[$parity-1] = $cyc_no;
      last;
    }
  }

  close OUTE;
  return;
} ##ReadOUTE()

############################################################################
sub ReadHF() {   #6/12/2006 5:06PM A.Kramida
############################################################################
  my $file = 'RCEINP.HF';
  if ( !-f $file ) {
    return;
  }
  open(HF, "<$file") || return;

  print "Reading HF parameter values...\n";
  my $parity = 0;
  my $param_count = 0;
  my $s;
  while ( (defined($s = <HF>)) && ($s !~ /PARAMETER/) ) {
    next;
  }

  while ( $s =~ /PARAMETER/ ) {
    $parity++;
    $param_count = 0;

    # Read parameter names, flags, and values
    while ( (defined($s = <HF>)) && ($s !~ /^    2/) ) {
      if ( $s =~ /(.{11})(.{4})(.{14})/ ) {
        my ($par_name, $flag, $value) = ($1, $2, $3);
        $par_name =~ s/^\s+|\s+$//g;
        $flag =~ s/^\s+|\s+$//g;
        $value =~ s/^\s+|\s+$//g;
        $param_count++;
        $params[$parity-1]->{'HF'} = [] unless defined $params[$parity-1]->{'sq'};
        $params[$parity-1]->{'HF'}->[$param_count-1] = [$par_name, $flag, $value];
      }

      next;
    }
    while ( (defined($s = <HF>)) && ($s !~ /PARAMETER/) ) {
      next;
    }
  }

  close HF;
  return;
} ##ReadHF()

############################################################################
sub get_electron_parity($) {   #03/09/2010 7:30AM
############################################################################
  my $ell = shift;
  return ($ell =~ /[pfhkmoruwy]/) ? -1 : 0;
} ##get_electron_parity($)

############################################################################
sub get_shell_parity($$) {   #03/09/2010 7:29AM
############################################################################
  my ($sh,$occup) = @_;
  my $par = ($occup % 2) ? get_electron_parity($sh) : 0;
  return $par;
} ##get_shell_parity($$)

############################################################################
sub get_parity_char($) {   #03/09/2010 7:37AM
############################################################################
  my $par = shift;
  return $par ? '*' : '';
} ##get_parity_char($)

############################################################################
sub complete_shells($) {    #04/17/2013 12:43PM
############################################################################
  my $sh = shift;
  my $n = 2*(2*$L_moment{$sh} + 1);
  return $n;
} ##complete_shells($)

############################################################################
sub get_lowest_complete_shell() {   #11/21/2013 11:42AM
############################################################################
  my $lowest_complete_shell = '';
  for ( my $par = 1; $par <= 2; $par++ ) {
    last if $lowest_complete_shell;
    foreach my $n_conf (keys %{$shells[$par-1]} ) {
      my ($sh,$occup) = @{$shells[$par-1]->{$n_conf}->[$shell_ord_seq{1} - 1]};
      my $L = $sh;
      $L =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
      if ( $occup == &complete_shells($L) ) {
        $lowest_complete_shell = $sh;
        last;
      }
    }
  }
  return $lowest_complete_shell;
} ##get_lowest_complete_shell()

############################################################################
sub init_vars() {    #11/22/2013 6:44AM
############################################################################
  # Initialize global variables
  @confs = ({},{});
  %confs36 = ();
  @shells = ({},{});
  @params = ({},{}); # Parameters for the two config. sets
  @sdx = ({},{});    # Standard deviations of params for the two config. sets
  @param_flags = ({},{});   # Parameter flags for the two config. sets
  @param_group_counts = ({},{});   # Number of linked parameters in each group for the two config. sets
  @param_scaling = (0.85,1.00,0.85,0.85,0.85);
  @max_cyc_no = (0,0);    # Max. cycle numner for the two config. sets
  @CI = ([],[]);
  @basis = ({},{});
  @energies = ({},{});
  @vectors = ({},{});
  @term_labels_RCG = ({},{});
  @term_labels_map = ({},{});
  @parities = ('','');
  @RCE_lev = ({},{});
  @map_RCG_RCE = ({},{});
  @Emax = ();
  @num_states_in_blocks = ();
  %shell_order_in36 = ();
  %shell_order_in36_back = ();
  %shell_ord_seq = ();
  %conf_nums = ();
  %L_moment = ('s' => 0, 'p' => 1, 'd' => 2, 'f' => 3, 'g' => 4, 'h' => 5,
                  'i' => 6, 'k' => 7, 'l' => 8, 'm' => 9, 'n' =>10, 'o' =>11,
                  'q' =>12, 'r' =>13, 't' =>14, 'u' =>15, 'v' =>16, 'w' =>17,
                  'x' =>18, 'y' =>19, 'z' =>20);
  %no_genealogy_shells =
    ( 's'=>1, 's.p2'=>1, 's.p4'=>1, 's.p' => 1, 'p.s' => 1, 's.s' => 1,
      's.d2'=>1, 's.d8'=>1, 's.s.d2'=>1, 's.s.d8'=>1,
      'p5' => 1, 'd9' => 1, 'f13' => 1, 'p' => 1, 'd' => 1, 'f' => 1, 'g' => 1, 'h' => 1, 'i' => 1
    );

  # Initialize the basis hashes
  $basis[0]->{'LS'} = {};
  $basis[0]->{'JJ'} = {};
  $basis[1]->{'LS'} = {};
  $basis[1]->{'JJ'} = {};

  $vectors[0]->{'LS'} = {};
  $vectors[0]->{'JJ'} = {};
  $vectors[1]->{'LS'} = {};
  $vectors[1]->{'JJ'} = {};

} ##init_vars()

############################################################################
sub read_RCG_options() {    #11/22/2013 6:48AM
############################################################################
  # Read RCG options
  my $s;
  while ( (defined($s = <OUTG11>)) && ($s !~ /0RCG MOD 11 +(\S+) COUPLING /) ) {
    next;
  }
  if ( $s !~ /0RCG MOD 11 +(\S+) COUPLING / ) {
    die "Coupling scheme line not found at the top of OUTG11 file"
  } else {
    my $first_coupling = $1;
    if ( $first_coupling ne 'LS' ) {
      die "First coupling must be LS for this program to run. Please re-run RCG with coupling = LS";
    }
  }

  while ( (defined($s = <OUTG11>)) && ($s !~ / PRINT=[0-9 ]{5}/) ) {
    next;
  }

  if ( $s !~ / PRINT=\d\d\d\d\d / ) {
    die "Options line not found at the top of OUTG11 file";
  }

  #IV=0   119 001111101 99.0 99.0  PRINT=00010 500 2 0 0 0
  #IV=0   119 001101101 99.0 99.0  PRINT=00010   1 2 0 0 0

  if ( $s =~ / (\d)(\d)(\d)(\d)(\d)(\d)(\d)\d(\d) \d\d[.]\d \d\d[.]\d  PRINT=[0-9 ]{3}([0-9 ])[0-9 ]{4}([0-9 ])/ ) {
    my ($print_LS, $print_JJ, $print_cpl3, $print_cpl4, $print_cpl5, $print_cpl6, $print_cpl7, $print_matr, $print_intermed_numbers, $print_matr1) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10);
    if ( ($print_LS+0 != 0) || ($print_JJ != 0)
       || ($print_cpl3 != 1) || ($print_cpl4 != 1) || ($print_cpl5 != 1) || ($print_cpl6 != 1) || ($print_cpl7 != 1)
       || ($print_intermed_numbers+0 == 0) || ($print_matr !=1) || ($print_matr1 != 0)
    ) {
      die "Please re-run RCG with zeros in cols 31, 32, and 71 and 1 in cols. 33, 34, 35, 36, 37, 39 and 69 on the 1st line of ING11.";
    }
  } else {
    die "Options line format error at the top of OUTG11 file";
  }
} ##read_rcg_options()

############################################################################
sub read_basis() {    #11/22/2013 7:08AM
############################################################################
  #my ($J,$J_prev) = @_;
  my ($J, $J_prev) = (-1,1000);

  &fill_last_shells();
  &fill_open_shells();
  &make_conf_templates();

  # Read the basis state definitions printed by CALCFC in OUTG11
  my $s;
  my $parity = 0;
  while ( $s = <OUTG11> ) {
    if ( $s =~ /STARTING CALCFC/ ) {
      # Switch to next parity
      $parity++;
    }
  #B III 3s              PARAMETER VALUES IN  1000.0 CM-1 (HR TIMES 0.80 1.00 0.80 0.80)     2  S 2  S 0  P 0  S 1  P 0  D 0  S 0  D 0
    if ( $s =~ /PARAMETER VALUES IN[^(]+\(\S+ TIMES ([.0-9]+) ([.0-9]+) ([.0-9]+) ([.0-9]+)\)/ ) {
      # Scaling factors for Fk(il,il), Zeta(il), Fk(il,jl), Gk(il,jl)
      @param_scaling = ($1,$2,$3,$4);
      last; # Stop and commend control to the following section
    }
    if ( $s !~ /^0J= *(\d+[.]\d) +CONFIGURATION/ ) {
      next;
    }
    if ( $s =~ /^0J= *(\d+[.]\d) +CONFIGURATION +(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)/ ) {
      $J = $1;
      $J_prev = $J;

      # Init the basis hashes for this J
      $basis[$parity-1]->{'LS'}->{$J} = {};
      $basis[$parity-1]->{'JJ'}->{$J} = {};

      # Read the rest of the config shell definitions
      while ( defined($s = <OUTG11>) && ($s =~ /^ +(\d+)* *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)/) ) {
        next;
      }
      $s = <OUTG11>;
    }

    # Process shell definition lines
    my $state_num = 0;
    my $n_shell = -1;
    my $last_shell = 0;
    my $L_format = -1;
    $s = <OUTG11> if ($s =~ /^\s*------$/);
    #if ( ($J eq '2.0') && ($parity == 1)  ) {
    #  $J = $J;
    #}
    while ( defined($s = <OUTG11>) && ($s =~ /^(.{4})(.{1,3})  \(.{3}(.{3,4})\)(.{3,4})( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}( {1,4}\(.{3}(.{3,4})\)(.{3,4})){0,1}/) ) {
      # Parse the LS parent terms
      $state_num++;
      my $lead_token = "$1$2";
      my $n_conf;
      my @p = ($3,$4,$6,$7,$9,$10,$12,$13,$15,$16,$18,$19,$21,$22,$24,$25);
      if ( $L_format == -1 ) {
        $L_format = ($s =~ /^.{68}\)/);
        #if ( $lead_token =~ /^ +([\d*]+) *(\d+) *$/ ) {
        #  $L_format = 1;
        #} else {
        #  $L_format = 0;
        #}
      }
      #if ( $L_format ) {
      #  $lead_token =~ /^.{4}(.+)$/;
      #  $n_conf = $1 + 0;
      #} else {
      #  $n_conf = $lead_token + 0;
      #}
      $n_conf = 0;
      if ( $lead_token =~ /^([ \d]{3}) +$/ ) {
        $n_conf = $1 + 0;
      } elsif ( $lead_token =~ /^([ \d]{4}| [*]{3})([ \d]{3})$/ ) {
        $n_conf = $2 + 0;
      }
      if ( !$n_conf ) {
        die "Unrecognized format of CALCFC output in OUTG11";
      }
      #if ( ($parity == 2) && ($J eq '5.0') && ($state_num == 1049) ) {
      #  $J = $J;
      #}
      #$L_format = &fill_LS_shell($parity,$state_num,$J,$2,$3,$4,$6,$7,$9,$10,$12,$13,$15,$16,$18,$19,$21,$22,$24,$25);
      #my $n_conf = ($2 ne ' ' ? $2 : $1) + 0;
      &fill_LS_desig($parity,$state_num,$J,$n_conf,@p);
      if ( !$L_format ) {
          # Now read the JJ classification
        $s = <OUTG11>;
        if ( $s =~ /^\+ +(\d+) +\( *(\d+) +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}/ ) {
          my $c_no  = $basis[$parity-1]->{'LS'}->{$J}->{$1}->{'cn'};
      #    &fill_JJ_shell($parity,$1,$c_no,$J,$3,$4,$5,$7,$8,$9,$11,$12,$13,$15,$16,$17,$19,$20,$21,$23,$24,$25,$27,$28,$29,$31,$32,$33);
          &fill_JJ_desig($parity,$1,$c_no,$J,$3,$4,$5,$7,$8,$9,$11,$12,$13,$15,$16,$17,$19,$20,$21,$23,$24,$25,$27,$28,$29,$31,$32,$33);
        }
      }

    }

    if ( $L_format ) {
      while ( defined($s = <OUTG11>) && ($s =~ /^ +(\d+) +\( *(\d+) +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}( +\( *\d+ +(\d+\S) *(\d+[.]\d)\) *(\d+[.]\d)){0,1}/) ) {
        my $c_no  = $basis[$parity-1]->{'LS'}->{$J}->{$1}->{'cn'};
        #if ( !defined($c_no) ) {
        #  $s = $s;
        #}
        &fill_JJ_desig($parity,$1,$c_no,$J,$3,$4,$5,$7,$8,$9,$11,$12,$13,$15,$16,$17,$19,$20,$21,$23,$24,$25,$27,$28,$29,$31,$32,$33);
      }
    }
  }
  return ($s);
} ##read_basis()

############################################################################
sub read_basis_labels($$) {    #11/22/2013 7:30AM
############################################################################
# Continue to read the OUTG11 file. Find and read the LS basis state labels printed by ENERGY
  my ($s,$fix_trailing) = @_;
  my $parity = 0;
  my ($J, $J_prev) = (-1,1000);
  my $n_conf = '';

  #my $start = 0;
  my $prev_bsn = 0;
  do {
    if ($s =~ /ENERGY MATRIX   \(    LS COUPLING\)       J= *(\d+\.\d) {7}(.{6})(.{12}) {5}CONFIG +(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)/ ) {
      $J = $1;
      if ( $J <= $J_prev ) {
        $parity++;
        if ( (!defined($basis[$parity-1]->{'LS'})) || (!defined($shells[$parity-1]->{'1'})) || (!defined($confs[$parity-1]->{'1'}))) {
          die "Preliminary quantum numbers were not found in OUTG11. Set '1' in column 69 of ING11 and re-run RCG.";
        }
      }
      $J_prev = $J;

      my $spectrum = $2;
      my $conf = $3;
      $n_conf = $4;
      $conf =~ s/^\s+|\s+$//g;
      $prev_bsn = 0;

      # Fill in all config names of this parity
      while ( defined($s = <OUTG11>) && ($s =~ /^.{55}$spectrum(.{12}) +(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+) +(\S) *(\d+)/) ) {
        $conf = $1;
        $n_conf = $2;
        $conf =~ s/^\s+|\s+$//g;
      }
      #$start = 1;
    }

    while ( defined($s = <OUTG11>) && ($s !~ /^ +\d+:/) ) {
      next;
    }

    #if ( ($parity == 2) && ($J eq '5.0') ) {
    #  $parity = 2;
    #}
    # Fill in the LS basis state labels while reading the energy matrix section of OUTG11
    if ( $s =~ /^ +(\d+):(.{12}) (\(.{8}) *(\d+|[*]{3}) /) {
      my ($c_num, $conf1, $term1, $basis_state_num) = ($1,$2,$3,$4);
      $conf1 =~ s/^\s+|\s+$//g;
      $term1 =~ s/^\s+|\s+$//g;
      if ( $basis_state_num eq '***' ) {
        $basis_state_num = $prev_bsn + 1;
      }
      #if ( $basis_state_num == 999 ) {
      #  $basis_state_num = $basis_state_num;
      #}
      $prev_bsn = $basis_state_num;
      $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'label'} = $term1;
      &fix_trailing_index($parity, $J, $basis_state_num, $term1) if $fix_trailing;
      if ( $conf1 ne $confs[$parity-1]->{$c_num} ) {
        my $c_name1 = $confs[$parity-1]->{$c_num};
        die "Error in parity $parity:\nConfig name mismatch for J=$J, state No. $basis_state_num;\nIn energy matrix conf name is $conf1; in header it is $c_name1";
      }
    }
    while ( defined($s = <OUTG11>) && ($s =~ /^ *(\d+):(.{12}) (\(.{8}) *(\d+|[*]{3}) /) ) {
      my ($c_num, $conf1, $term1, $basis_state_num) = ($1,$2,$3,$4);
      $conf1 =~ s/^\s+|\s+$//g;
      $term1 =~ s/^\s+|\s+$//g;
      if ( $basis_state_num eq '***' ) {
        $basis_state_num = $prev_bsn + 1;
      }
      $prev_bsn = $basis_state_num;
      $basis[$parity-1]->{'LS'}->{$J}->{$basis_state_num}->{'label'} = $term1;

      &fix_trailing_index($parity, $J, $basis_state_num, $term1) if $fix_trailing;

      if ( $conf1 ne $confs[$parity-1]->{$c_num} ) {
        my $c_name1 = $confs[$parity-1]->{$c_num};
        die "Error in parity $parity:\nConfig name mismatch for J=$J, state No. $basis_state_num;\nIn energy matrix conf name is $conf1; in header it is $c_name1";
      }
    }

    while ( defined($s = <OUTG11>) && ($s !~ /ENERGY MATRIX   \(    LS COUPLING\)/) ) {

      if ( $s =~ /PARAMETER VALUES/ ) {
        $s = &read_OUTG11_params(2,$s);
        last;
      }

      last if ($s =~ /ELEC DIP SPECTRUM/);

      my $cpl1 = '';
      if ( $s =~ /[01]  EIGENVALUES      \(J=/ ) {
        &read_energies($J,$parity);
        &read_vectors($J, $cpl1, $parity);
      } elsif ( $s =~ /EIGENVECTORS   \( *(\S+) COUPLING/ ) {
        $cpl1 = $1;
        &read_vectors($J, $cpl1,$parity);
      }
      next;
    }
  } while ( defined($s) && ($s !~ /ELEC DIP SPECTRUM/) );
  return $s;
} ##read_basis_labels($$)

############################################################################
sub fill_last_shells() {    #11/22/2013 11:49AM
############################################################################
  @last_shells = ([],[]);
  for ( my $parity = 1; $parity <= 2; $parity++ ) {
    foreach my $nc ( sort {$a<=>$b} keys %{$confs[$parity-1]} ) {
      my $last_shell = &get_last_shell($parity,$nc);  # 0,1,...
      my $last_shell_reordered = &get_last_shell_reordered($parity,$nc); # 0,1,...
      my $last_shell_summation_order = &get_last_shell_summation_order($parity,$nc);  # 0,1,...
      $last_shells[$parity-1]->[$nc-1] = [$last_shell, $last_shell_reordered, $last_shell_summation_order];
    }
  }
} ##fill_last_shells()

############################################################################
sub make_conf_templates() {   #12/04/2013 4:21PM
############################################################################
  @LS_templates = ([],[]);
  for ( my $parity = 1; $parity <= 2; $parity++ ) {
    foreach my $nc ( sort {$a<=>$b} keys %{$confs[$parity-1]} ) {
      my $template = &make_LS_template($parity,$nc);
      $LS_templates[$parity-1]->[$nc-1] = $template;
      $template = &make_JJ_template($parity,$nc);
      $JJ_templates[$parity-1]->[$nc-1] = $template;
    }
  }
} ##make_conf_templates()

############################################################################
sub make_LS_template(@) {  #4/127/2013 01:40PM A.Kramida
############################################################################
  my ($parity,$n_conf) = @_;

  #if ( ($n_conf == 4) && ($parity == 1) ) {
  #  $n_conf = $n_conf;
  #}
  my ($last_shell, $last_shell_reordered, $last_shell_summation_order) = @{$last_shells[$parity-1]->[$n_conf-1]};  # 0,1,...
  my $reordered = ($last_shell_reordered == $last_shell_summation_order ? 0 : 1);
  my $last_shell_in_this_conf = ($reordered ? $last_shell_reordered : $last_shell_summation_order);

  my $fill_str = '';
  my $acc_parity = 0;
  my $first_filled = -1;
  my $final_term = '';
  my ($LS_curr,$LS_accum, $LS_curr_reordered, $LS_accum_reordered) = ('1S','1S','1S','1S');
  my $LS_last = '1S';

  # Determine the accumulated parity after adding each shell in the order of summation given in OUTG11
  my @accum_par = ();
  my $acc_par = 0;
  #for ( my $ns = 0; $ns <= $last_shell; $ns++ ) {
  my $only_complete_shells = 1;
  for ( my $ns = 0; $ns <= $last_shell_summation_order; $ns++ ) {
    my ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    $acc_par++ if $sh_par_char;
    $accum_par[$ns] = ( $acc_par % 2 > 0 ) ? '*' : '';
    my ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    my $L = $sh;
    $L =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
    my $complete_shell = ($occup == &complete_shells($L));
    if ($occup && !$complete_shell) {
      $only_complete_shells = 0;
    }
  }

  my @shells_summation_order = ();
  my $shell_ind = 0;
  my $prev_accum_par_char = '';
  my $accum_par_char = '';
  my $final_term = '';
  my $desig = '';
  my $desig_no_shell_num = '';
  my $desig_has_genealogy = 0;
  my $nL_prev = 0;
  for ( my $ns = 0; $ns <= $last_shell_summation_order; $ns++ ) {
    # Get the new sequential order of shells from the hash
    my $n_shell_display_order = $ns;
    if ( $shell_ord_seq{$ns+1} ) {
      $n_shell_display_order = $shell_ord_seq{$ns+1} - 1;
    }
    $reordered = 1 unless ($ns == $n_shell_display_order);

    ($LS_curr,$LS_accum) = &get_LS_terms($ns);
    my ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$ns]};

    my $L = $sh;
    $L =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
    my $nL = $sh + 0; #$sh*10000 + $L_moment{$L};

    my $complete_shell = ($occup == &complete_shells($L));
    if ( !$occup || ( $complete_shell && !$open_shells{$ns} && ($sh eq $lowest_complete_shell) && !$reordered && 
      (!$only_complete_shells || ($sh eq $lowest_complete_shell)) && ($ns != $last_shell_in_this_conf) ) 
    ) {
      $nL_prev = $nL;
      next;
    }

    $accum_par_char = (($prev_accum_par_char && !$sh_par_char) || (!$prev_accum_par_char && $sh_par_char)) ? '*' : '';

    if ( $complete_shell ) {
      $LS_curr = '';
      $LS_accum = '' unless ( $ns == $last_shell_summation_order );
    } else {
      $LS_curr .= $sh_par_char;
      $LS_accum .= $accum_par_char;
    }
    if ( $ns == $last_shell_summation_order ) {
      $final_term = $LS_accum;
      $LS_accum = '';
    }
    $prev_accum_par_char = $accum_par_char;
    $occup = '' if $occup == 1;
    $desig .= ($desig ? '.' : '');
    my $added_desig = '';
    $desig .= $sh . $occup;
    $added_desig .= $sh . $occup;
    if ( !$complete_shell || ($nL > $nL_prev) || $open_shells{$ns} ) {
      $desig_no_shell_num .= ($desig_no_shell_num ? '.' : '');
      $desig_no_shell_num .= $L . $occup;
    }

    if ( ($shell_ind || ($ns < $last_shell) ) &&
      (&next_shell_requires_genealogy($parity,$n_conf,$ns,$last_shell_summation_order,$desig_no_shell_num)
        && !defined($no_genealogy_shells{$desig_no_shell_num})
        || ($reordered && (($ns == $last_shell) || ($ns == $last_shell_summation_order)) && $desig_has_genealogy)
        || (($ns == $last_shell) && !defined($no_genealogy_shells{$L . $occup}))
      )
    ) {
      $desig .= ".<" . $LS_curr . ">" if $LS_curr;
      $added_desig .= ".<" . $LS_curr . ">" if $LS_curr;
      $desig_has_genealogy = 1;
    }
    if ( $shell_ind && &next_shell_requires_accum_LS($parity,$n_conf,$ns,$last_shell_summation_order) ) {
      $desig .= ".(" . $LS_accum . ")" if $LS_accum;
      $added_desig .= ".(" . $LS_accum . ")" if $LS_accum;
      $desig_has_genealogy = 1;
    }
    my $shell_hash = {'LS_curr' => $LS_curr, 'LS_accum' => $LS_accum, 'name' => $sh . $occup,
      'code_name' => $L . $occup,'occup'=>$occup, 'sh'=>$sh, 'L' => $L, 'desig' => $desig, 'added_desig' => $added_desig};
    $shells_summation_order[$ns] = $shell_hash;
    $shell_ind++;
    $nL_prev = $nL;
  }

  if ( $reordered ) {
    for ( my $ns = 0; $ns <= 7; $ns++ ) {
      next unless ($shell_ord_seq{$ns+1});
      my $shell_ind = $shell_ord_seq{$ns+1} - 1;
      next unless defined $shells_summation_order[$shell_ind];
      my $shell_hash = $shells_summation_order[$shell_ind];
      my $added_desig = $shell_hash->{'added_desig'};
      $fill_str .= ($fill_str ? '.' : '') . $added_desig if $added_desig;
    }
  } else {
    $fill_str = $desig;
  }

  $fill_str .= "\t$final_term";

  # Convert the string template to an array
  my $template = [];
  while ( $fill_str =~ /^([^#]+){0,1}\#(\d+)(.*)$/ ) {
    # Strip the '#' from the term indexes
    push(@{$template}, $1, $2);
    $fill_str = $3;
  }
  # Store the remainder in the last element of the template
  push(@{$template}, $fill_str);
  return $template;
} ##make_LS_template(@)

############################################################################
sub get_LS_terms($) {    #12/04/2013 2:46PM
############################################################################
  my $ns = shift; # The shell number, 0..7
  # Shell terms and intermediate terms after summation with previous shells
  #my @fill = ();
  #for ( my $i = 0; $i <= 7; $i++ ) {
  #  $fill[$i] = [$p[2*$i],$p[2*$i+1]];
  #}
  my @term_numbers = ('#'.(2*$ns+1), '#'.(2*$ns+2));  # Term numbers are increased by 1 to make them non-zero
  return @term_numbers;
} ##get_LS_terms($)

############################################################################
sub make_JJ_template(@) {  #12/5/2013 8:33AM A.Kramida
############################################################################
  my ($parity,$n_conf) = @_;
  my ($last_shell, $last_shell_reordered, $last_shell_summation_order) = @{$last_shells[$parity-1]->[$n_conf-1]};  # 0,1,...
  #if ( ($n_conf == 6) && ($parity == 1) ) {
  #  $n_conf = $n_conf;
  #}

  my $fill_str = '';
  my $acc_parity = 0;
  my $acc_par_char = '';
  my $first_filled = -1;
  my $reordered = 0;
  my $J_last = '0';
  my $J_prev_to_last = '0';
  my $J_last_shell = '0';
  my ($nL_last, $nL_prev) = (0,0);
  my $dot = '';
  for (my $ns = 0; $ns <= 7; $ns++ ) {
    my ($n_shell, $sh, $occup, $sh_par_char) = ();
    if ($ns <= $last_shell) {
      # Get the new sequential order of shells from the hash
      $n_shell = $ns;
      if ( $shell_ord_seq{$ns+1} ) {
        $n_shell = $shell_ord_seq{$ns+1} - 1;
      }
      $reordered = 1 if ($ns != $n_shell);

      ($sh,$occup,$sh_par_char) = @{$shells[$parity-1]->{$n_conf}->[$n_shell]};
      if ( $sh_par_char ) {
        $acc_parity++;
      }
      if ( $acc_parity % 2 > 0 ) {
        $acc_par_char = '*';
      }
    }
    my $L = $sh;
    $L =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code

    #($sh,$occup,$sh_par_char)      = @{$shells[$parity-1]->{$n_conf}->[$n_shell]};
    my ($sh1,$occup1,$sh_par_char1) = @{$shells[$parity-1]->{$n_conf}->[$ns]};
    my $L1 = $sh1;
    $L1 =~ s/[0-9]//g; # Strip the principal quantum number, leaving only the orbital code
    if ( $occup1 && ($occup1 != &complete_shells($L1)) ) {
      $nL_prev = $nL_last;
      $nL_last = $sh1;
      $nL_last =~ s/[a-z]//g;
      $nL_last = $nL_last*10000 + $L_moment{$L1};

      #my $LS = $fill[$ns]->[0];
      #my $J1 = &get_J($fill[$ns]->[1]);
      #my $J2 = &get_J($fill[$ns]->[2]);
      my ($LS, $J1, $J2) = &get_JJ_terms($ns);
      $J_prev_to_last = $J_last;
      $J_last = $J2;  # Not re-ordered !  These are tracked in order to determine the final term.
      $J_last_shell = $J1;
    }
    if ( $ns > $last_shell ) {
      next;
    }
    if ( $occup > 0 ) {
      # Get the J value for the shell and the accumulated J value
      #my $LS_reordered = $fill[$n_shell]->[0];
      #my $J1_reordered = &get_J($fill[$n_shell]->[1]);
      #my $J2_reordered = &get_J($fill[$n_shell]->[2]);
      my ($LS_reordered, $J1_reordered, $J2_reordered) = &get_JJ_terms($n_shell);

      if ( $occup == &complete_shells($L) ) {
        if ( ($fill_str =~ /^(.+)\.\(([^()]+)\)\.\(([^()]+)\)$/) && !$reordered ) {
          # Move the previously given parent term and total J outside the last
          # complete shell
          $fill_str = "$1.$sh$occup.($2).($3)";
          $dot = '.';
        } elsif ( ($ns != $last_shell) && !$open_shells{$ns} ) {
          if ($sh ne $lowest_complete_shell ) {
            # Skip the 1s2 in 1s2.2s, but include the 2s2 in 2s2.3s. Give no intermediate term.
            $dot = '.' if $fill_str;
            $fill_str .= "$dot$sh$occup";
            $dot = '.';
          }
        } elsif ( ($sh ne $lowest_complete_shell) && !$open_shells{$ns} ) {
          $dot = '.' if $fill_str;
          if ( $nL_last == 0 ) {
            $fill_str .= "$dot$sh$occup.($LS_reordered$sh_par_char<$J1_reordered>).($J2_reordered)";
          } else {
            $nL_prev = $nL_last;
            $nL_last = $sh;
            $nL_last =~ s/[a-z]//g;
            $nL_last = $nL_last*10000 + $L_moment{$L};

            #my ($LS, $J1, $J2) = &get_JJ_terms($ns);
            $J_prev_to_last = $J_last;
            $J_last = $J2_reordered;
            $J_last_shell = $J1_reordered;
            $fill_str .= "$dot$sh$occup.($LS_reordered$sh_par_char<$J1_reordered>)";
          }
          $dot = '.';
        } else {
          $dot = '.' if $fill_str;
          $fill_str .= "$dot$sh$occup";
          $dot = '.';
        }

      } else {
        $first_filled = $n_shell unless $first_filled >= 0;

        $fill_str .= "$dot$sh";
        $fill_str .= $occup if ( $occup > 1 );
        $dot = '.';

        if ( ($LS_reordered ne '') && ($reordered || (($sh ne 's') || ($fill_str =~ /\.\(.+\)/))) ) {
          if ( $sh ne 's' ) {
            $fill_str .= ".($LS_reordered$sh_par_char<$J1_reordered>)";
          }
          if ( (($n_shell != $first_filled) && ($ns != $last_shell)) || ($reordered && $n_shell != $last_shell_reordered)) {
            $fill_str .= ".($J2_reordered)";
          }
        }
      }
    }
    if ( ($ns == $last_shell) && !$reordered) {
      # If, due to complete last shells, the previous parent term has
      # moved up tp the very end, detach it and convert into the final term
      if ($fill_str =~ /^([^()]+)\.\(([^()<>]+)<([^<>]+)>\)\.\(([^()]*)\)$/ ) {
        $fill_str = "$1.($2<$3>)";
      }
      if ( $fill_str =~ /^([^()]+)\.\(([^()<>]+)<([^<>]+)>\)$/ ) {
        $fill_str = "$1.($2<$3>)";
      }
      # If there are only two intermediate shell J's and no intermediate accumulating J's,
      # make the final JJ term out of the two intermediate J's
      if ( $fill_str =~ /^([^<> ]*)<([^<> ]*)>([^<> ]*)<([^<> ]*)>\)(\.\([^()<>]+\)){0,1}$/ ) {
        # Detach the last shell's LSJ if it is a singly occupied shell or a shell with one hole
        $fill_str =~ s/([spdfghiklmno]|p5|d9|f13)\.\([^.()]+\)$/$1/;
      }
      # Take the last intermediate accumulating J,
      # make the final JJ term out of it and the last shell's J
      my $fs1 = $fill_str;
      my @Js = ();
      my $nJ = 0;
      my $J11 = '';
      my $term = '';
      while ( $fs1 =~ /^(.+)\.\(([#0-9\/]+)\)\./ ) {
        $nJ++;
        $J11 = $2;
        $term = $J11 unless $term; # The regex finds the last intermediate J first
        $Js[$nJ - 1] = $J11;
        # Replace round parentheses with angular brackets
        $fs1 =~ s/\.\($J11\)\./.<$J11>./;
      }
      if ( $fill_str =~ /^(.+)<([^.<>()]+)>\)$/ ) {
        # Detach the last shell's LSJ if it is a singly occupied shell or a shell with one hole
        $fs1 =~ s/([spdfghiklmno]|p5|d9|f13)\.\([^.()]+\)$/$1/;
        $fill_str = $fs1;
      }

    }
  }
  # Remove the unnecessary intermediate states
  $fill_str =~ s/s\.\(#\d+<#\d+>\)\./s./g;

  # Append the final term to the state name
  my $final_term = ($nL_last > $nL_prev) ? "($J_prev_to_last,$J_last_shell)$acc_par_char"
                                         : "($J_last_shell,$J_prev_to_last)$acc_par_char";
  $fill_str .= "\t$final_term";

  # Convert the string template to an array
  my $template = [];
  while ( $fill_str =~ /^([^#]+){0,1}\#(\d+)(.*)$/ ) {
    # Strip the '#' from the term indexes
    push(@{$template}, $1, $2);
    $fill_str = $3;
  }
  # Store the remainder in the last element of the template
  push(@{$template}, $fill_str);
  return $template;
}  ##make_JJ_template

############################################################################
sub get_JJ_terms($) {    #12/5/2013 8:33AM A.Kramida
############################################################################
  my $ns = shift; # The shell number, 0..7
  # Shell LS term and J and intermediate J values after summation with previous shells
  #my @fill = ();
  #for ( my $i = 0; $i <= 7; $i++ ) {
  #  my $j = $i*3;
  #  $fill[$i] = [$p[$j],$p[$j+1],$p[$j+2]];
  #}
  my @term_numbers = ('#'.(3*$ns+1),'#'.(3*$ns+2),'#'.(3*$ns+3));  # Term numbers are increased by 1 to make them non-zero
  return @term_numbers;
} ##get_JJ_terms($)

############################################################################
sub fill_LS_desig(@) {  #12/5/2013 8:33AM A.Kramida
############################################################################
  my ($parity,$state_num,$J,$n_conf, @p) = @_;
  #if (($parity == 2) && ($state_num == 4) && ($J == 2)) {
  #  $J = $J;
  #}
  foreach (@p) {
    $_ =~ s/^\s+|\s+$//g; # strip leading and trailing spaces
  }
  $n_conf += 0;

  $basis[$parity-1]->{'LS'}->{$J}->{$state_num} = {};
  my $L_format = $p[11]; # Will be non-empty if > 4 shells

  my $template = $LS_templates[$parity-1]->[$n_conf-1];

  my $fill_str = '';
  my $is_ref = 1;
  my ($n_last,$L_last,$occ_last) = (0,0,0);
  for (my $i = 0; $i <= $#{$template}; $i++) {
    my $tok = $template->[$i];
    $is_ref = !$is_ref;
    $tok = ($tok ? $p[$tok-1] : '') if $is_ref;
    $fill_str .= $tok;
    if ( $tok =~ /(\d+)([a-z]+)(\d*)/ ) {
      ($n_last,$L_last,$occ_last) = ($1,$2,$3);
      $occ_last = 1 unless $occ_last;
    }
  }

  $fill_str =~ s/[<]/\(/g;
  $fill_str =~ s/[>]/\)/g;
  $fill_str =~ s/[*](\d+)/$1/g;

  if ( $fill_str =~ /^.+\t$/) {
    $fill_str .= '1S';
  }
  $fill_str =~ s/^([0-9]+p6\.[0-9]+d\.)\(2D\)\.([0-9]+s)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+d\.)\(2D\)\.([0-9]+d)\.\(2D\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+d\.)\(2D\)\.([0-9]+p)\.\(2P\*\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+d\.)\(2D\)\.([0-9]+f)\.\(2F\*\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+f\.)\(2F\*\)\.([0-9]+d)\.\(2D\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.)(4f)\.(4d)\t/\1\3.\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+p\.)\(2P\*\)\.([0-9]+d)\.\(2D\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+s\.)\(2S\)\.([0-9]+f)\.\(2F\*\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+s\.)\(2S\)\.([0-9]+d)\.\(2D\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+p6\.[0-9]+s\.)\(2S\)\.([0-9]+p)\.\(2P\*\)\t/\1\2\t/;
  $fill_str =~ s/^([0-9]+s\.[0-9]+p5\.)\(2P\*\)\.(\([13]P\*\))/\1\2/;
  $fill_str =~ s/([0-9]+d)\.\(2D\)\t/\1\t/;

  if ( $fill_str =~ /4f3.\(2D/ ) {
    $fill_str=$fill_str;
  }
  if ( $fill_str =~ /4f4.\(1D/ ) {
    $fill_str=$fill_str;
  }
  if ( $fill_str =~ /4f4\t3K/ ) {
    $fill_str=$fill_str;
  }
  $basis[$parity-1]->{'LS'}->{$J}->{$state_num}->{'cn'} = $n_conf;
  $basis[$parity-1]->{'LS'}->{$J}->{$state_num}->{'sh'} = $fill_str;
  $basis[$parity-1]->{'LS'}->{$J}->{$state_num}->{'lastsh'} = [$n_last,$L_last,$occ_last];

  return $L_format;
} ##fill_LS_desig(@)

############################################################################
sub fill_JJ_desig(@) {  #12/5/2013 8:33AM A.Kramida
############################################################################
  my ($parity,$state_num,$n_conf,$J,@p) = @_;

  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num} = {};

  my $fill_str = '';
  my $template = $JJ_templates[$parity-1]->[$n_conf-1];
  #if ( ($n_conf == 6) && ($parity == 1)  ) {
  #  $n_conf = $n_conf;
  #}
  my $is_ref = 1;
  my ($n_last,$L_last,$occ_last) = (0,0,0);
  for (my $i = 0; $i <= $#{$template}; $i++) {
    my $tok = $template->[$i];
    $is_ref = !$is_ref;
    if ($is_ref) {
      $tok = ($tok ? $p[$tok-1] : '') if $is_ref;
      if ( $tok =~ /^\d+\.\d+$/ ) {
        # This token is a floating-point J value - convert it to a string
        $tok = &get_J($tok);
      }
    }
    $fill_str .= $tok;
    if ( $tok =~ /(\d+)([a-z]+)(\d*)/ ) {
      ($n_last,$L_last,$occ_last) = ($1,$2,$3);
      $occ_last = 1 unless $occ_last;
    }
  }

  $fill_str =~ s/[<]/\(/g;
  $fill_str =~ s/[>]/\)/g;

  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num}->{'cn'} = $n_conf;
  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num}->{'sh'} = $fill_str;
  $basis[$parity-1]->{'JJ'}->{$J}->{$state_num}->{'lastsh'} = [$n_last,$L_last,$occ_last];

}  ##fill_JJ_desig

#############################################################################
#sub FindLev($$$$$) {   #03/03/2017 8:30AM A.Kramida
#############################################################################
#  my ($par, $J, $E, $nc, $term) = @_;
#  my $key = join("\t",$nc,$term);
#  my $num_RCG_e = $term_labels_map[$par-1]->{$J}->{$key};
#  my $E1 = $energies[$par-1]->{$J}->{$num_RCG_e}->[0];
#  if (abs($E1 - $E) > 0.03) {
#    die "Energy mismatch in OUTG11 between transition list ($E) and level list ($E1), parity $par, J=$J";
#  }
#  return $num_RCG_e;
#} ##FindLev($$$$$) {

############################################################################
sub FindLev($$$$$) {   #1/5/2005 4:03PM A.Kramida
############################################################################
  my ($par, $J, $E, $nc, $term) = @_;
  my ($E1, $min_diff, $n_lev) = (-10000, 1e10, 0);
  my ($num_lead_bs, $lead_c_no, $lead_term) = (0,0,'');

  my @states = @{$basis_hash[$par-1]->{$J}->[$nc-1]->{$term}};
  foreach my $st (@states) {
    my ($A, $num_RCG_e, $num_state) = @{$st};
    $E1 = $energies[$par-1]->{$J}->{$num_RCG_e}->[0];
    next if (abs($E1 - $E) > 0.03);
    if ( abs($E1-$E) < $min_diff ) {
      $min_diff = abs($E1-$E);
      $n_lev = $num_RCG_e;
    }
    my $max_RCG_comp = $#{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$num_RCG_e-1]};
    $max_RCG_comp = 5 if $max_RCG_comp > 5;
    for ( my $j = 0; $j <= 30; $j++ ) {
      if ( ($min_diff <= 0.001*$j) && (abs($E1-$E) <= 0.001*$j) ) {
        for ( my $m = 0; $m <= $max_RCG_comp; $m++ ) {
          my ($A,$num_bas) = @{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$num_RCG_e-1]->[$m]};
          $num_lead_bs = $num_bas;
          $lead_c_no  = $basis[$par-1]->{'LS'}->{$J}->{$num_lead_bs}->{'cn'};
          $lead_term  = $basis[$par-1]->{'LS'}->{$J}->{$num_lead_bs}->{'label'};
          if ( ($lead_c_no == $nc) && ($lead_term eq $term) ) {
            last;
          }
        }
        if ( $num_lead_bs && ($lead_c_no == $nc) && ($lead_term eq $term)  ) {
          $n_lev = $num_RCG_e;
          last;
        }
      }
      last if ($n_lev && $num_lead_bs);
    }
  }
  if ( $n_lev && !$num_lead_bs ) {
    my $max_RCG_comp = $#{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$n_lev-1]};
    #$max_RCG_comp = 2 if $max_RCG_comp > 2;
    for ( my $m = 0; $m <= $max_RCG_comp; $m++ ) {
      my ($A,$num_bas) = @{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$n_lev-1]->[$m]};
      $num_lead_bs = $num_bas;
      $lead_c_no  = $basis[$par-1]->{'LS'}->{$J}->{$num_lead_bs}->{'cn'};
      $lead_term  = $basis[$par-1]->{'LS'}->{$J}->{$num_lead_bs}->{'label'};
      if ( ($lead_c_no == $nc) && ($lead_term eq $term) ) {
        last;
      }
    }
  }
  $E1 = $energies[$par-1]->{$J}->{$n_lev}->[0];
  if ( $n_lev && (abs($E1-$E) > 0.03) ) {
    print "\nWarning: too large energy differnce between identified RCG levels, parity $par, J=$J: expected $E, found $E1";
  }
  return ($n_lev, $E1, $lead_c_no, $lead_term);
} ##FindLev($$$$$)

############################################################################
sub FindLev1($$$$$$;$$) {   #8/23/2018 5:01PM A.Kramida
############################################################################
  my ($n_lev, $par, $J, $E, $nc, $term, $outg11_line_no, $call_no) = @_;
  my ($num_lead_bs, $lead_c_no, $lead_term) = (0,0,'');
  my $E1 = $energies[$par-1]->{$J}->{$n_lev}->[0];
  if (abs($E1 - $E) > 0.03) {
    print "\nWarning: too large energy differnce between identified RCG levels, parity $par, J=$J: expected $E, found $E1";
  }
  if (($par < 1) || ($nc < 1)) {
    die "\nInvalid parity index $par or configuration index $nc";
  }
  if (ref($basis_hash[$par-1]->{$J}->[$nc-1]->{$term}) ne 'ARRAY') {
    die "\nUndefined basis hash entry. Parity $par, J=$J, conf index $nc, term $term, OUTG11 line no $outg11_line_no, call no $call_no";
  }
  my @states = @{$basis_hash[$par-1]->{$J}->[$nc-1]->{$term}};
  if ( $n_lev && !$num_lead_bs ) {
    my $max_RCG_comp = $#{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$n_lev-1]};
    #$max_RCG_comp = 2 if $max_RCG_comp > 2;
    for ( my $m = 0; $m <= $max_RCG_comp; $m++ ) {
      my ($A,$num_bas) = @{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$n_lev-1]->[$m]};
      $num_lead_bs = $num_bas;
      $lead_c_no  = $basis[$par-1]->{'LS'}->{$J}->{$num_lead_bs}->{'cn'};
      $lead_term  = $basis[$par-1]->{'LS'}->{$J}->{$num_lead_bs}->{'label'};
      if ( ($lead_c_no == $nc) && ($lead_term eq $term) ) {
        last;
      }
    }
  }
  return ($n_lev, $E1, $lead_c_no, $lead_term);
} ##FindLev1($$$$$$)

############################################################################
sub get_vectors_sorted() {    #12/10/2013 2:08PM
############################################################################
  return @vectors_sorted;
} ##get_vectors_sorted()

############################################################################
sub get_params() {    #12/10/2013 4:17PM
############################################################################
  return (\@params,\@CI);
} ##get_params()

############################################################################
sub get_leading_LS_term {   #12/11/2013 5:47PM A.Kramida
############################################################################
  my ($par,$J,$n_lev) = @_;
  # Take the leading term as the level designation
  my $par_code = $parities[$par-1];
  my $shells = '';
  #foreach my $num_bas (sort {abs($vectors[$par-1]->{$cpl}->{$J}->{$n_lev}->{$b})<=>
  #                          abs($vectors[$par-1]->{$cpl}->{$J}->{$n_lev}->{$a})}
  #                      keys %{$vectors[$par-1]->{$cpl}->{$J}->{$n_lev}}) {
  #  $shells = $basis[$par-1]->{$cpl}->{$J}->{$num_bas}->{'shells'};
  #  # Correct the final parity
  #  if ($par_code eq 'e') {
  #    $shells =~ s/[*]$//;
  #  } else {
  #    $shells .= '*' unless $shells =~ /[*]$/;
  #  }
  #  last;
  #}
  my ($A,$num_bas) = @{$vectors_sorted[$par-1]->{'LS'}->{$J}->[$n_lev-1]->[0]};
  my $cpl = 'LS';
  my $shells = $basis[$par-1]->{$cpl}->{$J}->{$num_bas}->{'sh'};
  # Correct the final parity
  if ($par_code eq 'e') {
    $shells =~ s/[*]$//;
  } else {
    $shells .= '*' unless $shells =~ /[*]$/;
  }
  return $shells;
} ##get_leading_term


1;
