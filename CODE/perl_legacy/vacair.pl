#!perl
use strict;

sub nair1($) {
#    Calculates the deviation of the refractive index of air from unity
#   (n-1) by a formula derived by E.R.Peck and K.Reeder:
#      JOSA 62 (1972), p.958-962.
#      Sigma is a wavenumber in cm-1.
  my $sigma = shift;
  $sigma *= $sigma;
# square of sigma is stored in the same variable;
#  it does not change the value of variable used in the function call,
#  because it is transferred to the function as a value-parameter,
#  not as a var-parameter
  my $n = 8060.51e-8 + 2480990/(132.274e8-$sigma) + 17455.7/(39.32957e8-$sigma);
  return $n;
}

sub Lair($) {
#    Calculates the wavelength in air from the known wavelength in vacuum
#   using the nair1 function.
#      Wavelengths are in angstroems
  my $Lvac = shift;
  my $L = $Lvac/(1 + &nair1(1e8/$Lvac));
  return $L;
}

sub Lvac($) {
#    Calculates the wavelength in vacuum from the known wavelength in air
#   using the nair1 function.
#      Wavelengths are in angstroems *)
  my $La = shift;
  my $tol = 1e-12;

  my ($Lv1,$Lv0) = ($La,$La);
  if ($La > 1500.0) {
    while ( 1 ) {
      $Lv1 = $La*(1 + &nair1(1e8/$Lv0));
      if (abs($Lv1 - $Lv0) <= $tol) {
        last;
      }
      $Lv0 = $Lv1;
    }
  } else {
    $Lv1 = $La;
  }
  return $Lv1;
}

1;