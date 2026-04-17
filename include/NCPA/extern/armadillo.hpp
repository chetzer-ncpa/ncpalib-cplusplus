#pragma once

#ifndef NCPA_ARMADILLO_INCLUDED
#  if __has_include( "armadillo" )
#    if __has_include( "fftw3.h" )
#      define ARMA_USE_FFTW3
#    endif
#    include <armadillo>
#    define NCPA_ARMADILLO_INCLUDED
#  else
static_assert( false, "Armadillo header library not accessible" );
#  endif
#endif
