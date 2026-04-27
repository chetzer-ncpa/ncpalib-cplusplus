#pragma once

#ifndef NCPA_ARMADILLO_INCLUDED
#  if __has_include( "armadillo" )
#    if __has_include( "fftw3.h" )
#      ifndef ARMA_USE_FFTW3
#        define ARMA_USE_FFTW3
#      endif
#      ifdef ARMA_DONT_USE_FFTW3
#        undef ARMA_DONT_USE_FFTW3
#      endif
#    endif
#    include <armadillo>
#    define NCPA_ARMADILLO_INCLUDED
#  else
static_assert( false, "Armadillo header library not accessible" );
#  endif
#endif
