#pragma once

#ifndef NCPA_ARMADILLO_INCLUDED

#  ifdef ARMA_USE_FFTW3
#    undef ARMA_USE_FFTW3
#  endif

#  define NCPA_ARMADILLO_INCLUDED
#  include <armadillo>

#else
static_assert( false, "Armadillo header library not accessible" );
#endif
