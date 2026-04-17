#pragma once

#ifndef NCPA_ARMADILLO_INCLUDED
#  if __has_include( "fftw3.h" )
#    define NCPA_ARMADILLO_INCLUDED
#    define ARMA_USE_FFTW3
#    include <armadillo>
#  else
static_assert(
    false, "FFTW3 header file not found, can't tell armadillo to use it!" );
#  endif
#endif
