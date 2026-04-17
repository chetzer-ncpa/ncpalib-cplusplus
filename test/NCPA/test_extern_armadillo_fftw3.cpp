#if __has_include( "armadillo" )
#  if __has_include( "fftw3.h" )
#    include "NCPA/extern/armadillo_with_fftw3.hpp"
#    include "NCPA/geographic.hpp"
#    include "NCPA/gtest.hpp"

#    include <gmock/gmock.h>
#    include <gtest/gtest.h>

using namespace std;

TEST( NCPAExternLibraryArmadilloWithFFTW3Test, ArmadilloIncludedCorrectly ) {
#    ifdef ARMA_INCLUDES
    bool armadillo_included = true;
#    else
    bool armadillo_included = false;
#    endif
    EXPECT_TRUE( armadillo_included );
}

TEST( NCPAExternLibraryArmadilloWithFFTW3Test, ArmaUseFFTW3SetCorrectly ) {
#    ifdef ARMA_USE_FFTW3
    bool use_fftw3 = true;
#    else
    bool use_fftw3 = false;
#    endif
    EXPECT_TRUE( use_fftw3 );
}

#  endif
#endif
