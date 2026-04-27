// only run this test if armadillo is actually available
#if __has_include( "armadillo" )
#  include "NCPA/extern/armadillo_no_fftw3.hpp"
#  include "NCPA/geographic.hpp"
#  include "NCPA/gtest.hpp"

#  include <gmock/gmock.h>
#  include <gtest/gtest.h>

using namespace std;

TEST( NCPAExternLibraryArmadilloWithoutFFTW3Test,
      ArmadilloIncludedCorrectly ) {
#  ifdef ARMA_INCLUDES
    bool armadillo_included = true;
#  else
    bool armadillo_included = false;
#  endif
    EXPECT_TRUE( armadillo_included );
}

TEST( NCPAExternLibraryArmadilloWithoutFFTW3Test,
      ArmaUseFFTW3UnsetCorrectly ) {
#  ifdef ARMA_USE_FFTW3
    bool use_fftw3_set = true;
#  else
    bool use_fftw3_set = false;
#  endif
    EXPECT_FALSE( use_fftw3_set );
}


#endif
