// only run this test if armadillo is actually available
#if __has_include( "armadillo" )
#  include "NCPA/extern/armadillo.hpp"
#  include "NCPA/geographic.hpp"
#  include "NCPA/gtest.hpp"

#  include <gmock/gmock.h>
#  include <gtest/gtest.h>

using namespace std;

TEST( NCPAExternLibraryArmadilloTest,
      ArmadilloIncludedCorrectlyWithoutFFTW3 ) {
#  ifdef ARMA_INCLUDES
    bool armadillo_included = true;
#  else
    bool armadillo_included = false;
#  endif
    EXPECT_TRUE( armadillo_included );
}
#endif
