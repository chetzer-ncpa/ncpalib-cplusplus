#include "NCPA/arrays.hpp"
#include "NCPA/constants.hpp"
#include "NCPA/dsp.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <utility>

using namespace std;
using namespace NCPA::dsp;
using namespace NCPA::arrays;
using namespace testing;

TEST( NCPADSPLibraryTest, HannWindowIsCorrect ) {
    size_t nwin = 11;
    HannWindow<double> hann( nwin );
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hann[ i ],
                          ( 0.5
                            - 0.5
                                  * cos( 2.0 * NCPA::constants::PI * (double)i
                                         / (double)( nwin - 1 ) ) ) );
    }
}

TEST( NCPADSPLibraryTest, HammingWindowIsCorrect ) {
    size_t nwin = 11;
    HammingWindow<double> hamm( nwin );
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hamm[ i ],
                          ( 0.54
                            - 0.46
                                  * cos( 2.0 * NCPA::constants::PI * (double)i
                                         / (double)( nwin - 1 ) ) ) );
    }
}

TEST( NCPADSPLibraryTest, HanningWindowIsHannWindow ) {
    size_t nwin = 11;
    HanningWindow<double> hann( nwin );
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hann[ i ],
                          ( 0.5
                            - 0.5
                                  * cos( 2.0 * NCPA::constants::PI * (double)i
                                         / (double)( nwin - 1 ) ) ) );
    }
}