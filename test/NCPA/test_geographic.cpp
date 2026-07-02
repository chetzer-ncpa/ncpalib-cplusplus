#include "NCPA/geographic.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <limits>
#include <numbers>
#include <utility>

using namespace std;
using namespace NCPA::geographic;
using namespace testing;

TEST( NCPAGeographicLibraryTest, ConstantEarthRadiusReturnsExpected ) {
    EXPECT_DOUBLE_EQ( earthradius(), 6371.009 );
}