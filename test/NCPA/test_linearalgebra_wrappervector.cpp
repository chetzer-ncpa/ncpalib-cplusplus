#include "NCPA/math.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

typedef double test_t;

#define _TEST_EQ_       EXPECT_DOUBLE_EQ
#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_DOUBLE_EQ
#define _TEST_TITLE_    NCPALinearAlgebraLibraryWrapperVectorTest

class _TEST_TITLE_ : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            v1.set( { 1.0, 2.0, 3.0, 4.0, 5.0 } );
            wrapper1 = WrapperVector<test_t>( v1 );
        }  // void TearDown() override {}

        // declare stuff here

        const test_t zero = NCPA::math::zero<test_t>(),
                     one  = NCPA::math::one<test_t>();
        WrapperVector<test_t> wrapper1, wrapper2;
        Vector<test_t> newvec;
        details::dense_vector<test_t> v1;
        size_t vecsize = 5;
};



TEST_F( _TEST_TITLE_, WrapperSizeWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
}

TEST_F( _TEST_TITLE_, WrapperCopyConstructorWorks ) {
    wrapper2 = WrapperVector<test_t>( wrapper1 );
    EXPECT_TRUE( wrapper1.equals( wrapper2 ) );
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_FALSE( wrapper1.equals( wrapper2 ) );
}

TEST_F( _TEST_TITLE_, WrapperSwapWorks ) {
    swap( wrapper1, wrapper2 );
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
}

TEST_F( _TEST_TITLE_, WrapperEqualityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_TRUE( wrapper1 == wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( _TEST_TITLE_, WrapperInequalityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_FALSE( wrapper1 != wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_FALSE( wrapper1 != wrapper2 );
}

TEST_F( _TEST_TITLE_, WrapperAssignmentOperatorWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
    wrapper2 = wrapper1;
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( _TEST_TITLE_, WrapperGetAllowsModification ) {
    wrapper1.get( 0 ) = 5.0;
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 5.0 );
}

TEST_F( _TEST_TITLE_, WrapperIndexingOperatorAllowsModification ) {
    EXPECT_DOUBLE_EQ( wrapper1[0], wrapper1.get( 0 ) );
    wrapper1[0] = 5.0;
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 5.0 );
}

TEST_F( _TEST_TITLE_, WrapperIndexOperatorAllowsRead ) {
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], (test_t)( i + 1 ) );
    }
}

TEST_F( _TEST_TITLE_, WrapperIndexOperatorAllowsWrite ) {
    wrapper1[ 0 ] = 5.0;
    for ( size_t i = 1; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], (test_t)( i + 1 ) );
    }
    EXPECT_DOUBLE_EQ( wrapper1[ 0 ], 5.0 );
}

TEST_F( _TEST_TITLE_, WrapperAsStdReturnsCorrectVector ) {
    EXPECT_ARRAY_DOUBLE_EQ( vecsize, wrapper1.as_std(), v1.as_std() );
}

TEST_F( _TEST_TITLE_, WrapperClearWorks ) {
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_THROW( { wrapper1[ vecsize ] = 2.0; }, std::logic_error );
}

TEST_F( _TEST_TITLE_, WrapperResizeWorks ) {
    wrapper1.resize( 8 );
    EXPECT_EQ( wrapper1.size(), 8 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], (test_t)( i + 1 ) );
    }
    for ( size_t i = 5; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], 0.0 );
    }
}

TEST_F( _TEST_TITLE_, WrapperAsArrayReturnsCorrectValues ) {
    test_t *d = NCPA::arrays::zeros<test_t>( vecsize );
    wrapper1.as_array( 5, d );
    EXPECT_ARRAY_DOUBLE_EQ( 5, wrapper1.as_std(), v1.as_std() );
}

TEST_F( _TEST_TITLE_, WrapperSetSetsSingleValue ) {
    wrapper1.set( 0, 6 );
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 6.0 );
}

TEST_F( _TEST_TITLE_, WrapperSetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    wrapper1.set( newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( _TEST_TITLE_, WrapperSetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    wrapper1.set( 9, newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( _TEST_TITLE_, WrapperScaleWorks ) {
    wrapper1.scale( 2.0 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( _TEST_TITLE_, WrapperScaleByVectorWorks ) {
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ),
                          (test_t)( i + 1 ) * (test_t)( i + 1 ) );
    }
}

TEST_F( _TEST_TITLE_, WrapperDotWorks ) {
    EXPECT_DOUBLE_EQ( wrapper1.dot( wrapper1 ),
                      (test_t)( 1 + 4 + 9 + 16 + 25 ) );
}

TEST_F( _TEST_TITLE_, WrapperAddWorksWithScalar ) {
    wrapper1.add( 3.5 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) + 3.5 );
    }
}

TEST_F( _TEST_TITLE_, WrapperAddWorksWithVector ) {
    wrapper1.add( wrapper1 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

// TEST_F( _TEST_TITLE_, WrapperAddThrowsExceptionOnMismatchedSizes ) {
//     wrapper2 = wrapper1;
//     wrapper1.resize( 8 );
//     EXPECT_THROW( { wrapper2 += wrapper1; }, std::invalid_argument );
// }

TEST_F( _TEST_TITLE_, WrapperPlusEqualOperatorWorks ) {
    wrapper1 += wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( _TEST_TITLE_, WrapperPlusEqualOperatorWorksWithScalar ) {
    wrapper1 += 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) + 3.5 );
    }
}


TEST_F( _TEST_TITLE_, WrapperMinusEqualOperatorWorks ) {
    wrapper1 -= wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( _TEST_TITLE_, WrapperMinusEqualOperatorWorksWithScalar ) {
    wrapper1 -= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) - 3.5 );
    }
}

TEST_F( _TEST_TITLE_, WrapperTimesEqualOperatorWorks ) {
    wrapper1 *= wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ),
                          (test_t)( i + 1 ) * (test_t)( i + 1 ) );
    }
}

TEST_F( _TEST_TITLE_, WrapperTimesEqualOperatorWorksWithScalar ) {
    wrapper1 *= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) * 3.5 );
    }
}

TEST_F( _TEST_TITLE_, WrapperDivideEqualOperatorWorksWithScalar ) {
    wrapper1 /= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) / 3.5 );
    }
}

TEST_F( _TEST_TITLE_, WrapperPlusOperatorWorks ) {
    cout << "Attempting to add vectors" << endl;
    newvec = wrapper1 + wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( newvec[i], 2.0*wrapper1[i] );
    }
}

TEST_F( _TEST_TITLE_, WrapperMinusOperatorWorks ) {
    newvec = wrapper1 - wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( newvec[i], 0.0 );
    }
}

TEST_F( _TEST_TITLE_, WrapperTimesOperatorWorks ) {
    newvec = wrapper1 * wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( newvec[i], wrapper1[i]*wrapper1[i] );
    }
}