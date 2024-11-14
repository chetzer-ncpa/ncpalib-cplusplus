#include "NCPA/arrays.hpp"
#include "NCPA/gtest.hpp"
#include "NCPA/linearalgebra.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <memory>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace testing;
using namespace std;
using namespace NCPA::linear;

typedef double test_t;

class NCPALinearAlgebraDenseVectorTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            svd      = std::vector<test_t>( { 1.0, 2.0, 3.0, 4.0, 5.0 } );
            dvec_std = details::dense_vector<test_t>( svd );

            wrapper1 = Vector<test_t>( dvec_std.clone() );
            vecsize = 5;
        }  // void TearDown() override {}

        // declare stuff here
        details::dense_vector<test_t> dvec_initlist { 1.0, 2.0, 3.0, 4.0,
                                                      5.0 };
        details::dense_vector<test_t> dvec_std, dvec_copy, dvec;
        std::vector<test_t> svd;

        Vector<test_t> wrapper1, wrapper2;
        size_t vecsize;
};

TEST_F( NCPALinearAlgebraDenseVectorTest, DefaultConstructorWorks ) {
    EXPECT_EQ( dvec.size(), 0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, InitializerListConstructorWorks ) {
    EXPECT_EQ( dvec_initlist.size(), vecsize );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_initlist.get( i ), (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, StdVectorConstructorWorks ) {
    EXPECT_EQ( dvec_std.size(), vecsize );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, CopyConstructorWorks ) {
    dvec = details::dense_vector<test_t>( dvec_std );
    EXPECT_EQ( dvec.size(), vecsize );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec.get( i ), (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, CopyConstructorMakesDeepCopy ) {
    dvec = details::dense_vector<test_t>( dvec_std );
    dvec_std.clear();
    EXPECT_EQ( dvec.size(), vecsize );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec.get( i ), (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, ClearClearsContents ) {
    EXPECT_EQ( dvec_std.size(), vecsize );
    dvec_std.clear();
    EXPECT_EQ( dvec_std.size(), 0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, SetSetsSingleValue ) {
    dvec_std.set( 0, 6 );
    EXPECT_DOUBLE_EQ( dvec_std.get( 0 ), 6.0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, SetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    dvec_std.set( newcontents );
    EXPECT_EQ( dvec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, SetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    dvec_std.set( 9, newcontents );
    EXPECT_EQ( dvec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, CountNonzeroIndicesIsCorrect ) {
    EXPECT_EQ( dvec_std.count_nonzero_indices(), vecsize );
    dvec_std.set( 0, 0.0 );
    EXPECT_EQ( dvec_std.count_nonzero_indices(), 4 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, NonzeroIndicesReportedCorrectly ) {
    std::vector<size_t> inds = NCPA::arrays::index_vector<size_t>( vecsize );
    EXPECT_ARRAY_EQ( 5, dvec_std.nonzero_indices(), inds );

}

TEST_F( NCPALinearAlgebraDenseVectorTest, ResizeResizesArray ) {
    dvec_std.resize( 6 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) );
    }
    EXPECT_DOUBLE_EQ( dvec_std.get( vecsize ), 0.0 );
    dvec_std.resize( 3 );
    for ( size_t i = 0; i < 3; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) );
    }
    EXPECT_ANY_THROW( { dvec_std.get( 4 ); } );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AsStdMatches ) {
    EXPECT_ARRAY_DOUBLE_EQ( 5, dvec_std.as_std(), svd );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, CloneClones ) {
    std::unique_ptr<details::abstract_vector<test_t>> vptr = dvec_std.clone();
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( vptr->get( i ), (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, CloneClonesDeep ) {
    std::unique_ptr<details::abstract_vector<test_t>> vptr = dvec_std.clone();
    dvec_std.clear();
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( vptr->get( i ), (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AsArrayMatches ) {
    test_t *darr = NCPA::arrays::zeros<test_t>( vecsize );
    dvec_std.as_array( vecsize, darr );
    EXPECT_ARRAY_DOUBLE_EQ( 5, darr, svd );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, ScaleWorks ) {
    dvec_std.scale( 2.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, ScaleByVectorWorks ) {
    dvec_std.scale( dvec_std );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ),
                          (test_t)( i + 1 ) * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, DotWorks ) {
    EXPECT_DOUBLE_EQ( dvec_std.dot( dvec_std ),
                      (test_t)( 1 + 4 + 9 + 16 + 25 ) );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AddWorksWithScalar ) {
    dvec_std.add( 3.5 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AddWorksWithVector ) {
    dvec_std.add( dvec_std );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AddWorksWithModifiedVector ) {
    dvec_std.add( dvec_std, -1.0 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AddThrowsExceptionOnMismatchedSizes ) {
    dvec_std.resize( 3 );
    EXPECT_THROW( { dvec_initlist.add( dvec_std ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, PlusEqualOperatorWorks ) {
    dvec_std += dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, PlusEqualOperatorWorksWithScalar ) {
    dvec_std += 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, MinusEqualOperatorWorks ) {
    dvec_std -= dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, MinusEqualOperatorWorksWithScalar ) {
    dvec_std -= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) - 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, TimesEqualOperatorWorks ) {
    dvec_std *= dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ),
                          (test_t)( i + 1 ) * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, TimesEqualOperatorWorksWithScalar ) {
    dvec_std *= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) * 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, DivideEqualOperatorWorksWithScalar ) {
    dvec_std /= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), (test_t)( i + 1 ) / 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, EqualityOperatorWorks ) {
    EXPECT_TRUE( dvec_std == dvec_initlist );
    dvec_std.set( 0, 5.0 );
    EXPECT_FALSE( dvec_std == dvec_initlist );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, InequalityOperatorWorks ) {
    EXPECT_FALSE( dvec_std != dvec_initlist );
    dvec_std.set( 0, 5.0 );
    EXPECT_TRUE( dvec_std != dvec_initlist );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, SwapFunctionWorks ) {
    dvec_std.set( 0, 5.0 );
    ASSERT_DOUBLE_EQ( dvec_initlist.get( 0 ), 1.0 );
    std::swap( dvec_std, dvec_initlist );
    ASSERT_DOUBLE_EQ( dvec_initlist.get( 0 ), 5.0 );
    ASSERT_DOUBLE_EQ( dvec_std.get( 0 ), 1.0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, AssignmentOperatorWorks ) {
    dvec_std *= 2.0;
    dvec      = dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( dvec.get( i ), (test_t)( i + 1 ) * 2.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperSizeWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperCopyConstructorWorks ) {
    wrapper2 = Vector<test_t>( wrapper1 );
    EXPECT_TRUE( wrapper1.equals( wrapper2 ) );
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_FALSE( wrapper1.equals( wrapper2 ) );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperSwapWorks ) {
    swap( wrapper1, wrapper2 );
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperEqualityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_TRUE( wrapper1 == wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperInequalityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_FALSE( wrapper1 != wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_TRUE( wrapper1 != wrapper2 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperAssignmentOperatorWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
    wrapper2 = wrapper1;
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperGetAllowsModification ) {
    wrapper1.get( 0 ) = 5.0;
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 5.0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperIndexOperatorAllowsRead ) {
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperIndexOperatorAllowsWrite ) {
    wrapper1[ 0 ] = 5.0;
    for ( size_t i = 1; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], (test_t)( i + 1 ) );
    }
    EXPECT_DOUBLE_EQ( wrapper1[ 0 ], 5.0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperAsStdReturnsCorrectVector ) {
    EXPECT_ARRAY_DOUBLE_EQ( 5, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperClearWorks ) {
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_THROW( { wrapper1[ vecsize ] = 2.0; }, std::logic_error );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperResizeWorks ) {
    wrapper1.resize( 8 );
    EXPECT_EQ( wrapper1.size(), 8 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], (test_t)( i + 1 ) );
    }
    for ( size_t i = 5; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperAsArrayReturnsCorrectValues ) {
    test_t *d = NCPA::arrays::zeros<test_t>( vecsize );
    wrapper1.as_array( 5, d );
    EXPECT_ARRAY_DOUBLE_EQ( 5, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperSetSetsSingleValue ) {
    wrapper1.set( 0, 6 );
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 6.0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperSetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    wrapper1.set( newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperSetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    wrapper1.set( 9, newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperScaleWorks ) {
    wrapper1.scale( 2.0 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperScaleByVectorWorks ) {
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ),
                          (test_t)( i + 1 ) * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperDotWorks ) {
    EXPECT_DOUBLE_EQ( wrapper1.dot( wrapper1 ),
                      (test_t)( 1 + 4 + 9 + 16 + 25 ) );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperAddWorksWithScalar ) {
    wrapper1.add( 3.5 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperAddWorksWithVector ) {
    wrapper1.add( wrapper1 );
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperAddThrowsExceptionOnMismatchedSizes ) {
    wrapper2 = wrapper1;
    wrapper1.resize( 3 );
    EXPECT_THROW( { wrapper2.add( wrapper1 ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperPlusEqualOperatorWorks ) {
    wrapper1 += wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperPlusEqualOperatorWorksWithScalar ) {
    wrapper1 += 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperNegativeOperatorWorks ) {
    wrapper2 = -wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper2.get( i ), -(test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperMinusEqualOperatorWorks ) {
    wrapper1 -= wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperMinusEqualOperatorWorksWithScalar ) {
    wrapper1 -= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) - 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperTimesEqualOperatorWorks ) {
    wrapper1 *= wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ),
                          (test_t)( i + 1 ) * (test_t)( i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperTimesEqualOperatorWorksWithScalar ) {
    wrapper1 *= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) * 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperDivideEqualOperatorWorksWithScalar ) {
    wrapper1 /= 3.5;
    for ( size_t i = 0; i < 5; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), (test_t)( i + 1 ) / 3.5 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperPlusOperatorWorks ) {
    wrapper2 = wrapper1 + wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], 2.0*wrapper1[i] );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperMinusOperatorWorks ) {
    wrapper2 = wrapper1 - wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, WrapperTimesOperatorWorks ) {
    wrapper2 = wrapper1 * wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], wrapper1[i]*wrapper1[i] );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, ZeroWorksWithSingleIndex ) {
    dvec_std.zero( 2 );
    EXPECT_DOUBLE_EQ( dvec_std.get( 2 ), 0.0 );
}

TEST_F( NCPALinearAlgebraDenseVectorTest, ZeroWorksWithVector ) {
    std::vector<size_t> zero_out;
    for (size_t i = 0; i < vecsize; i++) {
        zero_out.push_back( i );
    }
    dvec_std.zero( zero_out );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraDenseVectorTest, ZeroWorksWithInitList ) {
    dvec_std.zero( { 0, 1, 2, 3, 4 } );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( dvec_std.get( i ), 0.0 );
    }
}