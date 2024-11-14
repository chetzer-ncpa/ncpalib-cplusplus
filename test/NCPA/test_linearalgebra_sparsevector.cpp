#define __NCPALIB_DEBUG true

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

class NCPALinearAlgebraSparseVectorTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            svd = std::vector<double>(
                { 0.0, 0.0, 3.0, 0.0, 5.0, 6.0, 0.0, 0.0 } );
            svd_opposite = std::vector<double>(
                { 1.0, 2.0, 0.0, 4.0, 0.0, 0.0, 7.0, 8.0 } );
            vec_initlist = { 0.0, 0.0, 3.0, 0.0, 5.0, 6.0, 0.0, 0.0 };

            vec_std      = details::sparse_vector<double>( svd );
            vec_opposite = details::sparse_vector<double>( svd_opposite );

            wrapper1         = Vector<double>( vec_std.clone() );
            wrapper_opposite = Vector<double>( vec_opposite.clone() );

            svd_dot_svd = 0.0;
            for (auto it = svd.cbegin(); it != svd.cend(); ++it) {
                svd_dot_svd += (*it) * (*it);
            }
        }  // void TearDown() override {}

        // declare stuff here
        details::sparse_vector<double> vec_initlist, vec_std, vec_opposite, vec;
        std::vector<double> svd, svd_opposite;
        double svd_dot_svd;

        Vector<double> wrapper1, wrapper2, wrapper_opposite;
        size_t vecsize = 8;
};

TEST_F( NCPALinearAlgebraSparseVectorTest, DefaultConstructorWorks ) {
    EXPECT_EQ( vec.size(), 0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, StdvectorConstructorWorks ) {
    // details::sparse_vector<double> vec2( svd );
    // vec_std.set( svd );
    EXPECT_EQ( vec_std.size(), vecsize );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] );
    }
}


TEST_F( NCPALinearAlgebraSparseVectorTest, InitializerListConstructorWorks ) {
    EXPECT_EQ( vec_initlist.size(), vecsize );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_initlist.get( i ), svd[ i ] );
    }
}




TEST_F( NCPALinearAlgebraSparseVectorTest, CopyConstructorWorks ) {
    vec = details::sparse_vector<double>( vec_std );
    EXPECT_EQ( vec.size(), vecsize );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( vec.get( i ), svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, CopyConstructorMakesDeepCopy ) {
    vec = details::sparse_vector<double>( vec_std );
    vec_std.clear();
    EXPECT_EQ( vec.size(), vecsize );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( vec.get( i ), svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ClearClearsContents ) {
    EXPECT_EQ( vec_std.size(), vecsize );
    vec_std.clear();
    EXPECT_EQ( vec_std.size(), 0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, SetSetsSingleValue ) {
    vec_std.set( 0, 6 );
    EXPECT_DOUBLE_EQ( vec_std.get( 0 ), 6.0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, SetReplacesContents ) {
    std::vector<double> newcontents( 9 );
    vec_std.set( newcontents );
    EXPECT_EQ( vec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, SetWithArrayReplacesContents ) {
    double *newcontents = NCPA::arrays::zeros<double>( 9 );
    vec_std.set( 9, newcontents );
    EXPECT_EQ( vec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ResizeResizesArray ) {
    vec_std.resize( vecsize*2 );
    EXPECT_ARRAY_DOUBLE_EQ( vecsize, vec_std, svd );
    EXPECT_DOUBLE_EQ( vec_std.get( vecsize ), 0.0 );
    vec_std.resize( vecsize/2 );
    for ( size_t i = 0; i < vec_std.size(); i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] );
    }
    EXPECT_ANY_THROW( { vec_std.get( vecsize ); } );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AsStdMatches ) {
    EXPECT_ARRAY_DOUBLE_EQ( vecsize, vec_std.as_std(), svd );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, CloneClones ) {
    std::unique_ptr<details::abstract_vector<double>> vptr = vec_std.clone();
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vptr->get( i ), svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, CloneClonesDeep ) {
    std::unique_ptr<details::abstract_vector<double>> vptr = vec_std.clone();
    vec_std.clear();
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vptr->get( i ), svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AsArrayMatches ) {
    double *darr = NCPA::arrays::zeros<double>( vecsize );
    vec_std.as_array( vecsize, darr );
    EXPECT_ARRAY_DOUBLE_EQ( vecsize, darr, svd );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ScaleWorks ) {
    vec_std.scale( 2.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ScaleByVectorWorks ) {
    vec_std.scale( vec_std );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ),
                          svd[ i ] * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, DotWorks ) {
    EXPECT_DOUBLE_EQ( vec_std.dot( vec_std ), svd_dot_svd );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AddWorksWithScalar ) {
    vec_std.add( 3.5 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AddWorksWithVector ) {
    vec_std.add( vec_std );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AddWorksWithModifiedVector ) {
    vec_std.add( vec_std, -1.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AddThrowsExceptionOnMismatchedSizes ) {
    vec_std.resize( 3 );
    EXPECT_THROW( { vec_initlist.add( vec_std ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, PlusEqualOperatorWorks ) {
    vec_std += vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, PlusEqualOperatorWorksWithScalar ) {
    vec_std += 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, MinusEqualOperatorWorks ) {
    vec_std -= vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, MinusEqualOperatorWorksWithScalar ) {
    vec_std -= 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] - 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, TimesEqualOperatorWorks ) {
    vec_std *= vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ),
                          svd[ i ] * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, TimesEqualOperatorWorksWithScalar ) {
    vec_std *= 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] * 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, DivideEqualOperatorWorksWithScalar ) {
    vec_std /= 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), svd[ i ] / 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, EqualityOperatorWorks ) {
    EXPECT_TRUE( vec_std == vec_initlist );
    vec_std.set( 0, 5.0 );
    EXPECT_FALSE( vec_std == vec_initlist );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, InequalityOperatorWorks ) {
    EXPECT_FALSE( vec_std != vec_initlist );
    vec_std.set( 0, 5.0 );
    EXPECT_TRUE( vec_std != vec_initlist );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, SwapFunctionWorks ) {
    vec_std.set( 0, 5.0 );
    ASSERT_DOUBLE_EQ( vec_initlist.get( 0 ), 0.0 );
    std::swap( vec_std, vec_initlist );
    EXPECT_DOUBLE_EQ( vec_initlist.get( 0 ), 5.0 );
    EXPECT_DOUBLE_EQ( vec_std.get( 0 ), 0.0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, AssignmentOperatorWorks ) {
    vec_std *= 2.0;
    vec      = vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( vec.get( i ), svd[ i ] * 2.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, CountNonzeroIndicesIsCorrect ) {
    ASSERT_EQ( vec_std.count_nonzero_indices(), 3 );
    vec_std.set( 0, 1.0 );
    EXPECT_EQ( vec_std.count_nonzero_indices(), 4 );
    ASSERT_EQ( vec_opposite.count_nonzero_indices(), 5 );
    vec_opposite.zero( 0 );
    EXPECT_EQ( vec_opposite.count_nonzero_indices(), 4 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperSizeWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperCopyConstructorWorks ) {
    wrapper2 = Vector<double>( wrapper1 );
    EXPECT_TRUE( wrapper1.equals( wrapper2 ) );
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_FALSE( wrapper1.equals( wrapper2 ) );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperSwapWorks ) {
    swap( wrapper1, wrapper2 );
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperEqualityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_TRUE( wrapper1 == wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperInequalityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_FALSE( wrapper1 != wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_TRUE( wrapper1 != wrapper2 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperAssignmentOperatorWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
    wrapper2 = wrapper1;
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperGetAllowsModification ) {
    wrapper1.get( 0 ) = 5.0;
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 5.0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperIndexOperatorAllowsRead ) {
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperIndexOperatorAllowsWrite ) {
    wrapper1[ 0 ] = 5.0;
    for ( size_t i = 1; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], svd[ i ] );
    }
    EXPECT_DOUBLE_EQ( wrapper1[ 0 ], 5.0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperAsStdReturnsCorrectVector ) {
    EXPECT_ARRAY_DOUBLE_EQ( vecsize, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperClearWorks ) {
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_THROW( { wrapper1[ vecsize ] = 2.0; }, std::logic_error );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperResizeWorks ) {
    wrapper1.resize( 8 );
    EXPECT_EQ( wrapper1.size(), 8 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], svd[ i ] );
    }
    for ( size_t i = vecsize; i < 8; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1[ i ], 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperAsArrayReturnsCorrectValues ) {
    double *d = NCPA::arrays::zeros<double>( vecsize );
    wrapper1.as_array( vecsize, d );
    EXPECT_ARRAY_DOUBLE_EQ( vecsize, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperSetSetsSingleValue ) {
    wrapper1.set( 0, 6 );
    EXPECT_DOUBLE_EQ( wrapper1.get( 0 ), 6.0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperSetReplacesContents ) {
    std::vector<double> newcontents( 9 );
    wrapper1.set( newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperSetWithArrayReplacesContents ) {
    double *newcontents = NCPA::arrays::zeros<double>( 9 );
    wrapper1.set( 9, newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperScaleWorks ) {
    wrapper1.scale( 2.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperScaleByVectorWorks ) {
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ),
                          svd[ i ] * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperDotWorks ) {
    EXPECT_DOUBLE_EQ( wrapper1.dot( wrapper1 ), svd_dot_svd );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperAddWorksWithScalar ) {
    wrapper1.add( 3.5 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), svd[ i ] + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperAddWorksWithVector ) {
    wrapper1.add( wrapper1 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperAddThrowsExceptionOnMismatchedSizes ) {
    wrapper2 = wrapper1;
    wrapper1.resize( 3 );
    EXPECT_THROW( { wrapper2.add( wrapper1 ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperPlusEqualOperatorWorks ) {
    wrapper1 += wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperPlusEqualOperatorWorksWithScalar ) {
    wrapper1 += 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), svd[ i ] + 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperNegativeOperatorWorks ) {
    wrapper2 = -wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper2.get( i ), -svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperMinusEqualOperatorWorks ) {
    wrapper1 -= wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperMinusEqualOperatorWorksWithScalar ) {
    wrapper1 -= 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), svd[ i ] - 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperTimesEqualOperatorWorks ) {
    wrapper1 *= wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ),
                          svd[ i ] * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperTimesEqualOperatorWorksWithScalar ) {
    wrapper1 *= 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), svd[ i ] * 3.5 );
    }
}

// TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperDivideEqualOperatorWorks ) {
//     wrapper1 /= wrapper1;
//     for ( size_t i = 0; i < vecsize; i++ ) {
//         EXPECT_DOUBLE_EQ( wrapper1.get( i ), 1.0 );
//     }
// }

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperDivideEqualOperatorWorksWithScalar ) {
    wrapper1 /= 3.5;
    for ( size_t i = 0; i < vecsize; i++ ) {
        EXPECT_DOUBLE_EQ( wrapper1.get( i ), svd[ i ] / 3.5 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperPlusOperatorWorks ) {
    wrapper2 = wrapper1 + wrapper_opposite;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], (double)(i+1) );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperMinusOperatorWorks ) {
    wrapper2 = wrapper1 - wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, WrapperTimesOperatorWorks ) {
    wrapper2 = wrapper1 * wrapper1;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], wrapper1[i]*wrapper1[i] );
    }

    wrapper2 = wrapper1 * wrapper_opposite;
    for (auto i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( wrapper2[i], 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ZeroWorksWithSingleIndex ) {
    vec_std.zero( 2 );
    EXPECT_DOUBLE_EQ( vec_std.get( 2 ), 0.0 );
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ZeroWorksWithVector ) {
    std::vector<size_t> zero_out;
    for (size_t i = 0; i < vecsize; i++) {
        zero_out.push_back( i );
    }
    vec_std.zero( zero_out );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 0.0 );
    }
}

TEST_F( NCPALinearAlgebraSparseVectorTest, ZeroWorksWithInitList ) {
    vec_std.zero( { 0, 1, 2, 3, 4, 5, 6, 7 } );
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_DOUBLE_EQ( vec_std.get( i ), 0.0 );
    }
}