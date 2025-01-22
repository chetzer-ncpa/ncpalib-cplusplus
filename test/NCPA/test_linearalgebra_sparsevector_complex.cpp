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

#define _TEST_ARRAY_EQ_ EXPECT_ARRAY_COMPLEX_DOUBLE_EQ
#define _TEST_EQ_       EXPECT_COMPLEX_DOUBLE_EQ
typedef complex<double> test_t;

class NCPALinearAlgebraComplexSparseVectorTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            zero    = test_t( 0.0, 0.0 );
            testval = test_t( 6.0, -6.0 );

            svd          = std::vector<test_t>( {
                zero,
                zero,
                { 3.0, 3.0 },
                zero,
                { 5.0, 5.0 },
                { 6.0, 6.0 },
                zero,
                zero
            } );
            svd_opposite = std::vector<test_t>( {
                { 1, 1 },
                { 2, 2 },
                zero,
                { 4, 4 },
                zero,
                zero,
                { 7, 7 },
                { 8, 8 }
            } );
            vec_initlist = {
                zero,         zero,         { 3.0, 3.0 },
                                  zero,
                { 5.0, 5.0 },
                                  { 6.0, 6.0 },
                                  zero,         zero
            };

            vec_std      = sparse_vector<test_t>( svd );
            vec_opposite = sparse_vector<test_t>( svd_opposite );

            wrapper1         = Vector<test_t>( vec_std.clone() );
            wrapper_opposite = Vector<test_t>( vec_opposite.clone() );

            svd_dot_svd = zero;
            for ( auto it = svd.cbegin(); it != svd.cend(); ++it ) {
                svd_dot_svd += ( *it ) * ( *it );
            }
        }  // void TearDown() override {}

        // declare stuff here
        sparse_vector<test_t> vec_initlist, vec_std, vec_opposite,
            vec;
        std::vector<test_t> svd, svd_opposite;
        test_t svd_dot_svd, zero, testval;

        Vector<test_t> wrapper1, wrapper2, wrapper_opposite;
        size_t vecsize = 8;
};

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, DefaultConstructorWorks ) {
    EXPECT_EQ( vec.size(), 0 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, StdvectorConstructorWorks ) {
    // sparse_vector<test_t> vec2( svd );
    // vec_std.set( svd );
    EXPECT_EQ( vec_std.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, vec_std, svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        InitializerListConstructorWorks ) {
    EXPECT_EQ( vec_initlist.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, vec_initlist, svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, CopyConstructorWorks ) {
    vec = sparse_vector<test_t>( vec_std );
    EXPECT_EQ( vec.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, vec, svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        CopyConstructorMakesDeepCopy ) {
    vec = sparse_vector<test_t>( vec_std );
    vec_std.clear();
    EXPECT_EQ( vec.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, vec, svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ClearClearsContents ) {
    ASSERT_EQ( vec_std.size(), vecsize );
    vec_std.clear();
    EXPECT_EQ( vec_std.size(), 0 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, IsZeroWorksAsExpected ) {
    EXPECT_TRUE( vec.is_zero() );
    EXPECT_FALSE( vec_std.is_zero() );
    vec_std.zero();
    EXPECT_TRUE( vec_std.is_zero() );
    vec_std.clear();
    EXPECT_TRUE( vec_std.is_zero() );
    vec_std.resize(5).set( NCPA::math::zero<test_t>() );
    EXPECT_TRUE( vec_std.is_zero() );
    vec_std.set( 0, NCPA::math::one<test_t>() );
    EXPECT_FALSE( vec_std.is_zero() );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, SetSetsSingleValue ) {
    vec_std.set( 0, testval );
    _TEST_EQ_( vec_std.get( 0 ), testval );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, SetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    vec_std.set( newcontents );
    EXPECT_EQ( vec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( vec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        SetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    vec_std.set( 9, newcontents );
    EXPECT_EQ( vec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( vec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ResizeResizesArray ) {
    vec_std.resize( vecsize * 2 );
    _TEST_ARRAY_EQ_( vecsize, vec_std, svd );
    _TEST_EQ_( vec_std.get( vecsize ), zero );
    vec_std.resize( vecsize / 2 );
    for ( size_t i = 0; i < vec_std.size(); i++ ) {
        _TEST_EQ_( vec_std.get( i ), svd[ i ] );
    }
    EXPECT_ANY_THROW( { vec_std.get( vecsize ); } );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, AsStdMatches ) {
    _TEST_ARRAY_EQ_( vecsize, vec_std.as_std(), svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, CloneClones ) {
    std::unique_ptr<abstract_vector<test_t>> vptr = vec_std.clone();
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vptr->get( i ), svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, CloneClonesDeep ) {
    std::unique_ptr<abstract_vector<test_t>> vptr = vec_std.clone();
    vec_std.clear();
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vptr->get( i ), svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, AsArrayMatches ) {
    test_t *darr = NCPA::arrays::zeros<test_t>( vecsize );
    vec_std.as_array( vecsize, darr );
    _TEST_ARRAY_EQ_( vecsize, darr, svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ScaleWorks ) {
    vec_std.scale( 2.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ScaleByVectorWorks ) {
    vec_std.scale( vec_std );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] * svd[ i ] ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, DotWorks ) {
    _TEST_EQ_( vec_std.dot( vec_std ), svd_dot_svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, AddWorksWithScalar ) {
    vec_std.add( testval );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] + testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, AddWorksWithVector ) {
    vec_std.add( vec_std );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        AddWorksWithModifiedVector ) {
    vec_std.add( vec_std, -1.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        AddThrowsExceptionOnMismatchedSizes ) {
    vec_std.resize( 3 );
    EXPECT_THROW( { vec_initlist.add( vec_std ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, PlusEqualOperatorWorks ) {
    vec_std += vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        PlusEqualOperatorWorksWithScalar ) {
    vec_std += testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] + testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, MinusEqualOperatorWorks ) {
    vec_std -= vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        MinusEqualOperatorWorksWithScalar ) {
    vec_std -= testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] - testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, TimesEqualOperatorWorks ) {
    vec_std *= vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] * svd[ i ] ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        TimesEqualOperatorWorksWithScalar ) {
    vec_std *= testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] * testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        DivideEqualOperatorWorksWithScalar ) {
    vec_std /= testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), ( svd[ i ] / testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, EqualityOperatorWorks ) {
    EXPECT_TRUE( vec_std == vec_initlist );
    vec_std.set( 0, testval );
    EXPECT_FALSE( vec_std == vec_initlist );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, InequalityOperatorWorks ) {
    EXPECT_FALSE( vec_std != vec_initlist );
    vec_std.set( 0, testval );
    EXPECT_TRUE( vec_std != vec_initlist );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, SwapFunctionWorks ) {
    vec_std.set( 0, testval );
    ASSERT_COMPLEX_DOUBLE_EQ( vec_initlist.get( 0 ), zero );
    std::swap( vec_std, vec_initlist );
    _TEST_EQ_( vec_initlist.get( 0 ), testval );
    _TEST_EQ_( vec_std.get( 0 ), zero );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, AssignmentOperatorWorks ) {
    vec_std *= 2.0;
    vec      = vec_std;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec.get( i ), ( svd[ i ] * 2.0 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        CountNonzeroIndicesIsCorrect ) {
    ASSERT_EQ( vec_std.count_nonzero_indices(), 3 );
    vec_std.set( 0, 1.0 );
    EXPECT_EQ( vec_std.count_nonzero_indices(), 4 );
    ASSERT_EQ( vec_opposite.count_nonzero_indices(), 5 );
    vec_opposite.zero( 0 );
    EXPECT_EQ( vec_opposite.count_nonzero_indices(), 4 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperSizeWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperCopyConstructorWorks ) {
    wrapper2 = Vector<test_t>( wrapper1 );
    EXPECT_TRUE( wrapper1.equals( wrapper2 ) );
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_FALSE( wrapper1.equals( wrapper2 ) );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperSwapWorks ) {
    swap( wrapper1, wrapper2 );
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperEqualityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_TRUE( wrapper1 == wrapper2 );
    wrapper1.set( 0, testval );
    EXPECT_FALSE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperInequalityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_FALSE( wrapper1 != wrapper2 );
    wrapper1.set( 0, testval );
    EXPECT_TRUE( wrapper1 != wrapper2 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperAssignmentOperatorWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
    wrapper2 = wrapper1;
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperGetAllowsModification ) {
    wrapper1.get( 0 ) = testval;
    _TEST_EQ_( wrapper1.get( 0 ), testval );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperIndexingOperatorAllowsModification ) {
    _TEST_EQ_( wrapper1[0], wrapper1.get( 0 ) );
    wrapper1[0] = complex<double>( 5.0, -5.0 );
    _TEST_EQ_( wrapper1.get( 0 ), complex<double>( 5.0, -5.0 ) );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperIndexOperatorAllowsRead ) {
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1[ i ], svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperIndexOperatorAllowsWrite ) {
    wrapper1[ 0 ] = testval;
    for ( size_t i = 1; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1[ i ], svd[ i ] );
    }
    _TEST_EQ_( wrapper1[ 0 ], testval );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperAsStdReturnsCorrectVector ) {
    _TEST_ARRAY_EQ_( vecsize, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperClearWorks ) {
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_THROW( { wrapper1[ vecsize ] = 2.0; }, std::logic_error );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperResizeWorks ) {
    wrapper1.resize( 8 );
    EXPECT_EQ( wrapper1.size(), 8 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1[ i ], svd[ i ] );
    }
    for ( size_t i = vecsize; i < 8; i++ ) {
        _TEST_EQ_( wrapper1[ i ], zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperAsArrayReturnsCorrectValues ) {
    test_t *d = NCPA::arrays::zeros<test_t>( vecsize );
    wrapper1.as_array( vecsize, d );
    _TEST_ARRAY_EQ_( vecsize, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperSetSetsSingleValue ) {
    wrapper1.set( 0, testval );
    _TEST_EQ_( wrapper1.get( 0 ), testval );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperSetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    wrapper1.set( newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperSetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    wrapper1.set( 9, newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperScaleWorks ) {
    wrapper1.scale( 2.0 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperScaleByVectorWorks ) {
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] * svd[ i ] ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperDotWorks ) {
    _TEST_EQ_( wrapper1.dot( wrapper1 ), svd_dot_svd );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperAddWorksWithScalar ) {
    wrapper1.add( testval );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] + testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperAddWorksWithVector ) {
    wrapper1.add( wrapper1 );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperAddThrowsExceptionOnMismatchedSizes ) {
    wrapper2 = wrapper1;
    wrapper1.resize( 3 );
    EXPECT_THROW( { wrapper2.add( wrapper1 ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperPlusEqualOperatorWorks ) {
    wrapper1 += wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), 2.0 * svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperPlusEqualOperatorWorksWithScalar ) {
    wrapper1 += testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] + testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperNegativeOperatorWorks ) {
    wrapper2 = -wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2.get( i ), -svd[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperMinusEqualOperatorWorks ) {
    wrapper1 -= wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperMinusEqualOperatorWorksWithScalar ) {
    wrapper1 -= testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] - testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperTimesEqualOperatorWorks ) {
    wrapper1 *= wrapper1;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] * svd[ i ] ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperTimesEqualOperatorWorksWithScalar ) {
    wrapper1 *= testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] * testval ) );
    }
}

// TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
// WrapperDivideEqualOperatorWorks ) {
//     wrapper1 /= wrapper1;
//     for ( size_t i = 0; i < vecsize; i++ ) {
//         _TEST_EQ_( wrapper1.get( i ), 1.0 );
//     }
// }

TEST_F( NCPALinearAlgebraComplexSparseVectorTest,
        WrapperDivideEqualOperatorWorksWithScalar ) {
    wrapper1 /= testval;
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), ( svd[ i ] / testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperPlusOperatorWorks ) {
    wrapper2 = wrapper1 + wrapper_opposite;
    for (size_t i = 0; i < vecsize; i++) {
        EXPECT_COMPLEX_DOUBLE_EQ( wrapper2[ i ], test_t( i+1, i+1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperMinusOperatorWorks ) {
    wrapper2 = wrapper1 - wrapper1;
    for ( auto i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2[ i ], zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, WrapperTimesOperatorWorks ) {
    wrapper2 = wrapper1 * wrapper1;
    for ( auto i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2[ i ], ( wrapper1[ i ] * wrapper1[ i ] ) );
    }

    wrapper2 = wrapper1 * wrapper_opposite;
    for ( auto i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2[ i ], zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ZeroWorksWithSingleIndex ) {
    vec_std.zero( 2 );
    _TEST_EQ_( vec_std.get( 2 ), zero );
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ZeroWorksWithVector ) {
    std::vector<size_t> zero_out;
    for ( size_t i = 0; i < vecsize; i++ ) {
        zero_out.push_back( i );
    }
    vec_std.zero( zero_out );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexSparseVectorTest, ZeroWorksWithInitList ) {
    vec_std.zero( { 0, 1, 2, 3, 4, 5, 6, 7 } );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( vec_std.get( i ), zero );
    }
}
