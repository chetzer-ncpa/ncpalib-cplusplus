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

class NCPALinearAlgebraComplexDenseVectorTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
            svd      = std::vector<test_t>( {
                { 1.0, 1.0 },
                { 2.0, 2.0 },
                { 3.0, 3.0 },
                { 4.0, 4.0 },
                { 5.0, 5.0 }
            } );
            dvec_std = details::dense_vector<test_t>( svd );

            wrapper1 = Vector<test_t>( dvec_std.clone() );
            vecsize  = 5;
            testval  = test_t( 6.0, -6.0 );
            zero     = test_t( 0.0, 0.0 );
        }  // void TearDown() override {}

        // declare stuff here
        details::dense_vector<test_t> dvec_initlist {
            { 1.0, 1.0 },
            { 2.0, 2.0 },
            { 3.0, 3.0 },
            { 4.0, 4.0 },
            { 5.0, 5.0 }
        };
        details::dense_vector<test_t> dvec_std, dvec_copy, dvec;
        std::vector<test_t> svd;

        Vector<test_t> wrapper1, wrapper2;
        size_t vecsize;
        test_t testval, zero;
};

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, DefaultConstructorWorks ) {
    EXPECT_EQ( dvec.size(), 0 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, InitializerListConstructorWorks ) {
    EXPECT_EQ( dvec_initlist.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, dvec_initlist, svd );
    // for ( size_t i = 0; i < 5; i++ ) {
    //     _TEST_EQ_( dvec_initlist.get( i ), test_t( i+1,i+1) );
    // }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, StdVectorConstructorWorks ) {
    EXPECT_EQ( dvec_std.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, dvec_std, svd );
    // for ( size_t i = 0; i < 5; i++ ) {
    //     _TEST_EQ_( dvec_std.get( i ), test_t( i+1,i+1) );
    // }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, CopyConstructorWorks ) {
    dvec = details::dense_vector<test_t>( dvec_std );
    EXPECT_EQ( dvec.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, dvec_std, svd );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, CopyConstructorMakesDeepCopy ) {
    dvec = details::dense_vector<test_t>( dvec_std );
    dvec_std.clear();
    EXPECT_EQ( dvec.size(), vecsize );
    _TEST_ARRAY_EQ_( vecsize, dvec, svd );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ClearClearsContents ) {
    EXPECT_EQ( dvec_std.size(), vecsize );
    dvec_std.clear();
    EXPECT_EQ( dvec_std.size(), 0 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, SetSetsSingleValue ) {
    dvec_std.set( 0, testval );
    _TEST_EQ_( dvec_std.get( 0 ), testval );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, SetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    dvec_std.set( newcontents );
    EXPECT_EQ( dvec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, SetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    dvec_std.set( 9, newcontents );
    EXPECT_EQ( dvec_std.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, CountNonzeroIndicesIsCorrect ) {
    EXPECT_EQ( dvec_std.count_nonzero_indices(), vecsize );
    dvec_std.set( 0, zero );
    EXPECT_EQ( dvec_std.count_nonzero_indices(), 4 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, NonzeroIndicesReportedCorrectly ) {
    vector<size_t> inds = NCPA::arrays::index_vector<size_t>( vecsize );
    EXPECT_ARRAY_EQ( 5, dvec_std.nonzero_indices(), inds );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ResizeResizesArray ) {
    dvec_std.resize( 6 );
    _TEST_ARRAY_EQ_( vecsize, dvec_std, svd );
    _TEST_EQ_( dvec_std.get( vecsize ), zero );
    dvec_std.resize( 3 );
    _TEST_ARRAY_EQ_( 3, dvec_std, svd );
    EXPECT_ANY_THROW( { dvec_std.get( 4 ); } );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AsStdMatches ) {
    _TEST_ARRAY_EQ_( 5, dvec_std.as_std(), svd );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, CloneClones ) {
    std::unique_ptr<details::abstract_vector<test_t>> vptr = dvec_std.clone();
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( vptr->get( i ), test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, CloneClonesDeep ) {
    std::unique_ptr<details::abstract_vector<test_t>> vptr = dvec_std.clone();
    dvec_std.clear();
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( vptr->get( i ), test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AsArrayMatches ) {
    test_t *darr = NCPA::arrays::zeros<test_t>( vecsize );
    dvec_std.as_array( vecsize, darr );
    _TEST_ARRAY_EQ_( 5, darr, svd );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ScaleWorks ) {
    dvec_std.scale( 2.0 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), 2.0 * test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ScaleByVectorWorks ) {
    dvec_std.scale( dvec_std );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ),
                   ( test_t( i + 1, i + 1 ) * test_t( i + 1, i + 1 ) ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, DotWorks ) {
    test_t expected;
    for ( auto it = svd.begin(); it != svd.end(); ++it ) {
        expected += ( *it ) * ( *it );
    }
    _TEST_EQ_( dvec_std.dot( dvec_std ), expected );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AddWorksWithScalar ) {
    dvec_std.add( testval );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), ( test_t( i + 1, i + 1 ) + testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AddWorksWithVector ) {
    dvec_std.add( dvec_std );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), ( 2.0 * test_t( i + 1, i + 1 ) ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AddWorksWithModifiedVector ) {
    dvec_std.add( dvec_std, -1.0 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AddThrowsExceptionOnMismatchedSizes ) {
    dvec_std.resize( 3 );
    EXPECT_THROW( { dvec_initlist.add( dvec_std ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, PlusEqualOperatorWorks ) {
    dvec_std += dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), 2.0 * test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, PlusEqualOperatorWorksWithScalar ) {
    dvec_std += testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), ( test_t( i + 1, i + 1 ) + testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, MinusEqualOperatorWorks ) {
    dvec_std -= dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, MinusEqualOperatorWorksWithScalar ) {
    dvec_std -= testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), ( test_t( i + 1, i + 1 ) - testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, TimesEqualOperatorWorks ) {
    dvec_std *= dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ),
                   ( test_t( i + 1, i + 1 ) * test_t( i + 1, i + 1 ) ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, TimesEqualOperatorWorksWithScalar ) {
    dvec_std *= testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), ( test_t( i + 1, i + 1 ) * testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, DivideEqualOperatorWorksWithScalar ) {
    dvec_std /= testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), ( test_t( i + 1, i + 1 ) / testval ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, EqualityOperatorWorks ) {
    EXPECT_TRUE( dvec_std == dvec_initlist );
    dvec_std.set( 0, 5.0 );
    EXPECT_FALSE( dvec_std == dvec_initlist );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, InequalityOperatorWorks ) {
    EXPECT_FALSE( dvec_std != dvec_initlist );
    dvec_std.set( 0, 5.0 );
    EXPECT_TRUE( dvec_std != dvec_initlist );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, SwapFunctionWorks ) {
    dvec_std.set( 0, testval );
    _TEST_EQ_( dvec_initlist.get( 0 ), test_t( 1, 1 ) );
    std::swap( dvec_std, dvec_initlist );
    _TEST_EQ_( dvec_initlist.get( 0 ), testval );
    _TEST_EQ_( dvec_std.get( 0 ), test_t( 1, 1 ) );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, AssignmentOperatorWorks ) {
    dvec_std *= 2.0;
    dvec      = dvec_std;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( dvec.get( i ), ( test_t( i + 1, i + 1 ) * 2.0 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperSizeWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperCopyConstructorWorks ) {
    wrapper2 = Vector<test_t>( wrapper1 );
    EXPECT_TRUE( wrapper1.equals( wrapper2 ) );
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_FALSE( wrapper1.equals( wrapper2 ) );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperSwapWorks ) {
    swap( wrapper1, wrapper2 );
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_EQ( wrapper2.size(), vecsize );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperEqualityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_TRUE( wrapper1 == wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperInequalityOperatorWorks ) {
    wrapper2 = wrapper1;
    EXPECT_FALSE( wrapper1 != wrapper2 );
    wrapper1.set( 0, 5.0 );
    EXPECT_TRUE( wrapper1 != wrapper2 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperAssignmentOperatorWorks ) {
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_EQ( wrapper2.size(), 0 );
    EXPECT_FALSE( wrapper1 == wrapper2 );
    wrapper2 = wrapper1;
    EXPECT_EQ( wrapper2.size(), vecsize );
    EXPECT_EQ( wrapper1.size(), vecsize );
    EXPECT_TRUE( wrapper1 == wrapper2 );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperGetAllowsModification ) {
    wrapper1.get( 0 ) = testval;
    _TEST_EQ_( wrapper1.get( 0 ), testval );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperIndexOperatorAllowsRead ) {
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper1[ i ], test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperIndexOperatorAllowsWrite ) {
    wrapper1[ 0 ] = testval;
    for ( size_t i = 1; i < 5; i++ ) {
        _TEST_EQ_( wrapper1[ i ], test_t( i + 1, i + 1 ) );
    }
    _TEST_EQ_( wrapper1[ 0 ], testval );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperAsStdReturnsCorrectVector ) {
    _TEST_ARRAY_EQ_( 5, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperClearWorks ) {
    wrapper1.clear();
    EXPECT_EQ( wrapper1.size(), 0 );
    EXPECT_THROW( { wrapper1[ vecsize ] = 2.0; }, std::logic_error );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperResizeWorks ) {
    wrapper1.resize( 8 );
    EXPECT_EQ( wrapper1.size(), 8 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1[ i ], test_t( i + 1, i + 1 ) );
    }
    for ( size_t i = 5; i < 8; i++ ) {
        _TEST_EQ_( wrapper1[ i ], zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperAsArrayReturnsCorrectValues ) {
    test_t *d = NCPA::arrays::zeros<test_t>( vecsize );
    wrapper1.as_array( 5, d );
    _TEST_ARRAY_EQ_( 5, wrapper1.as_std(), svd );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperSetSetsSingleValue ) {
    wrapper1.set( 0, testval );
    _TEST_EQ_( wrapper1.get( 0 ), testval );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperSetReplacesContents ) {
    std::vector<test_t> newcontents( 9 );
    wrapper1.set( newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperSetWithArrayReplacesContents ) {
    test_t *newcontents = NCPA::arrays::zeros<test_t>( 9 );
    wrapper1.set( 9, newcontents );
    EXPECT_EQ( wrapper1.size(), 9 );
    for ( auto i = 0; i < 9; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperScaleWorks ) {
    wrapper1.scale( 2.0 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), 2.0 * test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperScaleByVectorWorks ) {
    wrapper1.scale( wrapper1 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ),
                   (test_t( i + 1, i + 1 ) * test_t( i + 1, i + 1 )) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperDotWorks ) {
    test_t expected;
    for ( auto it = svd.begin(); it != svd.end(); ++it ) {
        expected += ( *it ) * ( *it );
    }
    _TEST_EQ_( wrapper1.dot( wrapper1 ), expected );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperAddWorksWithScalar ) {
    wrapper1.add( testval );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), (test_t( i + 1, i + 1 ) + testval) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperAddWorksWithVector ) {
    wrapper1.add( wrapper1 );
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), 2.0 * test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperAddThrowsExceptionOnMismatchedSizes ) {
    wrapper2 = wrapper1;
    wrapper1.resize( 3 );
    EXPECT_THROW( { wrapper2.add( wrapper1 ); }, std::invalid_argument );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperPlusEqualOperatorWorks ) {
    wrapper1 += wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), 2.0 * test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperPlusEqualOperatorWorksWithScalar ) {
    wrapper1 += testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), (test_t( i + 1, i + 1 ) + testval) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperNegativeOperatorWorks ) {
    wrapper2 = -wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper2.get( i ), -test_t( i + 1, i + 1 ) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperMinusEqualOperatorWorks ) {
    wrapper1 -= wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperMinusEqualOperatorWorksWithScalar ) {
    wrapper1 -= testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), (test_t( i + 1, i + 1 ) - testval) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperTimesEqualOperatorWorks ) {
    wrapper1 *= wrapper1;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ),
                   (test_t( i + 1, i + 1 ) * test_t( i + 1, i + 1 )) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperTimesEqualOperatorWorksWithScalar ) {
    wrapper1 *= testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), (test_t( i + 1, i + 1 ) * testval) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperDivideEqualOperatorWorksWithScalar ) {
    wrapper1 /= testval;
    for ( size_t i = 0; i < 5; i++ ) {
        _TEST_EQ_( wrapper1.get( i ), (test_t( i + 1, i + 1 ) / testval) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperPlusOperatorWorks ) {
    wrapper2 = wrapper1 + wrapper1;
    for ( auto i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2[ i ], 2.0 * wrapper1[ i ] );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperMinusOperatorWorks ) {
    wrapper2 = wrapper1 - wrapper1;
    for ( auto i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2[ i ], zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, WrapperTimesOperatorWorks ) {
    wrapper2 = wrapper1 * wrapper1;
    for ( auto i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( wrapper2[ i ], (wrapper1[ i ] * wrapper1[ i ]) );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ZeroWorksWithSingleIndex ) {
    dvec_std.zero( 2 );
    _TEST_EQ_( dvec_std.get( 2 ), zero );
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ZeroWorksWithVector ) {
    std::vector<size_t> zero_out;
    for ( size_t i = 0; i < vecsize; i++ ) {
        zero_out.push_back( i );
    }
    dvec_std.zero( zero_out );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), zero );
    }
}

TEST_F( NCPALinearAlgebraComplexDenseVectorTest, ZeroWorksWithInitList ) {
    dvec_std.zero( { 0, 1, 2, 3, 4 } );
    for ( size_t i = 0; i < vecsize; i++ ) {
        _TEST_EQ_( dvec_std.get( i ), zero );
    }
}
