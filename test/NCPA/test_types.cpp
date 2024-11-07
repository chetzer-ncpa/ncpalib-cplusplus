#include "NCPA/types.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <list>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest-spi.h>

using namespace std;
using namespace testing;
using namespace NCPA::types;

class NCPAtypesTest : public ::testing::Test {
    protected:
        void SetUp() override {  // define stuff here
           
        }  // void TearDown() override {}

        // declare stuff here
        int i;
    unsigned int ui;
    long l;
    unsigned long ul;
    short s;
    unsigned short us;
    long long ll;
    unsigned long long ull;
    size_t st;
    float f;
    double d;
    long double ld;
    std::complex<double> cd;
    std::complex<float> cf;
    std::complex<long double> cld;

    std::string ss;
    std::list<double> sld;
    std::map<std::string, int> msi;
    std::map<double,double> mdd;
    std::vector<double> svd;
    std::vector<int> svi;
    std::vector<std::complex<double>> svcd;

    double ad[ 5 ] = { 1, 2, 3, 4, 5 };
double *pd = new double[ 5 ];
    int *pi;
    float *pf;
    std::string *pss;
    std::list<double> *psld;
    
};

template<typename T>
void test_deleteable(
    T obj, typename std::enable_if<is_deleteable<T>::value, int>::type ENABLER
           = 0 ) {}

template<typename T>
void test_deleteable(
    T obj, typename std::enable_if<!is_deleteable<T>::value, int>::type ENABLER
           = 0 ) {
    throw std::invalid_argument( "Does not satisfy is_deleteable" );
}

template<typename T>
void test_iterable(
    T obj,
    typename std::enable_if<is_iterable<T>::value, int>::type ENABLER = 0 ) {}

template<typename T>
void test_iterable(
    T obj,
    typename std::enable_if<!is_iterable<T>::value, int>::type ENABLER = 0 ) {
    throw std::invalid_argument( "Does not satisfy is_iterable" );
}

template<typename T>
void test_complex(
    T obj,
    typename std::enable_if<is_complex<T>::value, int>::type ENABLER = 0 ) {}

template<typename T>
void test_complex(
    T obj,
    typename std::enable_if<!is_complex<T>::value, int>::type ENABLER = 0 ) {
    throw std::invalid_argument( "Does not satisfy is_complex" );
}


template<typename T>
void test_numeric(
    T obj,
    typename std::enable_if<is_numeric<T>::value, int>::type ENABLER = 0 ) {}

template<typename T>
void test_numeric(
    T obj,
    typename std::enable_if<!is_numeric<T>::value, int>::type ENABLER = 0 ) {
    throw std::invalid_argument( "Does not satisfy is_numeric" );
}

template<typename T>
void test_dereferenceable(
    T obj,
    typename std::enable_if<is_dereferenceable<T>::value, int>::type ENABLER
    = 0 ) {}

template<typename T>
void test_dereferenceable(
    T obj,
    typename std::enable_if<!is_dereferenceable<T>::value, int>::type ENABLER
    = 0 ) {
    throw std::invalid_argument( "Does not satisfy is_dereferenceable" );
}

template<typename T, typename U>
void test_iterable_of(
    T obj, U prim,
    typename std::enable_if<!is_iterable_of<T, U>::value, int>::type ENABLER
    = 0 ) {
    throw std::invalid_argument( "Does not satisfy is_iterable_of" );
}

template<typename T, typename U>
void test_iterable_of(
    T obj, U prim,
    typename std::enable_if<is_iterable_of<T, U>::value, int>::type ENABLER
    = 0 ) {}

TEST_F( NCPAtypesTest, IsIterableAcceptsVector ) {
    EXPECT_NO_THROW( { test_iterable( svd ); } );
}

TEST_F( NCPAtypesTest, IsIterableAcceptsList ) {
     EXPECT_NO_THROW( { test_iterable( sld ); } );
}

TEST_F( NCPAtypesTest, IsIterableAcceptsMap ) {
    EXPECT_NO_THROW( { test_iterable( mdd ); } );
}

TEST_F( NCPAtypesTest, IsIterableDoesNotAcceptScalar ) {
    EXPECT_THROW( { test_iterable( d ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsIterableDoesNotAcceptArray ) {
    EXPECT_THROW( { test_iterable( ad ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsIterableDoesNotAcceptPointer ) {
    
    EXPECT_THROW( { test_iterable( pd ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsIterableOfAcceptsVectorOfSameType ) {
    
    EXPECT_NO_THROW( { test_iterable_of( svd, d ); } );
}

TEST_F( NCPAtypesTest, IsIterableOfAllowsFloatToDouble ) {
    EXPECT_NO_THROW( { test_iterable_of( svd, f ); } );
}

TEST_F( NCPAtypesTest, IsIterableOfAcceptsIntToDouble ) {
    EXPECT_NO_THROW( { test_iterable_of( svd, i ); } );
}

TEST_F( NCPAtypesTest, IsIterableOfAcceptsDoubleToInt ) {
    EXPECT_NO_THROW( { test_iterable_of( svi, d ); } );
}

TEST_F( NCPAtypesTest, IsIterableOfDoesNotAllowVectorOfInconvertibleType ) {
   EXPECT_THROW(
        { test_iterable_of( svd, ss ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsIterableOfDoesNotAllowComplexToDouble ) {
    EXPECT_THROW(
        { test_iterable_of( svd, cd ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsIterableOfAllowsDoubleToComplex ) {
    EXPECT_NO_THROW( { test_iterable_of( svcd, d ); } );
}

TEST_F( NCPAtypesTest, IsDeleteableAllowsObjects ) {
    EXPECT_NO_THROW( { test_deleteable( ss ); } );
}

TEST_F( NCPAtypesTest, IsDeleteableDoesNotAllowBaseTypes ) {
    EXPECT_THROW( { test_deleteable( d ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsDereferenceableAllowsPointers ) {
    EXPECT_NO_THROW( { test_dereferenceable( pi ); } );
    EXPECT_NO_THROW( { test_dereferenceable( pd ); } );
    EXPECT_NO_THROW( { test_dereferenceable( pf ); } );
    EXPECT_NO_THROW( { test_dereferenceable( pss ); } );
    EXPECT_NO_THROW( { test_dereferenceable( psld ); } );
}

TEST_F( NCPAtypesTest, IsDereferenceableDoesNotAllowFundamentalTypes ) {
    EXPECT_THROW( { test_dereferenceable( i ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( d ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( f ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsDereferenceableDoesNotAllowObjects ) {
    EXPECT_THROW( { test_dereferenceable( ss ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( sld ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( msi ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsDereferenceableAllowsIterators ) {
    EXPECT_NO_THROW( { test_dereferenceable( ss.begin() ); } );
    EXPECT_NO_THROW( { test_dereferenceable( sld.end() ); } );
    EXPECT_NO_THROW( { test_dereferenceable( msi.begin() ); } );
}

TEST_F( NCPAtypesTest, IsComplexAllowsComplexTypes ) {
    EXPECT_NO_THROW( { test_complex( cd ); } );
    EXPECT_NO_THROW( { test_complex( cf ); } );
    EXPECT_NO_THROW( { test_complex( cld ); } );
}

TEST_F( NCPAtypesTest, IsComplexDoesNotAllowArithmeticTypes ) {
    EXPECT_THROW( { test_complex( d ); }, std::invalid_argument );
    EXPECT_THROW( { test_complex( f ); }, std::invalid_argument );
    EXPECT_THROW( { test_complex( ld ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsComplexDoesNotAllowOtherClasses ) {
    EXPECT_THROW( { test_complex( ss ); }, std::invalid_argument );
    EXPECT_THROW( { test_complex( svd ); }, std::invalid_argument );
}

TEST_F( NCPAtypesTest, IsNumericAllowsIntegerTypes ) {
    EXPECT_NO_THROW( { test_numeric( i ); } );
    EXPECT_NO_THROW( { test_numeric( ui ); } );
    EXPECT_NO_THROW( { test_numeric( l ); } );
    EXPECT_NO_THROW( { test_numeric( ul ); } );
    EXPECT_NO_THROW( { test_numeric( s ); } );
    EXPECT_NO_THROW( { test_numeric( us ); } );
    EXPECT_NO_THROW( { test_numeric( ll); } );
    EXPECT_NO_THROW( { test_numeric( ull ); } );
    EXPECT_NO_THROW( { test_numeric( st ); } );
}

TEST_F( NCPAtypesTest, IsNumericAllowsFloatingPointTypes ) {
    EXPECT_NO_THROW( { test_numeric( f ); } );
    EXPECT_NO_THROW( { test_numeric( d ); } );
    EXPECT_NO_THROW( { test_numeric( ld ); } );
    EXPECT_NO_THROW( { test_numeric( cd ); } );
}