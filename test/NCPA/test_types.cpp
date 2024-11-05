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

TEST( NCPAtypesTest, IsIterableAcceptsVector ) {
    std::vector<double> v;
    EXPECT_NO_THROW( { test_iterable( v ); } );
}

TEST( NCPAtypesTest, IsIterableAcceptsList ) {
    list<double> v;
    EXPECT_NO_THROW( { test_iterable( v ); } );
}

TEST( NCPAtypesTest, IsIterableAcceptsMap ) {
    map<double, double> v;
    EXPECT_NO_THROW( { test_iterable( v ); } );
}

TEST( NCPAtypesTest, IsIterableDoesNotAcceptScalar ) {
    double v;
    EXPECT_THROW( { test_iterable( v ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsIterableDoesNotAcceptArray ) {
    double v[ 5 ] = { 1, 2, 3, 4, 5 };
    EXPECT_THROW( { test_iterable( v ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsIterableDoesNotAcceptPointer ) {
    double *v = new double[ 5 ];
    EXPECT_THROW( { test_iterable( v ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsIterableOfAcceptsVectorOfSameType ) {
    vector<double> Container;
    double Element;
    EXPECT_NO_THROW( { test_iterable_of( Container, Element ); } );
}

TEST( NCPAtypesTest, IsIterableOfAllowsFloatToDouble ) {
    vector<double> Container;
    float Element;
    EXPECT_NO_THROW( { test_iterable_of( Container, Element ); } );
}

TEST( NCPAtypesTest, IsIterableOfAcceptsIntToDouble ) {
    vector<double> Container;
    int Element;
    EXPECT_NO_THROW( { test_iterable_of( Container, Element ); } );
}

TEST( NCPAtypesTest, IsIterableOfAcceptsDoubleToInt ) {
    vector<int> Container;
    double Element;
    EXPECT_NO_THROW( { test_iterable_of( Container, Element ); } );
}

TEST( NCPAtypesTest, IsIterableOfDoesNotAllowVectorOfInconvertibleType ) {
    vector<double> Container;
    std::string Element;
    EXPECT_THROW(
        { test_iterable_of( Container, Element ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsIterableOfDoesNotAllowComplexToDouble ) {
    vector<double> Container;
    complex<double> Element;
    EXPECT_THROW(
        { test_iterable_of( Container, Element ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsIterableOfAllowsDoubleToComplex ) {
    vector<complex<double>> Container;
    double Element;
    EXPECT_NO_THROW( { test_iterable_of( Container, Element ); } );
}

TEST( NCPAtypesTest, IsDeleteableAllowsObjects ) {
    std::string s;
    EXPECT_NO_THROW( { test_deleteable( s ); } );
}

TEST( NCPAtypesTest, IsDeleteableDoesNotAllowBaseTypes ) {
    double d;
    EXPECT_THROW( { test_deleteable( d ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsDereferenceableAllowsPointers ) {
    int *i;
    double *d;
    float *f;
    std::string *s;
    std::list<double> *ld;
    EXPECT_NO_THROW( { test_dereferenceable( i ); } );
    EXPECT_NO_THROW( { test_dereferenceable( d ); } );
    EXPECT_NO_THROW( { test_dereferenceable( f ); } );
    EXPECT_NO_THROW( { test_dereferenceable( s ); } );
    EXPECT_NO_THROW( { test_dereferenceable( ld ); } );
}

TEST( NCPAtypesTest, IsDereferenceableDoesNotAllowFundamentalTypes ) {
    int i;
    double d;
    float f;
    EXPECT_THROW( { test_dereferenceable( i ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( d ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( f ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsDereferenceableDoesNotAllowObjects ) {
    std::string s;
    std::list<double> ld;
    std::map<std::string, int> m;
    EXPECT_THROW( { test_dereferenceable( s ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( ld ); }, std::invalid_argument );
    EXPECT_THROW( { test_dereferenceable( m ); }, std::invalid_argument );
}

TEST( NCPAtypesTest, IsDereferenceableAllowsIterators ) {
    std::string s;
    std::list<double> ld;
    std::map<std::string, int> m;
    EXPECT_NO_THROW( { test_dereferenceable( s.begin() ); } );
    EXPECT_NO_THROW( { test_dereferenceable( ld.end() ); } );
    EXPECT_NO_THROW( { test_dereferenceable( m.begin() ); } );
}
