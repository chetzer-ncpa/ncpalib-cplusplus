#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/defines.hpp"

#include <memory>

DECLARE_CLASS_AND_SWAP_2NAMESPACE( ValidationTest, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( ValidationTestSuite, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( NullaryValidationTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( TypedValidationTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( UnaryValidationTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( BinaryValidationTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( ListValidationTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( IsEqualToTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( IsNotEqualToTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( IsNotOneOfTest, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( IsOneOfTest, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( WasSetTest, NCPA, config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( IsBetweenTest, NCPA, config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( IsGreaterThanTest, NCPA,
                                                  config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( IsGreaterThanOrEqualToTest,
                                                  NCPA, config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( IsLessThanTest, NCPA,
                                                  config )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_2NAMESPACE( IsLessThanOrEqualToTest,
                                                  NCPA, config )

namespace NCPA {
    namespace config {
        typedef std::unique_ptr<ValidationTest> test_ptr_t;

        // convenience functions
        template<typename T>
        test_ptr_t IsLessThan( const T& val ) {
            return test_ptr_t( new IsLessThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsLessThan( const T& val, bool equalOK ) {
            return equalOK
                     ? test_ptr_t( new IsLessThanOrEqualToTest<T>( val ) )
                     : test_ptr_t( new IsLessThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsGreaterThan( const T& val ) {
            return test_ptr_t( new IsGreaterThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsGreaterThan( const T& val, bool equalOK ) {
            return equalOK
                     ? test_ptr_t( new IsGreaterThanOrEqualToTest<T>( val ) )
                     : test_ptr_t( new IsGreaterThanTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsEqualTo( const T& val ) {
            return test_ptr_t( new IsEqualToTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsNotEqualTo( const T& val ) {
            return test_ptr_t( new IsNotEqualToTest<T>( val ) );
        }

        template<typename T>
        test_ptr_t IsZero() {
            return IsEqualTo<T>( (T)0 );
        }

        template<typename T>
        test_ptr_t IsNotZero() {
            return IsNotEqualTo<T>( (T)0 );
        }

        template<typename T>
        test_ptr_t IsNegative() {
            return IsLessThan<T>( (T)0 );
        }

        template<typename T>
        test_ptr_t IsPositive() {
            return IsGreaterThan<T>( (T)0 );
        }

        inline test_ptr_t IsNotEmptyString() {
            return IsNotEqualTo<std::string>( "" );
        }

        inline test_ptr_t IsEmptyString() {
            return IsEqualTo<std::string>( "" );
        }
    }  // namespace config
}  // namespace NCPA
