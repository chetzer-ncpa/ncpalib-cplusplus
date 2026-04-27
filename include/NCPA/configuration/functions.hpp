#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/types.hpp"

#include <string>

namespace NCPA {
    namespace config {

        template<typename PARAMTYPE>
        TypedParameter<PARAMTYPE> *make_typed( Parameter& p ) {
            return dynamic_cast<TypedParameter<PARAMTYPE> *>( &p );
        }

        template<typename PARAMTYPE>
        ScalarParameter<PARAMTYPE> *make_scalar( Parameter& p ) {
            return dynamic_cast<ScalarParameter<PARAMTYPE> *>( &p );
        }

        template<typename PARAMTYPE>
        VectorParameter<PARAMTYPE> *make_vector( Parameter& p ) {
            return dynamic_cast<VectorParameter<PARAMTYPE> *>( &p );
        }

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
