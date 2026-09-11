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

        template<typename T>
        validation_function_t<T> _binary_comparison(
            const std::function<bool( T, T )>& cmp, const T& a,
            std::string msg ) {
            using std::to_string;
            validation_function_t<T> func
                = [ cmp, a, msg ]( T d ) -> validation_status_t {
                return cmp( d, a )
                         ? validation_status_t { test_result_t::PASSED, "" }
                         : validation_status_t { test_result_t::FAILED, msg };
            };
            return func;
        }

        template<typename T>
        TypedValidation<T> is_less_than( const T& val, bool equalok = false,
                                         std::string errmsg = "" ) {
            using std::to_string;
            if (errmsg.empty()) {
                std::string addendum = ( equalok ? "or equal to " : "" );
                errmsg = "Supplied value must be less than " + addendum
                       + to_string( val );
            }
            return std::move( TypedValidation<T>( _binary_comparison<T>(
                ( equalok ? []( T x, T y ) { return x <= y; }
                          : []( T x, T y ) { return x < y; } ),
                val, errmsg ) ) );
        }

        template<typename T>
        TypedValidation<T> is_equal_to( const T& val,
                                        std::string errmsg = "" ) {
            using std::to_string;
            if (errmsg.empty()) {
                errmsg = "Supplied value must be equal to " + to_string( val );
            }
            return std::move( TypedValidation<T>( _binary_comparison<T>(
                []( T x, T y ) { return x == y; }, val, errmsg ) ) );
        }

        template<typename T>
        TypedValidation<T> is_not_equal_to( const T& val,
                                            std::string errmsg = "" ) {
            using std::to_string;
            if (errmsg.empty()) {
                errmsg = "Supplied value must be equal to " + to_string( val );
            }
            return std::move( TypedValidation<T>( _binary_comparison<T>(
                []( T x, T y ) { return x != y; }, val, errmsg ) ) );
        }

        template<typename T>
        TypedValidation<T> is_greater_than( const T& val, bool equalok = false,
                                            std::string errmsg = "" ) {
            using std::to_string;
            if (errmsg.empty()) {
                errmsg = "Supplied value must be greater than "
                       + to_string( val );
            }
            return std::move( TypedValidation<T>( _binary_comparison<T>(
                ( equalok ? []( T x, T y ) { return x >= y; }
                          : []( T x, T y ) { return x > y; } ),
                val, errmsg ) ) );
        };

        template<typename T>
        TypedValidation<T> is_positive() {
            return std::move(
                is_greater_than( static_cast<T>( 0 ), false,
                                 "Supplied value must be positive." ) );
        }

        template<typename T>
        TypedValidation<T> is_negative() {
            return std::move(
                is_less_than( static_cast<T>( 0 ), false,
                              "Supplied value must be negative." ) );
        }

        template<typename T>
        TypedValidation<T> is_zero() {
            return std::move( is_equal_to( static_cast<T>( 0 ),
                                           "Supplied value must be zero" ) );
        }

        template<typename T>
        TypedValidation<T> is_nonzero() {
            return std::move( is_not_equal_to(
                static_cast<T>( 0 ), "Supplied value must be nonzero" ) );
        }

        inline validation_status_t merge( const validation_status_t& a, const validation_status_t& b ) {
            validation_status_t merged;
            if (a.result == test_result_t::FAILED || b.result == test_result_t::FAILED) {
                merged.result = test_result_t::FAILED;
            } else if (a.result == test_result_t::PENDING || b.result == test_result_t::PENDING) {
                merged.result = test_result_t::PENDING;
            } else if (a.result == test_result_t::NONE && b.result == test_result_t::NONE) { 
                merged.result = test_result_t::NONE;
            } else {
                merged.result = test_result_t::PASSED;
            }

            if (!a.message.empty()) {
                if (b.message.empty()) {
                    merged.message = a.message;
                } else {
                    merged.message = a.message + "\n" + b.message;
                }
            } else {
                if (!b.message.empty()) {
                    merged.message = b.message;
                }
            }
            
            return merged;
        }
    }  // namespace config
}  // namespace NCPA
