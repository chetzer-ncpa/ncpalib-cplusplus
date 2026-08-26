#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <functional>
#include <string>
#include <utility>

namespace NCPA {
    namespace config {

        // template<typename T>
        // using validation_function_t = std::function<validation_status_t( T )>;

        class Validation : public Cloneable<Validation> {
            public:
                Validation() {}

                virtual ~Validation() {}

                Validation( const Validation& other ) {}

                Validation( Validation&& other ) noexcept {
                    swap( *this, other );
                }

                friend void swap( Validation& a, Validation& b ) noexcept {}

                // template<typename T>
                // static const PredefinedValidation<T> predefined;

                // protected:
                //     validation_status_t _status;
        };

        

        // template<typename T>
        // validation_function_t<T>& _binary_comparison( const std::function<bool( T, T )>& cmp,
        //                                const T& a, std::string msg ) {
        //     using std::to_string;
        //     static validation_function_t<T> func
        //         = [ cmp, a, msg ]( T d ) -> validation_status_t {
        //         return cmp( d, a ) ? validation_status_t { test_result_t::PASSED, "" }
        //                            : validation_status_t { test_result_t::FAILED, msg };
        //     };
        //     return func;
        // }

        // template<typename T>
        // TypedValidation<T> is_less_than( const T& val, bool equalok = false,
        //                                  std::string errmsg = "" ) {
        //     using std::to_string;
        //     if (errmsg.empty()) {
        //         std::string addendum = ( equalok ? "or equal to " : "" );
        //         errmsg = "Supplied value must be less than " + addendum
        //                + to_string( val );
        //     }
        //     return std::move( TypedValidation<T>( _binary_comparison<T>(
        //         ( equalok ? []( T x, T y ) { return x <= y; }
        //                   : []( T x, T y ) { return x < y; } ),
        //         val, errmsg ) ) );
        // }

        // template<typename T>
        // TypedValidation<T> is_equal_to( const T& val,
        //                                 std::string errmsg = "" ) {
        //     using std::to_string;
        //     if (errmsg.empty()) {
        //         errmsg = "Supplied value must be equal to " + to_string( val );
        //     }
        //     return std::move( TypedValidation<T>( _binary_comparison<T>(
        //         []( T x, T y ) { return x == y; }, val, errmsg ) ) );
        // }

        // template<typename T>
        // TypedValidation<T> is_not_equal_to( const T& val,
        //                                     std::string errmsg = "" ) {
        //     using std::to_string;
        //     if (errmsg.empty()) {
        //         errmsg = "Supplied value must be equal to " + to_string( val );
        //     }
        //     return std::move( TypedValidation<T>( _binary_comparison<T>(
        //         []( T x, T y ) { return x != y; }, val, errmsg ) ) );
        // }

        // template<typename T>
        // TypedValidation<T> is_greater_than( const T& val, bool equalok = false,
        //                                     std::string errmsg = "" ) {
        //     using std::to_string;
        //     if (errmsg.empty()) {
        //         errmsg = "Supplied value must be greater than "
        //                + to_string( val );
        //     }
        //     return std::move( TypedValidation<T>( _binary_comparison<T>(
        //                 ( equalok ? []( T x, T y ) { std::cout << x << " >= " << y << std::endl; return x >= y; }
        //                           : []( T x, T y ) { std::cout << x << " > " << y << std::endl; return x > y; } ),
        //                 val, errmsg ) ) );
        // };

        // template<typename T>
        // TypedValidation<T> is_positive() {
        //     return std::move(
        //         is_greater_than( static_cast<T>( 0 ), false,
        //                          "Supplied value must be positive." ) );
        // }

        // template<typename T>
        // TypedValidation<T> is_negative() {
        //     return std::move(
        //         is_less_than( static_cast<T>( 0 ), false,
        //                       "Supplied value must be negative." ) );
        // }

        // template<typename T>
        // TypedValidation<T> is_zero() {
        //     return std::move( is_equal_to( static_cast<T>( 0 ),
        //                                    "Supplied value must be zero" ) );
        // }

        // template<typename T>
        // TypedValidation<T> is_nonzero() {
        //     return std::move( is_not_equal_to(
        //         static_cast<T>( 0 ), "Supplied value must be nonzero" ) );
        // }

    }  // namespace config
}  // namespace NCPA

// template<typename T>
// class PredefinedValidation {
//     public:
//         using validation_function_t = std::function<validation_status_t( T )>;

//         static validation_function_t& _binary_comparison(
//             const std::function<bool( T, T )>& cmp, const T& a,
//             std::string msg ) {
//             using std::to_string;
//             static validation_function_t func
//                 = [ cmp, a, msg ]( T d ) -> validation_status_t {
//                 return cmp( d, a )
//                          ? validation_status_t {
//                          test_result_t::PASSED,
//                                                  "" }
//                          : validation_status_t {
//                          test_result_t::FAILED,
//                                                  msg };
//             };
//             return func;
//         }

//         static TypedValidation<T> is_equal_to( const T& val,
//                                                std::string errmsg
//                                                = "" ) {
//             using std::to_string;
//             if (errmsg.empty()) {
//                 errmsg = "Supplied value must be equal to "
//                        + to_string( val );
//             }
//             return std::move( TypedValidation<T>(
//             _binary_comparison(
//                 []( T x, T y ) { return x == y; }, val, errmsg ) )
//                 );
//         }

//         static TypedValidation<T> is_not_equal_to( const T& val,
//                                                    std::string
//                                                    errmsg = "" ) {
//             using std::to_string;
//             if (errmsg.empty()) {
//                 errmsg = "Supplied value must be equal to "
//                        + to_string( val );
//             }
//             return std::move( TypedValidation<T>(
//             _binary_comparison(
//                 []( T x, T y ) { return x != y; }, val, errmsg ) )
//                 );
//         }

//         static TypedValidation<T> is_less_than( const T& val,
//                                                 bool equalok =
//                                                 false, std::string
//                                                 errmsg = "" ) {
//             using std::to_string;
//             if (errmsg.empty()) {
//                 std::string addendum
//                     = ( equalok ? "or equal to " : "" );
//                 errmsg = "Supplied value must be less than " +
//                 addendum
//                        + to_string( val );
//             }
//             return std::move( TypedValidation<T>(
//             _binary_comparison(
//                 ( equalok ? []( T x, T y ) { return x <= y; }
//                           : []( T x, T y ) { return x < y; } ),
//                 val, errmsg ) ) );
//         }

//         static TypedValidation<T> is_greater_than( const T& val,
//                                                    bool equalok
//                                                    = false,
//                                                    std::string
//                                                    errmsg = "" ) {
//             using std::to_string;
//             if (errmsg.empty()) {
//                 errmsg = "Supplied value must be greater than "
//                        + to_string( val );
//             }
//             return std::move( TypedValidation<T>(
//             _binary_comparison(
//                 ( equalok ? []( T x, T y ) { std::cout << x << " >=
//                 " << y << std::endl; return x >= y; }
//                           : []( T x, T y ) { std::cout << x << " > "
//                           << y << std::endl; return x > y; } ),
//                 val, errmsg ) ) );
//         }

//         static TypedValidation<T> is_positive() {
//             return std::move( is_greater_than(
//                 static_cast<T>( 0 ), false,
//                 "Supplied value must be positive." ) );
//         }

//         static TypedValidation<T> is_negative() {
//             return std::move(
//                 is_less_than( static_cast<T>( 0 ), false,
//                               "Supplied value must be negative." )
//                               );
//         }

//         static TypedValidation<T> is_zero() {
//             return std::move( is_equal_to(
//                 static_cast<T>( 0 ), "Supplied value must be zero" )
//                 );
//         }

//         static TypedValidation<T> is_nonzero() {
//             return std::move(
//                 is_not_equal_to( static_cast<T>( 0 ),
//                                  "Supplied value must be nonzero" )
//                                  );
//         }
// };
