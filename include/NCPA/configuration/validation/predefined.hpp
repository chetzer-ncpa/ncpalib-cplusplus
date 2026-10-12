#pragma once

#include "NCPA/configuration/validation/declarations.hpp"
#include "NCPA/configuration/validation/TypedValidation.hpp"
#include "NCPA/strings.hpp"

#include <functional>
#include <memory>

namespace NCPA {
    namespace config {
        namespace tests {

            // might end up using these
            // template<typename T>
            // using _unary_cmp_t = std::function<bool( T )>;
            // template<typename T>
            // using _binary_cmp_t = std::function<bool( T, T )>;

            // template<typename T>
            // static _binary_cmp_t<T> _gt_cmp( T a, T b ) {
            //     _binary_cmp_t<T> func = [a,b]() { return (a > b); };
            //     return func;
            // }

            template<typename T>
            static validation_function_t<T> _eq( T val,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ val, if_failed ]( T a ) {
                    return a == val ? ValidationStatus( test_result_t::PASSED )
                                    : ValidationStatus( test_result_t::FAILED,
                                                        if_failed );
                };
                return func;
            };

            template<typename T>
            static validation_function_t<T> _ge( T val,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ val, if_failed ]( T a ) {
                    return a >= val ? ValidationStatus( test_result_t::PASSED )
                                    : ValidationStatus( test_result_t::FAILED,
                                                        if_failed );
                };
                return func;
            };

            template<typename T>
            static validation_function_t<T> _gt( T val,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ val, if_failed ]( T a ) {
                    return a > val ? ValidationStatus( test_result_t::PASSED )
                                   : ValidationStatus( test_result_t::FAILED,
                                                       if_failed );
                };
                return func;
            };

            template<typename T>
            static validation_function_t<T> _in( const std::vector<T>& vals,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ vals, if_failed ]( T a ) {
                    for (auto& cmpval : vals) {
                        if (a == cmpval) {
                            return ValidationStatus( test_result_t::PASSED );
                        }
                    }
                    return ValidationStatus( test_result_t::FAILED,
                                             if_failed );
                };
            }

            template<typename T>
            static validation_function_t<T> _not_in(
                const std::vector<T>& vals,
                const std::string& if_failed = "" ) {
                validation_function_t<T> func = [ vals, if_failed ]( T a ) {
                    for (auto& cmpval : vals) {
                        if (a == cmpval) {
                            return ValidationStatus( test_result_t::FAILED,
                                                     if_failed );
                        }
                    }
                    return ValidationStatus( test_result_t::PASSED );
                };
            }

            template<typename T>
            static validation_function_t<T> _le( T val,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ val, if_failed ]( T a ) {
                    return a <= val ? ValidationStatus( test_result_t::PASSED )
                                    : ValidationStatus( test_result_t::FAILED,
                                                        if_failed );
                };
                return func;
            };

            template<typename T>
            static validation_function_t<T> _lt( T val,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ val, if_failed ]( T a ) {
                    return a < val ? ValidationStatus( test_result_t::PASSED )
                                   : ValidationStatus( test_result_t::FAILED,
                                                       if_failed );
                };
                return func;
            };

            template<typename T>
            static validation_function_t<T> _ne( T val,
                                                 const std::string& if_failed
                                                 = "" ) {
                validation_function_t<T> func = [ val, if_failed ]( T a ) {
                    return a != val ? ValidationStatus( test_result_t::PASSED )
                                    : ValidationStatus( test_result_t::FAILED,
                                                        if_failed );
                };
                return func;
            };

            template<typename T>
            static validation_ptr_t<T> eq( T val,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be > " << val << ".";
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _eq( val, oss.str() ) ) );
            }

            template<typename T>
            static validation_ptr_t<T> gt( T val,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be > " << val << ".";
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _gt( val, oss.str() ) ) );
            }

            template<typename T>
            static validation_ptr_t<T> ge( T val,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be >= " << val << ".";
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _ge( val, oss.str() ) ) );
            }

            template<typename T>
            static validation_ptr_t<T> in( const std::vector<T>& vals,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be one of { "
                        << NCPA::strings::join( vals, ", " ) << " }."
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _in( vals, oss.str() ) ) );
            }

            template<typename T>
            static validation_ptr_t<T> le( T val,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be < " << val << ".";
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _le( val, oss.str() ) ) );
            }

            template<typename T>
            static validation_ptr_t<T> lt( T val,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be < " << val << ".";
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _lt( val, oss.str() ) ) );
            }

            template<typename T>
            static validation_ptr_t<T> ne( T val,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must be > " << val << ".";
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _ne( val, oss.str() ) ) );
            }

             template<typename T>
            static validation_ptr_t<T> not_in( const std::vector<T>& vals,
                                           const std::string& msg = "" ) {
                std::ostringstream oss;
                if (msg.empty()) {
                    oss << "Value must not be any of { "
                        << NCPA::strings::join( vals, ", " ) << " }."
                } else {
                    oss << msg;
                }
                return validation_ptr_t<T>(
                    new TypedValidation<T>( _not_in( vals, oss.str() ) ) );
            }


        }  // namespace tests
    }  // namespace config
}  // namespace NCPA
