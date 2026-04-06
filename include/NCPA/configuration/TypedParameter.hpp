#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/Parameter.hpp"
#include "NCPA/types.hpp"
#include "NCPA/units.hpp"

#include <limits>
#include <type_traits>

namespace NCPA {
    namespace config {

        template<typename PARAMTYPE>
        class TypedParameter : public Parameter {
            public:
                using value_type = PARAMTYPE;

                TypedParameter() : Parameter() {}

                TypedParameter( const ValidationTest& newtest ) :
                    Parameter( newtest ) {}

                TypedParameter( const test_ptr_t& newtest ) :
                    Parameter( newtest ) {}

                TypedParameter( std::initializer_list<test_ptr_t> new_tests ) :
                    Parameter( new_tests ) {}

                virtual ~TypedParameter() {}

                TypedParameter( const TypedParameter<PARAMTYPE>& other ) :
                    Parameter( *this ) {}

                TypedParameter( TypedParameter<PARAMTYPE>&& other ) noexcept :
                    TypedParameter() {
                    ::swap( *this, other );
                }

                TypedParameter<PARAMTYPE>& operator=(
                    TypedParameter<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( TypedParameter<PARAMTYPE>& a,
                                      TypedParameter<PARAMTYPE>& b ) noexcept;

                virtual parameter_type_t type() const override {
                    return parameter_type<PARAMTYPE>();
                }

                virtual PARAMTYPE& value( size_t n = 0 ) {
                    return this->get( n );
                }

                virtual const PARAMTYPE& value( size_t n = 0 ) const {
                    return this->get( n );
                }

                virtual PARAMTYPE& get( size_t n = 0 )             = 0;
                virtual const PARAMTYPE& get( size_t n = 0 ) const = 0;
                virtual std::vector<PARAMTYPE> get_vector() const  = 0;

            protected:
                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             NCPA::types::has_to_string<T>::value, int>::type
                         = 0>
                std::string _as_string( size_t n = 0 ) const {
                    return to_string( this->get( n ) );
                }

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             !( NCPA::types::has_to_string<T>::value ),
                             int>::type = 0>
                std::string _as_string( size_t n = 0 ) const {
                    throw std::out_of_range(
                        "No as_string() conversion defined!" );
                }
        };
    }  // namespace config
}  // namespace NCPA

template<typename T>
void swap( NCPA::config::TypedParameter<T>& a,
           NCPA::config::TypedParameter<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::Parameter&>( a ),
            static_cast<NCPA::config::Parameter&>( b ) );
}
