#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/types.hpp"
#include "NCPA/units.hpp"

#include <limits>
#include <type_traits>

namespace NCPA {
    namespace config {

        template<typename PARAMTYPE>
        class TypedParameter : public BaseParameter {
            public:
                using value_type = PARAMTYPE;

                TypedParameter() : BaseParameter() {}

                TypedParameter( const ValidationTest& newtest ) :
                    BaseParameter( newtest ) {}

                TypedParameter( const test_ptr_t& newtest ) :
                    BaseParameter( newtest ) {}

                TypedParameter( std::initializer_list<test_ptr_t> new_tests ) :
                    BaseParameter( new_tests ) {}

                virtual ~TypedParameter() {}

                TypedParameter( const TypedParameter<PARAMTYPE>& other ) :
                    BaseParameter( *this ) {}

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

                // virtual PARAMTYPE& value( size_t n = 0 )             = 0;
                // virtual const PARAMTYPE& value( size_t n = 0 ) const = 0;
                virtual PARAMTYPE get( size_t n = 0 ) const = 0;

                // virtual ScalarWithUnits<PARAMTYPE> get_with_units(
                //     size_t n = 0 ) const {
                //     throw std::logic_error( "Parameter has no units defined" );
                // }

                virtual std::vector<PARAMTYPE> get_vector() const = 0;

                // virtual VectorWithUnits<PARAMTYPE> get_vector_with_units()
                //     const {
                //     throw std::logic_error( "Parameter has no units defined" );
                // }

                template<typename ASTYPE, typename std::enable_if<
                    std::is_convertible<PARAMTYPE,ASTYPE>::value, int
                >::type = 0>
                ASTYPE as() const {
                    return static_cast<ASTYPE>(this->get());
                }

                template<typename ASTYPE, typename std::enable_if<
                    std::is_convertible<PARAMTYPE,ASTYPE>::value, int
                >::type = 0>
                ScalarWithUnits<ASTYPE> as_with_units() const {
                    return ScalarWithUnits<ASTYPE>( static_cast<ASTYPE>(this->get()), this->get_units() );
                }

                template<typename ASTYPE, typename std::enable_if<
                    (!(std::is_convertible<PARAMTYPE,ASTYPE>::value)), int
                >::type = 0>
                ScalarWithUnits<ASTYPE> as_with_units() const {
                    return ScalarWithUnits<ASTYPE>( static_cast<ASTYPE>(this->as_double()), this->get_units() );
                }

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
    ::swap( static_cast<NCPA::config::BaseParameter&>( a ),
            static_cast<NCPA::config::BaseParameter&>( b ) );
}
