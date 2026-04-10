#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/VectorParameter.hpp"
#include "NCPA/units.hpp"

#include <type_traits>

namespace NCPA {
    namespace config {
        using namespace NCPA::units;

        template<typename PARAMTYPE>
        class VectorParameterWithUnits<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public VectorParameter<PARAMTYPE> {
            public:
                // default
                VectorParameterWithUnits() : VectorParameter<PARAMTYPE>() {}

                // scalar, no units
                VectorParameterWithUnits( PARAMTYPE defaultval ) :
                    VectorParameter<PARAMTYPE>(),
                    _uvalue { VectorWithUnits<PARAMTYPE>( 1, defaultval ) } {}

                // vector, no units
                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval ) :
                    VectorParameter<PARAMTYPE>(),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                // scalar, units separate
                VectorParameterWithUnits( PARAMTYPE defaultval,
                                          units_ptr_t u ) :
                    VectorParameter<PARAMTYPE>(),
                    _uvalue { VectorWithUnits<PARAMTYPE>(
                        std::vector<PARAMTYPE>( 1, defaultval ), u ) } {}

                // scalar with units
                VectorParameterWithUnits(
                    const ScalarWithUnits<PARAMTYPE>& defaultval ) :
                    VectorParameter<PARAMTYPE>(),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                // vector, units separate
                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u ) :
                    VectorParameter<PARAMTYPE>(),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval, u ) } {}

                // scalar with units
                VectorParameterWithUnits(
                    const VectorWithUnits<PARAMTYPE>& defaultval ) :
                    VectorParameter<PARAMTYPE>(), _uvalue { defaultval } {}

                VectorParameterWithUnits( const ValidationTest& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ) {}

                VectorParameterWithUnits( const test_ptr_t& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ) {}

                VectorParameterWithUnits(
                    std::initializer_list<test_ptr_t> new_tests ) :
                    VectorParameter<PARAMTYPE>( new_tests ) {}

                VectorParameterWithUnits( PARAMTYPE defaultval,
                                          const ValidationTest& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { std::vector<PARAMTYPE>( 1, defaultval ) } {}

                VectorParameterWithUnits( PARAMTYPE defaultval,
                                          const test_ptr_t& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { std::vector<PARAMTYPE>( 1, defaultval ) } {}

                VectorParameterWithUnits(
                    PARAMTYPE defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    VectorParameter<PARAMTYPE>( new_tests ),
                    _uvalue { std::vector<PARAMTYPE>( 1, defaultval ) } {}

                VectorParameterWithUnits( PARAMTYPE defaultval, units_ptr_t u,
                                          const ValidationTest& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { std::vector<PARAMTYPE>( 1, defaultval ), u } {}

                VectorParameterWithUnits( PARAMTYPE defaultval, units_ptr_t u,
                                          const test_ptr_t& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { std::vector<PARAMTYPE>( 1, defaultval ), u } {}

                VectorParameterWithUnits(
                    PARAMTYPE defaultval, units_ptr_t u,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    VectorParameter<PARAMTYPE>( new_tests ),
                    _uvalue { std::vector<PARAMTYPE>( 1, defaultval ), u } {}

                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    const ValidationTest& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    const test_ptr_t& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    VectorParameter<PARAMTYPE>( new_tests ),
                    _uvalue { defaultval, u } {}

                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const ValidationTest& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { defaultval, u } {}

                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const test_ptr_t& newtest ) :
                    VectorParameter<PARAMTYPE>( newtest ),
                    _uvalue { defaultval, u } {}

                VectorParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<test_ptr_t> new_tests ) :
                    VectorParameter<PARAMTYPE>( new_tests ),
                    _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const VectorParameterWithUnits<PARAMTYPE>& other ) :
                    VectorParameter<PARAMTYPE>( other ) {
                    _uvalue = other._uvalue;
                }

                VectorParameterWithUnits(
                    VectorParameterWithUnits<PARAMTYPE>&& other ) noexcept :
                    VectorParameter<PARAMTYPE>() {
                    ::swap( *this, other );
                }

                virtual ~VectorParameterWithUnits() {}

                VectorParameterWithUnits<PARAMTYPE>& operator=(
                    VectorParameterWithUnits<PARAMTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>(
                    VectorParameterWithUnits<PARAMTYPE>& a,
                    VectorParameterWithUnits<PARAMTYPE>& b ) noexcept;

                virtual param_ptr_t clone() const override {
                    return param_ptr_t(
                        new VectorParameterWithUnits<PARAMTYPE>( *this ) );
                }

                virtual std::string as_string() const override {
                    if (this->has_units()) {
                        return VectorParameter<PARAMTYPE>::as_string() + " "
                             + this->get_units()->to_string();
                    } else {
                        return VectorParameter<PARAMTYPE>::as_string();
                    }
                }

                virtual std::string as_string( size_t n ) const override {
                    if (this->has_units()) {
                        return VectorParameter<PARAMTYPE>::as_string( n ) + " "
                             + this->get_units()->to_string();
                    } else {
                        return VectorParameter<PARAMTYPE>::as_string( n );
                    }
                }

                virtual VectorParameterWithUnits<PARAMTYPE>& convert_units(
                    units_ptr_t u ) override {
                    _uvalue.convert_units( u );
                    return *this;
                }

                virtual PARAMTYPE get( size_t n = 0 ) const override {
                    return ( _uvalue.empty() ? this->_nullval()
                                             : _uvalue.at( n ) );
                }

                virtual units_ptr_t get_units() const override {
                    return _uvalue.get_units();
                }

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<T>::value, int>::type = 0>
                ScalarWithUnits<PARAMTYPE> get_with_units( size_t n
                                                           = 0 ) const {
                    return _uvalue.get_scalar( n );
                }

                virtual std::vector<PARAMTYPE> get_vector() const override {
                    return _uvalue.std();
                }

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<T>::value, int>::type = 0>
                VectorWithUnits<PARAMTYPE> get_vector_with_units() const {
                    return _uvalue;
                }

                virtual bool has_units() const {
                    return ( _uvalue.get_units() != nullptr );
                }

                virtual void set_units( units_ptr_t u ) {
                    _uvalue.set_units( u );
                }

                virtual size_t size() const override { return _uvalue.size(); }

            protected:
                VectorWithUnits<PARAMTYPE> _uvalue;
        };
    }  // namespace config
}  // namespace NCPA

template<typename PARAMTYPE>
void swap( NCPA::config::VectorParameterWithUnits<PARAMTYPE>& a,
           NCPA::config::VectorParameterWithUnits<PARAMTYPE>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::VectorParameter<PARAMTYPE>&>( a ),
            static_cast<NCPA::config::VectorParameter<PARAMTYPE>&>( b ) );
    swap( a._uvalue, b._uvalue );
}
