#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/parameters/BaseParameter.hpp"
#include "NCPA/configuration/parameters/VectorParameter.hpp"
#include "NCPA/units.hpp"

#include <type_traits>

namespace NCPA {
    namespace config {
        using namespace NCPA::units;

        template<typename PARAMTYPE>
        class VectorParameterWithUnits<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public Cloneable<VectorParameterWithUnits<PARAMTYPE>,
                               VectorParameter<PARAMTYPE>> {
                using cloneable_t
                    = Cloneable<VectorParameterWithUnits<PARAMTYPE>,
                                VectorParameter<PARAMTYPE>>;

            public:
                // default
                VectorParameterWithUnits( const std::string& tag = "" ) :
                    cloneable_t( tag ) {}

                // scalar, no units
                VectorParameterWithUnits( const std::string& tag,
                                          PARAMTYPE defaultval ) :
                    cloneable_t( tag ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( 1, defaultval ) } {}

                // vector, no units
                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval ) :
                    cloneable_t( tag ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const TypedValidation<PARAMTYPE>& v ) :
                    cloneable_t( tag, v ) {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const TypedValidation<PARAMTYPE> *v ) :
                    cloneable_t( tag, v ) {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    cloneable_t( tag, v ) {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    cloneable_t( tag, v ) {}

                // scalar, units separate
                VectorParameterWithUnits( const std::string& tag,
                                          PARAMTYPE defaultval,
                                          units_ptr_t u ) :
                    cloneable_t( tag ),
                    _uvalue { VectorWithUnits<PARAMTYPE>(
                        std::vector<PARAMTYPE>( 1, defaultval ), u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u, const TypedValidation<PARAMTYPE>& v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>(
                        std::vector<PARAMTYPE>( 1, defaultval ), u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u, const TypedValidation<PARAMTYPE> *v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>(
                        std::vector<PARAMTYPE>( 1, defaultval ), u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>(
                        std::vector<PARAMTYPE>( 1, defaultval ), u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>(
                        std::vector<PARAMTYPE>( 1, defaultval ), u ) } {}

                // scalar with units
                VectorParameterWithUnits(
                    const std::string& tag,
                    const ScalarWithUnits<PARAMTYPE>& defaultval ) :
                    cloneable_t( tag ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const ScalarWithUnits<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const ScalarWithUnits<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const ScalarWithUnits<PARAMTYPE>& defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const ScalarWithUnits<PARAMTYPE>& defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval ) } {}

                // vector, units separate
                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u ) :
                    cloneable_t( tag ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval, u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE>& v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval, u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE> *v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval, u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval, u ) } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    cloneable_t( tag, v ),
                    _uvalue { VectorWithUnits<PARAMTYPE>( defaultval, u ) } {}

                // vector with units
                VectorParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE>& defaultval ) :
                    cloneable_t( tag ), _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    cloneable_t( tag, v ), _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    cloneable_t( tag, v ), _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE>& defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    cloneable_t( tag, v ), _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE>& defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    cloneable_t( tag, v ), _uvalue { defaultval } {}

                VectorParameterWithUnits(
                    const VectorParameterWithUnits<PARAMTYPE>& other ) :
                    cloneable_t( other ), _uvalue( other._uvalue ) {}

                VectorParameterWithUnits(
                    VectorParameterWithUnits<PARAMTYPE>&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~VectorParameterWithUnits() = default;

                VectorParameterWithUnits<PARAMTYPE>& operator=(
                    VectorParameterWithUnits<PARAMTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap(
                    VectorParameterWithUnits<PARAMTYPE>& a,
                    VectorParameterWithUnits<PARAMTYPE>& b ) noexcept {
                    using std::swap;
                    swap(
                        static_cast<NCPA::config::VectorParameter<PARAMTYPE>&>(
                            a ),
                        static_cast<NCPA::config::VectorParameter<PARAMTYPE>&>(
                            b ) );
                    swap( a._uvalue, b._uvalue );
                }

                virtual std::string as_string( bool throw_on_error
                                               = false ) const override {
                    if (this->has_units()) {
                        return VectorParameter<PARAMTYPE>::as_string() + " "
                             + this->get_units()->to_string();
                    } else {
                        return VectorParameter<PARAMTYPE>::as_string();
                    }
                }

                virtual std::string as_string( size_t n = 0,
                                               bool throw_on_error
                                               = false ) const override {
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

// template<typename PARAMTYPE>
// void swap( NCPA::config::VectorParameterWithUnits<PARAMTYPE>& a,
//            NCPA::config::VectorParameterWithUnits<PARAMTYPE>& b ) noexcept {
//     using std::swap;
//     swap( static_cast<NCPA::config::VectorParameter<PARAMTYPE>&>( a ),
//             static_cast<NCPA::config::VectorParameter<PARAMTYPE>&>( b ) );
//     swap( a._uvalue, b._uvalue );
// }
