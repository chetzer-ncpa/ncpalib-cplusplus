#pragma once

#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/parameters/BaseParameter.hpp"
#include "NCPA/configuration/parameters/ScalarParameter.hpp"
#include "NCPA/units.hpp"

#include <type_traits>

namespace NCPA {
    namespace config {
        using namespace NCPA::units;

        template<typename PARAMTYPE>
        class ScalarParameterWithUnits<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            : public Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                               ScalarParameter<PARAMTYPE>> {
            public:
                using parent_t = Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                                           ScalarParameter<PARAMTYPE>>;
                using ScalarParameter<PARAMTYPE>::as_string;

                ScalarParameterWithUnits( const std::string& tag = "" ) :
                    parent_t( tag ) {}

                ScalarParameterWithUnits( const std::string& tag,
                                          PARAMTYPE defaultval ) :
                    parent_t( tag ), _uvalue { defaultval } {}

                ScalarParameterWithUnits( const std::string& tag,
                                          PARAMTYPE defaultval,
                                          units_ptr_t u ) :
                    parent_t( tag ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const ScalarWithUnits<PARAMTYPE>& defaultval ) :
                    ScalarParameterWithUnits<PARAMTYPE>(
                        tag, defaultval.get(), defaultval.get_units() ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval ) :
                    ScalarParameterWithUnits(
                        tag, defaultval.empty() ? 0 : defaultval[ 0 ] ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u ) :
                    ScalarParameterWithUnits(
                        tag, defaultval.empty() ? 0 : defaultval[ 0 ], u ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE>& defaultval ) :
                    ScalarParameterWithUnits(
                        tag, defaultval.empty()
                                 ? ScalarWithUnits<PARAMTYPE>(
                                       0, defaultval.get_units() )
                                 : defaultval.get_scalar( 0 ) ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& ptr ) :
                    parent_t( tag, ptr ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    parent_t( tag, tests ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, defaultval, v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, defaultval, v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    parent_t( tag, defaultval, v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    parent_t( tag, defaultval, v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u, const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, ScalarWithUnits<PARAMTYPE>( defaultval, u ),
                              v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u, const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, ScalarWithUnits<PARAMTYPE>( defaultval, u ),
                              v ) {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    parent_t( tag, v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag, PARAMTYPE defaultval,
                    units_ptr_t u,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    parent_t( tag, tests ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    parent_t( tag, tests ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    parent_t( tag, tests ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    parent_t( tag, v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    parent_t( tag, tests ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    parent_t( tag, v ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const std::string& tag,
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    parent_t( tag, tests ), _uvalue( defaultval ) {}

                ScalarParameterWithUnits(
                    const ScalarParameterWithUnits<PARAMTYPE>& other ) :
                    parent_t( other ) {
                    _uvalue = other._uvalue;
                }

                ScalarParameterWithUnits(
                    ScalarParameterWithUnits<PARAMTYPE>&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~ScalarParameterWithUnits() {}

                ScalarParameterWithUnits<PARAMTYPE>& operator=(
                    ScalarParameterWithUnits<PARAMTYPE> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap(
                    ScalarParameterWithUnits<PARAMTYPE>& a,
                    ScalarParameterWithUnits<PARAMTYPE>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<parent_t&>( a ),
                          static_cast<parent_t&>( b ) );
                    swap( a._uvalue, b._uvalue );
                }

                virtual std::string as_string( size_t n = 0,
                                               bool throw_on_error
                                               = false ) const override {
                    if (this->has_units()) {
                        return ScalarParameter<PARAMTYPE>::as_string( n ) + " "
                             + this->get_units()->to_string();
                    } else {
                        return ScalarParameter<PARAMTYPE>::as_string( n );
                    }
                }

                virtual bool has_units() const {
                    return ( _uvalue.get_units() != nullptr );
                }

                virtual void set_units( units_ptr_t u ) {
                    _uvalue.set_units( u );
                }

                virtual units_ptr_t get_units() const override {
                    return _uvalue.get_units();
                }

                virtual ScalarParameterWithUnits<PARAMTYPE>& convert_units(
                    units_ptr_t u ) override {
                    _uvalue.convert_units( u );
                    return *this;
                }

                virtual PARAMTYPE get( size_t n = 0 ) const override {
                    return _uvalue.get();
                }

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<T>::value, int>::type = 0>
                ScalarWithUnits<PARAMTYPE> get_with_units( size_t n
                                                           = 0 ) const {
                    return _uvalue;
                }

                virtual std::vector<PARAMTYPE> get_vector() const override {
                    return std::vector<PARAMTYPE> { _uvalue.get() };
                }

                template<typename T = PARAMTYPE,
                         typename std::enable_if<
                             std::is_floating_point<T>::value, int>::type = 0>
                VectorWithUnits<PARAMTYPE> get_vector_with_units() const {
                    return VectorWithUnits<PARAMTYPE>( _uvalue );
                }

            protected:
                ScalarWithUnits<PARAMTYPE> _uvalue;

                virtual void _check_has_units() const override {
                    if (_uvalue.get_units() == nullptr) {
                        throw std::logic_error(
                            "Parameter does not have units!" );
                    }
                }
        };
    }  // namespace config
}  // namespace NCPA
