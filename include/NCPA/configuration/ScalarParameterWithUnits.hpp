#pragma once

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ScalarParameter.hpp"
#include "NCPA/units.hpp"

#include <type_traits>

namespace NCPA {
    namespace config {
        using namespace NCPA::units;

        template<typename PARAMTYPE>
        class ScalarParameterWithUnits<
            PARAMTYPE, typename std::enable_if<
                           std::is_floating_point<PARAMTYPE>::value>::type>
            :  // public ScalarParameter<PARAMTYPE>,
               public Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                                ScalarParameter<PARAMTYPE>> {
            public:
                // using ScalarParameter<PARAMTYPE>::init;

                // virtual void init( const ScalarWithUnits<PARAMTYPE>& u ) {
                //     _uvalue = u;
                // }

                // virtual void init( const ScalarWithUnits<PARAMTYPE> *u ) {
                //     init( *u );
                // }

                // virtual void init( const std::vector<PARAMTYPE>& defaultval
                // ) {
                //     init(( defaultval.empty() ? 0 : defaultval.at( 0 ) ))
                // }

                // virtual void init( const std::vector<PARAMTYPE>& defaultval,
                //                    units_ptr_t u ) {
                //     init( ScalarWithUnits<PARAMTYPE>(
                //         ( defaultval.empty() ? 0 : defaultval.at( 0 ) ), u )
                //         );
                // }

                // virtual void init(
                //     const VectorWithUnits<PARAMTYPE>& defaultval ) {
                //     init( defaultval.empty()
                //               ? ScalarWithUnits<PARAMTYPE>(
                //                     0.0, defaultval.get_units() )
                //               : defaultval.get_scalar( 0 ) );
                // }

                ScalarParameterWithUnits() {}

                ScalarParameterWithUnits( PARAMTYPE defaultval ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval } {}

                ScalarParameterWithUnits( PARAMTYPE defaultval,
                                          units_ptr_t u ) :
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const ScalarWithUnits<PARAMTYPE>& defaultval ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval } {}

                ScalarParameterWithUnits(
                    const VectorWithUnits<PARAMTYPE>& defaultval ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval } {}

                ScalarParameterWithUnits(
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { v } {}

                ScalarParameterWithUnits(
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { v } {}

                ScalarParameterWithUnits(
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& ptr ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { ptr } {}

                ScalarParameterWithUnits(
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { tests } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval, v } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval, v } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval, v } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> { defaultval, v } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> {
                        ScalarWithUnits<PARAMTYPE>( defaultval, u ), v
                    } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>> {
                        ScalarWithUnits<PARAMTYPE>( defaultval, u ), v
                    } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval, units_ptr_t u,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    PARAMTYPE defaultval, units_ptr_t u,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( tests ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    ScalarWithUnits<PARAMTYPE> defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( tests ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( tests ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const std::vector<PARAMTYPE>& defaultval, units_ptr_t u,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( tests ),
                    _uvalue { ScalarWithUnits<PARAMTYPE>( defaultval, u ) } {}

                ScalarParameterWithUnits(
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    const TypedValidation<PARAMTYPE> *v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    const std::unique_ptr<TypedValidation<PARAMTYPE>>& v ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( v ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const VectorWithUnits<PARAMTYPE> defaultval,
                    std::initializer_list<TypedValidation<PARAMTYPE>> tests ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( tests ),
                    _uvalue { defaultval } {}

                ScalarParameterWithUnits(
                    const ScalarParameterWithUnits<PARAMTYPE>& other ) :
                    Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                              ScalarParameter<PARAMTYPE>>( other ) {
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
                    swap( static_cast<
                              Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                                        ScalarParameter<PARAMTYPE>>&>( a ),
                          static_cast<
                              Cloneable<ScalarParameterWithUnits<PARAMTYPE>,
                                        ScalarParameter<PARAMTYPE>>&>( b ) );
                    swap( a._uvalue, b._uvalue );
                }

                virtual std::string as_string( size_t n = 0 ) const override {
                    if (this->has_units()) {
                        return ScalarParameter<PARAMTYPE>::as_string( n ) + " "
                             + this->get_units()->to_string();
                    } else {
                        return ScalarParameter<PARAMTYPE>::as_string( n );
                    }
                }

                // NCPA_CLONE_METHOD( ScalarParameterWithUnits<PARAMTYPE>,
                //                    BaseParameter )

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

// template<typename PARAMTYPE>
// void swap( NCPA::config::ScalarParameterWithUnits<PARAMTYPE>& a,
//            NCPA::config::ScalarParameterWithUnits<PARAMTYPE>& b ) noexcept {
//     using std::swap;
//     swap( static_cast<NCPA::config::ScalarParameter<PARAMTYPE>&>( a ),
//             static_cast<NCPA::config::ScalarParameter<PARAMTYPE>&>( b ) );
//     swap( a._uvalue, b._uvalue );
// }
