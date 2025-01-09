#pragma once

#include "NCPA/parameters/declarations.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

DECLARE_SWAP_FUNCTION( NCPA::params::details::GreaterTest );
DECLARE_SWAP_FUNCTION( NCPA::params::details::LessTest );
DECLARE_SWAP_FUNCTION( NCPA::params::details::EqualTest );
DECLARE_SWAP_FUNCTION( NCPA::params::details::OneOfTest );
DECLARE_SWAP_FUNCTION( NCPA::params::Validator );
DECLARE_SWAP_FUNCTION( NCPA::params::RangeValidator );

namespace NCPA {
    namespace params {
        namespace details {
            template<typename TESTTYPE>
            class ValidationTest {
                public:
                    virtual ~ValidationTest() {}

                    virtual std::unique_ptr<ValidationTest<TESTTYPE>>
                        clone() = 0;
                    virtual bool test( TESTTYPE testval ) const = 0;
                    virtual std::string description() const     = 0;
            };

            template<typename TESTTYPE>
            class GreaterTest : public ValidationTest<TESTTYPE> {
                public:
                    GreaterTest( TESTTYPE minval, bool equalOK = false ) :
                        _min { minval }, _equalOK { equalOK } {}

                    GreaterTest( const GreaterTest<TESTTYPE>& other ) :
                        _min { other._min }, _equalOK { other._equalOK } {}

                    GreaterTest( GreaterTest<TESTTYPE>&& source ) noexcept {
                        ::swap( *this, source );
                    }

                    virtual ~GreaterTest() {}

                    friend void ::swap<TESTTYPE>(
                        GreaterTest<TESTTYPE>& a,
                        GreaterTest<TESTTYPE>& b ) noexcept;

                    GreaterTest<TESTTYPE>& operator=(
                        GreaterTest<TESTTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual std::unique_ptr<ValidationTest<TESTTYPE>>
                        clone() override {
                        return std::unique_ptr<ValidationTest<TESTTYPE>>(
                            new GreaterTest( *this ) );
                    }

                    virtual bool test( TESTTYPE testval ) const override {
                        if ( _equalOK ) {
                            return ( testval >= _min );
                        } else {
                            return ( testval > _min );
                        }
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "Value " << ( _equalOK ? ">= " : "> " ) << _min;
                        return oss.str();
                    }

                protected:
                    TESTTYPE _min;
                    bool _equalOK;
            };

            template<typename TESTTYPE>
            class LessTest : public ValidationTest<TESTTYPE> {
                public:
                    LessTest( TESTTYPE maxval, bool equalOK = false ) :
                        _max { maxval }, _equalOK { equalOK } {}

                    LessTest( const LessTest<TESTTYPE>& other ) :
                        _max { other._max }, _equalOK { other._equalOK } {}

                    LessTest( LessTest<TESTTYPE>&& source ) noexcept {
                        ::swap( *this, source );
                    }

                    virtual ~LessTest() {}

                    friend void ::swap<TESTTYPE>(
                        LessTest<TESTTYPE>& a,
                        LessTest<TESTTYPE>& b ) noexcept;

                    LessTest<TESTTYPE>& operator=( LessTest<TESTTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual bool test( TESTTYPE testval ) const override {
                        if ( _equalOK ) {
                            return ( testval <= _max );
                        } else {
                            return ( testval < _max );
                        }
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "Value " << ( _equalOK ? "<= " : "< " ) << _max;
                        return oss.str();
                    }

                protected:
                    TESTTYPE _max;
                    bool _equalOK;
            };

            template<typename TESTTYPE>
            class EqualTest : public ValidationTest<TESTTYPE> {
                public:
                    EqualTest( TESTTYPE compval, bool use_equals = true ) :
                        _val { compval }, _equals { use_equals } {}

                    EqualTest( const EqualTest<TESTTYPE>& other ) :
                        _val { other._val }, _equals { other._equals } {}

                    EqualTest( EqualTest<TESTTYPE>&& source ) noexcept {
                        ::swap( *this, source );
                    }

                    virtual ~EqualTest() {}

                    friend void ::swap<TESTTYPE>(
                        EqualTest<TESTTYPE>& a,
                        EqualTest<TESTTYPE>& b ) noexcept;

                    EqualTest<TESTTYPE>& operator=(
                        EqualTest<TESTTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual bool test( TESTTYPE testval ) const override {
                        return ( ( testval == _val ) == _equals );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "Value " << ( _equals ? "== " : "!= " ) << _val;
                        return oss.str();
                    }

                protected:
                    TESTTYPE _val;
                    bool _equals;
            };

            template<typename TESTTYPE>
            class OneOfTest : public ValidationTest<TESTTYPE> {
                public:
                    OneOfTest( std::vector<TESTTYPE> one_of,
                               bool return_if_in = true ) :
                        _testvec { one_of }, _return_if_in { return_if_in } {}

                    OneOfTest( const OneOfTest<TESTTYPE>& other ) :
                        _testvec { other._testvec },
                        _return_if_in { other._return_if_in } {}

                    OneOfTest( OneOfTest<TESTTYPE>&& source ) noexcept {
                        ::swap( *this, source );
                    }

                    virtual ~OneOfTest() {}

                    friend void ::swap<TESTTYPE>(
                        OneOfTest<TESTTYPE>& a,
                        OneOfTest<TESTTYPE>& b ) noexcept;

                    OneOfTest<TESTTYPE>& operator=(
                        OneOfTest<TESTTYPE> other ) {
                        ::swap( *this, other );
                        return *this;
                    }

                    virtual bool test( TESTTYPE testval ) const override {
                        for ( auto it = _testvec.cbegin();
                              it != _testvec.cend(); ++it ) {
                            if ( *it == testval ) {
                                return _return_if_in;
                            }
                        }
                        return !_return_if_in;
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "Value "
                            << ( _return_if_in ? "in { " : "not in { " );
                        for ( auto it = _testvec.cbegin();
                              it != _testvec.cend(); ++it ) {
                            if ( it != _testvec.cbegin() ) {
                                oss << ", ";
                            }
                            oss << *it;
                        }
                        oss << " }";
                        return oss.str();
                    }

                protected:
                    std::vector<TESTTYPE> _testvec;
                    bool _return_if_in;
            };

        }  // namespace details

        template<typename TESTTYPE>
        class Validator {
            public:
                Validator() {}

                Validator( const Validator<TESTTYPE>& other ) :
                    Validator(),
                    _tests { other._tests },
                    _errors { other._errors } {}

                Validator( Validator<TESTTYPE>&& source ) noexcept {
                    ::swap( *this, source );
                }

                virtual ~Validator() {}

                friend void ::swap<TESTTYPE>(
                    Validator<TESTTYPE>& a, Validator<TESTTYPE>& b ) noexcept;

                Validator<TESTTYPE>& operator=( Validator<TESTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual bool validate(
                    const Parameter<TESTTYPE>& param ) {
                    _errors.clear();
                    bool ok = true;
                    for ( auto it = _tests.cbegin(); it != _tests.cend();
                          ++it ) {
                        if ( !it->test( param.value() ) ) {
                            std::ostringstream oss;
                            oss << "Test failed for value " << param.value()
                                << ": " << it->description();
                            _errors.push_back( oss.str() );
                            ok = false;
                        }
                    }
                    return ok;
                }

                virtual Validator& add_test(
                    const details::ValidationTest<TESTTYPE>& test ) {
                    _tests.push_back( test );
                    return *this;
                }

            protected:
                std::vector<details::ValidationTest<TESTTYPE>> _tests;
                std::vector<std::string> _errors;
        };

        template<typename TESTTYPE>
        class RangeValidator : public Validator<TESTTYPE> {
            public:
                RangeValidator() : Validator<TESTTYPE>() {}

                RangeValidator( TESTTYPE minval, TESTTYPE maxval ) :
                    RangeValidator<TESTTYPE>() {
                    this->set_minimum( minval );
                    this->set_maximum( maxval );
                }

                RangeValidator( TESTTYPE minval ) :
                    RangeValidator<TESTTYPE>() {
                    this->set_minimum( minval );
                }

                RangeValidator( nullptr_t minval, TESTTYPE maxval ) :
                    RangeValidator<TESTTYPE>() {
                    this->set_maximum( maxval );
                }

                RangeValidator( TESTTYPE minval, bool mininclusive,
                                TESTTYPE maxval, bool maxinclusive ) :
                    RangeValidator<TESTTYPE>() {
                    this->set_minimum( minval, mininclusive );
                    this->set_maximum( maxval, maxinclusive );
                }

                RangeValidator( const RangeValidator<TESTTYPE>& other ) :
                    Validator<TESTTYPE>( other ) {}

                RangeValidator( RangeValidator<TESTTYPE>&& source ) noexcept {
                    ::swap( *this, source );
                }

                virtual ~RangeValidator() {}

                friend void ::swap<TESTTYPE>(
                    RangeValidator<TESTTYPE>& a,
                    RangeValidator<TESTTYPE>& b ) noexcept;

                RangeValidator<TESTTYPE>& operator=(
                    RangeValidator<TESTTYPE> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                RangeValidator<TESTTYPE>& set_minimum( TESTTYPE val,
                                                       bool inclusive
                                                       = true ) {
                    this->add_test( GreaterTest( val, inclusive ) );
                }

                RangeValidator<TESTTYPE>& set_maximum( TESTTYPE val,
                                                       bool inclusive
                                                       = true ) {
                    this->add_test( LessTest( val, inclusive ) );
                }
        };
    }  // namespace params
}  // namespace NCPA

template<typename T>
static void swap( NCPA::params::details::GreaterTest<T>& a,
                  NCPA::params::details::GreaterTest<T>& b ) noexcept {
    using std::swap;
    swap( a._min, b._min );
    swap( a._equalOK, b._equalOK );
}

template<typename T>
static void swap( NCPA::params::details::LessTest<T>& a,
                  NCPA::params::details::LessTest<T>& b ) noexcept {
    using std::swap;
    swap( a._max, b._max );
    swap( a._equalOK, b._equalOK );
}

template<typename T>
static void swap( NCPA::params::details::EqualTest<T>& a,
                  NCPA::params::details::EqualTest<T>& b ) noexcept {
    using std::swap;
    swap( a._val, b._val );
    swap( a._equals, b._equals );
}

template<typename T>
static void swap( NCPA::params::details::OneOfTest<T>& a,
                  NCPA::params::details::OneOfTest<T>& b ) noexcept {
    using std::swap;
    swap( a._testvec, b._testvec );
    swap( a._return_if_in, b._return_if_in );
}

template<typename T>
static void swap( NCPA::params::Validator<T>& a,
                  NCPA::params::Validator<T>& b ) noexcept {
    using std::swap;
    swap( a._tests, b._tests );
    swap( a._errors, b._errors );
}

template<typename T>
static void swap( NCPA::params::RangeValidator<T>& a,
                  NCPA::params::RangeValidator<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::params::Validator<T>&>( a ),
            static_cast<NCPA::params::Validator<T>&>( b ) );
}
