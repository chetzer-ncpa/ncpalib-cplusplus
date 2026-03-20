/**
 * configuration.hpp: Configuration and validation library
 *
 * The aim of this library is to generalize configuration parameters,
 * so that the use of type-specific
 *
 *
 * Example: declare a floating-point parameter with a validation
 * test that the value must be nonzero and less than 5.0:
 *      DoubleParameter d2( { IsNotZero<double>() } );
 *
 *
 * EXTENDING THE LIBRARY
 * To extend the validation test library, follow these steps:
 *
 *  1. Where it says DECLARE_NEW_TESTS_HERE, declare the new template
 * or class.  Match the existing syntax.  If you will need to declare multiple
 * template versions for different types (i.e. the program logic depends on the
 * type of the parameter), use the
 *          template<typename T, typename Enable = void>
 * syntax, otherwise use the
 *          template<typename T>
 * syntax.  If the test does not depend on type at all, declare the class with
 * no template.
 *
 *  2. Where it says DECLARE_SWAP_FUNCTIONS_HERE, declare a swap function for
 * the new template or class.  Match the existing syntax.
 *
 *  3. Where it says DEFINE_NEW_TESTS_HERE, add the definition(s) for the
 * template(s) declared in step 1.  The general syntax for a test whose logic
 * does not depend on the type of the parameter looks like this:

template<typename T>
class TESTNAME : public BASETESTNAME {
    public:
        TESTNAME( T value ) : BASETESTNAME( value ) {}

        // copy constructor, optional unless you've added member variables
        TESTNAME( const TESTNAME<T>& other ) {
            value1 = other.value1;
            value2 = other.value2;
            // etc
        }

        virtual ~TESTNAME() {}

        friend void ::swap<>( TESTNAME<T>& a,
                              TESTNAME<T>& b ) noexcept;

        virtual std::unique_ptr<ValidationTest> clone() const override {
            return std::unique_ptr<ValidationTest>( new TESTNAME<T>( *this ) );
        }

        virtual std::string description() const override {
            // output descriptive string
        }

        virtual ValidationTest& test(
            const _configuration_parameter *param ) override {

            T param_value = this->parameter_value( param );

            // test logic goes here.  A passed test should be indicated
            // by calling
            //      this->pass()
            // and a failed test should be indicated by calling
            //      this->fail()

            return *this;
        }
};

 * If the logic does depend on the parameter type, define multiple templates
 * with the class declarations changed to:

template<typename T>
class TESTNAME< T,
                typename std::enable_if<std::TYPETRAIT<T>::value>::type
              > : public BASETESTNAME {

 * In both above cases, TESTNAME is the name of the new test class and should
 * end in "Test".  BASETESTNAME should be one of the following:
 *      NullaryValidationTest: If the test does not involve any values (e.g.
 *                             was the parameter set or not)
 *      UnaryValidationTest<T>: If the test involves one value (e.g. is the
 *                              parameter greater than 1)
 *      BinaryValidationTest<T>: If the test involves two values (e.g. is the
 *                               parameter between -1 and 1)
 *
 *  4. (optional but recommended) Where it says
 * DEFINE_NEW_CONVENIENCE_FUNCTIONS_HERE, define a function following this
 * example:

template<typename T>
test_ptr_t IsEqualTo( const T& val ) {
    return test_ptr_t( new IsEqualToTest<T>( val ) );
}

 * If you choose not to do this, you will have to explicitly create the
 * unique_ptr instances to pass into the parameter validation methods.
 *
 *  5. Finally, where it says DEFINE_SWAP_FUNCTIONS_HERE, define the swap
 * function you declared in step 2.  Generally this will be of the form:

template<typename T>
void swap( NCPA::config::validation::TESTNAME<T>& a,
           NCPA::config::validation::TESTNAME<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::BASETESTNAME<T>&>( a ),
        static_cast<NCPA::config::validation::BASETESTNAME<T>&>( b ) );

    // if you added any new members, swap them here:
    swap( a.member1, b.member1 );
    swap( a.member2, b.member2 );
    // etc
}
 */


#pragma once

#include "NCPA/defines.hpp"
#include "NCPA/types.hpp"

#include <array>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace NCPA {
    namespace config {
        enum class test_status_t { NONE, PENDING, FAILED, PASSED };
        enum class parameter_form_t { UNDEF, SCALAR, VECTOR };
        enum class parameter_type_t {
            UNDEF,
            INTEGER,
            FLOAT,
            STRING,
            BOOLEAN,
            ENUM,
            COMPLEX,
            OTHER
        };
        template<typename KEYTYPE>
        class Configurable;
        // class ConfigurationParameter;
        // template<typename T>
        // class TypedParameter;
    }  // namespace config
}  // namespace NCPA

// class _configuration_parameter;
DECLARE_CLASS_AND_SWAP_2NAMESPACE( _configuration_parameter, NCPA, config )
// DECLARE_CLASS_AND_SWAP_2NAMESPACE( BooleanParameter, NCPA, config )
// DECLARE_CLASS_AND_SWAP_2NAMESPACE( DoubleParameter, NCPA, config )
// DECLARE_CLASS_AND_SWAP_2NAMESPACE( IntegerParameter, NCPA, config )
// DECLARE_CLASS_AND_SWAP_2NAMESPACE( StringParameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( ValidationTest, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( ValidationTestSuite, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( Parameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_3NAMESPACE( NullaryValidationTest, NCPA, config,
                                   validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( TypedValidationTest, NCPA, config,
                                      validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( UnaryValidationTest, NCPA, config,
                                      validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( BinaryValidationTest, NCPA, config,
                                      validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( ListValidationTest, NCPA, config,
                                      validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsEqualToTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsNotEqualToTest, NCPA, config,
                                      validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsNotOneOfTest, NCPA, config,
                                      validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsOneOfTest, NCPA, config, validation )
DECLARE_CLASS_AND_SWAP_3NAMESPACE( WasSetTest, NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsBetweenTest, NCPA, config,
                                                  validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsGreaterThanTest, NCPA,
                                                  config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsGreaterThanOrEqualToTest,
                                                  NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsLessThanTest, NCPA, config,
                                                  validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsLessThanOrEqualToTest,
                                                  NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( ConfigurationMap, NCPA, config )

template<typename DERIVEDTYPE, typename KEYTYPE>
void swap( NCPA::config::Configurable<KEYTYPE>& a,
           NCPA::config::Configurable<KEYTYPE>& b ) noexcept;

namespace NCPA {
    namespace config {
        typedef std::unique_ptr<ValidationTest> test_ptr_t;
    }  // namespace config
}  // namespace NCPA

namespace NCPA {
    namespace config {
        class ValidationTest {
            public:
                ValidationTest() : _status { test_status_t::PENDING } {}

                ValidationTest( const ValidationTest& other ) :
                    ValidationTest() {
                    _status = other._status;
                }

                ValidationTest( ValidationTest&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~ValidationTest() {}

                friend void ::swap( ValidationTest& a,
                                    ValidationTest& b ) noexcept;

                virtual ValidationTest& clear() {
                    _status = test_status_t::PENDING;
                    return *this;
                }

                virtual ValidationTest& fail() {
                    _status = test_status_t::FAILED;
                    return *this;
                }

                virtual bool failed() const {
                    return ( _status == test_status_t::FAILED );
                }

                virtual ValidationTest& pass() {
                    _status = test_status_t::PASSED;
                    return *this;
                }

                virtual bool passed() const {
                    return ( _status == test_status_t::PASSED );
                }

                virtual bool pending() const {
                    return ( _status == test_status_t::PENDING );
                }

                virtual test_status_t status() const { return _status; }

                virtual std::string status_message() const {
                    std::ostringstream oss;
                    oss << this->description() << ": ";
                    switch (this->status()) {
                        case test_status_t::PENDING:
                            oss << "PENDING";
                            break;
                        case test_status_t::PASSED:
                            oss << "PASSED";
                            break;
                        case test_status_t::FAILED:
                            oss << "FAILED";
                    }
                    return oss.str();
                }

                virtual ValidationTest& test(
                    const std::unique_ptr<_configuration_parameter>& param ) {
                    return this->test( param.get() );
                }

                virtual std::string description() const = 0;
                virtual ValidationTest& test(
                    const _configuration_parameter *param )           = 0;
                virtual std::unique_ptr<ValidationTest> clone() const = 0;

            protected:
                test_status_t _status;
        };

        class ValidationTestSuite {
            public:
                ValidationTestSuite() {}

                ValidationTestSuite( const ValidationTestSuite& other ) :
                    ValidationTestSuite() {
                    for (auto it = other._tests.begin();
                         it != other._tests.end(); ++it) {
                        this->append( **it );
                    }
                }

                ValidationTestSuite( ValidationTestSuite&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~ValidationTestSuite() {}

                ValidationTestSuite& operator=( ValidationTestSuite other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( ValidationTestSuite& a,
                                    ValidationTestSuite& b ) noexcept;

                virtual ValidationTestSuite& append(
                    const ValidationTest& newtest ) {
                    _tests.push_back( newtest.clone() );
                    return *this;
                }

                virtual ValidationTestSuite& append(
                    const ValidationTest *newtest ) {
                    _tests.push_back( newtest->clone() );
                    return *this;
                }

                virtual ValidationTestSuite& append(
                    const test_ptr_t& newtest ) {
                    _tests.push_back( newtest->clone() );
                    return *this;
                }

                virtual bool failed() const {
                    return ( this->status() == test_status_t::FAILED );
                }

                virtual std::vector<const ValidationTest *> failed_tests()
                    const {
                    std::vector<const ValidationTest *> failed_tests;
                    for (auto it = _tests.cbegin(); it != _tests.cend();
                         ++it) {
                        if (( *it )->failed()) {
                            failed_tests.push_back( it->get() );
                        }
                    }
                    return failed_tests;
                }

                virtual bool passed() const {
                    return ( this->status() == test_status_t::PASSED );
                }

                virtual bool pending() const {
                    return ( this->status() == test_status_t::PENDING );
                }

                virtual ValidationTestSuite& prepend(
                    const ValidationTest& newtest ) {
                    _tests.emplace( _tests.begin(), newtest.clone() );
                    return *this;
                }

                virtual ValidationTestSuite& prepend(
                    const test_ptr_t& newtest ) {
                    _tests.insert( _tests.begin(), newtest->clone() );
                    return *this;
                }

                virtual ValidationTestSuite& run_tests(
                    const _configuration_parameter *param,
                    bool short_circuit = false ) {
                    // _status = test_status_t::PENDING;
                    bool pass = true;
                    for (auto it = _tests.begin(); it != _tests.end(); ++it) {
                        pass = ( ( *it )->test( param ).passed() ) && pass;
                        if (short_circuit && !pass) {
                            return *this;
                        }
                    }
                    return *this;
                }

                virtual size_t size() const { return _tests.size(); }

                // if any tests have failed, return FAILED.  Otherwise, if any
                // tests are pending, return PENDING.  Otherwise return PASSED.
                virtual test_status_t status() const {
                    if (_tests.empty()) {
                        return test_status_t::PASSED;
                    }
                    bool pending = false;
                    for (auto it = _tests.cbegin(); it != _tests.cend();
                         ++it) {
                        switch (( *it )->status()) {
                            case test_status_t::FAILED:
                                return test_status_t::FAILED;
                                break;
                            case test_status_t::PENDING:
                                pending = true;
                                break;
                            case test_status_t::PASSED:
                                break;
                            default:
                                throw std::range_error(
                                    "ValidationTestSuite.status(): "
                                    "Unrecognized status type!" );
                        }
                    }
                    if (pending) {
                        return test_status_t::PENDING;
                    } else {
                        return test_status_t::PASSED;
                    }
                }

                virtual std::vector<std::string> status_messages(
                    bool failed_only = false ) {
                    std::vector<std::string> messages;
                    for (auto it = _tests.cbegin(); it != _tests.cend();
                         ++it) {
                        if (( !failed_only ) || ( *it )->failed()) {
                            messages.push_back( ( *it )->status_message() );
                        }
                    }
                    return messages;
                }

            private:
                // test_status_t _status;
                std::vector<std::unique_ptr<ValidationTest>> _tests;
        };

        class _configuration_parameter {
            public:
                _configuration_parameter() {}

                _configuration_parameter(
                    const _configuration_parameter& other ) {
                    _tests = other._tests;
                }

                virtual ~_configuration_parameter() {}

                friend void ::swap( _configuration_parameter& a,
                                    _configuration_parameter& b ) noexcept;

                virtual _configuration_parameter& append_test(
                    const ValidationTest& newtest ) {
                    _tests.append( newtest );
                    return *this;
                }

                virtual _configuration_parameter& append_test(
                    const test_ptr_t& newtest ) {
                    _tests.append( newtest );
                    return *this;
                }

                virtual _configuration_parameter& append_tests(
                    std::initializer_list<test_ptr_t> new_tests ) {
                    for (auto it = new_tests.begin(); it != new_tests.end();
                         ++it) {
                        this->append_test( *it );
                    }
                    return *this;
                }

                virtual const ValidationTestSuite& tests() const {
                    return _tests;
                }

                virtual std::vector<const ValidationTest *> failed_tests()
                    const {
                    return _tests.failed_tests();
                }

                virtual _configuration_parameter& prepend_test(
                    const ValidationTest& newtest ) {
                    _tests.prepend( newtest );
                    return *this;
                }

                virtual std::ostream& validation_report(
                    std::ostream& os, bool newline = true,
                    const std::string& prepend = "" ) const {
                    if (this->passed()) {
                        os << prepend << "All tests passed";
                    } else {
                        auto f = this->failed_tests();
                        for (auto it = f.begin(); it != f.end(); ++it) {
                            if (it != f.begin()) {
                                os << std::endl;
                            }
                            os << prepend
                               << "Failed: " << ( *it )->description();
                        }
                    }
                    if (newline) {
                        os << std::endl;
                    }
                    return os;
                }

                virtual _configuration_parameter& prepend_tests(
                    std::initializer_list<ValidationTest> new_tests ) {
                    for (auto it = new_tests.begin(); it != new_tests.end();
                         ++it) {
                        this->prepend_test( *it );
                    }
                    return *this;
                }

                virtual _configuration_parameter& validate( bool short_circuit
                                                            = false ) {
                    _tests.run_tests( this, short_circuit );
                    return *this;
                }

                virtual bool failed() const { return _tests.failed(); }

                virtual bool passed() const { return _tests.passed(); }

                virtual bool pending() const { return _tests.pending(); }

                virtual test_status_t validation_status() const {
                    return _tests.status();
                }

                virtual bool as_bool() const {
                    throw std::logic_error( "as_bool(): not implemented" );
                }

                virtual std::string as_string() const {
                    throw std::logic_error( "as_string(): not implemented" );
                }

                virtual int as_int() const {
                    throw std::logic_error( "as_int(): not implemented" );
                }

                virtual double as_double() const {
                    throw std::logic_error( "as_double(): not implemented" );
                }

                template<typename T>
                const T& get() const {
                    if (auto ptr
                        = dynamic_cast<const Parameter<T> *>( this )) {
                        return ptr->value();
                    } else {
                        throw std::out_of_range(
                            "Can't cast parameter to requested type!" );
                    }
                }

                template<typename T>
                T& get() {
                    if (auto ptr = dynamic_cast<Parameter<T> *>( this )) {
                        return ptr->value();
                    } else {
                        throw std::out_of_range(
                            "Can't cast parameter to requested type!" );
                    }
                }

                template<typename T>
                _configuration_parameter& set( T val ) {
                    if (auto ptr = dynamic_cast<Parameter<T> *>( this )) {
                        ptr->value() = val;
                    } else {
                        throw std::out_of_range(
                            "Can't cast parameter to requested type!" );
                    }
                    return *this;
                }

                virtual bool was_set() const = 0;
                virtual std::unique_ptr<_configuration_parameter> clone() const
                    = 0;
                virtual parameter_form_t form() const = 0;
                virtual parameter_type_t type() const = 0;

            protected:
                ValidationTestSuite _tests;
                // parameter_form_t _form;
                // parameter_type_t _type;


                // virtual parameter_form_t form() const { return _form; }

                // virtual parameter_type_t type() const { return _type; }

                // virtual _configuration_parameter& fromString(
                //     const std::string& in ) = 0;

                // // T is boolean specifically
                // template<typename T,
                //          typename std::enable_if<std::is_same<T,
                //          bool>::value,
                //                                  int>::type = 0>
                // T as() const {
                //     return this->as_bool();
                // }

                // // T is integral and not boolean
                // template<typename T,
                //          typename std::enable_if<
                //              ( std::is_integral<T>::value
                //                && !( std::is_same<T, bool>::value ) ),
                //              int>::type = 0>
                // T as() const {
                //     return static_cast<T>( this->as_int() );
                // }

                // // T is floating-point
                // template<typename T,
                //          typename std::enable_if<
                //              std::is_floating_point<T>::value, int>::type =
                //              0>
                // T as() const {
                //     return static_cast<T>( this->as_double() );
                // }

                // // T is string
                // template<typename T,
                //          typename std::enable_if<
                //              std::is_convertible<T, std::string>::value,
                //              int>::type = 0>
                // T as() const {
                //     return this->as_string();
                // }

                // // T is scalar but not fundamental or pointer
                // template<
                //     typename T,
                //     typename std::enable_if<
                //         ( std::is_scalar<T>::value
                //           && !( std::is_fundamental<T>::value
                //                 || std::is_pointer<T>::value
                //                 || std::is_same<T, std::string>::value ) ),
                //         int>::type = 0>
                // T as() const {
                //     return dynamic_cast<const Parameter<T> *>( this
                //     )->value();
                // }

                // // T is a vector


                // virtual bool as_bool() const {
                //     return ( this->as_int() != 0 );
                // }

                // virtual std::string as_string() const {
                //     throw std::logic_error( "as_string(): not implemented"
                //     );
                // }

                // virtual int as_int() const {
                //     throw std::logic_error( "as_int(): not implemented" );
                // }

                // virtual double as_double() const {
                //     throw std::logic_error( "as_double(): not implemented"
                //     );
                // }

                // template<typename T>
                // _configuration_parameter& set( const T& value ) {}
        };

        namespace validation {

            template<typename T>
            class TypedValidationTest : public ValidationTest {
                public:
                    TypedValidationTest() {}

                    TypedValidationTest(
                        const TypedValidationTest<T>& other ) {}

                    virtual ~TypedValidationTest() {}

                    friend void ::swap<>( TypedValidationTest<T>& a,
                                          TypedValidationTest<T>& b ) noexcept;

                    virtual const T& parameter_value(
                        const _configuration_parameter *param ) const {
                        return dynamic_cast<const Parameter<T> *>( param )
                            ->value();
                    }

                    virtual T value( size_t n ) const = 0;
            };

            class NullaryValidationTest : public ValidationTest {
                public:
                    NullaryValidationTest() : ValidationTest() {}

                    NullaryValidationTest(
                        const NullaryValidationTest& other ) {}

                    virtual ~NullaryValidationTest() {}

                    friend void ::swap( NullaryValidationTest& a,
                                        NullaryValidationTest& b ) noexcept;
            };

            template<typename T>
            class UnaryValidationTest : public TypedValidationTest<T> {
                public:
                    UnaryValidationTest( T value ) : TypedValidationTest<T>() {
                        _values.resize( 1 );
                        _values[ 0 ] = value;
                    }

                    UnaryValidationTest(
                        const UnaryValidationTest<T>& other ) {
                        _values = other._values;
                    }

                    virtual ~UnaryValidationTest() {}

                    friend void ::swap<>( UnaryValidationTest<T>& a,
                                          UnaryValidationTest<T>& b ) noexcept;

                    virtual T value( size_t n = 0 ) const override {
                        if (n > 0) {
                            throw std::out_of_range(
                                "Requested index out of range" );
                        }
                        return _values.at( 0 );
                    }

                    virtual const std::vector<T>& values() const {
                        return _values;
                    }

                protected:
                    std::vector<T> _values;
            };

            template<typename T>
            class BinaryValidationTest : public TypedValidationTest<T> {
                public:
                    BinaryValidationTest( T value1, T value2 ) :
                        TypedValidationTest<T>() {
                        _values.resize( 2 );
                        _values[ 0 ] = value1;
                        _values[ 1 ] = value2;
                    }

                    BinaryValidationTest(
                        const BinaryValidationTest<T>& other ) {
                        _values = other._values;
                    }

                    virtual ~BinaryValidationTest() {}

                    friend void ::swap<>(
                        BinaryValidationTest<T>& a,
                        BinaryValidationTest<T>& b ) noexcept;

                    virtual T value( size_t n ) const override {
                        if (n > 1) {
                            throw std::out_of_range(
                                "Requested index out of range" );
                        }
                        return _values.at( n );
                    }

                    virtual const std::vector<T>& values() const {
                        return _values;
                    }

                protected:
                    std::vector<T> _values;
                    // T _value1, _value2;
            };

            template<typename T>
            class ListValidationTest : public TypedValidationTest<T> {
                public:
                    ListValidationTest( std::vector<T> values ) :
                        TypedValidationTest<T>(), _values { values } {}

                    ListValidationTest( std::initializer_list<T> values ) :
                        TypedValidationTest<T>(), _values { values } {}

                    ListValidationTest( const ListValidationTest<T>& other ) {
                        _values = other._values;
                    }

                    virtual ~ListValidationTest() {}

                    friend void ::swap<>( ListValidationTest<T>& a,
                                          ListValidationTest<T>& b ) noexcept;

                    virtual T value( size_t n ) const override {
                        if (n >= _values.size()) {
                            throw std::out_of_range(
                                "Requested index out of range" );
                        }
                        return _values.at( n );
                    }

                    virtual const std::vector<T>& values() const {
                        return _values;
                    }

                protected:
                    std::vector<T> _values;
            };

            // specialization for integers
            template<typename T>
            class IsBetweenTest<
                T, typename std::enable_if<std::is_integral<T>::value>::type>
                : public BinaryValidationTest<T> {
                public:
                    IsBetweenTest( T value1, T value2, bool inclusive ) :
                        BinaryValidationTest<T>( value1, value2 ),
                        _inclusive { inclusive } {}

                    IsBetweenTest( const IsBetweenTest<T>& other ) {
                        _inclusive = other._inclusive;
                    }

                    virtual ~IsBetweenTest() {}

                    friend void ::swap<>( IsBetweenTest<T>& a,
                                          IsBetweenTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsBetweenTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is between " << this->value( 0 ) << " and "
                            << this->value( 1 );
                        if (_inclusive) {
                            oss << " (inclusive)";
                        }
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        int val = param->as_int();
                        if (_inclusive) {
                            if (val >= this->value( 0 )
                                && val <= this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        } else {
                            if (val > this->value( 0 )
                                && val < this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        }
                        return *this;
                    }

                private:
                    bool _inclusive;
            };

            // specialization for floating-point
            template<typename T>
            class IsBetweenTest<T, typename std::enable_if<
                                       std::is_floating_point<T>::value>::type>
                : public BinaryValidationTest<T> {
                public:
                    IsBetweenTest( T value1, T value2, bool inclusive ) :
                        BinaryValidationTest<T>( value1, value2 ),
                        _inclusive { inclusive } {}

                    IsBetweenTest( const IsBetweenTest<T>& other ) {
                        _inclusive = other._inclusive;
                    }

                    virtual ~IsBetweenTest() {}

                    friend void ::swap<>( IsBetweenTest<T>& a,
                                          IsBetweenTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsBetweenTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is between " << this->value( 0 ) << " and "
                            << this->value( 1 );
                        if (_inclusive) {
                            oss << " (inclusive)";
                        }
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        double val = param->as_double();
                        if (_inclusive) {
                            if (val >= this->value( 0 )
                                && val <= this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        } else {
                            if (val > this->value( 0 )
                                && val < this->value( 1 )) {
                                this->pass();
                            } else {
                                this->fail();
                            }
                        }
                        return *this;
                    }

                private:
                    bool _inclusive;
            };

            // specialization for integers
            template<typename T>
            class IsGreaterThanTest<
                T, typename std::enable_if<std::is_integral<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsGreaterThanTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    IsGreaterThanTest( const IsGreaterThanTest<T>& other ) {}

                    virtual ~IsGreaterThanTest() {}

                    friend void ::swap<>( IsGreaterThanTest<T>& a,
                                          IsGreaterThanTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsGreaterThanTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is greater than " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_int() > this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for floating-point
            template<typename T>
            class IsGreaterThanTest<
                T, typename std::enable_if<
                       std::is_floating_point<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsGreaterThanTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    IsGreaterThanTest( const IsGreaterThanTest<T>& other ) {}

                    virtual ~IsGreaterThanTest() {}

                    friend void ::swap<>( IsGreaterThanTest<T>& a,
                                          IsGreaterThanTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsGreaterThanTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is greater than " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_double() > this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for integers
            template<typename T>
            class IsGreaterThanOrEqualToTest<
                T, typename std::enable_if<std::is_integral<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsGreaterThanOrEqualToTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsGreaterThanOrEqualToTest() {}

                    friend void ::swap<>(
                        IsGreaterThanOrEqualToTest<T>& a,
                        IsGreaterThanOrEqualToTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsGreaterThanOrEqualToTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is greater than or equal to " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_int() >= this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for floating-point
            template<typename T>
            class IsGreaterThanOrEqualToTest<
                T, typename std::enable_if<
                       std::is_floating_point<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsGreaterThanOrEqualToTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsGreaterThanOrEqualToTest() {}

                    friend void ::swap<>(
                        IsGreaterThanOrEqualToTest<T>& a,
                        IsGreaterThanOrEqualToTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsGreaterThanOrEqualToTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is greater than or equal to " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_double() >= this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for integers
            template<typename T>
            class IsLessThanTest<
                T, typename std::enable_if<std::is_integral<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsLessThanTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsLessThanTest() {}

                    friend void ::swap<>( IsLessThanTest<T>& a,
                                          IsLessThanTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsLessThanTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is less than " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_int() < this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for floating-point
            template<typename T>
            class IsLessThanTest<T,
                                 typename std::enable_if<
                                     std::is_floating_point<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsLessThanTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsLessThanTest() {}

                    friend void ::swap<>( IsLessThanTest<T>& a,
                                          IsLessThanTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsLessThanTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is less than " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_double() < this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for integers
            template<typename T>
            class IsLessThanOrEqualToTest<
                T, typename std::enable_if<std::is_integral<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsLessThanOrEqualToTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsLessThanOrEqualToTest() {}

                    friend void ::swap<>(
                        IsLessThanOrEqualToTest<T>& a,
                        IsLessThanOrEqualToTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsLessThanOrEqualToTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is less than or equal to " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_int() <= this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // specialization for floating-point
            template<typename T>
            class IsLessThanOrEqualToTest<
                T, typename std::enable_if<
                       std::is_floating_point<T>::value>::type>
                : public UnaryValidationTest<T> {
                public:
                    IsLessThanOrEqualToTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsLessThanOrEqualToTest() {}

                    friend void ::swap<>(
                        IsLessThanOrEqualToTest<T>& a,
                        IsLessThanOrEqualToTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsLessThanOrEqualToTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is less than or equal to " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->as_double() <= this->value()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            template<typename T>
            class IsEqualToTest : public UnaryValidationTest<T> {
                public:
                    IsEqualToTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsEqualToTest() {}

                    friend void ::swap<>( IsEqualToTest<T>& a,
                                          IsEqualToTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsEqualToTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is equal to " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (this->value() == this->parameter_value( param )) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            template<typename T>
            class IsNotEqualToTest : public UnaryValidationTest<T> {
                public:
                    IsNotEqualToTest( T value ) :
                        UnaryValidationTest<T>( value ) {}

                    virtual ~IsNotEqualToTest() {}

                    friend void ::swap<>( IsNotEqualToTest<T>& a,
                                          IsNotEqualToTest<T>& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new IsNotEqualToTest<T>( *this ) );
                    }

                    virtual std::string description() const override {
                        std::ostringstream oss;
                        oss << "is not equal to " << this->value();
                        return oss.str();
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (this->value() != this->parameter_value( param )) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            class WasSetTest : public NullaryValidationTest {
                public:
                    WasSetTest() {}

                    virtual ~WasSetTest() {}

                    friend void ::swap( WasSetTest& a,
                                        WasSetTest& b ) noexcept;

                    virtual std::unique_ptr<ValidationTest> clone()
                        const override {
                        return std::unique_ptr<ValidationTest>(
                            new WasSetTest( *this ) );
                    }

                    virtual std::string description() const override {
                        return "was set";
                    }

                    virtual ValidationTest& test(
                        const _configuration_parameter *param ) override {
                        if (param->was_set()) {
                            this->pass();
                        } else {
                            this->fail();
                        }
                        return *this;
                    }
            };

            // DEFINE_NEW_TESTS_HERE

            // convenience functions
            template<typename T>
            test_ptr_t IsLessThan( const T& val ) {
                return test_ptr_t( new IsLessThanTest<T>( val ) );
            }

            template<typename T>
            test_ptr_t IsLessThan( const T& val, bool equalOK ) {
                return equalOK
                         ? test_ptr_t( new IsLessThanOrEqualToTest<T>( val ) )
                         : test_ptr_t( new IsLessThanTest<T>( val ) );
            }

            template<typename T>
            test_ptr_t IsGreaterThan( const T& val ) {
                return test_ptr_t( new IsGreaterThanTest<T>( val ) );
            }

            template<typename T>
            test_ptr_t IsGreaterThan( const T& val, bool equalOK ) {
                return equalOK ? test_ptr_t(
                                     new IsGreaterThanOrEqualToTest<T>( val ) )
                               : test_ptr_t( new IsGreaterThanTest<T>( val ) );
            }

            template<typename T>
            test_ptr_t IsEqualTo( const T& val ) {
                return test_ptr_t( new IsEqualToTest<T>( val ) );
            }

            template<typename T>
            test_ptr_t IsNotEqualTo( const T& val ) {
                return test_ptr_t( new IsNotEqualToTest<T>( val ) );
            }

            template<typename T>
            test_ptr_t IsZero() {
                return IsEqualTo<T>( (T)0 );
            }

            template<typename T>
            test_ptr_t IsNotZero() {
                return IsNotEqualTo<T>( (T)0 );
            }

            template<typename T>
            test_ptr_t IsNegative() {
                return IsLessThan<T>( (T)0 );
            }

            template<typename T>
            test_ptr_t IsPositive() {
                return IsGreaterThan<T>( (T)0 );
            }

            inline test_ptr_t IsNotEmptyString() {
                return IsNotEqualTo<std::string>( "" );
            }

            inline test_ptr_t IsEmptyString() {
                return IsEqualTo<std::string>( "" );
            }

            // DEFINE_NEW_CONVENIENCE_FUNCTIONS_HERE


        }  // namespace validation

        typedef std::unique_ptr<_configuration_parameter> param_ptr_t;

        template<typename T>
        class Parameter : public _configuration_parameter {
            public:
                Parameter() {}

                Parameter( const T& defaultval ) :
                    _value { defaultval }, _was_set { false } {}

                Parameter( const Parameter<T>& other ) :
                    _configuration_parameter( other ) {
                    _value   = other._value;
                    _was_set = other._was_set;
                }

                Parameter( std::initializer_list<test_ptr_t> tests ) :
                    Parameter<T>() {
                    this->append_tests( tests );
                }

                Parameter( const T& defaultval,
                           std::initializer_list<test_ptr_t> tests ) :
                    Parameter<T>( defaultval ) {
                    this->append_tests( tests );
                }

                virtual ~Parameter() {}

                virtual bool as_bool() const override {
                    return this->_as_bool();
                }

                virtual std::string as_string() const override {
                    return this->_as_string();
                    ;
                }

                virtual int as_int() const override { return this->_as_int(); }

                virtual double as_double() const override {
                    return this->_as_double();
                }

                T& value() { return _value; }

                const T& value() const { return _value; }

                void set( const T& newval ) {
                    _value   = newval;
                    _was_set = true;
                }

                virtual bool was_set() const override { return _was_set; }

                virtual param_ptr_t clone() const override {
                    return param_ptr_t( new Parameter<T>( *this ) );
                }

                virtual parameter_form_t form() const override {
                    return _form();
                }

                virtual parameter_type_t type() const { return _type(); }

            protected:
                T _value;
                bool _was_set = false;

                // normal case: is scalar
                template<typename U                         = T,
                         typename std::enable_if<std::is_scalar<U>::value,
                                                 int>::type = 0>
                parameter_form_t _form() const {
                    return parameter_form_t::SCALAR;
                }

                // special case: is complex
                template<typename U = T,
                         typename std::enable_if<
                             !( std::is_scalar<U>::value
                                || NCPA::types::is_complex<U>::value ),
                             int>::type = 0>
                parameter_form_t _form() const {
                    return parameter_form_t::VECTOR;
                }

                // otherwise vector
                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_scalar<U>::value )
                               && NCPA::types::is_complex<U>::value ),
                             int>::type = 0>
                parameter_form_t _form() const {
                    return parameter_form_t::SCALAR;
                }

                // general type detection
                template<typename U                         = T,
                         typename std::enable_if<std::is_same<U, bool>::value,
                                                 int>::type = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::BOOLEAN;
                }

                template<
                    typename U                         = T,
                    typename std::enable_if<( !( std::is_same<U, bool>::value )
                                              && std::is_integral<U>::value ),
                                            int>::type = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::INTEGER;
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( std::is_floating_point<U>::value ), int>::type
                         = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::INTEGER;
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<U>::value )
                               && std::is_convertible<U, std::string>::value ),
                             int>::type = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::STRING;
                }

                template<typename U                         = T,
                         typename std::enable_if<( std::is_enum<U>::value ),
                                                 int>::type = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::ENUM;
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( NCPA::types::is_complex<U>::value ), int>::type
                         = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::COMPLEX;
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<U>::value
                                  || NCPA::types::is_complex<U>::value
                                  || std::is_convertible<U, std::string>::value
                                  || std::is_enum<U>::value ) ),
                             int>::type = 0>
                parameter_type_t _type() const {
                    return parameter_type_t::COMPLEX;
                }

                // boolean conversions
                template<typename U = T,
                         typename std::enable_if<
                             std::is_convertible<U, bool>::value, int>::type
                         = 0>
                bool _as_bool() const {
                    return static_cast<bool>( this->value() );
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_convertible<U, bool>::value )
                               && std::is_convertible<U, std::string>::value ),
                             int>::type = 0>
                bool _as_bool() const {
                    return ( this->_as_string().size() > 0 );
                }

                template<
                    typename U = T,
                    typename std::enable_if<
                        !( std::is_convertible<U, bool>::value
                           || std::is_convertible<U, std::string>::value ),
                        int>::type = 0>
                bool _as_bool() const {
                    throw std::out_of_range( "Not convertible to bool" );
                }

                // string conversions
                template<typename U                         = T,
                         typename std::enable_if<std::is_same<U, bool>::value,
                                                 int>::type = 0>
                std::string _as_string() const {
                    return ( this->value() ? "true" : "false" );
                }

                template<typename U = T,
                         typename std::enable_if<
                             NCPA::types::is_complex<U>::value, int>::type = 0>
                std::string _as_string() const {
                    std::ostringstream oss;
                    double r = this->value().real();
                    double i = this->value().imag();
                    oss << "(" << r << "," << i << ")";
                    return oss.str();
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_same<U, bool>::value )
                               && std::is_fundamental<U>::value ),
                             int>::type = 0>
                std::string _as_string() const {
                    std::ostringstream oss;
                    oss << this->value();
                    return oss.str();
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_fundamental<U>::value )
                               && std::is_convertible<U, std::string>::value ),
                             int>::type = 0>
                std::string _as_string() const {
                    std::string s = this->value();
                    return s;
                }

                template<
                    typename U = T,
                    typename std::enable_if<
                        ( !( std::is_fundamental<U>::value
                             || std::is_convertible<U, std::string>::value ) ),
                        int>::type = 0>
                std::string _as_string() const {
                    throw std::out_of_range( "Not convertible to string" );
                }

                // integer conversions
                template<typename U = T,
                         typename std::enable_if<
                             ( std::is_integral<U>::value ), int>::type = 0>
                int _as_int() const {
                    return static_cast<int>( this->value() );
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( std::is_floating_point<U>::value ), int>::type
                         = 0>
                int _as_int() const {
                    return static_cast<int>(
                        this->value() + std::numeric_limits<U>::epsilon() );
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_arithmetic<U>::value )
                               && std::is_convertible<U, std::string>::value ),
                             int>::type = 0>
                int _as_int() const {
                    return static_cast<int>(
                        std::stoi( std::string( this->value() ) ) );
                }

                template<
                    typename U = T,
                    typename std::enable_if<
                        !( std::is_arithmetic<U>::value
                           || std::is_convertible<U, std::string>::value ),
                        int>::type = 0>
                int _as_int() const {
                    throw std::out_of_range( "Not convertible to int" );
                }

                // float conversions
                template<typename U = T,
                         typename std::enable_if<
                             std::is_convertible<U, double>::value, int>::type
                         = 0>
                double _as_double() const {
                    return static_cast<double>( this->value() );
                }

                template<typename U = T,
                         typename std::enable_if<
                             ( !( std::is_convertible<U, double>::value )
                               && std::is_convertible<U, std::string>::value ),
                             int>::type = 0>
                double _as_double() const {
                    return static_cast<double>(
                        std::stod( std::string( this->value() ) ) );
                }

                template<
                    typename U = T,
                    typename std::enable_if<
                        !( std::is_convertible<U, double>::value
                           || std::is_convertible<U, std::string>::value ),
                        int>::type = 0>
                double _as_double() const {
                    throw std::out_of_range( "Not convertible to int" );
                }
        };

        using DoubleParameter  = Parameter<double>;
        using IntegerParameter = Parameter<int>;
        using StringParameter  = Parameter<std::string>;
        using BooleanParameter = Parameter<bool>;

        // template<typename T>
        // class Parameter : public _configuration_parameter {
        //     public:
        //         Parameter() {}

        //         Parameter( const T& defaultval ) :
        //             _value { defaultval }, _was_set { false } {}

        //         Parameter( std::initializer_list<test_ptr_t> tests ) :
        //             Parameter<T>() {
        //             this->append_tests( tests );
        //         }

        //         Parameter( const T& defaultval,
        //                    std::initializer_list<test_ptr_t> tests ) :
        //             Parameter<T>( defaultval ) {
        //             this->append_tests( tests );
        //         }

        //         Parameter( const Parameter<T>& other ) :
        //             _configuration_parameter( other ) {
        //             _value   = other._value;
        //             _was_set = other._was_set;
        //         }

        //         virtual ~Parameter() {}

        //         virtual bool as_bool() const override {
        //             return ( this->as_int() != 0 );
        //         }

        //         virtual std::string as_string() const override {
        //             throw std::logic_error( "as_string(): not implemented"
        //             );
        //         }

        //         virtual int as_int() const override {
        //             throw std::logic_error( "as_int(): not implemented" );
        //         }

        //         virtual double as_double() override const {
        //             throw std::logic_error( "as_double(): not implemented"
        //             );
        //         }

        //         T& value() { return _value; }

        //         const T& value() const { return _value; }

        //         void set( const T& newval ) {
        //             _value   = newval;
        //             _was_set = true;
        //         }

        //         virtual bool was_set() const override { return _was_set; }

        //         virtual std::unique_ptr<_configuration_parameter> clone()
        //             const override {
        //             return std::unique_ptr<_configuration_parameter>(
        //                 new Parameter<T>( *this ) );
        //         }

        //     protected:
        //         T _value;
        //         bool _was_set = false;
        // };

        // class StringParameter : public Parameter<std::string> {
        //     public:
        //         StringParameter() : Parameter<std::string>() {}

        //         StringParameter( const std::string& defaultval ) :
        //             Parameter<std::string>( defaultval ) {}

        //         StringParameter( std::initializer_list<test_ptr_t> tests ) :
        //             Parameter<std::string>( tests ) {}

        //         StringParameter( const std::string& defaultval,
        //                          std::initializer_list<test_ptr_t> tests ) :
        //             Parameter<std::string>( defaultval, tests ) {}

        //         virtual ~StringParameter() {}

        //         virtual std::string as_string() const override {
        //             return this->value();
        //         }

        //         virtual int as_int() const override {
        //             return std::stoi( this->value() );
        //         }

        //         virtual double as_double() const override {
        //             return std::stod( this->value() );
        //         }

        //         virtual std::unique_ptr<_configuration_parameter> clone()
        //             const override {
        //             return std::unique_ptr<_configuration_parameter>(
        //                 new StringParameter( *this ) );
        //         }
        // };

        // class IntegerParameter : public Parameter<int> {
        //     public:
        //         IntegerParameter() : Parameter<int>() {}

        //         IntegerParameter( const int& defaultval ) :
        //             Parameter<int>( defaultval ) {}

        //         IntegerParameter( std::initializer_list<test_ptr_t> tests )
        //         :
        //             Parameter<int>( tests ) {}

        //         IntegerParameter( const int& defaultval,
        //                           std::initializer_list<test_ptr_t> tests )
        //                           :
        //             Parameter<int>( defaultval, tests ) {}

        //         virtual ~IntegerParameter() {}

        //         virtual std::string as_string() const override {
        //             std::ostringstream oss;
        //             oss << _value;
        //             return oss.str();
        //         }

        //         virtual int as_int() const override { return this->value();
        //         }

        //         virtual double as_double() const override {
        //             return (double)this->value();
        //         }

        //         virtual std::unique_ptr<_configuration_parameter> clone()
        //             const override {
        //             return std::unique_ptr<_configuration_parameter>(
        //                 new IntegerParameter( *this ) );
        //         }
        // };

        // class DoubleParameter : public Parameter<double> {
        //     public:
        //         DoubleParameter() : Parameter<double>() {}

        //         DoubleParameter( const double& defaultval ) :
        //             Parameter<double>( defaultval ) {}

        //         DoubleParameter( std::initializer_list<test_ptr_t> tests ) :
        //             Parameter<double>( tests ) {}

        //         DoubleParameter( const double& defaultval,
        //                          std::initializer_list<test_ptr_t> tests ) :
        //             Parameter<double>( defaultval, tests ) {}

        //         virtual ~DoubleParameter() {}

        //         virtual std::string as_string() const override {
        //             std::ostringstream oss;
        //             oss << _value;
        //             return oss.str();
        //         }

        //         virtual int as_int() const override {
        //             return (int)( this->value() );
        //         }

        //         virtual double as_double() const override {
        //             return this->value();
        //         }

        //         virtual bool as_bool() const override {
        //             return ( this->value() != 0.0 );
        //         }

        //         virtual std::unique_ptr<_configuration_parameter> clone()
        //             const override {
        //             return std::unique_ptr<_configuration_parameter>(
        //                 new DoubleParameter( *this ) );
        //         }
        // };

        // class BooleanParameter : public Parameter<bool> {
        //     public:
        //         BooleanParameter() : Parameter<bool>() {}

        //         BooleanParameter( const bool& defaultval ) :
        //             Parameter<bool>( defaultval ) {}

        //         BooleanParameter( std::initializer_list<test_ptr_t> tests )
        //         :
        //             Parameter<bool>( tests ) {}

        //         BooleanParameter( const bool& defaultval,
        //                           std::initializer_list<test_ptr_t> tests )
        //                           :
        //             Parameter<bool>( defaultval, tests ) {}

        //         virtual ~BooleanParameter() {}

        //         virtual std::string as_string() const override {
        //             return ( this->value() ? "true" : "false" );
        //         }

        //         virtual int as_int() const override {
        //             return ( this->value() ? 1 : 0 );
        //         }

        //         virtual double as_double() const override {
        //             return ( this->value() ? 1.0 : 0.0 );
        //         }

        //         virtual bool as_bool() const override { return
        //         this->value(); }

        //         virtual std::unique_ptr<_configuration_parameter> clone()
        //             const override {
        //             return std::unique_ptr<_configuration_parameter>(
        //                 new BooleanParameter( *this ) );
        //         }
        // };

        // example:
        // class PropagationModel :
        //      public Configurable< PropagationModel,
        //      model_param_t >{ ... }
        //
        // PropagationModel model;
        // model.add_parameter( model_param_t::FREQUENCY,
        //      DoubleParameter( { validation::IsPositive } ) );
        //
        //

        template<typename KEYTYPE>
        using param_pair_t = std::pair<KEYTYPE, param_ptr_t>;

        template<typename KEYTYPE>
        class ConfigurationMap
            : public std::unordered_map<KEYTYPE, param_ptr_t> {
            public:
                ConfigurationMap() {}

                ConfigurationMap( const ConfigurationMap& other ) {
                    for (auto it = other.begin(); it != other.end(); ++it) {
                        this->emplace(
                            std::make_pair( it->first, it->second->clone() ) );
                    }
                }

                ConfigurationMap( ConfigurationMap&& other ) noexcept {
                    ::swap( *this, other );
                }

                ConfigurationMap& operator=( ConfigurationMap other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ConfigurationMap<KEYTYPE>& a,
                                      ConfigurationMap<KEYTYPE>& b ) noexcept;

                virtual ~ConfigurationMap() {}
        };

        template<typename KEYTYPE>
        class Configurable {
            public:
                Configurable() { this->define_parameters(); }

                Configurable( const Configurable<KEYTYPE>& other ) :
                    Configurable<KEYTYPE>() {
                    _parameters = other._parameters;
                }

                virtual ~Configurable() {}

                void add_parameter( KEYTYPE key,
                                    const _configuration_parameter *param ) {
                    // _parameters.emplace(
                    //     param_pair_t<KEYTYPE>{ key, param->clone() } );
                    _parameters[ key ] = param->clone();
                }

                void add_parameter( KEYTYPE key, const param_ptr_t param ) {
                    // _parameters.emplace(
                    //     param_pair_t<KEYTYPE>{ key, param->clone() } );
                    _parameters[ key ] = param->clone();
                }

                void add_parameter( KEYTYPE key,
                                    const _configuration_parameter& param ) {
                    // _parameters.emplace(
                    //     param_pair_t<KEYTYPE>{ key, param.clone() } );
                    _parameters[ key ] = param.clone();
                }

                _configuration_parameter& parameter( KEYTYPE key ) {
                    return *( _parameters.at( key ).get() );
                }

                const _configuration_parameter& parameter(
                    KEYTYPE key ) const {
                    return *( _parameters.at( key ).get() );
                }

                void copy_parameter( const KEYTYPE& key,
                                     const param_ptr_t& ptr ) {
                    // _parameters[ key ] = ptr->clone();
                    // return static_cast<DERIVEDTYPE&>( *this );
                    this->copy_parameter( key, *ptr );
                }

                void copy_parameter( const KEYTYPE& key,
                                     const _configuration_parameter& ptr ) {
                    _parameters[ key ] = ptr.clone();
                }

                virtual void define_parameters() {}

                void validate_parameters() {
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        it->second->validate();
                    }
                }

                std::string validation_report() const {
                    std::ostringstream oss;
                    validation_report( oss, false );
                    return oss.str();
                }

                std::ostream& validation_report( std::ostream& os,
                                                 bool newline = true ) const {
                    bool firsttime = true;
                    for (auto it = _parameters.begin();
                         it != _parameters.end(); ++it) {
                        if (firsttime) {
                            firsttime = false;
                        } else {
                            os << std::endl;
                        }
                        os << it->first << ":" << std::endl;
                        it->second->validation_report( os, false, "  " );
                    }
                    if (newline) {
                        os << std::endl;
                    }
                    return os;
                }

                bool passed() const {
                    bool pass = true;
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        pass = pass && it->second->passed();
                    }
                    return pass;
                }

                bool failed() const { return !( this->passed() ); }

                std::vector<_configuration_parameter *> invalid() const {
                    std::vector<_configuration_parameter *> inv;
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        if (it->second->failed()) {
                            inv.push_back( it->second.get() );
                        }
                    }
                    return inv;
                }

                template<typename PARAMTYPE>
                void set( KEYTYPE key, PARAMTYPE value ) {
                    if (!has_parameter( key )) {
                        return this->add_parameter(
                            key,
                            param_ptr_t( new Parameter<PARAMTYPE>( value ) ) );
                    } else {
                        if (auto sub = dynamic_cast<Parameter<PARAMTYPE> *>(
                                &this->parameter( key ) )) {
                            sub->set( value );
                        } else {
                            throw std::logic_error(
                                "Can't cast parameter to requested type!" );
                        }
                    }
                }

                template<typename PARAMTYPE>
                const PARAMTYPE& get( KEYTYPE key ) const {
                    if (auto sub = dynamic_cast<const Parameter<PARAMTYPE> *>(
                            &this->parameter( key ) )) {
                        return sub->value();
                    } else {
                        throw std::logic_error(
                            "Can't cast parameter to requested type!" );
                    }
                }

                bool has_parameter( KEYTYPE key ) const {
                    return ( _parameters.find( key ) != _parameters.cend() );
                }

                void copy_parameters_from( const Configurable<KEYTYPE>& other,
                                bool create_if_missing = true ) {
                    for (auto it = other.parameters().cbegin();
                         it != other.parameters().cend(); ++it) {
                        if (create_if_missing
                            || this->has_parameter( it->first )) {
                            this->copy_parameter( it->first,
                                                  it->second->clone() );
                        }
                    }
                }

                void copy_parameters_to( Configurable<KEYTYPE>& other,
                              bool create_if_missing = true ) const {
                    other.copy_parameters_from( *this );
                }

                virtual ConfigurationMap<KEYTYPE>& parameters() {
                    return _parameters;
                }

                virtual const ConfigurationMap<KEYTYPE>& parameters() const {
                    return _parameters;
                }

            private:
                // std::unordered_map<KEYTYPE,
                //                    std::unique_ptr<_configuration_parameter>>
                //     _parameters;
                ConfigurationMap<KEYTYPE> _parameters;
        };

        // template<typename KEYTYPE>
        // class ConfigurationMap
        //     : public std::unordered_map<KEYTYPE, param_ptr_t> {
        //     public:
        //         ConfigurationMap() {}

        //         ConfigurationMap( const ConfigurationMap& other ) {
        //             for (auto it = other.begin(); it != other.end(); ++it) {
        //                 this->emplace(
        //                     std::make_pair( it->first, it->second->clone() )
        //                     );
        //             }
        //         }

        //         ConfigurationMap( ConfigurationMap&& other ) noexcept {
        //             ::swap( *this, other );
        //         }

        //         ConfigurationMap& operator=( ConfigurationMap other ) {
        //             ::swap( *this, other );
        //             return *this;
        //         }

        //         friend void ::swap<>( ConfigurationMap<KEYTYPE>& a,
        //                               ConfigurationMap<KEYTYPE>& b )
        //                               noexcept;

        //         virtual ~ConfigurationMap() {}
        // };

        // template<typename DERIVEDTYPE, typename KEYTYPE>
        // class Configurable {
        //     public:
        //         Configurable() {}

        //         Configurable(
        //             const Configurable<DERIVEDTYPE, KEYTYPE>& other ) :
        //             Configurable<DERIVEDTYPE, KEYTYPE>() {
        //             _parameters = other._parameters;
        //         }

        //         virtual ~Configurable() {}

        //         DERIVEDTYPE& add_parameter(
        //             KEYTYPE key, const _configuration_parameter *param ) {
        //             // _parameters.emplace(
        //             //     param_pair_t<KEYTYPE>{ key, param->clone() } );
        //             _parameters[ key ] = param->clone();
        //             return static_cast<DERIVEDTYPE&>( *this );
        //         }

        //         DERIVEDTYPE& add_parameter(
        //             KEYTYPE key, const _configuration_parameter& param ) {
        //             // _parameters.emplace(
        //             //     param_pair_t<KEYTYPE>{ key, param.clone() } );
        //             _parameters[ key ] = param.clone();
        //             return static_cast<DERIVEDTYPE&>( *this );
        //         }

        //         _configuration_parameter& parameter( KEYTYPE key ) {
        //             return *( _parameters.at( key ).get() );
        //         }

        //         const _configuration_parameter& parameter(
        //             KEYTYPE key ) const {
        //             return *( _parameters.at( key ).get() );
        //         }

        //         DERIVEDTYPE& validate_parameters() {
        //             for (auto it = _parameters.cbegin();
        //                  it != _parameters.cend(); ++it) {
        //                 it->second->validate();
        //             }
        //             return static_cast<DERIVEDTYPE&>( *this );
        //         }

        //         std::string validation_report() const {
        //             std::ostringstream oss;
        //             validation_report( oss, false );
        //             return oss.str();
        //         }

        //         std::ostream& validation_report( std::ostream& os,
        //                                          bool newline = true ) const
        //                                          {
        //             bool firsttime = true;
        //             for (auto it = _parameters.begin();
        //                  it != _parameters.end(); ++it) {
        //                 if (firsttime) {
        //                     firsttime = false;
        //                 } else {
        //                     os << std::endl;
        //                 }
        //                 os << it->first << ":" << std::endl;
        //                 it->second->validation_report( os, false, "  " );
        //             }
        //             if (newline) {
        //                 os << std::endl;
        //             }
        //             return os;
        //         }

        //         bool passed() const {
        //             bool pass = true;
        //             for (auto it = _parameters.cbegin();
        //                  it != _parameters.cend(); ++it) {
        //                 pass = pass && it->second->passed();
        //             }
        //             return pass;
        //         }

        //         bool failed() const { return !( this->passed() ); }

        //         std::vector<_configuration_parameter *> invalid() const {
        //             std::vector<_configuration_parameter *> inv;
        //             for (auto it = _parameters.cbegin();
        //                  it != _parameters.cend(); ++it) {
        //                 if (it->second->failed()) {
        //                     inv.push_back( it->second.get() );
        //                 }
        //             }
        //             return inv;
        //         }

        //         template<typename PARAMTYPE>
        //         DERIVEDTYPE& set( KEYTYPE key, PARAMTYPE value ) {
        //             if (auto sub = dynamic_cast<Parameter<PARAMTYPE> *>(
        //                     &this->parameter( key ) )) {
        //                 sub->set( value );
        //             } else {
        //                 throw std::logic_error(
        //                     "Can't cast parameter to requested type!" );
        //             }
        //             return static_cast<DERIVEDTYPE&>( *this );
        //         }

        //         template<typename PARAMTYPE>
        //         const PARAMTYPE& get( KEYTYPE key ) const {
        //             if (auto sub = dynamic_cast<const Parameter<PARAMTYPE>
        //             *>(
        //                     &this->parameter( key ) )) {
        //                 return sub->value();
        //             } else {
        //                 throw std::logic_error(
        //                     "Can't cast parameter to requested type!" );
        //             }
        //         }

        //     private:
        //         // std::unordered_map<KEYTYPE,
        //         // std::unique_ptr<_configuration_parameter>>
        //         //     _parameters;
        //         ConfigurationMap<KEYTYPE> _parameters;
        // };

    }  // namespace config
}  // namespace NCPA

inline void swap( NCPA::config::ValidationTest& a,
                  NCPA::config::ValidationTest& b ) noexcept {
    using std::swap;
    swap( a._status, b._status );
}

inline void swap( NCPA::config::_configuration_parameter& a,
                  NCPA::config::_configuration_parameter& b ) noexcept {}

template<typename T>
void swap( NCPA::config::Parameter<T>& a,
           NCPA::config::Parameter<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::_configuration_parameter&>( a ),
            static_cast<NCPA::config::_configuration_parameter&>( b ) );
}

inline void swap( NCPA::config::ValidationTestSuite& a,
                  NCPA::config::ValidationTestSuite& b ) noexcept {
    using std::swap;
    swap( a._tests, b._tests );
}

inline void swap(
    NCPA::config::validation::NullaryValidationTest& a,
    NCPA::config::validation::NullaryValidationTest& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::ValidationTest&>( a ),
            static_cast<NCPA::config::ValidationTest&>( b ) );
}

template<typename T>
void swap( NCPA::config::validation::TypedValidationTest<T>& a,
           NCPA::config::validation::TypedValidationTest<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::ValidationTest&>( a ),
            static_cast<NCPA::config::ValidationTest&>( b ) );
}

template<typename T>
void swap( NCPA::config::validation::UnaryValidationTest<T>& a,
           NCPA::config::validation::UnaryValidationTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::TypedValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::TypedValidationTest<T>&>( b ) );
    swap( a._values, b._values );
}

template<typename T>
void swap( NCPA::config::validation::BinaryValidationTest<T>& a,
           NCPA::config::validation::BinaryValidationTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::TypedValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::TypedValidationTest<T>&>( b ) );
    swap( a._values, b._values );
}

template<typename T>
void swap( NCPA::config::validation::ListValidationTest<T>& a,
           NCPA::config::validation::ListValidationTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::TypedValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::TypedValidationTest<T>&>( b ) );
    swap( a._values, b._values );
}

template<typename T>
void swap( NCPA::config::validation::IsBetweenTest<T>& a,
           NCPA::config::validation::IsBetweenTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::BinaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::BinaryValidationTest<T>&>( b ) );
    swap( a._inclusive, b._inclusive );
}

template<typename T>
void swap( NCPA::config::validation::IsGreaterThanTest<T>& a,
           NCPA::config::validation::IsGreaterThanTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( b ) );
}

template<typename T>
void swap(
    NCPA::config::validation::IsGreaterThanOrEqualToTest<T>& a,
    NCPA::config::validation::IsGreaterThanOrEqualToTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( b ) );
}

template<typename T>
void swap( NCPA::config::validation::IsLessThanTest<T>& a,
           NCPA::config::validation::IsLessThanTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( b ) );
}

template<typename T>
void swap( NCPA::config::validation::IsLessThanOrEqualToTest<T>& a,
           NCPA::config::validation::IsLessThanOrEqualToTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( b ) );
}

template<typename T>
void swap( NCPA::config::validation::IsEqualToTest<T>& a,
           NCPA::config::validation::IsEqualToTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( b ) );
}

template<typename T>
void swap( NCPA::config::validation::IsNotEqualToTest<T>& a,
           NCPA::config::validation::IsNotEqualToTest<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( a ),
        static_cast<NCPA::config::validation::UnaryValidationTest<T>&>( b ) );
}

inline void swap( NCPA::config::validation::WasSetTest& a,
                  NCPA::config::validation::WasSetTest& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::NullaryValidationTest&>( a ),
        static_cast<NCPA::config::validation::NullaryValidationTest&>( b ) );
}

template<typename T>
void swap( NCPA::config::ConfigurationMap<T>& a,
           NCPA::config::ConfigurationMap<T>& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<std::unordered_map<T, NCPA::config::param_ptr_t>&>( a ),
        static_cast<std::unordered_map<T, NCPA::config::param_ptr_t>&>( b ) );
}

template<typename KEYTYPE>
void swap( NCPA::config::Configurable<KEYTYPE>& a,
           NCPA::config::Configurable<KEYTYPE>& b ) noexcept {
    using std::swap;
    swap( a._parameters, b._parameters );
}

// DEFINE_SWAP_FUNCTIONS_HERE
