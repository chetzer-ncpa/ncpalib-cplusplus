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

#include <array>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
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
    }
}  // namespace NCPA

// class _configuration_parameter;
DECLARE_CLASS_AND_SWAP_2NAMESPACE( _configuration_parameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( BooleanParameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( DoubleParameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( IntegerParameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( StringParameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( ValidationTest, NCPA, config )
DECLARE_CLASS_AND_SWAP_2NAMESPACE( ValidationTestSuite, NCPA, config )
DECLARE_TEMPLATE_AND_SWAP_2NAMESPACE( _typed_parameter, NCPA, config )
DECLARE_CLASS_AND_SWAP_3NAMESPACE( NullaryValidationTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( TypedValidationTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( UnaryValidationTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( BinaryValidationTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( ListValidationTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsEqualToTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsNotEqualToTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsNotOneOfTest, NCPA, config, validation )
DECLARE_TEMPLATE_AND_SWAP_3NAMESPACE( IsOneOfTest, NCPA, config, validation )
DECLARE_CLASS_AND_SWAP_3NAMESPACE( WasSetTest, NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsBetweenTest, NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsGreaterThanTest, NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsGreaterThanOrEqualToTest, NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsLessThanTest, NCPA, config, validation )
DECLARE_TYPELIMITED_TEMPLATE_AND_SWAP_3NAMESPACE( IsLessThanOrEqualToTest, NCPA, config, validation )

namespace NCPA {
    namespace config {
        // class BooleanParameter;
        // class DoubleParameter;
        // class IntegerParameter;
        // class StringParameter;
        // class ValidationTest;
        // class ValidationTestSuite;

        typedef std::unique_ptr<ValidationTest> test_ptr_t;

        // template<typename T>
        // class _typed_parameter;

        // namespace validation {
        // class NullaryValidationTest;
        // template<typename T>
        // class TypedValidationTest;
        // template<typename T>
        // class UnaryValidationTest;
        // template<typename T>
        // class BinaryValidationTest;
        // template<typename T>
        // class ListValidationTest;
        // template<typename T, typename Enable = void>
        // class IsBetweenTest;
        // template<typename T, typename Enable = void>
        // class IsGreaterThanTest;
        // template<typename T, typename Enable = void>
        // class IsGreaterThanOrEqualToTest;
        // template<typename T, typename Enable = void>
        // class IsLessThanTest;
        // template<typename T, typename Enable = void>
        // class IsLessThanOrEqualToTest;
        // template<typename T>
        // class IsEqualToTest;
        // template<typename T>
        // class IsNotEqualToTest;
        // template<typename T>
        // class IsOneOfTest;
        // template<typename T>
        // class IsNotOneOfTest;
        // class WasSetTest;

        // DECLARE_NEW_TESTS_HERE

        // }  // namespace validation
    }  // namespace config
}  // namespace NCPA

// DECLARE_SWAP_FUNCTION_2NAMESPACE( _configuration_parameter, NCPA, config )
// DECLARE_SWAP_FUNCTION_2NAMESPACE( BooleanParameter, NCPA, config )
// DECLARE_SWAP_FUNCTION_2NAMESPACE( DoubleParameter, NCPA, config )
// DECLARE_SWAP_FUNCTION_2NAMESPACE( IntegerParameter, NCPA, config )
// DECLARE_SWAP_FUNCTION_2NAMESPACE( StringParameter, NCPA, config )
// DECLARE_SWAP_FUNCTION_2NAMESPACE( ValidationTest, NCPA, config )
// DECLARE_SWAP_FUNCTION_2NAMESPACE( ValidationTestSuite, NCPA, config )
// DECLARE_SWAP_TEMPLATE_2NAMESPACE( _typed_parameter, NCPA, config )
// DECLARE_SWAP_FUNCTION_3NAMESPACE( NullaryValidationTest, NCPA, config,
//                                   validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( TypedValidationTest, NCPA, config,
//                                   validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( UnaryValidationTest, NCPA, config,
//                                   validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( BinaryValidationTest, NCPA, config,
//                                   validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( ListValidationTest, NCPA, config,
//                                   validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( IsEqualToTest, NCPA, config, validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( IsNotEqualToTest, NCPA, config, validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( IsNotOneOfTest, NCPA, config, validation )
// DECLARE_SWAP_TEMPLATE_3NAMESPACE( IsOneOfTest, NCPA, config, validation )
// DECLARE_SWAP_FUNCTION_3NAMESPACE( WasSetTest, NCPA, config, validation )


// void swap( NCPA::config::_configuration_parameter& a,
//            NCPA::config::_configuration_parameter& b ) noexcept;

// void swap( NCPA::config::BooleanParameter& a,
//            NCPA::config::BooleanParameter& b ) noexcept;

// void swap( NCPA::config::DoubleParameter& a,
//            NCPA::config::DoubleParameter& b ) noexcept;

// void swap( NCPA::config::IntegerParameter& a,
//            NCPA::config::IntegerParameter& b ) noexcept;

// void swap( NCPA::config::StringParameter& a,
//            NCPA::config::StringParameter& b ) noexcept;

// void swap( NCPA::config::ValidationTest& a,
//            NCPA::config::ValidationTest& b ) noexcept;

// void swap( NCPA::config::ValidationTestSuite& a,
//            NCPA::config::ValidationTestSuite& b ) noexcept;


// void swap( NCPA::config::validation::NullaryValidationTest& a,
//            NCPA::config::validation::NullaryValidationTest& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::TypedValidationTest<T>& a,
//            NCPA::config::validation::TypedValidationTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::UnaryValidationTest<T>& a,
//            NCPA::config::validation::UnaryValidationTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::BinaryValidationTest<T>& a,
//            NCPA::config::validation::BinaryValidationTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::ListValidationTest<T>& a,
//            NCPA::config::validation::ListValidationTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::IsBetweenTest<T>& a,
//            NCPA::config::validation::IsBetweenTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::IsGreaterThanTest<T>& a,
//            NCPA::config::validation::IsGreaterThanTest<T>& b ) noexcept;
// template<typename T>
// void swap(
//     NCPA::config::validation::IsGreaterThanOrEqualToTest<T>& a,
//     NCPA::config::validation::IsGreaterThanOrEqualToTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::IsLessThanTest<T>& a,
//            NCPA::config::validation::IsLessThanTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::IsLessThanOrEqualToTest<T>& a,
//            NCPA::config::validation::IsLessThanOrEqualToTest<T>& b )
//            noexcept;
// template<typename T>
// void swap( NCPA::config::validation::IsEqualToTest<T>& a,
//            NCPA::config::validation::IsEqualToTest<T>& b ) noexcept;
// template<typename T>
// void swap( NCPA::config::validation::IsNotEqualToTest<T>& a,
//            NCPA::config::validation::IsNotEqualToTest<T>& b ) noexcept;

// void swap( NCPA::config::validation::WasSetTest& a,
//            NCPA::config::validation::WasSetTest& b ) noexcept;

// DECLARE_SWAP_FUNCTIONS_HERE

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
                    const _configuration_parameter *param )
                    = 0;
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
                        return test_status_t::PENDING;
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

                virtual bool as_bool() const {
                    return ( this->as_int() != 0 );
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

                virtual std::string as_string() const = 0;
                virtual int as_int() const            = 0;
                virtual double as_double() const      = 0;
                virtual bool was_set() const          = 0;
                virtual std::unique_ptr<_configuration_parameter> clone() const
                    = 0;

            protected:
                ValidationTestSuite _tests;
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
                        return dynamic_cast<const _typed_parameter<T> *>(
                                   param )
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

            test_ptr_t IsNotEmptyString() {
                return IsNotEqualTo<std::string>( "" );
            }

            test_ptr_t IsEmptyString() {
                return IsEqualTo<std::string>( "" );
            }

            // DEFINE_NEW_CONVENIENCE_FUNCTIONS_HERE


        }  // namespace validation

        template<typename T>
        class _typed_parameter : public _configuration_parameter {
            public:
                _typed_parameter() {}

                _typed_parameter( const T& defaultval ) :
                    _value { defaultval }, _was_set { false } {}

                _typed_parameter( std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<T>() {
                    this->append_tests( tests );
                }

                _typed_parameter( const T& defaultval,
                                  std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<T>( defaultval ) {
                    this->append_tests( tests );
                }

                _typed_parameter( const _typed_parameter<T>& other ) :
                    _configuration_parameter( other ) {
                    _value   = other._value;
                    _was_set = other._was_set;
                }

                virtual ~_typed_parameter() {}

                T& value() { return _value; }

                const T& value() const { return _value; }

                void set( const T& newval ) {
                    _value   = newval;
                    _was_set = true;
                }

                virtual bool was_set() const override { return _was_set; }

            protected:
                T _value;
                bool _was_set = false;
        };

        class StringParameter : public _typed_parameter<std::string> {
            public:
                StringParameter() : _typed_parameter<std::string>() {}

                StringParameter( const std::string& defaultval ) :
                    _typed_parameter<std::string>( defaultval ) {}

                StringParameter( std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<std::string>( tests ) {}

                StringParameter( const std::string& defaultval,
                                 std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<std::string>( defaultval, tests ) {}

                virtual ~StringParameter() {}

                virtual std::string as_string() const override {
                    return this->value();
                }

                virtual int as_int() const override {
                    return std::stoi( this->value() );
                }

                virtual double as_double() const override {
                    return std::stod( this->value() );
                }

                virtual std::unique_ptr<_configuration_parameter> clone()
                    const override {
                    return std::unique_ptr<_configuration_parameter>(
                        new StringParameter( *this ) );
                }
        };

        class IntegerParameter : public _typed_parameter<int> {
            public:
                IntegerParameter() : _typed_parameter<int>() {}

                IntegerParameter( const int& defaultval ) :
                    _typed_parameter<int>( defaultval ) {}

                IntegerParameter( std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<int>( tests ) {}

                IntegerParameter( const int& defaultval,
                                  std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<int>( defaultval, tests ) {}

                virtual ~IntegerParameter() {}

                virtual std::string as_string() const override {
                    std::ostringstream oss;
                    oss << _value;
                    return oss.str();
                }

                virtual int as_int() const override { return this->value(); }

                virtual double as_double() const override {
                    return (double)this->value();
                }

                virtual std::unique_ptr<_configuration_parameter> clone()
                    const override {
                    return std::unique_ptr<_configuration_parameter>(
                        new IntegerParameter( *this ) );
                }
        };

        class DoubleParameter : public _typed_parameter<double> {
            public:
                DoubleParameter() : _typed_parameter<double>() {}

                DoubleParameter( const double& defaultval ) :
                    _typed_parameter<double>( defaultval ) {}

                DoubleParameter( std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<double>( tests ) {}

                DoubleParameter( const double& defaultval,
                                 std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<double>( defaultval, tests ) {}

                virtual ~DoubleParameter() {}

                virtual std::string as_string() const override {
                    std::ostringstream oss;
                    oss << _value;
                    return oss.str();
                }

                virtual int as_int() const override {
                    return (int)( this->value() );
                }

                virtual double as_double() const override {
                    return this->value();
                }

                virtual bool as_bool() const override {
                    return ( this->value() != 0.0 );
                }

                virtual std::unique_ptr<_configuration_parameter> clone()
                    const override {
                    return std::unique_ptr<_configuration_parameter>(
                        new DoubleParameter( *this ) );
                }
        };

        class BooleanParameter : public _typed_parameter<bool> {
            public:
                BooleanParameter() : _typed_parameter<bool>() {}

                BooleanParameter( const bool& defaultval ) :
                    _typed_parameter<bool>( defaultval ) {}

                BooleanParameter( std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<bool>( tests ) {}

                BooleanParameter( const bool& defaultval,
                                  std::initializer_list<test_ptr_t> tests ) :
                    _typed_parameter<bool>( defaultval, tests ) {}

                virtual ~BooleanParameter() {}

                virtual std::string as_string() const override {
                    return ( this->value() ? "true" : "false" );
                }

                virtual int as_int() const override {
                    return ( this->value() ? 1 : 0 );
                }

                virtual double as_double() const override {
                    return ( this->value() ? 1.0 : 0.0 );
                }

                virtual bool as_bool() const override { return this->value(); }

                virtual std::unique_ptr<_configuration_parameter> clone()
                    const override {
                    return std::unique_ptr<_configuration_parameter>(
                        new BooleanParameter( *this ) );
                }
        };

        // example:
        // class PropagationModel :
        //      public Configurable< PropagationModel,
        //      model_param_t >{ ... }
        //
        // PropagationModel model;
        // model.add_parameter( model_param_t::FREQUENCY,
        //      DoubleParameter( { validation::IsPositive } ) );
        template<typename KEYTYPE>
        using param_pair_t
            = std::pair<KEYTYPE, std::unique_ptr<_configuration_parameter>>;

        template<typename DERIVEDTYPE, typename KEYTYPE>
        class Configurable {
            public:
                virtual ~Configurable() {}

                DERIVEDTYPE& add_parameter(
                    KEYTYPE key, const _configuration_parameter *param ) {
                    // _parameters.emplace(
                    //     param_pair_t<KEYTYPE>{ key, param->clone() } );
                    _parameters[ key ] = param->clone();
                    return static_cast<DERIVEDTYPE&>( *this );
                }

                DERIVEDTYPE& add_parameter(
                    KEYTYPE key, const _configuration_parameter& param ) {
                    // _parameters.emplace(
                    //     param_pair_t<KEYTYPE>{ key, param.clone() } );
                    _parameters[ key ] = param.clone();
                    return static_cast<DERIVEDTYPE&>( *this );
                }

                _configuration_parameter& parameter( KEYTYPE key ) {
                    return *( _parameters.at( key ).get() );
                }

                const _configuration_parameter& parameter(
                    KEYTYPE key ) const {
                    return *( _parameters.at( key ).get() );
                }

                DERIVEDTYPE& validate_parameters() {
                    for (auto it = _parameters.cbegin();
                         it != _parameters.cend(); ++it) {
                        it->second->validate();
                    }
                    return static_cast<DERIVEDTYPE&>( *this );
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


            private:
                std::unordered_map<KEYTYPE,
                                   std::unique_ptr<_configuration_parameter>>
                    _parameters;
        };

    }  // namespace config
}  // namespace NCPA

void swap( NCPA::config::ValidationTest& a,
           NCPA::config::ValidationTest& b ) noexcept {
    using std::swap;
    swap( a._status, b._status );
}

void swap( NCPA::config::_configuration_parameter& a,
           NCPA::config::_configuration_parameter& b ) noexcept {}

template<typename T>
void swap( NCPA::config::_typed_parameter<T>& a,
           NCPA::config::_typed_parameter<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::config::_configuration_parameter&>( a ),
            static_cast<NCPA::config::_configuration_parameter&>( b ) );
}

void swap( NCPA::config::ValidationTestSuite& a,
           NCPA::config::ValidationTestSuite& b ) noexcept {
    using std::swap;
    swap( a._tests, b._tests );
}

void swap( NCPA::config::validation::NullaryValidationTest& a,
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

void swap( NCPA::config::validation::WasSetTest& a,
           NCPA::config::validation::WasSetTest& b ) noexcept {
    using std::swap;
    ::swap(
        static_cast<NCPA::config::validation::NullaryValidationTest&>( a ),
        static_cast<NCPA::config::validation::NullaryValidationTest&>( b ) );
}

// DEFINE_SWAP_FUNCTIONS_HERE
