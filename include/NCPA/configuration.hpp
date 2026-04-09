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

#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/Configurable.hpp"
#include "NCPA/configuration/ConfigurationMap.hpp"
#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/ScalarParameter.hpp"
#include "NCPA/configuration/ScalarParameterWithUnits.hpp"
#include "NCPA/configuration/TypedParameter.hpp"
#include "NCPA/configuration/validation.hpp"
#include "NCPA/configuration/ValidationTest.hpp"
#include "NCPA/configuration/ValidationTestSuite.hpp"
#include "NCPA/configuration/VectorParameter.hpp"
#include "NCPA/configuration/VectorParameterWithUnits.hpp"

// Make sure this is included last
#include "NCPA/configuration/functions.hpp"
#include "NCPA/configuration/Parameter.hpp"
