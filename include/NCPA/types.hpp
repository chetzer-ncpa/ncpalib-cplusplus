/**
 * Type traits library, to extend <type_traits>.
 *
 * This library defines some custom type_traits for specializing templates,
 * etc.
 *
 * Examples:
 * ==================================
 * General type_traits/concepts usage
 * ==================================
 * To designate different overloaded templates using a single criterion, e.g.
 * for integer and floating-point values, the general idea is to use:
 *
 * (before C++20, using type_traits and enable_if)
 * template<typename T>
 * std::vector<T> random_numbers( size_t N, T minrange, T maxrange,
 *          typename std::enable_if<std::is_floating_point<T>::value,int>::type
 *              ENABLER = 0 ) { ... }
 *
 * template<typename T> std::vector<T>
 * random_numbers( size_t N, T minrange, T maxrange,
 *          typename std::enable_if<std::is_integral<T>::value, int>::type
 *              ENABLER = 0 ) { ... }
 *
 * ==========================
 * Examples from this library
 * ==========================
 * Consider an is_deleteable() trait.  This should return positive if you
 * can call delete on it.  In this case, it is implemented using two traits:
 * is_destructible<T> to check that it can be destroyed, and is_fundamental<T>
 * to see if it is a fundamental type (which can be destroyed but not deleted).
 * The declaration for the new trait thus is:
 *
 * template<typename T>
 * struct is_deleteable {
 *      static constexpr bool value =
 *          std::is_destructible<T>::value
 *          && (!std::is_fundamental<T>::value );
 * };
 *
 * To define a trait to see whether a particular method is implemented, we
 * leverage the void_t type, which maps any sequence of types to void.  First
 * define a default template trait that inherits from false_type, such that it
 * always evaluates as false.  Then, define a specialized template with a
 * second parameter that only evaluates validly if the condition is met, i.e.
 * the method in question can be validly called.  This is called SFINAE.  For
 * example, to test if the dereference (*) operator is valid, we do:
 *
 * template<typename T, typename = void>
 *      struct is_dereferenceable : std::false_type {};
 *
 * template<typename T>
 *      struct is_dereferenceable<T,
 *          details::void_t<decltype(*std::declval<T>() )>> : std::true_type
 * {};
 *
 * Now, consider the is_iterable<T>() trait.  This should return positive if
 * T has begin() and end() methods, and those methods return iterators that
 * support the dereference operator.  However, if it is implemented just like
 * this, it will fail to compile if given e.g. a fundamental type like double,
 * because you can't test for the methods of a fundamental type, it's invalid.
 * Thus we also add a check to make sure T is a class, but that doesn't fix the
 * problem because in the instantiation of the template instance, logical
 * operators don't short-circuit.  So we have to do it in two steps.  First, we
 * define a helper function that checks that the three iterator functions
 * return valid objects (iterators) if passed an argument that evaluates to
 * true, and another identical function that returns false if passed an
 * argument that evaluates to false.  Normally in this function we would check
 * that the function is a function, e.g. for the function size():
 * std::is_member_function_pointer<decltype(&T::size)>::value However, if
 * methods are overloaded this will fail, as in the case of the begin() and
 * end() methods.  We therefore have to do a little sleight of hand, turning
 * the type T into a value of type T, calling begin() on that value, and
 * checking the return type to ensure it's an object.  We also test the
 * validity of the dereference operator as above.  The final product looks like
 * this:
 *
 * // helper function:
 *  template<typename T>
 *  constexpr bool _hasIteratorFunctions( std::true_type ) {
 *     return std::is_object<decltype( std::declval<T>().begin() )>::value
 *            && std::is_object<decltype( std::declval<T>().end() )>::value
 *            && is_dereferenceable<
 *                        decltype( std::declval<T>().begin() )>::value;
 * }
 *
 * template<typename T>
 * constexpr bool _hasIteratorFunctions( std::false_type ) {
 *     return false;
 * }
 *
 * template<class T>
 * struct is_iterable {
 *      static constexpr bool value
 *          = _hasIteratorFunctions<T>( std::is_class<T> {}
 * );
 *
 *
 * =================
 * Examples of Usage
 * =================
 * template<typename T>
 * void test_dereferenceable( T obj,
 *    typename std::enable_if<is_dereferenceable<T>::value, int>::type ENABLER
 * = 0 )
 * {}
 *
 * template<typename T>
 * void test_dereferenceable( T obj,
 *    typename std::enable_if<!is_dereferenceable<T>::value, int>::type ENABLER
 * = 0 ) { throw std::invalid_argument( "Does not satisfy is_dereferenceable"
 * );
 * }
 */


#pragma once

#include <type_traits>
#include <utility>

// #define ENABLE_IF( CONDITION ) \
//     typename std::enable_if<CONDITION::value, int>::type ENABLER = 0
// #define ENABLE_AND( CONDITION1, CONDITION2 ) \
//     typename std::enable_if<CONDITION1::value, int>::type ENABLER1 = 0, \
//     typename std::enable_if<CONDITION2::value, int>::type ENABLER2 = 0
// #define NO_DEFAULT_ENABLE_IF( CONDITION ) \
//     typename std::enable_if<CONDITION::value, int>::type ENABLER
// #define ENABLE_IF_T( CONDITION, T ) \
//     typename std::enable_if<CONDITION<T>::value, int>::type ENABLER = 0
// #define ENABLE_IF_TU( CONDITION, T, U ) \
//     typename std::enable_if<CONDITION<T, U>::value, int>::type ENABLER = 0

// #define ENABLE_IF_REAL( T ) ENABLE_IF( std::is_floating_point<T> )

namespace NCPA {
    namespace types {
        namespace details {
            // In C++11 we have to define void_t ourselves.

            // Helper class template to turn any type or types into void.
            template<typename... Ts>
            struct make_void {
                    typedef void type;
            };

            // declares the templated type alias
            template<typename... Ts>
            using void_t = typename make_void<Ts...>::type;
        }  // namespace details

        

        // decltype( *std::declval<T>() ):
        // declval<T>.method() lets you get a dummy instance of t.method()
        // without having to actually create an object.  It can't be actually
        // evaluated but you can use decltype() to get the type that the method
        // returns.  In this case, we are getting the type returned by the *
        // operator.  If the type actually exists (i.e. the dereference
        // operator is implemented by T), then the void_t template will exist
        // and the template function will be used.  If it doesn't, then the
        // void_t template can't be built, and we will default to the first
        // template, which always evaluates false.
        template<typename T, typename = void>
        struct is_dereferenceable : std::false_type {};

        template<typename T>
        struct is_dereferenceable<
            T, details::void_t<decltype( *std::declval<T>() )>>
            : std::true_type {};

        template<typename T, typename = void, typename = void>
        struct is_complex : std::false_type {};

        template<typename T>
        struct is_complex<T, details::void_t<decltype( std::declval<T>().real() )>,
        details::void_t<decltype( std::declval<T>().imag() )>>
            : std::true_type {};

        // Tester for deleteable: destructible and not a fundamental type.
        // Good example of a simple compound trait.
        template<typename T>
        struct is_deleteable {
                static constexpr bool value
                    = std::is_destructible<T>::value
                   && ( !std::is_fundamental<T>::value )
                   && ( !is_complex<T>::value );
        };

        namespace details {
            // // Testers for complex: is a class, and has real() and imag()
            // // methods
            // template<typename T>
            // constexpr bool _hasComplexFunctions( std::true_type ) {
            //     return std::is_object<
            //                decltype( std::declval<T>().real() )>::value
            //         && std::is_object<
            //                decltype( std::declval<T>().imag() )>::value;
            // }

            // template<typename T>
            // constexpr bool _hasComplexFunctions( std::false_type ) {
            //     return false;
            // }

            // Testers for iterables: is a class, and has begin(), end(), and
            // operator*() methods
            template<typename T>
            constexpr bool _hasIteratorFunctions( std::true_type ) {
                return std::is_object<
                           decltype( std::declval<T>().begin() )>::value
                    && std::is_object<
                           decltype( std::declval<T>().end() )>::value
                    && is_dereferenceable<
                           decltype( std::declval<T>().begin() )>::value;
            }

            template<typename T>
            constexpr bool _hasIteratorFunctions( std::false_type ) {
                return false;
            }
        }  // namespace details

        // // Complex type
        // template<typename T>
        // struct is_complex {
        //         static constexpr bool value
        //             = details::_hasComplexFunctions<T>( std::is_class<T> {} );
        // };

        // Numeric type: is_integral() or is_floating_point()
        template<typename T>
        struct is_numeric {
                static constexpr bool value
                    = std::is_arithmetic<T>::value || is_complex<T>::value;
        };

        

        
        template<class T>
        struct is_iterable {
                static constexpr bool value
                    = details::_hasIteratorFunctions<T>( std::is_class<T> {} );
        };

        template<typename T, typename U>
        struct is_iterable_of {
                static constexpr bool value
                    = ( is_iterable<T>::value
                        && std::is_convertible<
                            U, typename T::value_type>::value );
        };

        

    }  // namespace types
}  // namespace NCPA
