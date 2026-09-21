/**
 * NCPA pointers library
 * @version 1.0.0
 * @author Claus Hetzer
 * @date 2026-04-14
 *
 * Implements a generic base for cloneable members.
 */
#pragma once

#include <memory>
#include <type_traits>
#include <utility>

// #define NCPA_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE )               \
//     virtual std::unique_ptr<PARENTTYPE> clone() const override {     \
//         return std::unique_ptr<PARENTTYPE>( new THISTYPE( *this ) ); \
//     }

// #define NCPA_FRESH_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE )           \
//     virtual std::unique_ptr<PARENTTYPE> fresh_clone() const override { \
//         return std::unique_ptr<PARENTTYPE>( new THISTYPE() );          \
//     }

// #define NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )  \
//     NCPA_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE ) \
//     NCPA_FRESH_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE )

// #define DECLARE_NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )       \
//     virtual std::unique_ptr<PARENTTYPE> clone() const override; \
//     virtual std::unique_ptr<PARENTTYPE> fresh_clone() const override;

// #define DEFINE_NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )             \
//     std::unique_ptr<PARENTTYPE> THISTYPE::clone() const {            \
//         return std::unique_ptr<PARENTTYPE>( new THISTYPE( *this ) ); \
//     }                                                                \
//     std::unique_ptr<PARENTTYPE> THISTYPE::fresh_clone() const {      \
//         return std::unique_ptr<PARENTTYPE>( new THISTYPE() );        \
//     }

namespace NCPA {

    // template<typename BASE>
    // class Cloneable {
    //     public:
    //         Cloneable()                             = default;
    //         virtual ~Cloneable()                    = default;
    //         Cloneable( const Cloneable<BASE>& )     = default;
    //         Cloneable( Cloneable<BASE>&& ) noexcept = default;

    //         friend void swap( Cloneable<BASE>& a,
    //                           Cloneable<BASE>& b ) noexcept {}

    //         virtual std::unique_ptr<BASE> clone() const         = 0;
    //         virtual std::unique_ptr<BASE> default_clone() const = 0;
    // };

    template<typename BASE>
    class CloneBase {
        public:
            CloneBase()                             = default;
            virtual ~CloneBase()                    = default;
            CloneBase( const CloneBase<BASE>& )     = default;
            CloneBase( CloneBase<BASE>&& ) noexcept = default;

            virtual std::unique_ptr<BASE> clone() const         = 0;
            virtual std::unique_ptr<BASE> default_clone() const = 0;
    };

    template<typename DERIVED, typename BASE>
    class Cloneable : public virtual CloneBase<BASE> {
        public:
            Cloneable()                                      = default;
            virtual ~Cloneable()                             = default;
            Cloneable( const Cloneable<DERIVED, BASE>& )     = default;
            Cloneable( Cloneable<DERIVED, BASE>&& ) noexcept = default;

            friend void swap( Cloneable<DERIVED, BASE>& a,
                              Cloneable<DERIVED, BASE>& b ) noexcept {}

            virtual std::unique_ptr<BASE> clone() const override {
                return std::unique_ptr<BASE>(
                    new DERIVED( *static_cast<const DERIVED *>( this ) ) );
            }

            // template<typename U = DERIVED,
            //          typename std::enable_if<
            //              std::is_default_constructible<U>::value, int>::type
            //              ENABLER = 0>
            std::unique_ptr<BASE> default_clone() const override {
                static_assert(
                    std::is_default_constructible<DERIVED>::value,
                    "Cloneable methods must have a default constructor" );
                return std::unique_ptr<BASE>( new DERIVED() );
            }
    };

    template<typename U, typename Base>
    class can_clone_to {
        private:
            template<typename X>
            static auto test( int ) -> typename std::is_convertible<
                decltype( std::declval<const X>().clone() ),
                std::unique_ptr<Base>>::type;

            template<typename>
            static std::false_type test( ... );

        public:
            static const bool value = decltype( test<U>( 0 ) )::value;
    };
}  // namespace NCPA
