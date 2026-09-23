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

namespace NCPA {
    template<typename BASE>
    class CloneBase {
        public:
            using BaseType = BASE;

            CloneBase()                             = default;
            virtual ~CloneBase()                    = default;
            CloneBase( const CloneBase<BASE>& )     = default;
            CloneBase( CloneBase<BASE>&& ) noexcept = default;

            virtual std::unique_ptr<BASE> clone() const         = 0;
            virtual std::unique_ptr<BASE> default_clone() const = 0;

            friend void swap( CloneBase<BASE>& a, CloneBase<BASE>& b ) noexcept {}
    };

    template<typename BASE>
    class ConcreteCloneBase : public CloneBase<BASE> {
        public:
            ConcreteCloneBase()                                     = default;
            virtual ~ConcreteCloneBase()                            = default;
            ConcreteCloneBase( const ConcreteCloneBase<BASE>& )     = default;
            ConcreteCloneBase( ConcreteCloneBase<BASE>&& ) noexcept = default;

            std::unique_ptr<BASE> clone() const override {
                return std::unique_ptr<BASE>(
                    new BASE( static_cast<const BASE&>( *this ) ) );
            }

            std::unique_ptr<BASE> default_clone() const override {
                return std::unique_ptr<BASE>( new BASE() );
            }

            friend void swap( ConcreteCloneBase<BASE>& a,
                       ConcreteCloneBase<BASE>& b ) noexcept {
                using std::swap;
                swap( static_cast<CloneBase<BASE>&>( a ),
                      static_cast<CloneBase<BASE>&>( b ) );
            }
    };

    template<typename DERIVED, typename BASE_OR_INTERFACE>
    class Cloneable : public BASE_OR_INTERFACE {
        public:
            using cloneable_t = Cloneable<DERIVED, BASE_OR_INTERFACE>;

            // Cloneable() {}
            template<typename... Args>
            explicit Cloneable( Args&&...args ) :
                BASE_OR_INTERFACE( std::forward<Args>( args )... ) {}

            virtual ~Cloneable() = default;

            Cloneable( const Cloneable<DERIVED, BASE_OR_INTERFACE>& other ) :
                BASE_OR_INTERFACE( other ) {}

            Cloneable(
                Cloneable<DERIVED, BASE_OR_INTERFACE>&& other ) noexcept {
                swap( *this, other );
            }

            friend void swap(
                Cloneable<DERIVED, BASE_OR_INTERFACE>& a,
                Cloneable<DERIVED, BASE_OR_INTERFACE>& b ) noexcept {
                using std::swap;
                swap( static_cast<BASE_OR_INTERFACE&>( a ),
                      static_cast<BASE_OR_INTERFACE&>( b ) );
            }

            std::unique_ptr<typename BASE_OR_INTERFACE::BaseType> clone()
                const override {
                return std::unique_ptr<typename BASE_OR_INTERFACE::BaseType>(
                    new DERIVED( *static_cast<const DERIVED *>( this ) ) );
            }

            std::unique_ptr<typename BASE_OR_INTERFACE::BaseType>
                default_clone() const override {
                static_assert(
                    std::is_default_constructible<DERIVED>::value,
                    "Cloneable classes must have a default constructor" );
                return std::unique_ptr<typename BASE_OR_INTERFACE::BaseType>(
                    new DERIVED() );
                // return std::unique_ptr<BASE>( new DERIVED() );
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
