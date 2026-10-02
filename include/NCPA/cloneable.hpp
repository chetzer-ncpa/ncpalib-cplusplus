/**
 * NCPA Cloneable interface
 * @version 1.0.0
 * @author Claus Hetzer
 * @date 2026-09-29
 *
 * Implements a generic base for cloneable members.
 *
 * Conceptually, this library is used when you have a pure virtual base class
 * and a number of derived classes that you want to be polymorphic with that
 * class.  The base class should inherit from CloneBase thus:
 *
 * class AbstractBase : public CloneBase<AbstractBase> {};
 *
 * If it so happens that the base class is not pure virtual, use
 * ConcreteCloneBase instead of CloneBase.
 *
 * Concrete derived classes should then inherit from Cloneable (pure virtual
 * derived classes do not need to do anything special):
 *
 * class Derived : public Cloneable<Derived,Parent> {}
 *
 * where Parent is the immediate parent class, not the overall base class. This
 * allows for all necessary constructor arguments, methods, etc. to be
 * inherited properly.  Derived classes will have two methods available through
 * the interface:
 *
 * std::unique_ptr<AbstractBase> clone() const;
 * std::unique_ptr<AbstractBase> default_clone() const;
 *
 * The clone() method returns a unique_ptr to the base class containing a new
 * instance of the derived class (using its copy constructor). The
 * default_clone does the same except it uses the default constructor to create
 * the copy.  All Cloneable and CloneBase classes must be default-constructible
 * to enable this.
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

            friend void swap( CloneBase<BASE>& a,
                              CloneBase<BASE>& b ) noexcept {}
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
