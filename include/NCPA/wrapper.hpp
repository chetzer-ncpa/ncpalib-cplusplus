#pragma once

#include <memory>
#include <type_traits>

namespace NCPA {
    template<typename T>
    class Wrapper {
            static_assert( std::has_virtual_destructor<T>::value,
                           "Base must have a virtual destructor" );

        public:
            typedef std::unique_ptr<T> ptr_t;

        private:
            ptr_t _ptr;

        public:
            // Construct from base-type clone
            explicit Wrapper( ptr_t&& p ) : _ptr( std::move( p ) ) {}

            // Construct from any derived type U
            template<typename U>
            Wrapper( std::unique_ptr<U>&& p,
                     typename std::enable_if<
                         std::is_base_of<T, U>::value>::type * = 0 ) :
                _ptr( std::move( p ) ) {}

            template<typename U>
            Wrapper( const U& p,
                     typename std::enable_if<
                         std::is_base_of<T, U>::value
                         && can_clone_to<U, T>::value>::type * = 0 ) :
                _ptr( p.clone() ) {}

            template<typename U>
            Wrapper( U&& p, typename std::enable_if<
                                std::is_base_of<T, U>::value
                                && can_clone_to<U, T>::value>::type * = 0 ) :
                _ptr( p.clone() ) {}

            // 5. Copy constructor (polymorphic clone)
            Wrapper( const Wrapper& other ) :
                _ptr( other._ptr ? other._ptr->clone() : nullptr ) {}

            // 6. Move constructor
            Wrapper( Wrapper&& other ) noexcept :
                _ptr( std::move( other._ptr ) ) {}

            // 7. Copy assignment
            Wrapper& operator=( const Wrapper& other ) {
                if (this != &other) {
                    _ptr = other._ptr ? other._ptr->clone() : nullptr;
                }
                return *this;
            }

            // 8. Move assignment
            Wrapper& operator=( Wrapper&& other ) noexcept {
                _ptr = std::move( other._ptr );
                return *this;
            }

            friend void swap( Wrapper<T>& a, Wrapper<T>& b ) noexcept {
                using std::swap;
                swap( a._ptr, b._ptr );
            }

            // T *wrap( const T& target ) {
            //     _ptr = std::unique_ptr<T>( new T( target ) );
            //     return _ptr.get();
            // }

            // T *wrap( T&& target ) {
            //     _ptr = std::unique_ptr<T>( new T( std::move( target ) ) );
            //     return _ptr.get();
            // }

            T *operator->() { return _ptr.get(); }

            const T *operator->() const { return _ptr.get(); }
    };

}  // namespace NCPA
