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

#define NCPA_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE )               \
    virtual std::unique_ptr<PARENTTYPE> clone() const override {     \
        return std::unique_ptr<PARENTTYPE>( new THISTYPE( *this ) ); \
    }

#define NCPA_FRESH_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE )           \
    virtual std::unique_ptr<PARENTTYPE> fresh_clone() const override { \
        return std::unique_ptr<PARENTTYPE>( new THISTYPE() );          \
    }

#define NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )  \
    NCPA_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE ) \
    NCPA_FRESH_CLONE_METHOD_ONLY( THISTYPE, PARENTTYPE )

// #define NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )                      \
//     virtual std::unique_ptr<PARENTTYPE> clone() const override {       \
//         return std::unique_ptr<PARENTTYPE>( new THISTYPE( *this ) );   \
//     }                                                                  \
//     virtual std::unique_ptr<PARENTTYPE> fresh_clone() const override { \
//         return std::unique_ptr<PARENTTYPE>( new THISTYPE() );          \
//     }

#define DECLARE_NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )       \
    virtual std::unique_ptr<PARENTTYPE> clone() const override; \
    virtual std::unique_ptr<PARENTTYPE> fresh_clone() const override;

#define DEFINE_NCPA_CLONE_METHOD( THISTYPE, PARENTTYPE )             \
    std::unique_ptr<PARENTTYPE> THISTYPE::clone() const {            \
        return std::unique_ptr<PARENTTYPE>( new THISTYPE( *this ) ); \
    }                                                                \
    std::unique_ptr<PARENTTYPE> THISTYPE::fresh_clone() const {      \
        return std::unique_ptr<PARENTTYPE>( new THISTYPE() );        \
    }

namespace NCPA {
    template<typename BASE>
    class Cloneable {
        public:
            virtual ~Cloneable() {}

            virtual std::unique_ptr<BASE> clone() const       = 0;
            virtual std::unique_ptr<BASE> fresh_clone() const = 0;
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
