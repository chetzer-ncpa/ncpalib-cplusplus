/**
 * Concepts library, to extend <concepts>.
 *
 * This library defines some custom concepts for specializing templates.
 *
 * Examples:
 * ==================================
 * General concepts usage
 * ==================================
 * To designate different overloaded templates using a single criterion, e.g.
 * for integer and floating-point values, the general idea is to use:
 *
 * template<typename T>
 * requires std::floating_point<T>
 * std::vector<T> random_numbers( size_t N, T minrange, T maxrange )
 * { ... }
 *
 * template<typename T>
 * requires std::integral<T>
 * std::vector<T> random_numbers( size_t N, T minrange, T maxrange )
 * { ... }
 *
 *
 * ========
 * Examples
 * ========
 *
 */

#pragma once


#if __cplusplus >= 202002L
namespace NCPA::concepts {
#  include <concepts>
#  include <ranges>

        template<typename T>
        concept deleteable
            = std::destructible<T> && ( !std::is_fundamental_v<T> );

        template<typename T>
        concept iterable = std::ranges::range<T>;

        // template <std::ranges::range T, typename U>
        // concept EmbeddedTypeMatches = requires(T t) {
        //     { t.x } -> std::same_as<U>;
        // };

        template<typename R, typename V>
        concept iterable_of
            = std::ranges::range<R>
           && std::convertible_to<V, std::ranges::range_value_t<R>>;
}
#endif