#pragma once

#include <concepts>
#include <ranges>

namespace NCPA::concepts {
    template<typename T>
    concept deleteable = std::destructible<T> && ( !std::is_fundamental_v<T> );

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

}  // namespace NCPA::concepts
