#pragma once

#ifdef HAVE_NLOHMANN_JSON_HPP
#  undef HAVE_NLOHMANN_JSON_HPP
#endif

#if __has_include( "nlohmann/json.hpp" )
#  include "nlohmann/json.hpp"
#  define HAVE_NLOHMANN_JSON_HPP true
#else
#  define HAVE_NLOHMANN_JSON_HPP false
#endif
