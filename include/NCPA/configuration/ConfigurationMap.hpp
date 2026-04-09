#pragma once


#include "NCPA/configuration/declarations.hpp"
#include "NCPA/configuration/BaseParameter.hpp"
#include "NCPA/configuration/TypedParameter.hpp"

#include <unordered_map>

namespace NCPA {
    namespace config {
        template<typename KEYTYPE>
        class ConfigurationMap
            : public std::unordered_map<KEYTYPE, param_ptr_t> {
            public:
                ConfigurationMap() {}

                ConfigurationMap( const ConfigurationMap& other ) {
                    for (auto it = other.begin(); it != other.end(); ++it) {
                        this->emplace(
                            std::make_pair( it->first, it->second->clone() ) );
                    }
                }

                ConfigurationMap( ConfigurationMap&& other ) noexcept {
                    ::swap( *this, other );
                }

                ConfigurationMap& operator=( ConfigurationMap other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap<>( ConfigurationMap<KEYTYPE>& a,
                                      ConfigurationMap<KEYTYPE>& b ) noexcept;

                virtual ~ConfigurationMap() {}
        };
    }  // namespace config
}  // namespace NCPA

template<typename KEYTYPE>
void swap( NCPA::config::ConfigurationMap<KEYTYPE>& a,
           NCPA::config::ConfigurationMap<KEYTYPE>& b ) noexcept {
    using std::swap;
    swap( static_cast<std::unordered_map<KEYTYPE, NCPA::config::param_ptr_t>&>(
              a ),
          static_cast<std::unordered_map<KEYTYPE, NCPA::config::param_ptr_t>&>(
              b ) );
}
