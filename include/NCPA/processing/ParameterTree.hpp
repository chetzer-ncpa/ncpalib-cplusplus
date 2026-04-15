#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/Parameter.hpp"

#include <memory>
#include <string>
#include <vector>

namespace NCPA {
    namespace processing {
        class _parameter_tree_node;
    }
}  // namespace NCPA

void swap( NCPA::processing::_parameter_tree_node& a,
           NCPA::processing::_parameter_tree_node& b ) noexcept;
void swap( NCPA::processing::ParameterTree& a,
           NCPA::processing::ParameterTree& b ) noexcept;

namespace NCPA {
    namespace processing {

        class _parameter_tree_node {
            public:
                _parameter_tree_node() : _key { "" } {}

                _parameter_tree_node( const std::string& key,
                                      const parameter_ptr_t& param ) :
                    _parameter_tree_node( key, *param ) {}

                _parameter_tree_node( const std::string& key,
                                      const Parameter& param ) :
                    _key { key }, _parameter { param.clone() } {}

                _parameter_tree_node( const _parameter_tree_node& other ) :
                    _key { other._key },
                    _parameter { other._parameter->clone() } {}

                _parameter_tree_node( _parameter_tree_node&& other ) noexcept {
                    ::swap( *this, other );
                }

                _parameter_tree_node& operator=( _parameter_tree_node other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( _parameter_tree_node& a,
                                    _parameter_tree_node& b ) noexcept;

                std::string key() const { return _key; }

                parameter_ptr_t& parameter() { return _parameter; }

                const parameter_ptr_t& parameter() const { return _parameter; }

                _parameter_tree_node& set( const Parameter& param ) {
                    _parameter = param.clone();
                    return *this;
                }

                _parameter_tree_node& set( const parameter_ptr_t& param ) {
                    return this->set( *param );
                }

            private:
                std::string _key;
                parameter_ptr_t _parameter;
        };

        class ParameterTree {
            public:
                ParameterTree() {}

                ParameterTree( const ParameterTree& other ) :
                    _nodes { other._nodes } {}

                ParameterTree( ParameterTree&& other ) noexcept {
                    ::swap( *this, other );
                }

                virtual ~ParameterTree() {}

                ParameterTree& operator=( ParameterTree other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( ParameterTree& a,
                                  ParameterTree& b ) noexcept;

                virtual ParameterTree& add( const std::string& tag,
                                            const Parameter& param ) {
                    for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
                        if (it->key() == tag) {
                            it->set( param );
                            return *this;
                        }
                    }
                    _nodes.emplace_back( tag, param );
                    return *this;
                }

                virtual ParameterTree& add( const std::string& tag,
                                            const parameter_ptr_t& param ) {
                    return this->add( tag, *param );
                }

                std::vector<parameter_ptr_t> get(
                    const std::string& key ) const {
                    std::vector<parameter_ptr_t> params;
                    for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
                        if (it->key() == key) {
                            params.emplace_back( it->parameter()->clone() );
                        }
                    }
                    return params;
                }

            private:
                std::vector<_parameter_tree_node> _nodes;
        };
    }  // namespace processing
}  // namespace NCPA

inline void swap( NCPA::processing::_parameter_tree_node& a,
           NCPA::processing::_parameter_tree_node& b ) noexcept {
    using std::swap;
    swap( a._key, b._key );
    swap( a._parameter, b._parameter );
}

inline void swap( NCPA::processing::ParameterTree& a,
           NCPA::processing::ParameterTree& b ) noexcept {
    using std::swap;
    swap( a._nodes, b._nodes );
}
