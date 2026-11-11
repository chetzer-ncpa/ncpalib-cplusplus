#pragma once
/*





*/


#include "NCPA/arrays.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/exceptions.hpp"
#include "NCPA/cloneable.hpp"
#include "NCPA/interpolation.hpp"
#include "NCPA/math.hpp"
#include "NCPA/units.hpp"

#include <array>
#include <memory>
#include <vector>

namespace NCPA {
    namespace atmos {
        using NCPA::math::axis_t;
        using namespace NCPA::interpolation;

        namespace details {
            class _atmospheric_property
                : public CloneBase<_atmospheric_property> {
                public:
                    _atmospheric_property() = default;

                    _atmospheric_property( const std::string& key ) :
                        _key { key } {}

                    _atmospheric_property( const std::string& key,
                                           std::vector<axis_t> axes ) :
                        _key { key }, _axes { std::move( axes ) } {}

                    _atmospheric_property( const std::string& key,
                                           std::vector<std::string> axes ) :
                        _key { key } {
                        for (auto ax : axes) {
                            _axes.push_back( NCPA::math::axis( ax ) );
                        }
                    }

                    _atmospheric_property(
                        const std::string& key,
                        std::initializer_list<axis_t> axes ) :
                        _key { key }, _axes { axes } {}

                    _atmospheric_property(
                        const std::string& key,
                        std::initializer_list<std::string> axes ) :
                        _key { key } {
                        std::transform( axes.begin(), axes.end(),
                                        _axes.begin(), []( std::string c ) {
                                            return NCPA::math::axis( c );
                                        } );
                    }

                    virtual ~_atmospheric_property() = default;

                    virtual double _get()                      = 0;
                    virtual double _get( std::vector<double> ) = 0;
                    virtual double _get_first_derivative( axis_t ax1,
                                                          std::vector<double> )
                        = 0;
                    virtual double _get_second_derivative(
                        axis_t ax1, axis_t ax2, std::vector<double> ) = 0;

                    std::vector<axis_t>& axes() { return _axes; }

                    const std::vector<axis_t>& axes() const { return _axes; }

                    size_t dimensions() const { return _axes.size(); }

                    template<typename... Coords>
                    double get( Coords... coords ) {
                        std::vector<double> idx { static_cast<double>(
                            indices )... };
                        if (idx.empty()) {
                            return _get();
                        } else if (idx.size() == this->dimensions()) {
                            return _get( idx );
                        } else {
                            throw WrongNumberOfDimensions();
                        }
                    }

                    template<typename... Coords>
                    double getd( Coords... coords ) {
                        std::vector<double> idx { static_cast<double>(
                            indices )... };
                        if (idx.size() == this->dimensions()) {
                            return _get_first_derivative( _axes.front(), idx );
                        } else {
                            throw WrongNumberOfDimensions();
                        }
                    }

                    template<typename... Coords>
                    double getd( axis_t ax1, Coords... coords ) {
                        std::vector<double> idx { static_cast<double>(
                            indices )... };
                        if (idx.size() == this->dimensions()) {
                            return _get_first_derivative( ax1, idx );
                        } else {
                            throw WrongNumberOfDimensions();
                        }
                    }

                    template<typename... Coords>
                    double getdd( Coords... coords ) {
                        std::vector<double> idx { static_cast<double>(
                            indices )... };
                        if (idx.size() == this->dimensions()) {
                            return _get_second_derivative(
                                _axes.front(), _axes.front(), idx );
                        } else {
                            throw WrongNumberOfDimensions();
                        }
                    }

                    template<typename... Coords>
                    double getdd( axis_t ax1, axis_t ax2, Coords... coords ) {
                        std::vector<double> idx { static_cast<double>(
                            indices )... };
                        if (idx.size() == this->dimensions()) {
                            return _get_second_derivative( ax1, ax2, idx );
                        } else {
                            throw WrongNumberOfDimensions();
                        }
                    }

                    NCPA::arrays::SmartVector<double>& values() {
                        return _values;
                    }

                    const NCPA::arrays::SmartVector<double>& values() const {
                        return _values;
                    }

                protected:
                    std::string _key;
                    std::vector<axis_t> _axes;
                    NCPA::arrays::SmartVector<double> _values;

                    size_t _ax2dim( axis_t ax ) {
                        for (size_t dim = 0; dim < this->dimensions(); ++dim) {
                            if (_axes.at( dim ) == ax) {
                                return dim;
                            }
                        }
                        throw InvalidAxis();
                    }

                    size_t _ax2dim( std::string ax ) {
                        return _ax2dim( NCPA::math::axis( ax ) );
                    }
            };

            class _1d_atmospheric_property
                : public Cloneable<_atmospheric_property,
                                   _1d_atmospheric_property> {
                    typedef Cloneable<_atmospheric_property,
                                      _1d_atmospheric_property>
                        parent_t;

                public:
                    _1d_atmospheric_property()          = default;
                    virtual ~_1d_atmospheric_property() = default;

                    template<typename AxesType>
                    _1d_atmospheric_property(
                        const std::string& key, AxesType axes,
                        interpolator_1d_type_t interptype
                        = NCPA_ATMOSPHERE_DEFAULT_1D_INTERPOLATOR ) :
                        parent_t( key, axes ) {
                        if (axes.size() != 1) {
                            throw WrongNumberOfDimensions();
                        }
                        _interp = InterpolatorFactory<double, double>::build(
                            interptype );
                    }

                    virtual double _get() override {
                        throw WrongNumberOfDimensions(
                            "Not a scalar property" );
                    }

                    virtual double _get( std::vector<double> ) = 0;
                    virtual double _get_first_derivative( axis_t ax1,
                                                          std::vector<double> )
                        = 0;
                    virtual double _get_second_derivative(
                        axis_t ax1, axis_t ax2, std::vector<double> ) = 0;

                protected:
                    Interpolator1D<double, double> _interp;
            };
        }  // namespace details


    }  // namespace atmos
}  // namespace NCPA
