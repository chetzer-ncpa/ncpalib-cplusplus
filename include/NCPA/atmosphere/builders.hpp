#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_3d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty3D.hpp"
#include "NCPA/atmosphere/readers.hpp"
#include "NCPA/strings.hpp"
#include "NCPA/units.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace NCPA {
    namespace atmos {
        class AtmosphereFactory {
            public:
                static bool can_build( atmospheric_property_1d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_1d_t::TUPLE:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphericProperty1D build(
                    atmospheric_property_1d_t proptype ) {
                    // AtmosphericProperty1D prop;

                    switch (proptype) {
                        case atmospheric_property_1d_t::TUPLE:
                            return AtmosphericProperty1D( _atm_prop_1d_ptr_t(
                                new tuple_atmospheric_property_1d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D atmospheric property type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmospheric_property_2d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_2d_t::STRATIFIED:
                        case atmospheric_property_2d_t::GRID:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphericProperty2D build(
                    atmospheric_property_2d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_2d_t::STRATIFIED:
                            return AtmosphericProperty2D( _atm_prop_2d_ptr_t(
                                new stratified_atmospheric_property_2d() ) );
                            break;
                        case atmospheric_property_2d_t::GRID:
                            return AtmosphericProperty2D( _atm_prop_2d_ptr_t(
                                new grid_atmospheric_property_2d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D atmospheric property type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmospheric_property_3d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_3d_t::STRATIFIED:
                        case atmospheric_property_3d_t::GRID:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphericProperty3D build(
                    atmospheric_property_3d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_3d_t::STRATIFIED:
                            return AtmosphericProperty3D( _atm_prop_3d_ptr_t(
                                new stratified_atmospheric_property_3d() ) );
                            break;
                        case atmospheric_property_3d_t::GRID:
                            return AtmosphericProperty3D( _atm_prop_3d_ptr_t(
                                new grid_atmospheric_property_3d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 3-D atmospheric property type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmosphere_1d_t proptype ) {
                    switch (proptype) {
                        case atmosphere_1d_t::TUPLE:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static Atmosphere1D build( atmosphere_1d_t atmostype ) {
                    Atmosphere1D atmos;
                    switch (atmostype) {
                        case atmosphere_1d_t::TUPLE:
                            return Atmosphere1D(
                                _atm_1d_ptr_t( new tuple_atmosphere_1d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmosphere_2d_t proptype ) {
                    switch (proptype) {
                        case atmosphere_2d_t::GRID:
                        case atmosphere_2d_t::STRATIFIED:
                        case atmosphere_2d_t::PIECEWISE_STRATIFIED:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static Atmosphere2D build( atmosphere_2d_t atmostype ) {
                    Atmosphere1D atmos;
                    switch (atmostype) {
                        case atmosphere_2d_t::STRATIFIED:
                            return Atmosphere2D( _atm_2d_ptr_t(
                                new stratified_atmosphere_2d() ) );
                            break;
                        case atmosphere_2d_t::PIECEWISE_STRATIFIED:
                            return Atmosphere2D( _atm_2d_ptr_t(
                                new piecewise_stratified_atmosphere_2d() ) );
                            break;
                        case atmosphere_2d_t::GRID:
                            return Atmosphere2D(
                                _atm_2d_ptr_t( new grid_atmosphere_2d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmosphere_3d_t proptype ) {
                    switch (proptype) {
                        case atmosphere_3d_t::GRID:
                        case atmosphere_3d_t::STRATIFIED:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static Atmosphere3D build( atmosphere_3d_t atmostype ) {
                    switch (atmostype) {
                        case atmosphere_3d_t::STRATIFIED:
                            return Atmosphere3D( _atm_3d_ptr_t(
                                new stratified_atmosphere_3d() ) );
                            break;
                        case atmosphere_3d_t::GRID:
                            return Atmosphere3D(
                                _atm_3d_ptr_t( new grid_atmosphere_3d() ) );
                            break;

                        default:
                            throw std::range_error(
                                "Requested 3-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( reader_1d_t proptype ) {
                    switch (proptype) {
                        case reader_1d_t::NCPAPROP:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphereReader1D build( reader_1d_t atmostype ) {
                    switch (atmostype) {
                        case reader_1d_t::NCPAPROP:
                            return AtmosphereReader1D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_1d>(
                                    new ncpaprop_atmosphere_reader_1d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( reader_2d_t proptype ) {
                    switch (proptype) {
                        case reader_2d_t::NCPAPROP:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphereReader2D build( reader_2d_t atmostype ) {
                    switch (atmostype) {
                        case reader_2d_t::NCPAPROP:
                            return AtmosphereReader2D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_2d>(
                                    new ncpaprop_atmosphere_reader_2d() ) );
                            break;
                        case reader_2d_t::NCPAPROP_PIECEWISE_STRATIFIED:
                            return AtmosphereReader2D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_2d>(
                                    new ncpaprop_atmosphere_reader_stratified_2d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( reader_3d_t proptype ) {
                    switch (proptype) {
                        case reader_3d_t::NCPAPROP:
                        case reader_3d_t::NCPAPROP_STRATIFIED:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphereReader3D build( reader_3d_t atmostype ) {
                    switch (atmostype) {
                        case reader_3d_t::NCPAPROP:
                            return AtmosphereReader3D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_3d>(
                                    new ncpaprop_atmosphere_reader_3d() ) );
                            break;
                        case reader_3d_t::NCPAPROP_STRATIFIED:
                            return AtmosphereReader3D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_3d>(
                                    new ncpaprop_atmosphere_reader_stratified_3d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 3-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                // static Atmosphere1D build( const std::string& filename,
                //                            const std::string& headerfilename
                //                            = "" ) {
                //     std::string hfile( headerfilename );
                //     if ( headerfilename.size() == 0 ) {
                //         hfile = filename;
                //     }
                //     std::ifstream header_in( hfile );
                //     // read_header_from_stream( header_in );
                //     // header_in.close();

                //     std::ifstream in( filename );
                //     Atmosphere1D atm
                //         = read_values_from_stream( in, header_in );
                //     in.close();
                //     header_in.close();
                //     // headerlines.clear();
                //     return atm;
                // }

                // static std::vector<std::string> read_header_from_stream(
                //     std::istream& in ) {
                //     std::string line;
                //     std::vector<std::string> headerlines;
                //     if ( !in.good() ) {
                //         throw std::runtime_error(
                //             "Atmosphere1D - Header input stream not in good
                //             " "state" );
                //     }

                //     std::getline( in, line );
                //     while ( in.good() ) {
                //         // lines will either be comments (# ), field
                //         // descriptions (#% ), or field contents
                //         line = NCPA::strings::deblank( line );
                //         if ( line[ 0 ] == '#' ) {
                //             // check second character
                //             if ( line.size() > 1 && line[ 1 ] == '%' ) {
                //                 headerlines.push_back( line.substr( 2 ) );
                //             }  // otherwise it's a regular comment and can
                //             be
                //                // ignored
                //         }

                //         std::getline( in, line );
                //     }
                //     return headerlines;
                // }

                // static Atmosphere1D read_values_from_stream(
                //     std::istream& in, std::istream& header_in ) {
                //     std::vector<std::string> headerlines
                //         = read_header_from_stream( header_in );

                //     if ( !in.good() ) {
                //         throw std::runtime_error(
                //             "Atmosphere1D - Input stream not in good state"
                //             );
                //     }

                //     std::string line;
                //     std::vector<std::string> atmlines, scalarlines;
                //     std::ostringstream oss;  // for exceptions
                //     size_t i;                // repeated index variable
                //     tuple_atmosphere_1d atm;

                //     std::getline( in, line );
                //     while ( in.good() ) {
                //         // lines will either be comments (# ), field
                //         // descriptions (#% ), or field contents
                //         line = NCPA::strings::deblank( line );
                //         if ( line[ 0 ] == '#' ) {
                //             // check second character
                //             // if (line.size() > 1 && line[ 1 ] == '%') {
                //             // 	headerlines.push_back( line.substr( 2 ) );
                //             // } // otherwise it's a regular comment and can
                //             be
                //             // ignored

                //         } else if ( line.size() == 0 ) {
                //             // skip empty lines
                //         } else {
                //             atmlines.push_back( line );
                //         }

                //         std::getline( in, line );
                //     }
                //     // in.close();
                //     // cout << "Found " << headerlines.size() << " header
                //     // lines" << endl; cout << "Found " << atmlines.size()
                //     << "
                //     // data lines" << endl;

                //     // parse them out
                //     size_t nfields = headerlines.size();
                //     if ( nfields == 0 ) {
                //         throw std::runtime_error(
                //             "Atmosphere1D - No descriptive fields found." );
                //     }

                //     // hold contents
                //     std::vector<std::string> keys, fields;
                //     std::vector<size_t> column_numbers;
                //     std::vector<double> values;
                //     std::vector<units_ptr_t> units;

                //     for ( auto it = headerlines.cbegin();
                //           it != headerlines.cend(); ++it ) {
                //         fields = NCPA::strings::split(
                //             NCPA::strings::deblank( *it ), " ," );
                //         if ( fields.size() != 3 && fields.size() != 4 ) {
                //             oss << "Atmosphere1D - Error parsing descriptive
                //             "
                //                    "line:"
                //                 << std::endl
                //                 << line << std::endl
                //                 << "Must be formatted as:" << std::endl
                //                 << "column,key,units[,value]" << std::endl
                //                 << "Use column=0 and specify value for
                //                 scalar "
                //                    "quantities."
                //                 << std::endl;
                //             throw std::invalid_argument( oss.str() );
                //         }

                //         // process fields
                //         // field line needs to be parseable as an integer
                //         size_t col;
                //         try {
                //             col = (size_t)std::stoi( fields[ 0 ] );
                //         } catch ( std::invalid_argument& e ) {
                //             oss << "Atmosphere1D - Error parsing descriptive
                //             "
                //                    "line:"
                //                 << std::endl
                //                 << line << std::endl
                //                 << "First field not parseable as an integer"
                //                 << std::endl;
                //             throw std::invalid_argument( oss.str() );
                //         }

                //         double tempval = 0.0;
                //         units_ptr_t tempunits
                //             = NCPA::units::Units::from_string(
                //                 NCPA::strings::deblank( fields[ 2 ] ) );
                //         if ( tempunits == nullptr ) {
                //             throw std::invalid_argument(
                //                 "Unrecognized units string: " + fields[ 2 ]
                //                 );
                //         }
                //         if ( fields.size() == 4 ) {
                //             try {
                //                 tempval = std::stod(
                //                     NCPA::strings::deblank( fields[ 3 ] ) );
                //             } catch ( std::out_of_range& e ) {
                //                 oss << "Error converting data line:"
                //                     << std::endl
                //                     << *it << std::endl
                //                     << "Value " << fields[ 3 ]
                //                     << " is out of range for double "
                //                        "precision.";
                //                 throw std::runtime_error( oss.str() );
                //             }
                //         }

                //         // add to header vectors
                //         column_numbers.push_back( col );
                //         keys.push_back(
                //             NCPA::strings::deblank( fields[ 1 ] ) );
                //         units.push_back( tempunits );
                //         values.push_back( tempval );  // this will be
                //         ignored
                //                                       // for vector
                //                                       quantities
                //     }

                //     // check for uniqueness in keys
                //     std::vector<std::string> tempkeys( keys );
                //     std::sort( tempkeys.begin(), tempkeys.end() );
                //     std::vector<std::string>::iterator uit
                //         = std::unique( tempkeys.begin(), tempkeys.end() );
                //     if ( uit != tempkeys.end() ) {
                //         throw std::invalid_argument(
                //             "Atmosphere1D: key values are not unique" );
                //     }


                //     int nlines           = atmlines.size();
                //     size_t ncols         = 0;
                //     units_ptr_t depunits = nullptr;
                //     for ( i = 0; i < column_numbers.size(); i++ ) {
                //         ncols = column_numbers[ i ] > ncols
                //                   ? column_numbers[ i ]
                //                   : ncols;
                //         if ( column_numbers[ i ] == 1 ) {
                //             depunits = units[ i ];
                //         }
                //     }

                //     std::vector<double *> columns( ncols );
                //     for ( i = 0; i < ncols; i++ ) {
                //         double *col  = new double[ nlines ];
                //         columns[ i ] = col;
                //     }


                //     // step through the data lines
                //     size_t row = 0;
                //     std::vector<std::string>::const_iterator it;
                //     for ( it = atmlines.cbegin(); it != atmlines.cend();
                //           ++it ) {
                //         fields.clear();
                //         fields = NCPA::strings::split(
                //             NCPA::strings::deblank( *it ), " \t," );
                //         if ( fields.size() != ncols ) {
                //             oss << "Error parsing data line:" << std::endl
                //                 << *it << std::endl
                //                 << ncols << " columns expected, "
                //                 << fields.size() << " columns found.";
                //             for ( i = 0; i < ncols; i++ ) {
                //                 delete[] columns[ i ];
                //             }
                //             throw std::invalid_argument( oss.str() );
                //         }
                //         for ( i = 0; i < ncols; i++ ) {
                //             try {
                //                 double *thiscol = columns[ i ];
                //                 thiscol[ row ]  = std::stod( fields[ i ] );
                //             } catch ( std::out_of_range& e ) {
                //                 oss << "Error converting data line:"
                //                     << std::endl
                //                     << *it << std::endl
                //                     << "Value " << fields[ i ]
                //                     << " is out of range for double "
                //                        "precision.";
                //                 for ( i = 0; i < ncols; i++ ) {
                //                     delete[] columns[ i ];
                //                 }
                //                 throw std::runtime_error( oss.str() );
                //             } catch ( std::invalid_argument& e ) {
                //                 oss << "Error parsing data line:" <<
                //                 std::endl
                //                     << *it << std::endl
                //                     << "Can't parse field " << fields[ i ]
                //                     << " as a double";
                //                 for ( i = 0; i < ncols; i++ ) {
                //                     delete[] columns[ i ];
                //                 }
                //                 throw std::runtime_error( oss.str() );
                //             }
                //         }
                //         row++;
                //     }

                //     atm.set_axis(
                //         vector_u_t( nlines, columns[ 0 ], depunits ) );

                //     for ( i = 0; i < keys.size(); i++ ) {
                //         if ( column_numbers[ i ] == 0 ) {
                //             atm.add_property(
                //                 keys[ i ],
                //                 scalar_u_t( values[ i ], *units[ i ] ) );
                //         } else if ( column_numbers[ i ] > 1 ) {
                //             AtmosphericProperty1D prop
                //                 = AtmosphereFactory::build(
                //                     atmospheric_property_1d_t::TUPLE );
                //             prop.set(
                //                 atm.get_axis_vector(),
                //                 vector_u_t( nlines,
                //                             columns[ column_numbers[ i ] - 1
                //                             ], units[ i ] ) );
                //             atm.add_property( keys[ i ], prop );
                //         }
                //     }

                //     for ( i = 0; i < ncols; i++ ) {
                //         delete[] columns[ i ];
                //     }

                //     return Atmosphere1D(
                //         std::unique_ptr<abstract_atmosphere_1d>(
                //             new tuple_atmosphere_1d( atm ) ) );
                // }

            private:
                std::vector<std::string> headerlines;
        };
    }  // namespace atmos
}  // namespace NCPA
