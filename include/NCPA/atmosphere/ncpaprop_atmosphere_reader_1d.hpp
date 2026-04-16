#pragma once

#ifndef FILE_SEPARATOR
#  ifdef _WIN32
#    define FILE_SEPARATOR '\\'
#  else
#    define FILE_SEPARATOR '/'
#  endif
#endif

#include "NCPA/atmosphere/abstract_atmosphere_reader_1d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_1d&,
                  NCPA::atmos::ncpaprop_atmosphere_reader_1d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class ncpaprop_atmosphere_reader_1d
            : public _abstract_atmosphere_reader_1d {
            public:
                ncpaprop_atmosphere_reader_1d() :
                    _abstract_atmosphere_reader_1d() {}

                ncpaprop_atmosphere_reader_1d(
                    const ncpaprop_atmosphere_reader_1d& other ) :
                    ncpaprop_atmosphere_reader_1d() {}

                ncpaprop_atmosphere_reader_1d(
                    ncpaprop_atmosphere_reader_1d&& source ) noexcept :
                    ncpaprop_atmosphere_reader_1d() {
                    ::swap( *this, source );
                }

                virtual ~ncpaprop_atmosphere_reader_1d() {}

                ncpaprop_atmosphere_reader_1d& operator=(
                    ncpaprop_atmosphere_reader_1d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    ncpaprop_atmosphere_reader_1d& a,
                    ncpaprop_atmosphere_reader_1d& b ) noexcept;

                virtual std::unique_ptr<_abstract_atmosphere_reader_1d> clone()
                    const override {
                    return std::unique_ptr<_abstract_atmosphere_reader_1d>(
                        new ncpaprop_atmosphere_reader_1d( *this ) );
                }

                virtual Atmosphere1D read(
                    const std::string& filename ) override {
                    std::ifstream fs( filename, std::ios_base::in );
                    // std::cout << "Reading from single filename: " <<
                    // filename << std::endl;
                    return this->read( fs );
                }

                virtual Atmosphere1D read(
                    const std::vector<std::string>& filenames ) override {
                    // std::cout << "Reading from vector of filenames" <<
                    // std::endl;
                    std::ifstream fs1, fs2;
                    std::vector<std::istream *> strs;
                    switch (filenames.size()) {
                        case 0:
                            throw std::logic_error(
                                "ncpaprop_atmosphere_reader_1d.read(): No "
                                "filenames provided!" );
                            break;
                        case 1:
                            fs1.open( filenames[ 0 ], std::ios_base::in );
                            return this->read( fs1 );
                            break;
                        default:
                            fs1.open( filenames[ 0 ], std::ios_base::in );
                            fs2.open( filenames[ 1 ], std::ios_base::in );
                            strs = { &fs1, &fs2 };
                            return this->read( strs );
                            break;
                    }
                }

                virtual Atmosphere1D read( std::istream& in1 ) override {
                    // std::cout << "Reading from single stream" << std::endl;
                    return _read_from_streams( &in1 );
                }

                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams ) override {
                    // std::cout << "Reading from vector of streams" <<
                    // std::endl;
                    switch (streams.size()) {
                        case 0:
                            throw std::logic_error(
                                "ncpaprop_atmosphere_reader_1d.read(): No "
                                "streams provided!" );
                            break;
                        case 1:
                            return _read_from_streams( streams[ 0 ] );
                            break;
                        default:
                            return _read_from_streams( streams[ 0 ],
                                                       streams[ 1 ] );
                            break;
                    }
                }

            protected:
                Atmosphere1D _read_from_streams( std::istream *in,
                                                 std::istream *header_in
                                                 = nullptr ) {
                    std::vector<std::string> headerlines;
                    if (header_in != nullptr) {
                        headerlines = _read_header_lines( *header_in );
                    }

                    if (!in->good()) {
                        throw std::runtime_error(
                            "ncpaprop_atmosphere_reader_1d._read_from_streams("
                            "): Input stream not in good state" );
                    }

                    std::string line;
                    std::vector<std::string> atmlines, scalarlines;
                    std::ostringstream oss;  // for exceptions
                    size_t i;                // repeated index variable
                    tuple_atmosphere_1d atm;

                    std::getline( *in, line );
                    while (in->good()) {
                        // lines will either be comments (# ), field
                        // descriptions (#% ), or field contents
                        line = NCPA::strings::deblank( line );
                        if (line[ 0 ] == '#') {
                            // check second character
                            if (header_in == nullptr && line.size() > 1
                                && line[ 1 ] == '%') {
                                headerlines.push_back( line.substr( 2 ) );
                            }
                            // otherwise it's a regular comment and can be
                            // ignored

                        } else if (line.size() == 0) {
                            // skip empty lines
                        } else {
                            atmlines.push_back( line );
                        }

                        std::getline( *in, line );
                    }

                    // parse them out
                    size_t nfields = headerlines.size();
                    if (nfields == 0) {
                        throw std::runtime_error(
                            "ncpaprop_atmosphere_reader_1d._read_from_streams("
                            "): No descriptive fields found." );
                    }

                    // hold contents
                    std::vector<std::string> keys, fields;
                    std::vector<size_t> column_numbers;
                    std::vector<double> values;
                    std::vector<units_ptr_t> units;

                    for (auto it = headerlines.cbegin();
                         it != headerlines.cend(); ++it) {
                        fields = NCPA::strings::split(
                            NCPA::strings::deblank( *it ), " ," );
                        if (fields.size() != 3 && fields.size() != 4) {
                            oss << "ncpaprop_atmosphere_reader_1d._read_from_"
                                   "streams(): Error parsing descriptive line:"
                                << std::endl
                                << line << std::endl
                                << "Must be formatted as:" << std::endl
                                << "column,key,units[,value]" << std::endl
                                << "Use column=0 and specify value for scalar "
                                   "quantities."
                                << std::endl;
                            throw std::invalid_argument( oss.str() );
                        }

                        // process fields
                        // field line needs to be parseable as an integer
                        size_t col;
                        try {
                            col = (size_t)std::stoi( fields[ 0 ] );
                        } catch (std::invalid_argument& e) {
                            oss << "ncpaprop_atmosphere_reader_1d._read_from_"
                                   "streams(): Error parsing descriptive "
                                   "line:"
                                << std::endl
                                << line << std::endl
                                << "First field not parseable as an integer"
                                << std::endl;
                            throw std::invalid_argument( oss.str() );
                        }

                        double tempval = 0.0;
                        units_ptr_t tempunits
                            = NCPA::units::Units::from_string(
                                NCPA::strings::deblank( fields[ 2 ] ) );
                        if (tempunits == nullptr) {
                            throw std::invalid_argument(
                                "ncpaprop_atmosphere_reader_1d._read_from_"
                                "streams(): Unrecognized units string: "
                                + fields[ 2 ] );
                        }
                        if (fields.size() == 4) {
                            try {
                                tempval = std::stod(
                                    NCPA::strings::deblank( fields[ 3 ] ) );
                            } catch (std::out_of_range& e) {
                                oss << "ncpaprop_atmosphere_reader_1d._read_"
                                       "from_"
                                       "streams(): Error converting data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Value " << fields[ 3 ]
                                    << " is out of range for double "
                                       "precision.";
                                throw std::runtime_error( oss.str() );
                            }
                        }

                        // add to header vectors
                        column_numbers.push_back( col );
                        keys.push_back(
                            NCPA::strings::deblank( fields[ 1 ] ) );
                        units.push_back( tempunits );
                        values.push_back( tempval );  // this will be ignored
                                                      // for vector quantities
                    }

                    // check for uniqueness in keys
                    std::vector<std::string> tempkeys( keys );
                    std::sort( tempkeys.begin(), tempkeys.end() );
                    std::vector<std::string>::iterator uit
                        = std::unique( tempkeys.begin(), tempkeys.end() );
                    if (uit != tempkeys.end()) {
                        throw std::invalid_argument(
                            "ncpaprop_atmosphere_reader_1d._read_from_"
                            "streams(): key values are not unique" );
                    }


                    int nlines           = atmlines.size();
                    size_t ncols         = 0;
                    units_ptr_t depunits = nullptr;
                    for (i = 0; i < column_numbers.size(); i++) {
                        ncols = column_numbers[ i ] > ncols
                                  ? column_numbers[ i ]
                                  : ncols;
                        if (column_numbers[ i ] == 1) {
                            depunits = units[ i ];
                        }
                    }

                    std::vector<double *> columns( ncols );
                    for (i = 0; i < ncols; i++) {
                        double *col  = new double[ nlines ];
                        columns[ i ] = col;
                    }


                    // step through the data lines
                    size_t row = 0;
                    std::vector<std::string>::const_iterator it;
                    for (it = atmlines.cbegin(); it != atmlines.cend(); ++it) {
                        fields.clear();
                        fields = NCPA::strings::split(
                            NCPA::strings::deblank( *it ), " \t," );
                        if (fields.size() != ncols) {
                            oss << "ncpaprop_atmosphere_reader_1d._read_from_"
                                   "streams(): Error parsing data line:"
                                << std::endl
                                << *it << std::endl
                                << ncols << " columns expected, "
                                << fields.size() << " columns found.";
                            for (i = 0; i < ncols; i++) {
                                delete[] columns[ i ];
                            }
                            throw std::invalid_argument( oss.str() );
                        }
                        for (i = 0; i < ncols; i++) {
                            try {
                                double *thiscol = columns[ i ];
                                thiscol[ row ]  = std::stod( fields[ i ] );
                            } catch (std::out_of_range& e) {
                                oss << "ncpaprop_atmosphere_reader_1d._read_"
                                       "from_"
                                       "streams(): Error converting data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Value " << fields[ i ]
                                    << " is out of range for double "
                                       "precision.";
                                for (i = 0; i < ncols; i++) {
                                    delete[] columns[ i ];
                                }
                                throw std::runtime_error( oss.str() );
                            } catch (std::invalid_argument& e) {
                                oss << "ncpaprop_atmosphere_reader_1d._read_"
                                       "from_"
                                       "streams(): Error parsing data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Can't parse field " << fields[ i ]
                                    << " as a double";
                                for (i = 0; i < ncols; i++) {
                                    delete[] columns[ i ];
                                }
                                throw std::runtime_error( oss.str() );
                            }
                        }
                        row++;
                    }

                    atm.set_axis(
                        vector_u_t( nlines, columns[ 0 ], depunits ) );

                    for (i = 0; i < keys.size(); i++) {
                        if (column_numbers[ i ] == 0) {
                            atm.add_property(
                                keys[ i ],
                                scalar_u_t( values[ i ], *units[ i ] ) );
                        } else if (column_numbers[ i ] > 1) {
                            // AtmosphericProperty1D prop
                            //     = AtmosphereFactory::build(
                            //         atmospheric_property_1d_t::TUPLE );
                            tuple_atmospheric_property_1d prop;
                            prop.set(
                                atm.get_axis_vector(),
                                vector_u_t( nlines,
                                            columns[ column_numbers[ i ] - 1 ],
                                            units[ i ] ) );
                            atm.add_property( keys[ i ], prop );
                        }
                    }

                    for (i = 0; i < ncols; i++) {
                        delete[] columns[ i ];
                    }

                    // std::cout << "Read successful" << std::endl;

                    return Atmosphere1D(
                        std::unique_ptr<abstract_atmosphere_1d>(
                            new tuple_atmosphere_1d( atm ) ) );
                }

                std::vector<std::string> _read_header_lines(
                    std::istream& in ) {
                    std::string line;
                    std::vector<std::string> headerlines;
                    if (!in.good()) {
                        throw std::runtime_error(
                            "ncpaprop_atmosphere_reader_1d._read_header_lines("
                            "): Header input stream not in good "
                            "state" );
                    }

                    std::getline( in, line );
                    while (in.good()) {
                        // lines will either be comments (# ), field
                        // descriptions (#% ), or field contents
                        line = NCPA::strings::deblank( line );
                        if (line[ 0 ] == '#') {
                            // check second character
                            if (line.size() > 1 && line[ 1 ] == '%') {
                                headerlines.push_back( line.substr( 2 ) );
                            }  // otherwise it's a regular comment and can be
                               // ignored
                        }

                        std::getline( in, line );
                    }
                    return headerlines;
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_1d& a,
                  NCPA::atmos::ncpaprop_atmosphere_reader_1d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_1d&>( a ),
            dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_1d&>( b ) );
}
