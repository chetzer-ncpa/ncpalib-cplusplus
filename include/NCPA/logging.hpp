#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace NCPA {
    namespace logging {
        enum class log_level_t {
            NONE = 0,
            CRITICAL,
            ERROR,
            WARNING,
            INFO,
            DEBUG
        };

        namespace details {
            class _logger;
            class _file_logger;
        }  // namespace details
        class Logger;
    }  // namespace logging
}  // namespace NCPA

// template<typename T>
// NCPA::logging::Logger& operator<<(NCPA::logging::Logger& log, T const &
// msg);

namespace NCPA {
    namespace logging {
        namespace details {
            constexpr char _default_logger_name[]    = "default";
            constexpr log_level_t _default_log_level = log_level_t::WARNING;

            class _logger {
                public:
                    _logger( log_level_t level, std::ostream& stream ) :
                        _level { level },
                        _current_level { log_level_t::INFO },
                        _stream { &stream } {}

                    _logger( log_level_t level ) :
                        _logger( level, std::cout ) {}

                    _logger() : _logger( log_level_t::INFO, std::cout ) {}

                    virtual ~_logger() {}

                    virtual void flush() { stream().flush(); }

                    template<typename T>
                    _logger& log( log_level_t level, T message ) {
                        if ( level <= _level ) {
                            stream() << message;
                        }
                        return *this;
                    }

                    virtual _logger& set_level( log_level_t newlevel ) {
                        _level = newlevel;
                        return *this;
                    }

                    virtual _logger& set_stream( std::ostream& newstream ) {
                        _stream = &newstream;
                        return *this;
                    }

                    virtual std::ostream& stream() { return *_stream; }

                    template<typename T>
                    _logger operator<<( T msg ) {
                        // return log( _current_level, msg );
                        if ( _current_level <= _level ) {
                            stream() << msg;
                        }
                        return *this;
                    }

                    _logger operator<<(
                        NCPA::logging::log_level_t new_level ) {
                        _current_level = new_level;
                        return *this;
                    }

                private:
                    std::ostream *_stream;
                    NCPA::logging::log_level_t _level, _current_level;
            };

            class _file_logger : public _logger {
                public:
                    _file_logger( const std::string& filename,
                                  log_level_t level = _default_log_level ) :
                        _logger( level, _fstream ), _filename { filename } {
                        _open_stream();
                    }

                    virtual ~_file_logger() { _close_stream(); }

                    virtual void flush() override {
                        _close_stream();
                        _open_stream();
                    }

                private:
                    std::string _filename;
                    std::ofstream _fstream;

                    void _close_stream() { _fstream.close(); }

                    void _open_stream() {
                        _fstream.open( _filename, std::ios_base::out
                                                      | std::ios_base::app );
                        if (!_fstream.is_open()) {
                            throw std::runtime_error("Failed to open file " + _filename);
                        }
                    }
            };
        }  // namespace details

        class Logger {
            public:
                static Logger& logger() {
                    static Logger instance;
                    return instance;
                }

                void init( log_level_t loglevel
                           = details::_default_log_level ) {
                    if ( _loggers.count( _default ) == 0 ) {
                        add_logger( _default, loglevel );
                    }
                }

                Logger& add_logger( const std::string& loggername,
                                    log_level_t loglevel
                                    = details::_default_log_level ) {
                    assert_not_set( loggername );
                    _loggers[ loggername ] = details::_logger( loglevel );
                    return *this;
                }

                Logger& set_level( const std::string& loggername,
                                   log_level_t newlevel ) {
                    assert_set( loggername );
                    _loggers.at( loggername ).set_level( newlevel );
                    return *this;
                }

                Logger& set_level( log_level_t newlevel ) {
                    init();
                    return set_level( _default, newlevel );
                }

                Logger& set_output( const std::string& loggername,
                                    std::ostream& stream ) {
                    assert_set( loggername );
                    _loggers.at( loggername ).set_stream( stream );
                    return *this;
                }

                Logger& set_output( std::ostream& stream ) {
                    init();
                    return set_output( _default, stream );
                }

                Logger& set_output( const std::string& loggername,
                                    const std::string& filename ) {
                    _loggers[ loggername ] = details::_file_logger( filename );
                    return *this;
                }

                Logger& set_output( const std::string& filename ) {
                    return set_output( _default, filename );
                }

                template<typename T>
                friend Logger& operator<<( Logger& log, T const& msg ) {
                    return log.pass_on( msg );
                }

                friend Logger& operator<<( Logger& log,
                                       std::ostream& ( *f )(std::ostream&)) {
                    std::ostringstream s;
                    s << f;
                    return log << s.str();
                }

                template<typename T>
                Logger& pass_on( T msg ) {
                    init();
                    for ( auto it = _loggers.begin(); it != _loggers.end();
                          ++it ) {
                        it->second << msg;
                    }
                    return *this;
                }

                // template<typename T>
                // Logger& operator<<( T msg ) {
                //     init();
                //     for ( auto it = _loggers.begin(); it != _loggers.end();
                //           ++it ) {
                //         it->second << msg;
                //     }
                //     return *this;
                // }

                void assert_set( const std::string& loggername ) {
                    if ( _loggers.count( loggername ) == 0 ) {
                        throw std::range_error( "No logger named " + loggername
                                                + " has been set" );
                    }
                }

                void assert_not_set( const std::string& loggername ) {
                    if ( _loggers.count( loggername ) == 1 ) {
                        throw std::range_error( "A logger named " + loggername
                                                + " has already been set" );
                    }
                }

                Logger& flush() {
                    for ( auto it = _loggers.begin(); it != _loggers.end();
                          ++it ) {
                        it->second.flush();
                    }
                    return *this;
                }


            private:
                Logger() : _default { details::_default_logger_name } {
                    init();
                }

                std::unordered_map<std::string, details::_logger> _loggers;
                std::string _default;
        };

        Logger logger = Logger::logger();
    }  // namespace logging
}  // namespace NCPA
