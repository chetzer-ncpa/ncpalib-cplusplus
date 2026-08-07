#pragma once

#include <fstream>
#include <iostream>
#include <streambuf>
#include <vector>

namespace NCPA {
    namespace stream {
        class MultiStreamBuffer : public std::streambuf {
            public:
                MultiStreamBuffer() {}

                void add( std::ostream& os ) {
                    _buffers.push_back( os.rdbuf() );
                }

            protected:
                virtual int overflow( int c ) override {
                    if (c == EOF) {
                        return EOF;
                    }
                    for (auto *buf : _buffers) {
                        if (buf->sputc( c ) == EOF) {
                            return EOF;
                        }
                    }
                    return c;
                }

                virtual int sync() override {
                    for (auto *buf : _buffers) {
                        if (buf->pubsync() == -1) {
                            return -1;
                        }
                    }
                    return 0;
                }

            private:
                std::vector<std::streambuf *> _buffers;
        };

        class MultiOStream : public std::ostream {
            public:
                MultiOStream() : std::ostream( &_buffer ) {}

                void add( std::ostream& os ) { _buffer.add( os ); }

            private:
                MultiStreamBuffer _buffer;
        };

    }  // namespace stream
}  // namespace NCPA
