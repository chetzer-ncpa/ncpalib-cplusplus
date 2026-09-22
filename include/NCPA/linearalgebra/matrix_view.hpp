#include "NCPA/linearalgebra/abstract_matrix.hpp"
#include "NCPA/linearalgebra/abstract_vector.hpp"
#include "NCPA/linearalgebra/declarations.hpp"

#include <memory>

namespace NCPA {
    namespace linear {
        template<typename ELEMENTTYPE>
        class matrix_view : public abstract_matrix<ELEMENTTYPE> {
            public:
                matrix_view( const abstract_matrix<ELEMENTTYPE>& m ) {
                    _ptr = m.clone();
                }

                matrix_view(
                    const std::unique_ptr<abstract_matrix<ELEMENTTYPE>>& m ) {
                    _ptr = m->clone();
                }

                virtual ~matrix_view() {}

                friend void swap( matrix_view<ELEMENTTYPE>& a,
                                  matrix_view<ELEMENTTYPE>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<abstract_matrix<ELEMENTTYPE>&>( a ),
                          static_cast<abstract_matrix<ELEMENTTYPE>&>( b ) );
                }

                abstract_matrix<ELEMENTTYPE> *internal() { return _ptr.get(); }

                const abstract_matrix<ELEMENTTYPE> *internal() const {
                    return _ptr.get();
                }

            private:
                std::unique_ptr<abstract_matrix<ELEMENTTYPE>> _ptr;
        };
    }  // namespace linear
}  // namespace NCPA
