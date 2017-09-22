#include "dictionary.h"

dictionary::dictionary(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M):
					mesh(_mesh), phi_basis(_phi_basis), m(_m), M(_M)
{

}

dictionary::~dictionary()
{
}

