#ifndef PATCH_H
#define PATCH_H

#include "include.h"

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

/// Mesh dictionary learning and sparse coding namespace
namespace mdict {

/// 
class patch
{
	private:
		vector<index_t> vertices;
		mat T;
		mat xyz;
};

} // mdict

#endif // PATCH_H

