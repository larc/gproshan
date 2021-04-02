#ifndef CHE_PTX_H
#define CHE_PTX_H

#include "mesh/che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_ptx : public che
{
	public:
		che_ptx(const std::string & file);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_PTX_H

