#ifndef CHE_OBJ_H
#define CHE_OBJ_H

#include "mesh/che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_obj : public che
{
	public:
		che_obj(const std::string & file);

		static void write_file(const che * mesh, const std::string & file, const bool & color = false, const bool & pointcloud = false);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_OBJ_H

