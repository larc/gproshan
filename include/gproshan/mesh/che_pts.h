#ifndef CHE_PTS_H
#define CHE_PTS_H

#include "mesh/che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_pts : public che
{
	public:
		che_pts(const std::string & file);

		static void write_file(const che * mesh, const std::string & file);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_PTS_H

