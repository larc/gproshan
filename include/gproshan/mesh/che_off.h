#ifndef CHE_OFF_H
#define CHE_OFF_H

#include "mesh/che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_off : public che
{
	public:
		che_off(const std::string & file);

		enum type { OFF, NOFF, COFF, NCOFF };
		static void write_file(const che * mesh, const std::string & file, const che_off::type & off = OFF, const bool & pointcloud = false);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_OFF_H

