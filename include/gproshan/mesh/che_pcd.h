#ifndef CHE_PCD_H
#define CHE_PCD_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class che_pcd : public che
{
	public:
		che_pcd(const std::string & file);

		static void write_file(const che * mesh, const std::string & file, const bool & color = false);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_PCD_H

