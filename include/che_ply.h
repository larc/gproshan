#ifndef CHE_PLY_H
#define CHE_PLY_H

#include "che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_ply : public che
{
	public:
		che_ply(const std::string & file);
		che_ply(const che_ply & mesh);
		virtual ~che_ply() = default;

		static void write_file(const che * mesh, const std::string & file);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_PLY_H

