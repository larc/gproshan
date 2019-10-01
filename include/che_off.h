#ifndef CHE_OFF_H
#define CHE_OFF_H

#include "che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_off : public che
{
	public:
		che_off(const std::string & file);
		che_off(const che_off & mesh);
		virtual ~che_off();

		static void write_file(const che * mesh, const std::string & file);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_OFF_H

