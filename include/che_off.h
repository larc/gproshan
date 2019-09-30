#ifndef CHE_OFF_H
#define CHE_OFF_H

#include "che.h"

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

#endif // CHE_OFF_H

