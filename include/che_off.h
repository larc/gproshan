#ifndef CHE_OFF_H
#define CHE_OFF_H

#include "che.h"

class che_off : public che
{
	public:
		che_off(const string & file);
		virtual ~che_off();

		static void write_file(const che * mesh, const string & file);

	private:
		void read_file(const string & file);
};

#endif // CHE_OFF_H

