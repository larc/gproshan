#ifndef CHE_OBJ_H
#define CHE_OBJ_H

#include "che.h"

class che_obj : public che
{
	public:
		che_obj(const string & file);
		virtual ~che_obj() = default;

		static void write_file(const che * mesh, const string & file);

	private:
		void read_file(const string & file);
};

#endif // CHE_OBJ_H

