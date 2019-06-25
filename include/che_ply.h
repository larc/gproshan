#ifndef CHE_PLY_H
#define CHE_PLY_H

#include "che.h"

class che_ply : public che
{
	public:
		che_ply(const string & file);
		virtual ~che_ply() = default;

		static void write_file(const che * mesh, const string & file);

	private:
		void read_file(const string & file);
};

#endif // CHE_PLY_H

