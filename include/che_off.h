#ifndef CHE_OFF_H
#define CHE_OFF_H

#include "che.h"

class che_off : public che
{
	public:
		che_off(const size_t & n_v = 0, const size_t & n_f = 0);
		che_off(const vertex * vertices, const size_t & n_v, const index_t * faces, const size_t & n_f);
		che_off(const string & file);
		virtual ~che_off();
		void write_file(const string & file) const;

	private:
		void read_file(const string & file);
};

#endif // CHE_OFF_H

