#ifndef CHE_IMG_H
#define CHE_IMG_H

#include "che.h"

class che_img : public che
{
	public:
		che_img(const size_t & n_v = 0, const size_t & n_f = 0);
		che_img(const real_t * img, const size_t & n_rows, const size_t & n_cols);
		che_img(const string & file);
		virtual ~che_img();
		void write_file(const string & file) const;

	private:
		void read_file(const string & file);
};

#endif // CHE_IMG_H
