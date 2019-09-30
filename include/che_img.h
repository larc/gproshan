#ifndef CHE_IMG_H
#define CHE_IMG_H

#include "che.h"

class che_img : public che
{
	public:
		che_img(const std::string & file);
		che_img(const che_img & mesh);
		virtual ~che_img();
		void write_file(const std::string & file) const;

	private:
		void read_file(const std::string & file);
};

#endif // CHE_IMG_H
