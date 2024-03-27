#ifndef CHE_IMG_H
#define CHE_IMG_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class che_img : public che
{
	public:
		che_img(const std::string & file);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_IMG_H

