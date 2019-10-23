#ifndef CHE_IMG_H
#define CHE_IMG_H

#include "che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_img : public che
{
	public:
		che_img(const std::string & file);
		che_img(const che_img & mesh);
		virtual ~che_img();
		void write_file(const std::string & file) const;
		//void save_img(const std::string & file, size_t tam) const;
		static void save_img(const che * mesh, const std::string & file, size_t tam);

	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_IMG_H

