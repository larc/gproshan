#ifndef CHE_PTX_H
#define CHE_PTX_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class che_ptx : public che
{
	public:
		che_ptx(const std::string & file);

		static void write_file(const che * mesh, const std::string & file, const size_t & n_rows, const size_t & n_cols);


	private:
		void read_file(const std::string & file);
};


} // namespace gproshan

#endif // CHE_PTX_H

