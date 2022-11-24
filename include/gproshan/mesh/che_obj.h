#ifndef CHE_OBJ_H
#define CHE_OBJ_H

#include <gproshan/mesh/che.h>

#include <unordered_set>


// geometry processing and shape analysis framework
namespace gproshan {


class che_obj : public che
{
	public:
		che_obj(const std::string & file);

		static void write_file(const che * mesh, const std::string & file, const bool & color = false, const bool & pointcloud = false);

	private:
		void read_file(const std::string & file);

	public:
		struct parser
		{
			std::vector<vertex> vertices;
			std::vector<vec2> vtexcoords;
			std::vector<vertex> vnormals;
			std::vector<rgb_t> vcolors;
			std::vector<uvec3> trigs;
			std::vector<std::pair<std::string, index_t> > objects;
			std::unordered_set<std::string> mtllibs;

			parser(const std::string & file);
		};
};


} // namespace gproshan

#endif // CHE_OBJ_H

