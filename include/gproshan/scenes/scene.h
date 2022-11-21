#ifndef SCENE_H
#define SCENE_H

#include <gproshan/mesh/che.h>

#include <vector>
#include <string>
#include <unordered_map>


// geometry processing and shape analysis framework
namespace gproshan {


class scene: public che
{
	public:
		struct texture
		{
			vec3 * data = nullptr;
			size_t rows = 0;
			size_t cols = 0;
		};

		struct material
		{
			vec3 Ka = {0.2, 0.2, 0.2};
			vec3 Kd = {0.8, 0.8, 0.8};
			vec3 Ks = {1, 1, 1};
			real_t d = 1;	// Tr = 0, opposite
			real_t Ns = 0;
			index_t illum = 1;
			index_t map_Ka = NIL;
			index_t map_Kd = NIL;
		};

		struct object
		{
			index_t begin = 0;
			index_t end = 0;
			index_t idm = NIL;
		};

	protected:
		std::unordered_map<std::string, index_t> material_id;
		std::vector<std::string> material_name;
		std::vector<material> materials;

		std::unordered_map<std::string, index_t> texture_id;
		std::vector<std::string> texture_name;
		std::vector<texture> textures;

		vertex * texcoords = nullptr;

	public:
		scene(const std::string & file);
		~scene();
		void read_file(const std::string & file);
		bool load_obj(const std::string & file);
		bool load_mtl(const std::string & file);
		bool load_texture(const std::string & file);
};


} // namespace gproshan

#endif // SCENE_H

