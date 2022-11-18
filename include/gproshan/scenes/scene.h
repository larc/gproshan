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
		};

	protected:
		std::unordered_map<std::string, index_t> material_id;
		std::vector<std::string> material_name;
		std::vector<material> materials;

		std::unordered_map<std::string, index_t> texture_id;
		std::vector<std::string> texture_name;
		std::vector<material> textures;

	public:
		bool read_mtl(const std::string & file);
};


} // namespace gproshan

#endif // SCENE_H

