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
			unsigned char * data = nullptr;
			size_t width = 0;
			size_t height = 0;
			size_t spectrum = 0;
		};

		struct material
		{
			vec3 Ka = 1;
			vec3 Kd = 0.8;
			vec3 Ks = 0.2;
			float d = 1;	// Tr = 0, opposite
			float Ns = 10;
			float Ni = 0;
			int illum = 1;
			int map_Ka = -1;
			int map_Kd = -1;
			int map_Ks = -1;
			int map_d = -1;
			int map_bump = -1;
		};

		struct object
		{
			index_t begin = 0;
			index_t material_id = NIL;
		};

	public:
		std::unordered_map<std::string, index_t> material_id;
		std::vector<std::string> material_name;
		std::vector<material> materials;

		std::vector<std::string> texture_name;
		std::vector<texture> textures;

		std::vector<object> objects;

		index_t * trig_mat = nullptr;
		vec2 * texcoords = nullptr;
		bool load_scene = true;

	public:
		scene(const std::string & file);
		~scene();
		bool is_scene() const;
		bool is_pointcloud() const;
		void read_file(const std::string & file);
		bool load_obj(const std::string & file);
		bool load_mtl(const std::string & file);
		bool load_texture(const std::string & file);
};

struct scene_data
{
	scene::material * materials = nullptr;
	scene::texture * textures = nullptr;
	index_t * trig_mat = nullptr;
	vec2 * texcoords = nullptr;
};


} // namespace gproshan

#endif // SCENE_H
