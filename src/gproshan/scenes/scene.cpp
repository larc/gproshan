#include "gproshan/scenes/scene.h"

#include "gproshan/mesh/che_obj.h"

#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


scene::scene(const std::string & file)
{
	init(file);
}

scene::~scene()
{
	for(texture & tex: textures)
		delete tex.data;
	delete [] texcoords;
}

bool scene::is_scene() const
{
	return true;
}

bool scene::is_pointcloud() const
{
	return false;
}

void scene::read_file(const std::string & file)
{
	load_obj(file);
}

bool scene::load_obj(const std::string & file)
{
	const che_obj::parser p(file);

	const std::string path = file.substr(0, file.rfind('/') + 1);
	for(auto & m: p.mtllibs)
		if(!load_mtl(path + m))
			return false;

	alloc(p.trigs.size(), 0);
	texcoords = new vec2[n_vertices];

	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices; ++i)
	{
		const index_t & v = p.trigs[i].x();
		const index_t & t = p.trigs[i].y();
		GT[i] = p.vertices[v];
		VC[i] = p.vcolors[v];
		texcoords[i] = t != NIL ? p.vtexcoords[t] : vec2{-1, -1};
	}

	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices; ++i)
	{
		const index_t & trig = 3 * (i / 3);
		const index_t & n = p.trigs[i].z();
		VN[i] = n != NIL ? p.vnormals[n] : normalize((GT[trig + 1] - GT[trig]) * (GT[trig + 2] - GT[trig]));
	}

	for(auto & obj: p.objects)
		objects.push_back({obj.second, material_id[obj.first]});

	gproshan_log_var(objects.size());

	return true;
}

bool scene::load_mtl(const std::string & file)
{
	gproshan_error_var(file);

	FILE * fp = fopen(file.c_str(), "r");
	if(!fp) return false;

	std::unordered_map<std::string, index_t> texture_id;

	char line[256], str[64];
	while(fgets(line, sizeof(line), fp))
	{
		sscanf(line, "%s", str);
		switch(str[0])
		{
			case 'n':	// newmtl
			{
				sscanf(line, "%*s %s", str);
				material_id[str] = materials.size();
				material_name.push_back(str);
				materials.emplace_back();
				break;
			}
			case 'K':	// Ka, Kd, Ks
			{
				vec3 & rgb = str[1] == 'a' ? materials.back().Ka
							: str[1] == 'd' ? materials.back().Kd
							: materials.back().Ks;
				sscanf(line, "%*s %f %f %f", &rgb.x(), &rgb.y(), &rgb.z());
				break;
			}
			case 'd':	// d
			{
				real_t & d = materials.back().d;
				sscanf(line, "%*s %f", &d);
				break;
			}
			case 'T':	// Tr
			{
				real_t & d = materials.back().d;
				sscanf(line, "%*s %f", &d);
				d = 1 - d;
				break;
			}
			case 'N':	// Ns
			{
				real_t & N = str[1] == 's' ? materials.back().Ns : materials.back().Ni;
				sscanf(line, "%*s %f", &N);
				break;
			}
			case 'i':	// illum
			{
				index_t & illum = materials.back().illum;
				sscanf(line, "%*s %u", &illum);
				break;
			}
			case 'm':	// map_Ka, map_kd
			{
				index_t & m = str[5] == 'a' ? materials.back().map_Ka : materials.back().map_Kd;
				sscanf(line, "%*s %s", str);
				if(str[0] == '-') continue;		// ignoring map textures with options
				if(texture_id.find(str) == texture_id.end())
				{
					texture_id[str] = textures.size();
					texture_name.push_back(str);
				}
				m = texture_id[str];
				break;
			}
		}
	}

	fclose(fp);

	const std::string path = file.substr(0, file.rfind('/') + 1);
	for(auto & tex: texture_name)
	{
		for(char & c: tex)
			if(c == '\\') c = '/';

		if(!load_texture(path + tex))
			return false;
	}

	for(index_t i = 0; i < materials.size(); ++i)
	{
		const material & m = materials[i];
		gproshan_log_var(material_name[i]);
		gproshan_log_var(m.Ka);
		gproshan_log_var(m.Kd);
		gproshan_log_var(m.Ks);
		gproshan_log_var(m.d);
		gproshan_log_var(m.Ns);
		gproshan_log_var(m.illum);
		if(m.map_Ka != NIL)	gproshan_log_var(texture_name[m.map_Ka]);
		if(m.map_Kd != NIL)	gproshan_log_var(texture_name[m.map_Kd]);
	}

	gproshan_log_var(materials.size());
	gproshan_log_var(textures.size());

	return true;
}

bool scene::load_texture(const std::string & file)
{
	CImg<float> img(file.c_str());

	textures.emplace_back();
	texture & tex = textures.back();
	tex.rows = img.height();
	tex.cols = img.width();
	tex.data = new vec3[tex.rows * tex.cols];
	memcpy((float *) tex.data, img.data(), sizeof(vec3) * tex.rows * tex.cols);

	return true;
}


} // namespace gproshan

