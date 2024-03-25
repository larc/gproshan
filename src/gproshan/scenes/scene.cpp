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

	delete [] trig_mat;
	delete [] texcoords;
}

bool scene::is_scene() const
{
	return load_scene && size(objects) > 1;
}

bool scene::is_pointcloud() const
{
	return false;
}

void scene::read_file(const std::string & file)
{
	load_scene = load_obj(file);
}

bool scene::load_obj(const std::string & file)
{
	const che_obj::parser p(file);

	const std::string path = file.substr(0, file.rfind('/') + 1);
	for(auto & m: p.mtllibs)
		if(!load_mtl(path + m))
			return false;

	alloc(size(p.trigs), size(p.trigs));

	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices; ++i)
	{
		const index_t v = p.trigs[i].x();
		GT[i] = p.vertices[v];
		VC[i] = p.vcolors[v];
		VT[i] = i;
	}

	if(size(p.vtexcoords))
	{
		texcoords = new vec2[n_vertices];

		#pragma omp parallel for
		for(index_t i = 0; i < n_vertices; ++i)
		{
			const index_t t = p.trigs[i].y();
			texcoords[i] = t != NIL ? p.vtexcoords[t] : vec2{-1, -1};
		}
	}

	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices; ++i)
	{
		const index_t trig = 3 * (i / 3);
		const index_t n = p.trigs[i].z();
		VN[i] = n != NIL ? p.vnormals[n] : normalize(cross(GT[trig + 1] - GT[trig], GT[trig + 2] - GT[trig]));
	}

	for(auto & obj: p.objects)
		objects.push_back({obj.second, material_id[obj.first]});

	gproshan_log_var(size(objects));

	trig_mat = new index_t[n_vertices / 3];
	memset(trig_mat, -1, sizeof(index_t) * n_vertices / 3);

	#pragma omp parallel for
	for(index_t i = 0; i < size(objects) - 1; ++i)
	{
		const object & obj = objects[i];
		const index_t n = objects[i + 1].begin;

		for(index_t t = obj.begin; t < n; t += 3)
			trig_mat[t / 3] = obj.material_id;
	}

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
				material_id[str] = size(materials);
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
				float & d = materials.back().d;
				sscanf(line, "%*s %f", &d);
				break;
			}
			case 'T':	// Tr
			{
				float & d = materials.back().d;
				sscanf(line, "%*s %f", &d);
				d = 1 - d;
				break;
			}
			case 'N':	// Ns
			{
				float & N = str[1] == 's' ? materials.back().Ns : materials.back().Ni;
				sscanf(line, "%*s %f", &N);
				break;
			}
			case 'i':	// illum
			{
				int & illum = materials.back().illum;
				sscanf(line, "%*s %u", &illum);
				break;
			}
			case 'm':	// map_Ka, map_kd
			{
				int & m = str[4] == 'K' && str[5] == 'a' ? materials.back().map_Ka
						: str[4] == 'K' && str[5] == 'd' ? materials.back().map_Kd
						: str[4] == 'K' && str[5] == 's' ? materials.back().map_Ks
						: str[4] == 'd' ? materials.back().map_d
						: materials.back().map_bump;
				sscanf(line, "%*s %s", str);
				if(str[0] == '-') continue;		// ignoring map textures with options
				if(!texture_id.contains(str))
				{
					texture_id[str] = size(texture_name);
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

/*
	for(index_t i = 0; i < size(materials); ++i)
	{
		const material & m = materials[i];
		gproshan_log_var(material_name[i]);
		gproshan_log_var(m.Ka);
		gproshan_log_var(m.Kd);
		gproshan_log_var(m.Ks);
		gproshan_log_var(m.d);
		gproshan_log_var(m.Ns);
		gproshan_log_var(m.illum);
		if(m.map_Ka > -1)	gproshan_log_var(texture_name[m.map_Ka]);
		if(m.map_Kd > -1)	gproshan_log_var(texture_name[m.map_Kd]);
		if(m.map_Ks > -1)	gproshan_log_var(texture_name[m.map_Ks]);
		if(m.map_d > -1)	gproshan_log_var(texture_name[m.map_d]);
		if(m.map_bump > -1)	gproshan_log_var(texture_name[m.map_bump]);
	}
*/

	gproshan_log_var(size(materials));
	gproshan_log_var(size(textures));

	return true;
}

bool scene::load_texture(const std::string & file)
{
	try
	{
		CImg<unsigned char> img(file.c_str());
		img.mirror('y');

		textures.emplace_back();
		texture & tex = textures.back();
		tex.width = img.width();
		tex.height = img.height();
		tex.spectrum = img.spectrum();
		tex.data = new unsigned char[tex.width * tex.height * tex.spectrum];
		img.permute_axes("cxyz");
		memcpy(tex.data, img.data(), tex.width * tex.height * tex.spectrum);
	}
	catch(CImgException & e)
	{
		gproshan_error_var(e.what());
		return false;
	}

	return true;
}


} // namespace gproshan
