#include "gproshan/scenes/scene.h"


// geometry processing and shape analysis framework
namespace gproshan {


bool scene::read_mtl(const std::string & file)
{
	gproshan_error_var(file);

	FILE * fp = fopen(file.c_str(), "r");
	if(!fp) return false;

	char line[256], str[64];
	while(fgets(line, sizeof(line), fp))
	{
		switch(line[0])
		{
			case 'n':	// newmtl
			{
				sscanf(line, "%*s %s", str);
				material_id[str] = materials.size();
				material_name.push_back(str);
				materials.push_back({});
				break;
			}
			case 'K':	// Ka, Kd, Ks
			{
				vec3 & rgb = line[1] == 'a' ? materials.back().Ka
							: line[1] == 'd' ? materials.back().Kd
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
				real_t & d = materials.back().Ns;
				sscanf(line, "%*s %f", &d);
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
				index_t & m = line[5] == 'a' ? materials.back().map_Ka : materials.back().map_Kd;
				sscanf(line, "%*s %s", str);
				if(texture_id.find(str) == texture_id.end())
				{
					texture_id[str] = textures.size();
					texture_name.push_back(str);
					textures.push_back({});
				}
				m = texture_id[str];
				break;
			}
		}
	}

	fclose(fp);

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

	return true;
}


} // namespace gproshan

