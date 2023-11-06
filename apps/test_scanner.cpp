#include <gproshan/mesh/che_ply.h>
#include <gproshan/mesh/che_ptx.h>
#include <gproshan/scenes/scanner.h>
#include <gproshan/raytracing/embree.h>

#include <algorithm>
#include <cstdlib>


int main(int argc, char* argv[])
{
	if(argc < 2 || argc > 7)
	{
		std::cerr << "Correct usage: ./apps/test_scanner n_rows n_cols pc_radius ptx_folder jpg_folder" << std::endl;
		return -1;
 	}

	size_t n_rows = atoi(argv[2]);
	size_t n_cols = atoi(argv[3]);

	float pc_radius = atof(argv[4]);

	char * ptx_folder = argv[5];
	char * jpg_folder = argv[6];

	gproshan_log_var(ptx_folder);
	gproshan_log_var(jpg_folder);

	gproshan::che_ply * mesh_ply = new gproshan::che_ply(argv[1]);

	const gproshan::mat4 & model_mat = mesh_ply->normalize_box();

	gproshan::rt::raytracing * rt_embree = new gproshan::rt::embree({mesh_ply}, {model_mat}, false, pc_radius);

	gproshan::che * ptx_mesh = gproshan::scanner_ptx_jpg(rt_embree, n_rows, n_cols, {0, 0, 0}, jpg_folder + mesh_ply->name());

	std::string ptx_filename = ptx_folder + mesh_ply->name();
	gproshan::che_ptx::write_file(ptx_mesh, ptx_filename, n_rows, n_cols);

	delete mesh_ply;
	delete rt_embree;
	delete ptx_mesh;

	return 0;
}

