#include <gproshan/mesh/che_ply.h>
#include <gproshan/mesh/che_ptx.h>
#include <gproshan/scenes/scanner.h>
#include <gproshan/raytracing/rt_embree.h>

#include <algorithm>
#include <cstdlib>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>


void translate(glm::mat4 & model_mat, const gproshan::vertex & p)
{
	model_mat = glm::translate(model_mat, glm_vec3(p));
}

void scale(glm::mat4 & model_mat, const gproshan::real_t & s)
{
	model_mat = glm::scale(model_mat, {s, s, s});
}

glm::mat4 normalize_coordinates(gproshan::che_ply * mesh)
{
	glm::mat4 model_mat = glm::mat4(1);

	gproshan::vertex pmin = INFINITY;
	gproshan::vertex pmax = 0;

	for(gproshan::index_t v = 0; v < mesh->n_vertices; ++v)
	{
		const gproshan::vertex & p = mesh->point(v);

		pmin.x = std::min(pmin.x, p.x);
		pmin.y = std::min(pmin.y, p.y);
		pmin.z = std::min(pmin.z, p.z);

		pmax.x = std::max(pmax.x, p.x);
		pmax.y = std::max(pmax.y, p.y);
		pmax.z = std::max(pmax.z, p.z);
	}

	scale(model_mat, 2.0 / std::max({pmax.x - pmin.x, pmax.y - pmin.y, pmax.z - pmin.z}));
	translate(model_mat, - (pmax + pmin) / 2);

	return model_mat;
}

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

	gproshan::rt::raytracing * rt_embree;

	glm::mat4 model_mat = normalize_coordinates(mesh_ply);

	rt_embree = new gproshan::rt::embree({mesh_ply},{model_mat}, false, pc_radius);

	gproshan::che * ptx_mesh = scanner_ptx(mesh_ply, rt_embree, n_rows, n_cols, {0, 0, 0}, jpg_folder);

	std::string ptx_filename = ptx_folder + mesh_ply->name();
	gproshan::che_ptx::write_file(ptx_mesh, ptx_filename, n_rows, n_cols);

	return 0;
}

