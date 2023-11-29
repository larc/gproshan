#include <gproshan/features/key_points.h>

#include <cassert>
#include <cstring>
#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan {


key_points::key_points(che * mesh, const real_t & percent)
{
	compute_kps_areas(mesh, percent);
}

key_points::operator const std::vector<index_t> & () const
{
	return kps;
}

/// Efficient approach for interest points detection in non-rigid shapes
/// Cristian Jose Lopez Del Alamo; Luciano Arnaldo Romero Calla; Lizeth Joseline Fuentes Perez
/// DOI: 10.1109/CLEI.2015.7359459
void key_points::compute_kps_areas(che * mesh, const real_t & percent)
{
	std::vector<std::pair<real_t, index_t> > face_areas(mesh->n_trigs);

	#pragma omp parallel for
	for(index_t f = 0; f < mesh->n_trigs; ++f)
		face_areas[f] = { mesh->area_trig(f), f };

	std::ranges::sort(face_areas);

	is_kp.assign(mesh->n_vertices, false);
	kps.reserve(mesh->n_vertices);

	for(index_t he, t = 0; t < mesh->n_trigs; ++t)
	{
		he = che::mtrig * face_areas[t].second;
		for(index_t i = 0; i < che::mtrig; ++i)
		{
			const index_t & v = mesh->halfedge(he);
			if(!is_kp[v])
			{
				kps.push_back(v);
				is_kp[v] = true;
			}
			he = he_next(he);
		}
	}

	kps.resize(percent * size(kps));
	is_kp.assign(mesh->n_vertices, false);

	#pragma omp parallel for
	for(const index_t & v: kps)
		is_kp[v] = true;
}


} // namespace gproshan

