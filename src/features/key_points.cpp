#include "features/key_points.h"

#include <cassert>
#include <cstring>
#include <algorithm>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


key_points::key_points(che * mesh, const real_t & percent)
{
	compute_kps_areas(mesh, percent);
}

const index_t & key_points::operator[](const index_t & i) const
{
	assert(i < n_vertices);
	return kps[i];
}

const bool & key_points::operator()(const index_t & i) const
{
	assert(i < n_vertices);
	return is_kp[i];
}

const size_t & key_points::size() const
{
	return kps.size();
}

/// 
void key_points::compute_kps_areas(che * mesh, const real_t & percent)
{
	std::pair<real_t, index_t> face_areas;

	#pragma omp parallel for
	for(index_t t = 0; t < mesh->n_faces; ++t)
	{
//		face_areas[t].first = mesh->area_trig(t);
//		face_areas[t].second = t;
	}

//	sort(face_areas, face_areas + n_faces);

	is_kp.assign(false, mesh->n_vertices);
	kps.reserve(mesh->n_vertices);
	
	for(index_t t = 0; t < mesh->n_faces; ++t)
	{
		index_t he = che::mtrig;// * face_areas[t].second;
		for(index_t i = 0; i < che::mtrig; ++i)
		{
			const index_t & v = mesh->vt(he);
			if(!is_kp[v])
			{
				kps.push_back(v);
				is_kp[v] = true;
			}
			he = next(he);
		}
	}

	is_kp.assign(false, mesh->n_vertices);
	
	#pragma omp parallel for
	for(const index_t & v: kps)
		is_kp[v] = true;
}


} // namespace gproshan

