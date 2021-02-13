#include "features/key_points.h"

#include <cassert>
#include <cstring>
#include <algorithm>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


key_points::key_points(che * mesh, const real_t & percent)
{
	n_faces = mesh->n_faces;
	n_vertices = mesh->n_vertices;

	face_areas = new real_idx_t[n_faces];
	
	kps = new index_t[n_vertices];
	is_kp = new bool[n_vertices];

	n_kps = percent * n_vertices;
	compute_kps(mesh);
}

key_points::~key_points()
{
	delete [] face_areas;
	delete [] kps;
	delete [] is_kp;
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
	return n_kps;
}

void key_points::compute_kps(che * mesh)
{
	// compute faces areas

	#pragma omp parallel for
	for(index_t t = 0; t < n_faces; t++)
	{
		face_areas[t].first = mesh->area_trig(t);
		face_areas[t].second = t;
	}

	sort(face_areas, face_areas + n_faces);
	
	// compute kps
	memset(is_kp, 0, sizeof(bool) * n_vertices);

	index_t he, k = 0;
	for(index_t t = 0; t < n_faces; t++)
	{
		he = che::mtrig * face_areas[t].second;
		for(index_t i = 0; i < che::mtrig; i++)
		{
			const index_t & v = mesh->vt(he);
			if(!is_kp[v])
			{
				kps[k++] = v;
				is_kp[v] = 1;
			}
			he = next(he);
		}
	}
	
	// compute kps
	memset(is_kp, 0, sizeof(bool) * n_vertices);

	#pragma omp parallel for
	for(index_t i = 0; i < n_kps; i++)
		is_kp[kps[i]] = 1;
}


} // namespace gproshan

