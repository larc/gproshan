#include <gproshan/features/key_components.h>

#include <gproshan/geodesics/geodesics.h>

#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


key_components::key_components(che * mesh, const std::vector<index_t> & kps, const real_t & r): radio(r)
{
	n_vertices = mesh->n_vertices;

	comp = new index_t[n_vertices];
	comp_size = new size_t[n_vertices];

	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices; ++i)
	{
		comp[i] = i;
		comp_size[i] = 1;
	}

	n_comp = 0;
	compute_kcs(mesh, kps);
}

key_components::~key_components()
{
	delete [] comp;
	delete [] comp_size;
}

index_t key_components::operator()(const index_t i)
{
	assert(i < n_vertices);

	if(comp[i] == NIL) return NIL;
	return comp_idx[find(i)];
}

key_components::operator const size_t & () const
{
	return n_comp;
}

void key_components::compute_kcs(che * mesh, const std::vector<index_t> & kps)
{
	geodesics fm(mesh, kps);

	radio *= fm.radio();
	for(index_t i = 0; i < n_vertices && fm[fm(i)] <= radio; ++i)
		for(const index_t he: mesh->star(fm(i)))
			join(fm(i), mesh->halfedge(he_next(he)));

	for(index_t i = 0; i < n_vertices; ++i)
		if(comp[i] == i && comp_size[i] > 1)
			comp_idx[i] = n_comp++;
		else if(comp[i] == i) comp[i] = NIL;
}

index_t key_components::find(const index_t x)
{
	if(comp[x] == x) return x;
	return comp[x] = find(comp[x]);
}

bool key_components::join(index_t x, index_t y)
{
	x = find(x);
	y = find(y);

	if(x == y) return 0;

	comp_size[x] += comp_size[y];
	comp[y] = x;

	return 1;
}


} // namespace gproshan

