#include "patch.h"

#include "dictionary.h"

/// Mesh dictionary learning and sparse coding namespace
namespace mdict {

size_t patch::expected_nv = 3 * dictionary::F * (dictionary::F + 1);

void patch::init(che * mesh, const index_t & v, index_t * _level)
{
	index_t * level = _level ? _level : new index_t[mesh->n_vertices()];
	
	if(vertices.size()) clear();

	vertices.reserve(expected_nv << 1);
	memset(level, -1, sizeof(index_t) * mesh->n_vertices());
	
	link_t link;
	level[v] = 0;
	vertices.push_back(v);
	for(index_t i = 0; i < vertices.size(); i++)
	{
		const index_t & v = vertices[i];
		if(level[v] == dictionary::F)
			break;
		
		mesh->link(link, v);
		for(const index_t & he: link)
		{
			const index_t & u = mesh->vt(he);
			if(level[u] == NIL)
			{
				vertices.push_back(u);
				level[u] = level[v] + 1;
			}
		}

		link.clear();	
	}

	if(_level) delete [] level;
}

patch::operator const vector<index_t> & () const
{
	return vertices;
}

void patch::clear()
{
	vertices.clear();
}

} // mdict

