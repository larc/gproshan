#include "keycomponents.h"
#include "test.h"
#include <cstring>
#include <sstream>
#include <queue>

keycomponents::keycomponents(keypoints & kps, off & shape, size_t k_, percent_t percent)
{
	n_vertices = shape.get_nvertices();
	components = new size_t * [n_vertices];
	memset(components, 0, sizeof(size_t*)*n_vertices);

	marcados = new size_t[n_vertices];
	memset(marcados, 0, sizeof(size_t)*n_vertices);

	n_keypoints = n_vertices * percent;
	v_keypoints = kps.get_keypoints();

	k = k_;
	n_components = 0;
	is_keypoint = NULL;

	m_index = new queue<size_t*>[n_keypoints + 1];
	components_index = 0;
	//agrupar_kpoints(shape);
	calcular_krings(shape);
	contar_components();
}

//FM key components
keycomponents::keycomponents(keypoints & kps, off & shape, fastmarching & fm, distance_t radio_, percent_t percent)
{
	n_vertices = shape.get_nvertices();
	components_index = new size_t[n_vertices + 1];
	for(size_t i = 0; i <= n_vertices; i++)
		components_index[i] = i;

	marcados = NULL;

	n_keypoints = n_vertices * percent;
	v_keypoints = kps.get_keypoints();
	is_keypoint = new bool[n_vertices];
	memset(is_keypoint, 0, sizeof(bool) * n_vertices);

	for(index_t i = 0; i < n_keypoints; i++)
		is_keypoint[v_keypoints[i]] = 1;

	radio = radio_;
	n_components = 0;

	m_index = 0;
	components = 0;

	calcular_rings_fm(shape, fm);
//	contar_components();
}

keycomponents::~keycomponents()
{
	if(is_keypoint) delete [] is_keypoint;
	if(components) delete [] components;
	if(marcados) delete [] marcados;
	if(components_index) delete [] components_index;
	if(!m_index) return;
	size_t * q;
	for(size_t i = 0; i <= n_keypoints; i++)
		while(!m_index[i].empty())
		{
			q = m_index[i].front();
			m_index[i].pop();
			delete q;
		}
	delete [] m_index;
}

void keycomponents::print(ostream & os)
{
	for(size_t i = 0; i < n_vertices; i++)
		if(components[i])
			os<<*components[i]<<" ";
		else
			os<<0<<" ";

	os<<endl;
}

void keycomponents::print_fm(ostream & os)
{
	for(size_t i = 0; i < n_vertices; i++)
		os<<components_index[i]<<" ";

	os<<endl;
}


size_t keycomponents::get_ncomponents()
{
	return n_components;
}

void keycomponents::calcular_krings(off & shape)
{
	size_t * index;
	size_t nc = 0;

	for(size_t i = 0; i < n_keypoints; i++)
	{
		if(!components[v_keypoints[i]])
		{
			index = new size_t(++nc);
			m_index[nc].push(index);
		}
		else
			index = components[v_keypoints[i]];

		memset(marcados, 0, sizeof(size_t)*n_vertices);
		calcular_krings(shape, v_keypoints[i], index, k);
		/*
		stringstream ss;
		ss<<"../TEST/keycomponents/"<<i;
		ofstream os(ss.str());
		print(os);
		os.close();
		*/
	}

	nc++;

	vertex * centers = new vertex[nc];
	size_t * nro_x_comp = new size_t[nc];
	vertex_t * stds = new vertex_t[nc];

	for(size_t i = 1; i < nc; i++)
	{
		centers[i] = vertex(0,0,0);
		nro_x_comp[i] = 0;
		stds[i] = 0;
	}

	for(size_t i = 0; i < n_keypoints; i++)
	{
		centers[*components[v_keypoints[i]]] += shape(v_keypoints[i]);
		nro_x_comp[*components[v_keypoints[i]]]++;
	}

	for(size_t i = 0; i < nc; i++)
		centers[i] /= nro_x_comp[i];

	vertex_t dist;

	for(size_t i = 0; i < n_keypoints; i++)
	{
		dist = *(shape(v_keypoints[i]) - centers[*components[v_keypoints[i]]]);
		if(dist > stds[*components[v_keypoints[i]]])
			stds[*components[v_keypoints[i]]] = dist;
	}

	for(size_t i = 0; i < n_vertices; i++)
		if(components[i] && (*(shape(i) - centers[*components[i]])) > NSTD*stds[*components[i]])
			components[i] = 0;

	/*
	for(size_t i = 0; i < n_keypoints; i++)
	{
		dist = *(shape(v_keypoints[i]) - centers[*components[v_keypoints[i]]]);
		stds[*components[v_keypoints[i]]] += dist*dist;
	}

	for(size_t i = 0; i < nc; i++)
	{
		stds[i] /= sqrt(stds[i]);
		stds[i] /= nro_x_comp[i];
	}

	for(size_t i = 0; i < n_points; i++)
		if(components[i] && (*(shape(i) - centers[*components[i]])) > NSTD*stds[*components[i]])
			components[i] = 0;
	*/

	delete [] centers;
	delete [] nro_x_comp;
	delete [] stds;
}
//Funcion maestra
void keycomponents::calcular_krings(off & shape, size_t p, size_t * index, size_t kr)
{
	if(marcados[p] >= kr) return;

	size_t ir;
	marcados[p] = kr;
	if(!components[p])
		components[p] = index;
	else if(*components[p] != *index)
	{
		ir = *components[p];

		while(!m_index[ir].empty())
		{
			*m_index[ir].front() = *index;
			m_index[*index].push(m_index[ir].front());
			m_index[ir].pop();
		}
	}

	for(auto it: shape.get_rings(p))
		calcular_krings(shape, it, index, kr-1);
}

//with fast marching
void keycomponents::calcular_rings_fm(off & shape, fastmarching & fm)
{
	for(size_t i = 0; i < fm.get_npesos(); i++)
		for(size_t j: shape.get_rings(fm(i)))
			if(fm[j] <= radio)
				join(components_index, fm(i), j);

	size_t * contar = new size_t[n_vertices];
	size_t * contar_keypoints = new size_t[n_vertices];
	memset(contar, 0, n_vertices*sizeof(size_t));
	memset(contar_keypoints, 0, n_vertices*sizeof(size_t));

	for(size_t i = 0; i < n_vertices; i++)
	{
		contar[find(components_index, i)]++;
		if(is_keypoint[i])
			contar_keypoints[find(components_index, i)]++;
	}


	for(size_t i = 0; i < n_vertices; i++)
		if(find(components_index, i) == i)
			components_index[i] = 0;

	for(size_t i = 0; i < n_vertices; i++)
	{
		if(contar_keypoints[i] > 1)
			cout<< contar_keypoints[i];
		if(contar[i] > 1)
			cout << " / " << contar[i] << endl;
	}

	n_components = 0;
	for(size_t i = 0; i < n_vertices; i++)
		if(contar[i] > 1)
			n_components++;

	delete [] contar;
}

void keycomponents::contar_components()
{
	map<size_t, bool> nro_comp;
	for(size_t i = 0; i < n_vertices; i++)
		if(components[i])
			nro_comp[*components[i]] = 1;

	n_components = nro_comp.size();
}

void keycomponents::new_keypoints(off & shape, keypoints & kps, ostream & os, percent_t p)
{
	cout<<__FUNCTION__<<endl;
	//map<size_t, vector<size_t> > set_components;
	map<size_t, priority_queue<pair<area_t, size_t> > > set_components;

	for(size_t i = 0; i < n_vertices; i++)
		if(components[i])
			set_components[*components[i]].push(make_pair(kps.get_mean_area(i), i));

	size_t n;
	for(auto c: set_components)
	{
		n = c.second.size()*p;
		for(size_t i = 0; i < n; i++)
		{
			os<<c.second.top().second<<" ";
			c.second.pop();
		}
	}
	/*
	vertex * means = new vertex[set_components.size()];

	size_t ic, min_i;
	vertex_t min;
	ic = 0;
	for(auto c: set_components)
	{
		for(auto p: c.second)
			means[ic] += shape(p);

		means[ic] /= c.second.size();

		min = INFINITY;
		for(auto p: c.second)
		{
			vertex_t dt = *(means[ic] - shape(p));
			if(dt < min)
			{
				min = dt;
				min_i = p;
			}
		}
		os<<min_i<<" ";

		ic++;
	}
	*/

	os<<endl;
}

void keycomponents::agrupar_kpoints(off & shape)
{
	size_t * index;
	size_t nc = 0;
	bool * keypoints;
	keypoints = new bool[n_vertices];
	memset(keypoints, 0, sizeof(bool)*n_vertices);

	for(size_t i = 0; i < n_keypoints; i++)
		keypoints[v_keypoints[i]] = 1;

	for(size_t i = 0; i < n_keypoints; i++)
	{
		if(!components[v_keypoints[i]])
		{
			index = new size_t(++nc);
			m_index[nc].push(index);
		}
		else
			index = components[v_keypoints[i]];

		memset(marcados, 0, sizeof(size_t)*n_vertices);
		agrupar_kpoints(shape, v_keypoints[i], index, k, keypoints);
	}

	delete [] keypoints;
}

void keycomponents::agrupar_kpoints(off & shape, size_t p, size_t * index, size_t kr, bool * keypoints)
{
	if(marcados[p] >= kr) return;

	size_t ir;
	marcados[p] = kr;
	if( keypoints[p] && !components[p] )
		components[p] = index;
	else if( keypoints[p] && *components[p] != *index)
	{
		ir = *components[p];

		while(!m_index[ir].empty())
		{
			*m_index[ir].front() = *index;
			m_index[*index].push(m_index[ir].front());
			m_index[ir].pop();
		}
	}

	for(auto it: shape.get_rings(p))
		agrupar_kpoints(shape, it, index, kr-1, keypoints);
}

void join(size_t * n, size_t a, size_t b)
{
	size_t fa = find(n, a);
	size_t fb = find(n, b);
	if(fa == fb) return;
	if(fa > fb)	n[fb] = fa;
	else n[fa] = fb;
}

size_t find(size_t * n, size_t a)
{
	if(a == n[a]) return a;
	return n[a] = find(n, n[a]);
}
