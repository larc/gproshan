#include "off.h"

off::off(string file)
{
	rings = 0;
	ring_faces =0;

	if(file == "")
	{
		vertices = 0;
		faces = 0;
		n_vertices = n_faces = value = 0;
		return;
	}

	ifstream is(file);
	is>>*this;
	is.close();

	set_rings();
}

off::off(size_t v, size_t f)
{
	vertices = 0;
	faces = 0;
	n_vertices = v;
	n_faces = f;
	value = 0;
	rings = 0;
	ring_faces = 0;

	if(n_vertices && n_faces)
	{
		vertices = new vertex[n_vertices];
		faces = new face[n_faces];
	}
}

off::off(vector<vertex> & v_vertices, vector<face> & v_faces)
{
	rings = 0;
	value = 0;
	ring_faces = 0;
	n_vertices = v_vertices.size();
	n_faces = v_faces.size();
	if(n_vertices && n_faces)
	{
		vertices = new vertex[n_vertices];
		faces = new face[n_faces];

		for(size_t i = 0; i < n_vertices; i++)
			vertices[i] = v_vertices[i];

		for(size_t i = 0; i < n_faces; i++)
			faces[i] = v_faces[i];

		set_rings();
	}
}

off::~off()
{
	if(vertices)
		delete [] vertices;
	if(faces)
		delete [] faces;
	if(rings)
		delete [] rings;
	if(ring_faces)
		delete [] ring_faces;
}

size_t off::get_nvertices()
{
	return n_vertices;
}

size_t off::get_nfaces()
{
	return n_faces;
}

vertex & off::operator()(size_t i)
{
	return vertices[i];
}

face & off::operator[](size_t i)
{
	return faces[i];
}

ring_t & off::get_rings(size_t i)
{
	return rings[i];
}

face_t & off::get_faces(size_t i)
{
	return ring_faces[i];
}

void off::save(string file)
{
	ofstream os(file);
	os<<*this<<endl;
	os.close();
}

void off::generate_grid(size_t s, size_t f, real_t r)
{
	srand(time(NULL));

	size_t base = s/f;

	n_vertices = base + 1;
	n_vertices *= n_vertices;
	n_faces = base*base*2;

	vertices = new vertex[n_vertices];
	faces = new face[n_faces];

	s = base*f;

	real_t x = 0, y = 0;
	size_t i, j, k, rnd;

	i = 0;
	while(x <= s)
	{
		y = 0;
		while(y <= s)
		{
			vertices[i++] = vertex(x*r, y*r, 0);
			y += f;
		}
		x += f;
	}

	k = n_vertices - base - 1;
	i = 0;
	j = 0;
	rnd = 0;
	while(i < k)
	{
		rnd = rand() % 2;

		if(i % (base + 1))
		{
			faces[j][0] = i - 1;
			faces[j][1] = i + base;
			faces[j][2] = rnd ? i : i + base + 1;
			j++;
			faces[j][0] = rnd ? i : i - 1;
			faces[j][1] = rnd ? i + base : i + base + 1;
			faces[j][2] = rnd ? i + base + 1 : i;
			j++;
		}
		i++;
	}
}

void off::set_rings()
{
	if(rings) delete [] rings;
	rings = new ring_t[n_vertices];

	if(ring_faces) delete [] ring_faces;
	ring_faces = new face_t [n_vertices];

	for(size_t j = 0; j < n_faces; j++)
	{
		rings[faces[j][0]].insert(faces[j][1]);
		rings[faces[j][0]].insert(faces[j][2]);

		rings[faces[j][1]].insert(faces[j][0]);
		rings[faces[j][1]].insert(faces[j][2]);

		rings[faces[j][2]].insert(faces[j][0]);
		rings[faces[j][2]].insert(faces[j][1]);

		ring_faces[faces[j][0]].push_back(j);
		ring_faces[faces[j][1]].push_back(j);
		ring_faces[faces[j][2]].push_back(j);
	}
}

ostream & operator<<(ostream & os, off & o)
{
	os<<"OFF"<<endl;
	os<<o.n_vertices<<" "<<o.n_faces<<" "<<o.value<<endl;
	for(size_t i = 0; i < o.n_vertices; i++)
		os<<o.vertices[i]<<endl;
	for(size_t i = 0; i < o.n_faces; i++)
		os<<o.faces[i]<<endl;
	return os;
}

istream & operator>>(istream & is, off & o)
{
	string soff;
	is>>soff;
	is>>o.n_vertices>>o.n_faces>>o.value;
	o.vertices = new vertex[o.n_vertices];
	o.faces = new face[o.n_faces];
	for(size_t i = 0; i < o.n_vertices; i++)
		is>>o.vertices[i];
	for(size_t i = 0; i < o.n_faces; i++)
		is>>o.faces[i];
	return is;
}
