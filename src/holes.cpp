#include "holes.h"
#include <cstring>
#include <fstream>

holes::holes(keypoints & kps, off & shape)
{
	n_vertices = shape.get_nvertices();
	n_faces = shape.get_nfaces();

	holes_vertices = new bool[n_vertices];
	holes_faces = new bool[n_faces];
	neighbors_bound = new size_t[n_vertices];
	deleted_faces = new bool[n_faces];
	mapea = new size_t[n_vertices];
	tangentes = new curv_t[n_vertices];

	memset(holes_vertices, 0, sizeof(bool)*n_vertices);
	memset(holes_faces, 0, sizeof(bool)*n_faces);
	memset(neighbors_bound, 255, sizeof(size_t)*n_vertices);
	memset(deleted_faces, 0, sizeof(bool)*n_faces);
	memset(mapea, 255, sizeof(size_t)*n_vertices);

	mesh_holes = 0;
	curv_holes = 0;

	calcular_holes(kps, shape);
}

holes::~holes()
{
	delete [] holes_vertices;
	delete [] holes_faces;
	delete [] neighbors_bound;
	delete [] deleted_faces;
	delete [] mapea;
	delete [] tangentes;
	if(mesh_holes) delete [] mesh_holes;
	if(curv_holes) delete [] curv_holes;
}

void holes::save_off(string file)
{
	shape_holes.save(file);
}

void holes::print_holes(ostream & os)
{
	for(size_t i = 0; i < n_vertices; i++)
		if(holes_vertices[i]) os<<i<<endl;
}

size_t holes::getNroHoles()
{
	return index_holes.size();
}

void holes::calcular_holes(keypoints & kps, off & shape)
{
	for(size_t i = 0; i < n_vertices; i++)
		holes_vertices[i] = shape.get_rings(i).size() != kps.get_nro_faces(i);

	ofstream os("../TEST/holes/holes");
	print_holes(os);
	os.close();
	calcular_bounds(shape);
	find_holes();

	mesh_holes = new off[index_holes.size()];
	curv_holes = new pair<curv_t, curv_t>[index_holes.size()];
	pair<curv_t, size_t> * curv_index = new pair<curv_t, size_t>[index_holes.size()];
	for(size_t i = 0; i < index_holes.size(); i++)
	{
		repair_hole(shape, i);
		curv_index[i].first = curv_holes[i].second;
		curv_index[i].second = i;
		printf("%10lu\t\t%10.6f\t\t%10.6f\t\t%10.6f\n", i, curv_index[i].first, curv_holes[i].first, curv_holes[i].second);
	}
	sort(curv_index, curv_index + index_holes.size());
	for(size_t i = 0; i < index_holes.size(); i++)
		printf("%lu ", curv_index[i].second);

	delete [] curv_index;

	memset(holes_vertices, 0, sizeof(bool)*n_vertices);
	for(size_t i = 0; i < n_faces; i++)
		if(!deleted_faces[i])
		{
			holes_vertices[shape[i][0]] = 1;
			holes_vertices[shape[i][1]] = 1;
			holes_vertices[shape[i][2]] = 1;
		}

	vector<vertex> vertices;
	vector<face> faces;
	sub_mesh(shape, vertices, faces);

	for(size_t i = 0; i < index_holes.size(); i++)
		add_mesh(i, vertices, faces);

	shape_holes = *(new off(vertices,faces));
}

void holes::calcular_bounds(off & shape)
{
	size_t n_fv, a, b, tmp;
	for(size_t i = 0; i < n_faces; i++)
	{
		n_fv = holes_vertices[shape[i][0]] + holes_vertices[shape[i][1]] + holes_vertices[shape[i][2]];

		if(n_fv == 2)
		{
			n_fv = holes_vertices[shape[i][0]] + holes_vertices[shape[i][1]];
			if(n_fv != 2)
			{
				n_fv = holes_vertices[shape[i][1]] + holes_vertices[shape[i][2]];
				if(n_fv != 2)
				{
					a = shape[i][0];
					b = shape[i][2];
				}
				else
				{
					a = shape[i][1];
					b = shape[i][2];
				}
			}
			else
			{
				a = shape[i][0];
				b = shape[i][1];
			}

			if(neighbors_bound[a] == INFINITY)
				neighbors_bound[a] = b;
			else if(neighbors_bound[b] == INFINITY)
				neighbors_bound[b] = a;
			else
			{
				tmp = neighbors_bound[a];

				while(tmp != INFINITY)
				{
					if(b == tmp) break;

					neighbors_bound[a] = b;
					b = a;
					a = tmp;
					tmp = neighbors_bound[a];

				}

				neighbors_bound[a] = b;
			}
		}
		else if(n_fv == 3) //ignorar faces con 2 aristas en el hoyo
		{
			deleted_faces[i] = 1;
		}
	}
}

void holes::find_holes()
{
	size_t a, i = 0;
	char buffer[128];
	memset(holes_vertices, 0, sizeof(bool)*n_vertices);

	while(i < n_vertices)
	{
		while(i < n_vertices && (holes_vertices[i] || neighbors_bound[i] == INFINITY)) i++;
		if(i == n_vertices) break;

		a = i;
		index_holes.push_back(a);
		holes_vertices[a] = 1;

		sprintf(buffer,"../TEST/holes/holes_%lu", index_holes.size()-1);
		ofstream os(buffer);
		os<<a<<endl;
		i = neighbors_bound[a];

		while(a != i)
		{
			os<<i<<endl;
			holes_vertices[i] = 1;
			i = neighbors_bound[i];
		}
		os.close();
		i++;
	}
}

void holes::repair_hole(off & shape, size_t index)
{
	if(index >= index_holes.size()) return;

	curv_holes[index].first = 0;
	curv_holes[index].second = 0;

	vertex centro;
	size_t n = index_holes[index];
	size_t j, i = n;
	size_t tam = 1;

	memset(holes_vertices, 0, sizeof(bool)*n_vertices);
	holes_vertices[i] = 1;
	centro += shape(i);

	i = neighbors_bound[i];
	while(i != n)
	{
		centro += shape(i);
		tam++;
		holes_vertices[i] = 1;
		i = neighbors_bound[i];
	}
	centro /= tam;

	map<vertex, size_t> m_vertices;
	vector<vertex> vertices;
	vector<face> faces;

	j = 0;
	i = n;

	vertices.push_back(shape(i));
	m_vertices[shape(i)] = j++;

	i = neighbors_bound[i];
	while(i != n)
	{
		vertices.push_back(shape(i));
		m_vertices[shape(i)] = j++;

		i = neighbors_bound[i];
	}

	vertices.push_back(centro);
	m_vertices[centro] = j;

	face tface;
	for(i = 0; i < tam; i++)
	{
		tface[0] = i;
		tface[1] = (i + 1) % j;
		tface[2] = j;

		faces.push_back(tface);
	}

	divide_faces(2, tam, vertices, faces);
	//divide_faces(2, tam, vertices, m_vertices, faces);
	//divide_faces(2, vertices, m_vertices, faces);

	mesh_holes[index] = *(new off(vertices, faces));
	off & shape_hole = mesh_holes[index];
	off & shole = *sub_mesh(shape, 1);

	char buffer[128];
	sprintf(buffer,"../TEST/holes/sub_mesh.%lu.off", index);
	shole.save(buffer);
	sprintf(buffer,"../TEST/holes/hole_repair.%lu.off", index);
	shape_hole.save(buffer);

	MatrixXd p(3, shole.get_nvertices());
	for(i = 0; i < p.cols(); i++)
	{
		p(0,i) = shole(i).x;
		p(1,i) = shole(i).y;
		p(2,i) = shole(i).z;
	}

	MatrixXd m(3,1);
	MatrixXd V = pca(p, m);

	MatrixXd h(3, shape_hole.get_nvertices());
	for(i = 0; i < h.cols(); i++)
	{
		h(0,i) = shape_hole(i).x;
		h(1,i) = shape_hole(i).y;
		h(2,i) = shape_hole(i).z;

		h.col(i) -= m;
	}

	MatrixXd P = V.transpose()*p;
	MatrixXd H = V.transpose()*h;

	//ANALYSIS CURVATURE
	analysis_curvature_2(index, P, H, shole, shape_hole);

	//BIHARMONIC SPLINES INTERPOLANT
	biharmonic_interp_2(P, H);

	H = V*H;

	for(i = 0; i < H.cols(); i++)
	{
		H.col(i) += m;
		shape_hole(i).x = H(0,i);
		shape_hole(i).y = H(1,i);
		shape_hole(i).z = H(2,i);
	}

	sprintf(buffer,"../TEST/holes/hole_repair_2.%lu.off", index);
	shape_hole.save(buffer);
}

void holes::add_mesh(size_t index, vector<vertex> & vertices, vector<face> & faces)
{
	if(index >= index_holes.size()) return;
	off & shape_hole = mesh_holes[index];

	size_t n = index_holes[index];
	size_t i = n;
	size_t j = 0;
	size_t * map_vertex = new size_t[shape_hole.get_nvertices()];

	map_vertex[j++] = mapea[i];
	i = neighbors_bound[i];
	while(i != n)
	{
		map_vertex[j++] = mapea[i];
		//cout<<mapea[i]<<endl;
		i = neighbors_bound[i];
	}

	size_t nv = vertices.size();

	while(j < shape_hole.get_nvertices())
	{
		vertices.push_back(shape_hole(j));
		//cout<<nv<<endl;
		map_vertex[j++] = nv++;
	}

	face tface;
	for(i = 0; i < shape_hole.get_nfaces(); i++)
	{
		tface[2] = map_vertex[shape_hole[i][0]];
		tface[1] = map_vertex[shape_hole[i][1]];
		tface[0] = map_vertex[shape_hole[i][2]];

		faces.push_back(tface);
	}
}

void holes::divide_faces(size_t k, size_t nf, vector<vertex> & vertices, vector<face> & faces)
{
	if(!k) return;
	size_t nv = vertices.size();

	vertex c;
	face tface;
	size_t i0, i1, i2;
	size_t mod = nv + nf;
	size_t in = nv;

	for(size_t i = 0; i < nf; i++)
	{
		i0 = faces[i][0];
		i1 = faces[i][1];
		i2 = faces[i][2]; //centro

		c = (vertices[i0] + vertices[i1] + vertices[i2])/3;

		vertices.push_back(c);

		faces[i][0] = nv;
		faces[i][1] = (nv + 1) == mod ? in : nv + 1;

		tface[0] = i0;
		tface[1] = i1;
		tface[2] = nv;
		faces.push_back(tface);

		tface[0] = nv;
		tface[1] = i1;
		tface[2] = (nv + 1) == mod ? in : nv + 1;
		faces.push_back(tface);
		nv++;
	}

	divide_faces(k - 1, nf, vertices, faces);
}

void holes::divide_faces(size_t k, size_t nf, vector<vertex> & vertices, map<vertex, size_t> & m_vertices, vector<face> & faces)
{
	if(!k) return;
	size_t nv = vertices.size();

	vertex x, y;
	size_t i0, i1, i2;
	size_t ix, iy;

	face tface;

	for(size_t i = 0; i < nf; i++)
	{
		i0 = faces[i][0];
		i1 = faces[i][1];
		i2 = faces[i][2]; //centro

		x = (vertices[i0] + vertices[i2])/2;
		y = (vertices[i1] + vertices[i2])/2;

		if(m_vertices.find(x) == m_vertices.end())
		{
			vertices.push_back(x);
			m_vertices[x] = nv++;
		}
		if(m_vertices.find(y) == m_vertices.end())
		{
			vertices.push_back(y);
			m_vertices[y] = nv++;
		}

		ix = m_vertices[x];
		iy = m_vertices[y];

		faces[i][0] = ix;
		faces[i][1] = iy;

		tface[0] = i0;
		tface[1] = iy;
		tface[2] = ix;
		faces.push_back(tface);

		tface[0] = i0;
		tface[1] = i1;
		tface[2] = iy;
		faces.push_back(tface);
	}

	divide_faces(k-1, nf, vertices, m_vertices, faces);
}

void holes::divide_faces(size_t k, vector<vertex> & vertices, map<vertex, size_t> & m_vertices, vector<face> & faces)
{
	if(!k) return;
	size_t nv = vertices.size();
	size_t nf = faces.size();

	vertex a, b, c;
	size_t ia, ib, ic;
	size_t i0, i1, i2;
	face tface;
	for(size_t i = 0; i < nf; i++)
	{
		i0 = faces[i][0];
		i1 = faces[i][1];
		i2 = faces[i][2];

		a = (vertices[i0] + vertices[i1])/2;
		b = (vertices[i0] + vertices[i2])/2;
		c = (vertices[i1] + vertices[i2])/2;

		if(m_vertices.find(a) == m_vertices.end())
		{
			vertices.push_back(a);
			m_vertices[a] = nv++;
		}
		if(m_vertices.find(b) == m_vertices.end())
		{
			vertices.push_back(b);
			m_vertices[b] = nv++;
		}
		if(m_vertices.find(c) == m_vertices.end())
		{
			vertices.push_back(c);
			m_vertices[c] = nv++;
		}
		ia = m_vertices[a];
		ib = m_vertices[b];
		ic = m_vertices[c];

		faces[i][1] = ia;
		faces[i][2] = ib;

		tface[0] = i1;
		tface[1] = ic;
		tface[2] = ia;
		faces.push_back(tface);

		tface[0] = ia;
		tface[1] = ic;
		tface[2] = ib;
		faces.push_back(tface);

		tface[0] = i2;
		tface[1] = ib;
		tface[2] = ic;
		faces.push_back(tface);
	}

	divide_faces(k-1, vertices, m_vertices, faces);
}

void holes::biharmonic_interp_2(MatrixXd & P, MatrixXd & H)
{
	size_t n = P.cols();
	vertex_t x;

	MatrixXd A(n, n);
	MatrixXd pi(2,1);
	MatrixXd pj(2,1);
	for(size_t i = 0; i < n; i++)
	{
		pi(0,0) = P(0,i);
		pi(1,0) = P(1,i);
		for(size_t j = 0; j < n; j++)
		{
			pj(0,0) = P(0,j);
			pj(1,0) = P(1,j);
			x = (pi - pj).norm();
			A(i,j) = x*x*(log(x + 1e-18) - 1);
		}
	}

	FullPivLU<MatrixXd> dec(A);
	MatrixXd alpha = dec.solve(P.row(2).transpose());

	for(size_t i = 0; i < H.cols(); i++)
	{
		H(2,i) = 0;
		pi(0,0) = H(0,i);
		pi(1,0) = H(1,i);
		for(size_t j = 0; j < n; j++)
		{
			pj(0,0) = P(0,j);
			pj(1,0) = P(1,j);
			x = (pi - pj).norm();
			x = x*x*(log(x + 1e-18) - 1);
			H(2,i) += alpha(j,0)*x;
		}
	}
}

void holes::analysis_curvature_2(size_t index, MatrixXd & P, MatrixXd & H, off & p_mesh, off & h_mesh)
{
	size_t n = P.cols();
	vertex_t x;

	MatrixXd A(n, n);
	MatrixXd pi(2,1);
	MatrixXd pj(2,1);
	for(size_t i = 0; i < n; i++)
	{
		pi(0,0) = P(0,i);
		pi(1,0) = P(1,i);
		for(size_t j = 0; j < n; j++)
		{
			pj(0,0) = P(0,j);
			pj(1,0) = P(1,j);
			x = (pi - pj).norm();
			A(i,j) = x*(2*log(x + 1e-18) - 1);
		}
	}

	FullPivLU<MatrixXd> dec(A);
	MatrixXd alpha = dec.solve(P.row(2).transpose());

	//SUB MESH P
	curv_holes[index].first = analysis_curvature_2(alpha, P, P, p_mesh);

	//GENERADA H
	curv_holes[index].second = analysis_curvature_2(alpha, P, H, h_mesh);
}

curv_t holes::analysis_curvature_2(MatrixXd & alpha, MatrixXd & P, MatrixXd & H, off & mesh)
{
	MatrixXd pi(2,1);
	MatrixXd pj(2,1);
	vertex_t x;

	for(size_t i = 0; i < H.cols(); i++)
	{
		tangentes[i] = 0;
		pi(0,0) = H(0,i);
		pi(1,0) = H(1,i);
		for(size_t j = 0; j < P.cols(); j++)
		{
			pj(0,0) = P(0,j);
			pj(1,0) = P(1,j);
			x = (pi - pj).norm();
			x = x*(2*log(x + 1e-18) - 1);
			tangentes[i] += abs(alpha(j,0)*x);
		}
	}

	curv_t stag, curv = 0;
	for(size_t i = 0; i < mesh.get_nvertices(); i++)
	{
		stag = 0;
		for(auto j: mesh.get_rings(i))
			stag += abs(tangentes[j] - tangentes[i]);

		stag /= mesh.get_rings(i).size();
		curv += stag;
	}

	curv /= mesh.get_nvertices();
	return curv;
}

MatrixXd holes::pca(MatrixXd & p, MatrixXd & m)
{
	size_t n = p.cols();

	m<<0,0,0;

	for(size_t i = 0; i < n; i++)
		m += p.col(i);
	m /= n;

	MatrixXd r(3,n);
	for(size_t i = 0; i < n; i++)
		r.col(i) = p.col(i) - m;

	MatrixXd c = (r*r.transpose())/n;

	SelfAdjointEigenSolver<MatrixXd> eigensolver(c);

	MatrixXd v = eigensolver.eigenvectors();
	MatrixXd V = v;
	V.col(0) = v.col(2);
	V.col(2) = v.col(0);

	p = r;
	return V;
}

void holes::sub_mesh(off & shape, vector<vertex> & vertices, vector<face> & faces, size_t k)
{
	//Necesita que holes_vertices tenga los vertices marcados
	size_t i;
	while(k--)
	{
		for(i = 0; i < n_faces; i++)
			if(!deleted_faces[i])
				holes_faces[i] = holes_vertices[shape[i][0]] || holes_vertices[shape[i][1]] || holes_vertices[shape[i][2]];

		for(i = 0; i < n_faces; i++)
		{
			if(holes_faces[i])
			{
				holes_vertices[shape[i][0]] = 1;
				holes_vertices[shape[i][1]] = 1;
				holes_vertices[shape[i][2]] = 1;
			}
		}
	}

	size_t nv = 0;
	for(i = 0; i < n_vertices; i++)
		if(holes_vertices[i])
			mapea[i] = nv++;

	for(i = 0; i < n_vertices; i++)
		if(holes_vertices[i])
			vertices.push_back(shape(i));

	face tface;
	for(i = 0; i < n_faces; i++)
		if(holes_faces[i])
		{
			tface[0] = mapea[shape[i][0]];
			tface[1] = mapea[shape[i][1]];
			tface[2] = mapea[shape[i][2]];

			faces.push_back(tface);
		}
}

off * holes::sub_mesh(off & shape, size_t k)
{
	vector<vertex> vertices;
	vector<face> faces;

	sub_mesh(shape, vertices, faces, k);

	return new off(vertices, faces);
}
