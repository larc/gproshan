#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "off.h"
#include "keypoints.h"
#include "keycomponents.h"
#include "holes.h"
#include "test.h"
#include "testkeycomponents.h"
#include "testkeypoints.h"
#include "fastmarching.h"
#include "laplacian.h"
#include "che_off.h"
#include "dijkstra.h"
#include "geodesics.h"
#include "geodesics_ptp.h"
#include "test_geodesics_ptp.h"
#include "fairing_taubin.h"
#include "fairing_spectral.h"
#include "sampling.h"
#include "mdict/d_mesh.h"
#include "mdict/d_basis_dct.h"
#include "mdict/d_image_denoising.h"
#include "app_viewer.h"
#include "viewer/viewer.h"

using namespace std;

void generate_grid(size_t s, size_t f, string name);
void generate_grid_cylinder(const vertex_t & radio, vertex_t d_angle, size_t rings);
void main_testkeypoints();
double main_testkeycomponents();
void main_testholes();
void sampling_terrain(string file, size_t s, size_t K, string output = "patches_index", string planes = "normals");
void main_test_holes();
void generate_grid_obtuse(const size_t & nr, const size_t & nc, const vertex_t & d = 1);

int main(int nargs, const char ** args)
{
//	viewer_main(nargs, args);
	main_test_geodesics_ptp(nargs, args);

//	generate_grid_obtuse(81, 10);
//	generate_grid_cylinder(100, 1, 1000);
//	generate_grid(1000, 1, "tmp/grilla1m.off");
//	testkeycomponents prueba(50, 0.15);
//	prueba.one_test_fm("0001.null.0.off","");
//	if(nargs > 1) test_image_denoising(args[1]);


	//sampling_terrain("terreno.off", 4, 64);
//	sampling_shape(nargs, args);
//	if(nargs == 2) sampling_shape(args[1]);
//	distance_t radio = 0.04;

//	main_test_holes();
//	main_testkeycomponents();
//	main_testkeypoints();
	return 0;
}

void main_test_holes()
{
	string filename;
	float time;
	size_t h_size = 20;
	index_t histograma[20];

	while(cin >> filename)
	{
		debug(filename)
		filename = PATH_MDATA + filename + ".off";
		che * shape = new che_off(filename);

		size_t n_vertices, old_n_vertices = shape->n_vertices();
		size_t n_holes = shape->n_borders();
		debug(old_n_vertices)

		vector<index_t> * border_vertices;
		che ** holes;
		TIC(time) tie(border_vertices, holes) = fill_all_holes_meshes(shape); TOC(time)

		memset(histograma, 0, sizeof(histograma));
		vertex_t q, quality = 0;
		size_t n_faces = 0;
		for(index_t b = 0; b < n_holes; b++)
		{
			n_faces += holes[b]->n_faces();
			for(index_t t = 0; t < holes[b]->n_faces(); t++)
			{
				q = holes[b]->pdetriq(t);
				quality += q > 0.6;
				histograma[(index_t) (q * h_size)]++;
			}
		}

		ofstream os(PATH_TEST + "holes/" + shape->name() + ".hist");
		for(index_t h = 0; h < h_size; h++)
			os << h * 1.0 / h_size << " " << histograma[h] << endl;
		os.close();

		printf("%15s & %10lu & %10lu & %10.3f ", shape->name().c_str(), old_n_vertices, shape->n_vertices() - old_n_vertices, quality * 100 / n_faces);
		printf("& %10.3f ", time);

		for(index_t k = 1; k <=3; k++)
		{
			TIC(time) poisson(shape, old_n_vertices, k); TOC(time)
			printf("& %10.3f ", time);
		}

		n_vertices = old_n_vertices;
		index_t k = 2;
		TIC(time)
		for(index_t h = 0; h < n_holes; h++)
			if(holes[h])
			{
				old_n_vertices = n_vertices;
				biharmonic_interp_2(shape, old_n_vertices, n_vertices += holes[h]->n_vertices() - border_vertices[h].size(), border_vertices[h], k);
				delete holes[h];
			}
		TOC(time)
		printf("& %10.3f ", time);

		printf("\\\\\\hline\n");
		delete shape;
	}
}

void sampling_terrain(string file, size_t s, size_t K, string output, string planes)
{
	file = PATH_MDATA + file;
	che * shape_che = new che_off(file);
	size_t n = sqrt(shape_che->n_vertices());
	cout<<n<<endl;

	output = PATH_TEST + output;
	ofstream os(output);

	planes = PATH_TEST + planes;
	ofstream pos(planes);

	size_t v;
	for(size_t i = 0; i < n; i += s)
		for(size_t j = 0; j < n; j += s)
		{
			v = i * n + j;
			pos<<shape_che->normal(v)<<endl;

			vector<index_t> source;
			source.push_back(v);
			geodesics fm(shape_che, source, geodesics::FM, K);

			for(index_t k = 0; k < K; k++)
				os<<" "<<fm[k];
			os<<endl;
		}

	os.close();
	pos.close();
	delete shape_che;
}

void generate_grid(size_t s, size_t f, string name)
{
	off shape;
	shape.generate_grid(s, f, 0.01);
	shape.save(name);
}

void main_testkeypoints()
{
	char buffer[50];
	string comand;

	for(percent_t p = 0.01; p <= 0.50001; p += 0.01)
	{
		comand = "";
		sprintf(buffer, "tkps-%0.2f", p);
		string tmp(buffer);
		comand = comand + "mkdir " + PATH_KPS + tmp;
		if(system(comand.c_str())) cout<<"OK"<<endl;
		testkeypoints tkps(p);
		tkps.execute("database", buffer);
	}
}

double main_testkeycomponents()
{
	clock_t start = clock();
	for(int i = 15; i <= 15; i += 1) //Profundidad de k-rings
	{
		#pragma omp parallel for num_threads(8)

		for(int j = 1; j <= 50; j++) //Cantidad de kpoints a ser usados
		{
			double pj = j;
			pj /= 100;
			char buffer[50];
			string comand="";
			sprintf(buffer, "%d-%0.2f", i, pj);
			string tmp (buffer);
			comand = comand + "mkdir " + PATH_KCS + tmp;
			system (comand.c_str());
			comand="";
			comand = comand + "mkdir " + PATH_KPS + tmp;
			system (comand.c_str());
			testkeycomponents prueba(i, pj);
			cout<<buffer<<endl;
			prueba.execute("database", buffer);
		}
	}
	double time = clock() - start;
	return time /= 1E6;
//	testkeycomponents prueba(2, 0.50);
//	prueba.one_test("0001.null.0.off","");
}

void main_testholes()
{
	string file = "0001.holes.3.off";
	string filename = getFileName(file);

	off shape(PATH_DATA + file);
	//off shape("../DATA/nefertiti-entire.off");

	percent_t percent = 0.1;
	size_t nkps = shape.get_nvertices()*percent;

	keypoints kps(shape);

	ofstream oskp(PATH_KPS + filename + "kp");
	kps.print(oskp, nkps);
	oskp.close();

	holes sholes(kps,shape);
	sholes.save_off("../TEST/holes/mesh_refinada.off");
	cout<<"NRO HOLES: "<<sholes.getNroHoles()<<endl;
}

void generate_grid_obtuse(const size_t & nr, const size_t & nc, const vertex_t & d)
{
	size_t nvg = nr * nc;
	size_t nvi = ((nvg - nc) / (nc * 2)) * (nc - 1);
	size_t nv = nvg + nvi;
	vertex_t dd = 8 * d;

	vector<vertex> vertices(nv);
	for(index_t v = 0, i = 0; i < nr; i++)
	for(index_t j = 0; j < nc; j++, v++)
	{
		vertices[v].x = d * i;
		vertices[v].y = dd * j;
	}

	vector<index_t> faces(P * nvi * 6);
	index_t f = 0;
	for(index_t i = 0, v = nvg; v < nv; v++)
	{
		i = !((v - nvg) % (nc - 1)) ? ((v - nvg) / (nc - 1)) * 2 * nc : i + 1;
		debug(i)
		vertices[v] = vertices[i + nc];
		vertices[v].y += dd / 2;

		faces[f++] = v;
		faces[f++] = i + nc;
		faces[f++] = i;

		faces[f++] = v;
		faces[f++] = i;
		faces[f++] = i + 1;

		faces[f++] = v;
		faces[f++] = i + 1;
		faces[f++] = i + nc + 1;

		faces[f++] = v;
		faces[f++] = i + nc + 1;
		faces[f++] = i + 2 * nc + 1;

		faces[f++] = v;
		faces[f++] = i + 2 * nc + 1;
		faces[f++] = i + 2 * nc;

		faces[f++] = v;
		faces[f++] = i + 2 *nc;
		faces[f++] = i + nc;
		debug(i + 2 * nc + 1 < nv)
	}

	che_off mesh(vertices.data(), vertices.size(), faces.data(), faces.size() / P);
	mesh.write_file("tmp/grid_obtuse.off");
}

void generate_grid_cylinder(const vertex_t & radio, vertex_t d_angle, size_t rings)
{
	size_t n_vrings = 360 / d_angle;
	d_angle *= M_PI / 180;

	size_t n_vertices = 2 + rings * n_vrings;
	size_t n_faces = 2 * n_vrings * rings;

	index_t * faces = new index_t[n_faces * 3];
	vertex * vertices = new vertex[n_vertices];

	vertex_t hight = d_angle * radio;
	vertex_t z = 0;
	index_t v = 1;
	vertex_t a = 0;
	for(index_t r = 0; r < rings; r++, z += hight)
	{
		a = 0;
		for(index_t vr = 0; vr < n_vrings; vr++, v++)
		{
			vertices[v].x = radio * sin(a);
			vertices[v].y = radio * cos(a);
			vertices[v].z = z;
			a += d_angle;
		}
	}

	vertices[v++].z = z - hight;
	assert(v == n_vertices);

	v--;
	index_t f = 0;
	for(index_t vr = 0; vr < n_vrings; vr++)
	{
		faces[f++] = 0;
		faces[f++] = (vr % n_vrings) + 1;
		faces[f++] = ((vr + 1) % n_vrings) + 1;

		faces[f++] = v;
		faces[f++] = v - (vr % n_vrings) - 1;
		faces[f++] = v - ((vr + 1) % n_vrings) - 1;
	}

	v = 1;
	n_faces *= P;
	while(f < n_faces)
	{
		if(v % n_vrings)
		{
			faces[f++] = v;
			faces[f++] = v + n_vrings + 1;
			faces[f++] = v + 1;
			faces[f++] = v;
			faces[f++] = v + n_vrings;
			faces[f++] = v + n_vrings + 1;
		}
		else
		{
			faces[f++] = v;
			faces[f++] = v + 1;
			faces[f++] = v - n_vrings + 1;
			faces[f++] = v;
			faces[f++] = v + n_vrings;
			faces[f++] = v + 1;
		}
		v++;
	}


	debug(n_vrings)
	debug(n_vertices)

	che_off mesh(vertices, n_vertices, faces, f / P);
	mesh.write_file("tmp/grid_cylinder.off");

	delete [] vertices;
	delete [] faces;
}

