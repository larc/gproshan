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
#include "fairing_taubin.h"
#include "fairing_spectral.h"
#include "sampling.h"
#include "d_mesh_apps.h"
#include "d_mesh.h"
#include "d_image_denoising.h"
#include "app_viewer.h"

//DDG modified
#include "viewer/Viewer.h"

using namespace DDG;
using namespace std;

void generate_grid(size_t s, size_t f, string name);
void generate_grid_cylinder(const vertex_t & radio, vertex_t d_angle, size_t rings);
void main_testkeypoints();
double main_testkeycomponents();
void main_testholes();
void sampling_terrain(string file, size_t s, size_t K, string output = "patches_index", string planes = "normals");
float test_fastmarching(string file, size_t n_test = 10);
void main_test_fastmarching();
void main_test_holes();
void generate_grid_obtuse(const size_t & nr, const size_t & nc, const vertex_t & d = 1);
void main_solve_arma();

int main(int nargs, char ** args)
{
//	generate_grid_obtuse(81, 10);
//	generate_grid_cylinder(100, 1, 1000);
/*	mat A = {{2, 2, 0, 1}, {7, 3, 0, 1}, {4, 7, 0, 1}};
	A = A.t();
	vec x = {4.5, 2.5, 0, 1};
	vec a;
	a = solve(A, x);
	debug(a)
	debug(a > 0)
	debug(sum(a))
*/	
//	generate_grid(1000, 1, "tmp/grilla1m.off");
//	testkeycomponents prueba(50, 0.15);
//	prueba.one_test_fm("0001.null.0.off","");
//	if(nargs > 1) test_image_denoising(args[1]);
	
//	main_test_fastmarching();
	
	//sampling_terrain("terreno.off", 4, 64);
//	sampling_shape(nargs, args);
//	if(nargs == 2) sampling_shape(args[1]);
//	main_solve_arma();
//	distance_t radio = 0.04;

//	params_t params = { & K };

//	main_test_holes();
	viewer_main(nargs, args);
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

		printf("%15s & %10lu & %10lu & %10.3f ", shape->name().c_str(), old_n_vertices, shape->n_vertices() - old_n_vertices, quality  * 100 / n_faces);
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

void main_test_fastmarching()
{
	string filename;
	float time_s;

	while(cin >> filename)
	{
		debug(filename)
		time_s = test_fastmarching(filename, 1);
		
		filename = PATH_MDATA + filename + ".off";
		che * shape = new che_off(filename);
		vector<index_t> sources = {0};
/*
		float time;
		farthest_point_sampling_gpu(sources, time, shape, 1000, 0);
		debug(time)	
		

		map<index_t, index_t> deg;
		for(index_t d, v = 0; v < shape->n_vertices(); v++)
		{
			d = 0;
			for_star(he, shape, v) d++;
			if(shape->ot_evt(v) == NIL) d++;
			deg[d]++;
		}
	
		ofstream os("deg_hist");
		for(auto pp: deg)
			os << pp.first << " " << pp.second << endl;
		os.close();
		debug(system(("mv deg_hist " + PATH_TEST + "fastmarching/" + shape->name() + ".deg").c_str()))
*/
		
		index_t * rings = new index_t[shape->n_vertices()];
		index_t * sorted = new index_t[shape->n_vertices()];
		vector<index_t> limites;
		shape->sort_by_rings(rings, sorted, limites, {0});
		
		index_t * dist_rings = new index_t[limites.size() - 1];
		
		ofstream os("rings_dist");
		for(index_t i = 1; i < limites.size(); i++)
			os << i - 1 << " " << (dist_rings[i - 1] = limites[i] - limites[i - 1]) << endl;
		os.close();
	
		sort(dist_rings, dist_rings + limites.size() - 1);
		
		os.open("rings_dist_sort");
		for(index_t i = 1; i < limites.size(); i++)
			os << i - 1 << " " << dist_rings[i - 1] << endl;
		os.close();
		
		debug(system(("mv iter_error " + PATH_TEST + "fastmarching/fm_error_double_" + shape->name() + ".iter").c_str()))
		debug(system(("mv rings_dist " + PATH_TEST + "fastmarching/" + shape->name() + ".rings").c_str()))
		debug(system(("mv rings_dist_sort " + PATH_TEST + "fastmarching/" + shape->name() + ".s_rings").c_str()))
		
		delete shape;
		delete [] rings;
		delete [] sorted;
	}
}

distance_t * load_exact_geodesics(string file, size_t n, che * mesh, const index_t & s)
{
	ifstream is(PATH_TEST + ("exact_geodesics/" + file));
	distance_t * distances_e = new distance_t[n];
	
	if(is.good())
	{
		for(index_t i = 0; i < n; i++)
			is >> distances_e[i];
		is.close();
		return distances_e;
	}
	else
	{
		for(index_t i = 0; i < n; i++)
			distances_e[i] = *(mesh->gt(s) - mesh->gt(i));
	}

	return distances_e;
}

float test_fastmarching(string filename, size_t n_test)
{
	string file = filename + ".off";
//	filename = filename.substr(filename.find_last_of('/') + 1, filename.find_last_of('.') - filename.find_last_of('/') - 1);
	printf("%20s &", ("\\verb|" + filename + '|').c_str());
	file = PATH_MDATA + file;
	float time_c, time_p, time_s, time;

	che * shape = new che_off(file);
	vector<index_t> source = { 0 };
	printf("%12ld &", shape->n_vertices());
	
	distance_t * distances_e = load_exact_geodesics(filename, shape->n_vertices(), shape, source.front());
	
	index_t * rings = new index_t[shape->n_vertices()];
	index_t * sorted = new index_t[shape->n_vertices()];
	vector<index_t> limites;
	shape->sort_by_rings(rings, sorted, limites, source);
		
	distance_t * distances_c;
	
	time_c = 0;
	for(index_t t = 1; t < n_test; t++)
	{
		distances_c = parallel_fastmarching(shape, source.data(), source.size(), time, limites, sorted);
		time_c += time;
		delete [] distances_c;
	}

	//test iter error
	distances_c = parallel_fastmarching(shape, source.data(), source.size(), time, limites, sorted, 0, NULL, 1, distances_e);
	delete [] distances_c;
	// end test iter error
	
	distances_c = parallel_fastmarching(shape, source.data(), source.size(), time, limites, sorted);
	time_c += time;
	time_c /= n_test;

	time_s = 0;
	for(index_t t = 1; t < n_test; t++)
	{
		TIC(time) geodesics geodesic(shape, source); TOC(time)
		time_s += time;
	}
	
	TIC(time) geodesics geodesic(shape, source); TOC(time)
	time_s += time;
	time_s /= n_test;

	distance_t * distances_s = geodesic.distances;

	printf("%18.3f &%18.3f &%18.3f &", time_s, time_c, time_s/time_c);

	distance_t error_c = 0;
	distance_t error_s = 0;
	for(index_t v = 1; v < shape->n_vertices(); v++)
	{
		if(distances_e)
		{
			error_c += abs(distances_c[v] - distances_e[v]) / distances_e[v];
			error_s += abs(distances_s[v] - distances_e[v]) / distances_e[v];
		}
		else
		{
			distance_t d = *(shape->gt(0) - shape->gt(v));
			error_c += abs(distances_c[v] - d) / d;
			error_s += abs(distances_s[v] - d) / d;
		}
	}

	error_c /= shape->n_vertices() - 1;
	error_s /= shape->n_vertices() - 1;

	error_c *= 100;
	error_s *= 100;

	printf("%18.3e & ", error_s);
//	printf("%10s", error_s < error_c ? "\\color{blue}" : "\\color{red}");
	printf("%18.3e \\\\\\hline\n", error_c);
//	printf("%10s %18.3e \\\\\\hline\n", error_s < error_c ? "\\color{blue}" : "\\color{red}", error_s - error_c);

	delete shape;
	delete [] rings;
	delete [] sorted;
	delete [] distances_c;
	delete [] distances_e;

	return time_s;
}

void sampling_terrain(string file, size_t s, size_t K, string output, string planes)
{
	file = PATH_MDATA + file;
	che * shape_che = new che_off(file);
	size_t n =  sqrt(shape_che->n_vertices());
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
			geodesics fm(shape_che, source, K);

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
			comand =  comand + "mkdir " + PATH_KCS + tmp;
			system (comand.c_str());
			comand="";
			comand =  comand + "mkdir " + PATH_KPS + tmp;
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
		i = !((v - nvg) % (nc - 1)) ?  ((v - nvg) / (nc - 1)) * 2 * nc : i + 1;
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

void main_solve_arma()
{
	vertex v;
	mat A = { {2,7,4}, {2, 3, 7}, {0,0,0}, {1,1,1} };
	vec x = {4, 4, 0, 1};
	vec ans(&v[0], 3, false, true);
	cout<< v <<endl;
	ans = solve(A,x);
	cout<< v <<endl;
}

