#include "dictionary.h"

#include "sampling.h"
#include "mdict.h"
#include "che_poisson.h"
#include "che_fill_hole.h"

#include "viewer/viewer.h"

#include <cassert>
#include <CImg.h>
#include <fstream>

#ifndef CGAL_PATCH_DEFS
	#define CGAL_PATCH_DEFS
	#define CGAL_EIGEN3_ENABLED
	#define CGAL_USE_BOOST_PROGRAM_OPTIONS
	#define CGAL_USE_GMP
	#define DCGAL_USE_MPFR
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>


using namespace cimg_library;


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {

typedef real_t DFT;
typedef CGAL::Simple_cartesian<DFT> Data_Kernel;
typedef Data_Kernel::Point_3 DPoint;
typedef Data_Kernel::Vector_3 DVector;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;



size_t dictionary::L = 12;
size_t dictionary::K = 10;
size_t dictionary::T = 5;

dictionary::dictionary(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const real_t & _f, const bool & _learn, const bool & _d_plot):
					mesh(_mesh), phi_basis(_phi_basis), m(_m), M(_M), f(_f), learn(_learn), d_plot(_d_plot)
{
	A.eye(phi_basis->dim, m);
	dist = new real_t[mesh->n_vertices()]; 
}

dictionary::~dictionary()
{
	patch_t::del_index = true;
}

void dictionary::learning()
{
	gproshan_log(MDICT);

	string f_dict = tmp_file_path(mesh->name_size() + '_' + to_string(phi_basis->dim) + '_' + to_string(m) + '_' + to_string(f) + '_' + to_string(L) + ".dict");

	if(learn)
	{
		gproshan_log_var(f_dict);

		if(!A.load(f_dict))
		{
			A.eye(phi_basis->dim, m);
			A = normalise(A);
			gproshan_debug_var(phi_basis->radio);
			gproshan_debug_var(m);
			phi_basis->plot_atoms(A);
			KSVD(A, patches, L, K);
			A.save(f_dict);
		}
	}
	else A.eye(phi_basis->dim, m);
	gproshan_debug_var(phi_basis->radio);
	assert(A.n_rows == phi_basis->dim);
	assert(A.n_cols == m);
	if(d_plot)
	{
		phi_basis->plot_basis();
		phi_basis->plot_atoms(A);
	}
}

void dictionary::sparse_coding()
{
	gproshan_log(MDICT);
	
	vector<locval_t> locval;
	alpha = OMP_all(patches, phi_basis, A, L);
}

void dictionary::init_sampling()
{
	gproshan_log(MDICT);

	n_vertices = mesh->n_vertices();

	// load sampling
	if(M == 0)
	{
		M = mesh->n_vertices();
		phi_basis->radio = mesh->mean_edge();
	}
	else
	{
		sampling.reserve(M);
		if(!load_sampling(sampling, phi_basis->radio, mesh, M))
			cerr << "Failed to load sampling" << endl;
	}

	s_radio = phi_basis->radio;
	phi_basis->radio *= f;

	gproshan_debug_var(s_radio);
	gproshan_debug_var(phi_basis->radio);
}

void dictionary::load_curvatures(a_vec & curvatures)
{
	string f_curv = tmp_file_path(mesh->name_size()  + ".curv");
	string f_norm = tmp_file_path(mesh->name_size()  + ".n");

	if(! curvatures.load(f_curv))
	{
		curvatures.zeros(mesh->n_vertices());
		//real_t *mean_curvature = new real_t[mesh->n_vertices()];
		vector<index_t> points;

		map<size_t, char> non_rep;
		map<size_t, char>::iterator it;
		size_t d_fitting = 2;
		size_t d_monge = 2;
		size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;

		real_t min = INFINITY;
		real_t max = -INFINITY;
		a_mat normals;
		normals.zeros(3, mesh->n_vertices());

		for(index_t v = 0; v < mesh->n_vertices(); v++)
		{
			link_t linkv;
			mesh->link(linkv, v);
			for(const index_t & he: linkv)
			{
				link_t linku;
				const index_t & u = mesh->vt(he);
				mesh->link(linku, u);
				for(const index_t & he: linku)
				{
					it = non_rep.find(u);
					if(it == non_rep.end()) 
						points.push_back(u);
					
				}
			}
			assert(points.size() > min_points);
			vector<DPoint> in_points;
			in_points.reserve(points.size());
			for(const index_t & u: points)
				in_points.push_back(DPoint(mesh->gt(u).x, mesh->gt(u).y, mesh->gt(u).z));
			
			My_Monge_form monge_form;
			My_Monge_via_jet_fitting monge_fit;
			monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);

			vertex normal = mesh->normal(v);
			monge_form.comply_wrt_given_normal(DVector(normal.x, normal.y, normal.z));
			curvatures(v) = ( monge_form.principal_curvatures(0) + monge_form.principal_curvatures(0) ) / 2;


			normals(0, v) = monge_form.normal_direction()[0];
			normals(1, v) = monge_form.normal_direction()[1];
			normals(2, v) = monge_form.normal_direction()[2];
			//gproshan_debug_var(mean_curvature[v]);
			points.clear();
			non_rep.clear();
			
		}
		curvatures.save(f_curv);
		normals.save(f_norm);
	}
	gproshan_debug(curvatures ready);

}

void dictionary::load_features(vector<index_t> & v_feat, size_t & featsize)
{
	string f_feat = tmp_file_path(mesh->name()  + ".int");
	ifstream inp;
    inp.open(f_feat.c_str(), ifstream::in);
  
	size_t tam;
	index_t tmp;

	gproshan_debug_var(f_feat);
    if(inp.fail()){
		inp.clear(ios::failbit);
		// call the function using system
		//g++ -O3 *.cpp -lgsl -lCGAL -o harris3d
		//cmake -DCMAKE_BUILD_TYPE=Debug ..
		
		string command = "../../Harris3D-Cpp/harris3d " +   tmp_file_path(mesh->name()) + ".off" +  " ../tmp/example.prop"; 
		gproshan_debug_var(command);
		system(command.c_str()); 
		gproshan_debug(created);
		inp.close();
		inp.open(f_feat.c_str(), ifstream::in);
	}

	gproshan_debug(exists);
	inp>>featsize;
	//v_feat.resize(tam);
	for(int i=0; i<featsize; i++)
	{
		inp>>tmp;
		v_feat.push_back(tmp);
	}
	inp>>tam;
	for(int i=0; i<tam; i++)
	{
		inp>>tmp;
		v_feat.push_back(tmp);
	}

	inp.close();
}

void dictionary::init_patches(const bool & reset, const fmask_t & mask)
{
	gproshan_log(MDICT);

	if(reset)
	{
		patches.resize(M);
		patches_map.resize(n_vertices);

		#pragma omp parallel
		{
			index_t * toplevel = new index_t[n_vertices];

			#pragma omp for 
			for(index_t s = 0; s < M; s++)
			{
				index_t v = sample(s);
				patches[s].init(mesh, v, dictionary::T, phi_basis->radio, toplevel);
			}
			

			delete [] toplevel;
		}

		#ifndef NDEBUG
			size_t patch_avg_size = 0;
			size_t patch_min_size = NIL;
			size_t patch_max_size = 0;

			#pragma omp parallel for reduction(+: patch_avg_size)
			for(index_t s = 0; s < M; s++)
				patch_avg_size += patches[s].vertices.size();
			#pragma omp parallel for reduction(min: patch_min_size)
			for(index_t s = 0; s < M; s++)
				patch_min_size = min(patches[s].vertices.size(), patch_min_size);
			#pragma omp parallel for reduction(max: patch_max_size)
			for(index_t s = 0; s < M; s++)
				patch_max_size = max(patches[s].vertices.size(), patch_max_size);

			patch_avg_size /= M;
			gproshan_debug_var(patch_avg_size);
			gproshan_debug_var(patch_min_size);
			gproshan_debug_var(patch_max_size);
		#endif
	}

	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz(mesh, patches_map, s, mask);

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim);
		phi_basis->discrete(p.phi, p.xyz);
		p.phi = normalise(p.phi);
	}

/*	
#ifndef NDEBUG
	CImgList<real_t> imlist;
	for(index_t s = 0; s < M; s++)
		patches[s].save(phi_basis->radio, 16, imlist);
	imlist.save_ffmpeg_external("tmp/patches.mpg", 5);
#endif	

*/

	/*Saving Patches*/
/*
	ofstream os(tmp_file_path("patch-mat"));
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];
		p.save_z(os);
	}
	os.close();
	// DRAW NORMALS DEBUG
	for(index_t s = 0; s < M; s++)
	{
		viewer::vectors.push_back({patches[s].x(0), patches[s].x(1), patches[s].x(2)});
		a_vec r = patches[s].x + 0.02 * patches[s].normal();
		viewer::vectors.push_back({r(0), r(1), r(2)});
	}
	*/
}

real_t dictionary::mesh_reconstruction(const fmask_t & mask)
{
	gproshan_log(MDICT);

	assert(n_vertices == mesh->n_vertices());
	return mdict::mesh_reconstruction(mesh, M, phi_basis->get_radio(), patches, patches_map, A, alpha, dist);
}
void dictionary::update_alphas(a_mat & alpha, size_t threshold)
{
	size_t np_new = M - threshold;
	bool patches_covered[np_new];
	memset(patches_covered, 0, sizeof(patches_covered));
	size_t count = 0;

	// Choose the border patches using the threshold
	while(count < threshold)
	{
		#pragma omp parallel for
		for(index_t s = threshold; s < M; s++)
		{	

			if(!patches_covered[s-threshold])
			{	
				a_vec sum;
				sum.zeros();
				size_t c = 0;
				// Here updating alphas, we need a structure between patches and neighboor patches
				//We can simulate that structure by using patches map
				for(auto p: patches_map[s])
				{
					if(p.first < threshold || patches_covered[p.first-threshold])
					{
						sum += alpha.col(p.first);
					}	
					sum /= c;

				}
				alpha.col(s) = sum;
				patches_covered[s-threshold] = 1;
				count++;
			}	

		}
	}
	
	// update alphas of choosed patches
	// update the threshold
	// repeat until threshold reachs all patches
}

index_t dictionary::sample(const index_t & s)
{
	assert(s < M);
	if(sampling.size()) return sampling[s];
	return s;
}

const real_t & dictionary::operator[](const index_t & i) const
{
	assert(i < mesh->n_vertices());
	return dist[i];
}

void dictionary::draw_patches(index_t i)
{
	gproshan_debug_var(patches[i].vertices[0]);
	phi_basis->plot_patch(A*alpha.col(i),patches[i].xyz, patches[i].vertices[0]);
}

void dictionary::save_alpha(string file)
{
	alpha.save(file);
}

} // namespace gproshan::mdict

