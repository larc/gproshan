#include "mdict/msparse_coding.h"

#include "geodesics/sampling.h"
#include "mdict/mdict.h"
#include "mesh/che_off.h"
#include "mesh/che_poisson.h"
#include "mesh/che_fill_hole.h"

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



size_t msparse_coding::L = 12;
size_t msparse_coding::K = 10;
size_t msparse_coding::T = 5;

msparse_coding::msparse_coding(che *const & _mesh, basis *const & _phi_basis, const params & p): mesh(_mesh), phi_basis(_phi_basis), m_params(p)
{
	A.eye(phi_basis->dim(), m_params.n_atoms);
	dist = new real_t[mesh->n_vertices()]; 
}

msparse_coding::~msparse_coding()
{
}

void msparse_coding::learning()
{
	gproshan_log(MDICT);

	string f_dict = tmp_file_path(mesh->name_size() + '_' + to_string(phi_basis->dim()) + '_' + to_string(m_params.n_atoms) + '_' + to_string(m_params.f) + '_' + to_string(L) + ".dict");

	if(m_params.learn)
	{
		gproshan_log_var(f_dict);

		if(!A.load(f_dict))
		{
			//A.eye(phi_basis->dim(), m);
			//initialize with some alpha
			//random
			//arma::uvec r_ind = arma::randi<arma::uvec>(m, arma::distr_param(0, M));
			//A = alpha.cols(r_ind);
			a_mat R, E, U, V;
			a_vec s;
			svd(U, s, V, alpha);
			gproshan_debug(svd done!);	
			A = U.cols(0, m_params.n_atoms);
			gproshan_debug(svd done!);	
			A = normalise(A);
			gproshan_debug_var(phi_basis->radio());
			gproshan_debug_var(m_params.n_atoms);
			//
			
			phi_basis->plot_atoms(A);
			KSVD(A, patches, L, K);
			phi_basis->plot_atoms(A);
			A.save(f_dict);
		}
	}
	else A.eye(phi_basis->dim(), m_params.n_atoms);
	gproshan_debug_var(phi_basis->radio());
	assert(A.n_rows == phi_basis->dim());
	assert(A.n_cols == m_params.n_atoms);
	if(m_params.plot)
	{
		phi_basis->plot_basis();
		phi_basis->plot_atoms(A);
	}
}

void msparse_coding::sparse_coding()
{
	gproshan_log(MDICT);
	
	vector<locval_t> locval;
	alpha = OMP_all(patches, phi_basis, A, L);
}

void msparse_coding::init_sampling()
{
	gproshan_log(MDICT);

	n_vertices = mesh->n_vertices();

	// load sampling
	if(m_params.n_patches == 0)
	{
		m_params.n_patches = mesh->n_vertices();
		phi_basis->radio() = mesh->mean_edge();
	}
	else
	{
		sampling.reserve(m_params.n_patches);
		if(!load_sampling(sampling, phi_basis->radio(), mesh, m_params.n_patches))
			cerr << "Failed to load sampling" << endl;
	}

	s_radio = phi_basis->radio();
	phi_basis->radio() *= m_params.f;

	gproshan_debug_var(s_radio);
	gproshan_debug_var(phi_basis->radio());
}

void msparse_coding::load_features(vector<index_t> & v_feat, size_t & featsize)
{
	string f_feat = tmp_file_path(mesh->name() + ".int");
	ifstream inp;
	inp.open(f_feat.c_str(), ifstream::in);
	
	size_t tam;
	index_t tmp;

	gproshan_debug_var(f_feat);
	if(inp.fail())
	{
		inp.clear(ios::failbit);
		// call the function using system
		//g++ -O3 *.cpp -lgsl -lCGAL -o harris3d
		//cmake -DCMAKE_BUILD_TYPE=Debug ..
		
		string command = "../../harris3d/harris3d " + mesh->filename() + " " + f_feat + " " + tmp_file_path("example.prop"); 
		gproshan_debug_var(command);
		system(command.c_str()); 
		gproshan_debug(created);
		inp.close();
		inp.open(f_feat.c_str(), ifstream::in);
	}

	gproshan_debug(exists);
	inp >> featsize;
	for(index_t i = 0; i < featsize; i++)
	{
		inp>>tmp;
		v_feat.push_back(tmp);
	}

	inp >> tam;
	for(index_t i = 0; i < tam; i++)
	{
		inp >> tmp;
		v_feat.push_back(tmp);
	}

	inp.close();
}

void msparse_coding::init_patches(const bool & reset, const fmask_t & mask)
{
	gproshan_log(MDICT);

	if(reset)
	{
		patches.resize(m_params.n_patches);
		patches_map.resize(n_vertices);

		#pragma omp parallel
		{
			index_t * toplevel = new index_t[n_vertices];

			#pragma omp for 
			for(index_t s = 0; s < m_params.n_patches; s++)
			{
				index_t v = sample(s);
				patches[s].init(mesh, v, msparse_coding::T, phi_basis->radio(), toplevel);
			}
			

			delete [] toplevel;
		}

		#ifndef NDEBUG
			size_t patch_avg_size = 0;
			size_t patch_min_size = NIL;
			size_t patch_max_size = 0;

			#pragma omp parallel for reduction(+: patch_avg_size)
			for(index_t s = 0; s < m_params.n_patches; s++)
				patch_avg_size += patches[s].vertices.size();
			#pragma omp parallel for reduction(min: patch_min_size)
			for(index_t s = 0; s < m_params.n_patches; s++)
				patch_min_size = min(patches[s].vertices.size(), patch_min_size);
			#pragma omp parallel for reduction(max: patch_max_size)
			for(index_t s = 0; s < m_params.n_patches; s++)
				patch_max_size = max(patches[s].vertices.size(), patch_max_size);

			patch_avg_size /= m_params.n_patches;
			gproshan_debug_var(patch_avg_size);
			gproshan_debug_var(patch_min_size);
			gproshan_debug_var(patch_max_size);
		#endif
	}

	for(index_t s = 0; s < m_params.n_patches; s++)
		patches[s].reset_xyz(mesh, patches_map, s, mask);

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
		p.phi = normalise(p.phi);
	}

/*	
#ifndef NDEBUG
	CImgList<real_t> imlist;
	for(index_t s = 0; s < m_params.n_patches; s++)
		patches[s].save(phi_basis->radio(), 16, imlist);
	imlist.save_ffmpeg_external("tmp/patches.mpg", 5);
#endif	

*/

	/*Saving Patches*/
/*
	ofstream os(tmp_file_path("patch-mat"));
	for(index_t s = 0; s < m_params.n_patches; s++)
	{
		patch & p = patches[s];
		p.save_z(os);
	}
	os.close();
	// DRAW NORMALS DEBUG
	for(index_t s = 0; s < m_params.n_patches; s++)
	{
		viewer::vectors.push_back({patches[s].x(0), patches[s].x(1), patches[s].x(2)});
		a_vec r = patches[s].x + 0.02 * patches[s].normal();
		viewer::vectors.push_back({r(0), r(1), r(2)});
	}
	*/
}

real_t msparse_coding::mesh_reconstruction(const fmask_t & mask)
{
	gproshan_log(MDICT);

	assert(n_vertices == mesh->n_vertices());

	a_mat V(3, mesh->n_vertices(), arma::fill::zeros);

	patches_error.resize(m_params.n_patches);

	#pragma omp parallel for
	for(index_t p = 0; p < m_params.n_patches; p++)
	{
		patch & rp = patches[p];
		
		a_vec x = rp.phi * A * alpha.col(p);
			
		patches_error[p] = { accu(abs(x - rp.xyz.row(2).t())) / rp.vertices.size(), p };

		rp.xyz.row(2) = x.t();
	}

	sort(patches_error.begin(), patches_error.end());
	
	fprintf(stderr, "error %16s%16s\n", "best", "worst");
	for(index_t i = 0; i < 10; i++)
	{
		const index_t & best = patches_error[i].second;
		const index_t & worst = patches_error[m_params.n_patches - i - 1].second;
		
		fprintf(stderr, "%5d:%8u>%8u%8u>%8u\n", i, best, draw_patches(best), worst, draw_patches(worst));
	}
	
	#pragma omp parallel for
	for(index_t p = 0; p < m_params.n_patches; p++)
	{
		patch & rp = patches[p];
		rp.iscale_xyz(phi_basis->radio());
		rp.itransform();
	}

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		// simple means vertex
		if(patches_map[v].size() && (!mask || mask(v)))
		{
			a_vec mv = arma::zeros(3);		
			for(auto p: patches_map[v])
				mv += patches[p.first].xyz.col(p.second);

			V.col(v) = mv / patches_map[v].size();
		}
		else
		{
			V(0, v) = mesh->gt(v).x;
			V(1, v) = mesh->gt(v).y;
			V(2, v) = mesh->gt(v).z;
		}
	}

	// ------------------------------------------------------------------------

	vertex * new_vertices = (vertex *) V.memptr();

	real_t error = 0;
	real_t max_error = -1;

	#pragma omp parallel for reduction(+: error) reduction(max: max_error)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		dist[v] = *(new_vertices[v] - mesh->get_vertex(v));
		error += dist[v];
		max_error = max(max_error, dist[v]);
	}
		
	error /= mesh->n_vertices();

	gproshan_debug_var(mesh->n_vertices());
	gproshan_debug_var(error);
	gproshan_debug_var(max_error);

	mesh->set_vertices(new_vertices, mesh->n_vertices());
	che_off::write_file(mesh, "../tmp/recon_mesh");
	
	return max_error;
}

void msparse_coding::update_alphas(a_mat & alpha, size_t threshold)
{
	size_t np_new = m_params.n_patches - threshold;
	bool patches_covered[np_new];
	memset(patches_covered, 0, sizeof(patches_covered));
	size_t count = 0;

	// Choose the border patches using the threshold
	while(count < threshold)
	{
		#pragma omp parallel for
		for(index_t s = threshold; s < m_params.n_patches; s++)
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

index_t msparse_coding::sample(const index_t & s)
{
	assert(s < m_params.n_patches);
	if(sampling.size()) return sampling[s];
	return s;
}

const real_t & msparse_coding::operator[](const index_t & i) const
{
	assert(i < mesh->n_vertices());
	return dist[i];
}

const index_t & msparse_coding::draw_patches(const index_t & p)
{
	phi_basis->plot_patch(A * alpha.col(p), patches[p].xyz, patches[p].vertices[0]);
	return patches[p].vertices[0];
}

void msparse_coding::save_alpha(string file)
{
	alpha.save(file);
}


} // namespace gproshan::mdict

