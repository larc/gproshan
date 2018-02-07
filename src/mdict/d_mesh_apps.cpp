#include "d_mesh_apps.h"

#include "d_dict_learning.h"
#include "d_basis_cosine.h"
#include "d_basis_dct.h"
#include "sampling.h"
#include "che_fill_hole.h"
#include "che_poisson.h"
#include "viewer/viewer.h"

#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
#include <queue>
#include <set>
#include <functional>
#include <cassert>

// mesh dictionary learning and sparse coding namespace
namespace mdict {

void dictionary_learning_process(che * mesh, vector<index_t> & points, const size_t & freq, size_t & rt, const size_t & m, size_t & M, const distance_t & f, const index_t & pf, const bool & op_dict)
{
	size_t K = freq * rt;
	patch_t::del_index = false;
	bool all_points = M == 0;
	
	distance_t max_dist, radio;

	if(!all_points)
	{
		if(!load_sampling(points, radio, mesh, M)) return;
		max_dist = radio;
		assert(M == points.size());
	}
	else
	{
		M = mesh->n_vertices();
		max_dist = radio = 3 * mesh->mean_edge();
		debug_me(all vertices)
	}
	debug(M)
	
	auto p_vertex = [&](const index_t & s) -> index_t
	{
		if(all_points) return s;
		return points[s];
	};

	radio *= f;

	debug(radio)
	debug(max_dist)
	
	size_t min_nvp = 36; //MInimo numero de vertices por patche
	
	/* basis initialization */

	basis * phi_basis = new basis_dct(rt, radio);
//	basis * phi_basis = new basis_cosine(radio, rt, freq);
	phi_basis->plot_basis();	

	/***********************************************************************************************/
	
	size_t n_vertices = mesh->n_vertices();
		
/*	vector<index_t> * borders = fill_all_holes(mesh);
	if(!borders)
	{
		d_message(inpainting fail)
//		reconstruction = false;
		return;
	}

	delete [] borders;
	M = mesh->n_vertices();
	
	double time;

 	poisson(mesh, n_vertices, 2);*/
	
	double time;
	n_vertices = mesh->n_vertices();
	vector<patch_t> patches(M);
	vector<patches_map_t> patches_map(n_vertices);
	
	// Params ---------------------------------------------------------------------------------------
	// Gaussian
/*	
	vec cx(K); 
	vec cy(K); 
	vertex_t sigma = 2.5 * radio / sqrt(K);
	get_centers_gaussian(cx, cy, radio, K);

	params_t params = { & cx, & cy, & sigma };
*/
	// Cossine 


	// Init patches -------------------------------------------------------------------------------


	auto init_patches = [&]()
	{
		debug_me(init_patches)
		#pragma omp parallel for
		for(index_t s = 0; s < M; s++)
		{
			index_t v = p_vertex(s);

			patch_t & p = patches[s];
			geodesics fm(mesh, {v}, geodesics::FM, NIL, radio);
			p.n = fm.n_sorted_index();

			p.indexes = new index_t[p.n];
			fm.copy_sorted_index(p.indexes, p.n);
		}

		debug_me(init_patches)
		// problem with patches_map to be parallel
		for(index_t s = 0; s < M; s++)
		{
			patch_t & p = patches[s];
			p.reset_xyz(mesh, patches_map, s);
		}
		debug_me(init_patches)

		#pragma omp parallel for
		for(index_t s = 0; s < M; s++)
		{
			patch_t & p = patches[s];
			if(p.n > min_nvp)
			{
				jet_fit_directions(p);
			}
			else
				principal_curvatures(p, mesh);
			p.transform();
			p.phi.set_size(p.n, K);
			phi_basis->discrete(p.phi, p.xyz);
		}
		debug_me(init_patches)

		// ----------------------------------------------------------
		size_t patch_mean_size = 0;
		#pragma omp parallel for reduction(+: patch_mean_size)
		for(index_t s = 0; s < M; s++)
			patch_mean_size += patches[s].n;
		
		patch_mean_size /= M;
		debug(patch_mean_size)
		// ----------------------------------------------------------

	};
	
	d_message(init_patches)
	TIC(time)
	init_patches();
	TOC(time)
	debug(time)

	// Dictionary learning ------------------------------------------------------------------------

	size_t L = 10;	
	mat alpha(m, M, fill::zeros);
	mat A(K, m);
	A.eye();

	debug(K)
	debug(m)
	if(op_dict)
	{
		string fmesh_dict = "tmp/" + mesh->name_size() + '_' + to_string(K) + '_' + to_string(m) + ".a_dict";
	
		debug(fmesh_dict)
	
		if(!A.load(fmesh_dict))
		{
			A.eye(K, m);
			//A.random(K, m);
			
			d_message(Dictionary learning...)	
			TIC(time)
			KSVDT(A, patches, M, L);
			TOC(time)
			debug(time)
			A.save(fmesh_dict);
		}
	}
	
	debug(A.n_rows)
	debug(A.n_cols)
	phi_basis->plot_atoms(A);

	auto run_omp_all = [&]()
	{
		// Compute alphas -------------------------------------------------------------------------
		
		alpha.set_size(m, M);

		d_message(Compute alphas...)	
		TIC(time)
		OMP_all_patches_ksvt(alpha, A, patches, M, L);
		TOC(time)
		debug(time)
	};

	bool reconstruction = true;
		
	// Declare process functions ------------------------------------------------------------------
	
	auto denoising = [&]()
	{
		d_message(denoising)
		run_omp_all();
	};

	auto super_resolution = [&]()
	{
		run_omp_all();
		
		d_message(super_resolution)
		
		mesh->multiplicate_vertices();
		mesh->multiplicate_vertices();
	
		n_vertices = mesh->n_vertices();
		
		patches_map.clear();
		patches.resize(n_vertices);

		init_patches();
	};

	auto inpainting = [&]()
	{
		d_message(inpainting)
		
		vector<index_t> * borders = fill_all_holes(mesh);
		if(!borders)
		{
			d_message(inpainting fail)
			reconstruction = false;
			return;
		}

		delete [] borders;

 		poisson(mesh, n_vertices, 2);
		debug("before inpaiting")	
		size_t M_;
		if(all_points) M_ = mesh->n_vertices();
		else
		{
			parallel_farthest_point_sampling(points, mesh, mesh->n_vertices(), max_dist);
			M_ = points.size();
		}
		
		debug(M_)
		debug(M)

		patches.resize(M_);
		patches_map.resize(mesh->n_vertices());

		size_t size;
		for(index_t v, s = M; s < M_; s++)
		{
			v = p_vertex(s);
			patch_t & p = patches[s];

			geodesics fm(mesh, {v}, geodesics::FM, NIL, radio );
			size = fm.n_sorted_index();
			p.indexes = new index_t[size];
			fm.copy_sorted_index(p.indexes, size);
			
			p.n = size;
			p.reset_xyz(mesh, patches_map, s);

			if(p.n > min_nvp)
			{
				jet_fit_directions(p);
			}
			else
				principal_curvatures(p, mesh);
		//	principal_curvatures(p, mesh);
			p.transform();
			p.phi.set_size(p.n, K);
			phi_basis->discrete(p.phi, p.xyz);
		}

		M = M_;
	
		run_omp_all();
		
		n_vertices = mesh->n_vertices();
	};

	auto iterative_inpainting = [&]()
	{
		run_omp_all();
		
		d_message(iterative inpainting)
		
		size_t n_borders = mesh->n_borders();
		if(!n_borders) return;

		vector<index_t> * border_vertices = fill_all_holes(mesh);

		//with poisson
		poisson(mesh, n_vertices, 2);
		debug(mesh->n_vertices())

		//Contains index of patches which are borders
		set<index_t> border_patches;	
		
		for(index_t nb = 0; nb < n_borders; nb++)
		for(index_t b: border_vertices[nb])
		for(auto p_aux: patches_map[b])
			border_patches.insert(p_aux.first);

		delete [] border_vertices;
		
		debug(border_patches.size())
		
		index_t v;
		size_t size;
		vec tmp_alpha(m);
		tmp_alpha.zeros();

		patches_map.resize(mesh->n_vertices());

		for(index_t bp: border_patches)
		{
			//updating patch
			patch_t & p = patches[bp];
				
			v = p[0];
			geodesics fm(mesh, {v}, geodesics::FM, NIL, radio);

			delete [] p.indexes;
			p.indexes = new index_t[fm.n_sorted_index()];
		
			fm.copy_sorted_index(p.indexes, fm.n_sorted_index());
			size = fm.n_sorted_index();

			p.n = size;
			p.reset_xyz(mesh, patches_map, bp);

		//	principal_curvatures(p, mesh);
			if(p.n > min_nvp)
			{
				jet_fit_directions(p);
			}
			else
				principal_curvatures(p, mesh);

			p.transform();

			p.phi.set_size(p.n, K);
			phi_basis->discrete(p.phi, p.xyz);
		}
		
//
		// Second iteration
		size_t M_;
		if(all_points) M_ = mesh->n_vertices();
		else
		{
			debug(max_dist)	
			float time_g;
			//farthest_point_sampling_gpu(points, time_g, mesh, n_vertices, max_dist);
			debug(time_g)

			M_ = points.size();
		}

		patches.resize(M_);
		alpha.set_size(m, M_);
		patches_map.resize(mesh->n_vertices());

		debug(M_)	
		debug(M)

		index_t * levels = new index_t[M_];
		memset(levels, 0, sizeof(index_t) * M);
		memset(levels + M, 255, sizeof(index_t) * (M_ - M));


//
		// Including all patches to the mapping...
		debug_me(it inpaiting)

		index_t count = M;
		index_t level = 0;
		auto is_level = [&](const index_t & v, const index_t & level) -> bool
		{
			for(auto & m: patches_map[v])
				if(levels[m.first] < level)
					return true;

			return false;
		};
	
		while(count < M_)
		{
			level++;
			debug(level)
			for(index_t s = M; s < M_; s++)
			{
				v = p_vertex(s);
				patch_t & p = patches[s];

				if(levels[s] == NIL && is_level(v, level))
				{	
					levels[s] = level;
					count++;
					
					geodesics fm(mesh, {v}, geodesics::FM, NIL, radio );
					p.n = fm.n_sorted_index();
					p.indexes = new index_t[p.n];
					fm.copy_sorted_index(p.indexes, p.n);
						
					p.reset_xyz(mesh, patches_map, s);

			//		jet_fit_directions(p);
					if(p.n > min_nvp)
					{
						jet_fit_directions(p);
					}
					else
						principal_curvatures(p, mesh);
				//	principal_curvatures(p, mesh);
						
					p.transform();

					p.phi.set_size(p.n, K);
					phi_basis->discrete(p.phi, p.xyz);
				
				//	OMP_patch(alpha, A, s, p, L);	
					
					map<index_t, index_t> patches_alphas;

					for(index_t i = 0; i < p.n; i++)
					for(auto & pi: patches_map[p[i]])
				//		if(pi.first != s) patches_alphas[pi.first]++;
						if(levels[pi.first] < level) patches_alphas[pi.first]++;
					
					alpha.col(s).zeros();
					distance_t sum = 0;
					for(auto & pa: patches_alphas)
					{
						alpha.col(s) += pa.second * alpha.col(pa.first);
						sum += pa.second;
					}

					assert(sum > 0);
					alpha.col(s) /= sum;
				}
			}
			debug(count)
		}
		
		debug(count)
		debug(patches.size() - M)
		
		M = M_;

		n_vertices = mesh->n_vertices();

		delete [] levels;
	};


	auto non_local_inpainting = [&]()
	{
		run_omp_all();
		
		d_message(non local inpainting)
		
		size_t n_borders = mesh->n_borders();
		if(!n_borders) return;

		vector<index_t> * border_vertices = fill_all_holes(mesh);

		//with poisson
		poisson(mesh, n_vertices, 2);
		debug(mesh->n_vertices())

		//Contains index of patches which are borders
		set<index_t> border_patches;	
		
		for(index_t nb = 0; nb < n_borders; nb++)
		for(index_t b: border_vertices[nb])
		for(auto p_aux: patches_map[b])
			border_patches.insert(p_aux.first);

		delete [] border_vertices;
		
		debug(border_patches.size())
		
		index_t v;
		size_t size;
		vec tmp_alpha(m);
		tmp_alpha.zeros();

		patches_map.resize(mesh->n_vertices());

		for(index_t bp: border_patches)
		{
			//updating patch
			patch_t & p = patches[bp];
				
			v = p[0];
			geodesics fm(mesh, {v}, geodesics::FM, NIL, radio);

			delete [] p.indexes;
			p.indexes = new index_t[fm.n_sorted_index()];
		
			fm.copy_sorted_index(p.indexes, fm.n_sorted_index());
			size = fm.n_sorted_index();

			p.n = size;
			p.reset_xyz(mesh, patches_map, bp);

			if(p.n > min_nvp)
			{
				jet_fit_directions(p);
			}
			else
				principal_curvatures(p, mesh);
		//	jet_fit_directions(p);
		//	principal_curvatures(p, mesh);
			p.transform();

			p.phi.set_size(p.n, K);
			phi_basis->discrete(p.phi, p.xyz);
		}
		
//
		// Second iteration
		debug(max_dist)	
		float time_g;
		//farthest_point_sampling_gpu(points, time_g, mesh, n_vertices, max_dist);

		debug(time_g)
		
		patches.resize(points.size());
		alpha.set_size(m, points.size());

		debug(patches.size())	
		debug(points.size())

//		return;
		index_t * levels = new index_t[patches.size()];
		memset(levels, 0, sizeof(index_t) * M);
		memset(levels + M, 255, sizeof(index_t) * (patches.size() - M));


//
		// Including all patches to the mapping...
		debug_me(it inpaiting)


		distance_t d_min;
		index_t arg_min;

		for(index_t s = M; s < points.size(); s++)
		{
			v = p_vertex(s);		
			patch_t & p = patches[s];			
					
			geodesics fm(mesh, {v}, geodesics::FM, NIL, radio );
			p.n = fm.n_sorted_index();
			p.indexes = new index_t[p.n];
			fm.copy_sorted_index(p.indexes, p.n);
											
			p.reset_xyz(mesh, patches_map, s);

			if(p.n > min_nvp)
			{
				jet_fit_directions(p);
			}
			else
				principal_curvatures(p, mesh);
	//		jet_fit_directions(p);
		//	principal_curvatures(p, mesh);
											
			p.transform();

			p.phi.set_size(p.n, K);
			phi_basis->discrete(p.phi, p.xyz);
										
			OMP_patch(alpha, A, s, p, L);	
		
			arg_min = s;
			d_min = INFINITY;

			for(index_t i = 0; i < M; i++)
				if( norm(alpha.col(s) - alpha.col(i)) < d_min )
				{
					d_min = norm(alpha.col(s) - alpha.col(i));
					arg_min = i;
				}

	//		if(	norm (alpha.col(s) - alpha.col(arg_min)) != 0)
	//			debug("cambi alpha")
			alpha.col(s) = alpha.col(arg_min);
			
		}
	
		debug(patches.size() - M)
		
		M = patches.size();
		n_vertices = mesh->n_vertices();

		delete [] levels;
	};

	vector<function<void(void)> > processs_app = {
													denoising, 
													super_resolution,
													inpainting,
													iterative_inpainting,
													non_local_inpainting
												};

	
	// Call process functions --------------------------------------------------------------------
	
	if(pf < processs_app.size())
		processs_app[pf]();


	// Reconstruction -------------------------------------------------------------------------------

	if(reconstruction)
	{
		assert(n_vertices == mesh->n_vertices());

		d_message(Mesh reconstruction...)
		TIC(time)
		mesh_reconstruction(mesh, M, patches, patches_map, A, alpha);
		TOC(time)
		debug(time)
	}
	
	// Free vectors ---------------------------------------------------------------------------------
	
	points.clear();

	patch_t::del_index = true;
	debug_me(inpainting)
}

size_t sort_first_valid_vertex(index_t * indexes, const size_t & size, const vector<patches_map_t> & patches_map)
{
	index_t i = 0, f = size - 1;
	index_t tmp;

	while(i < f)
	{
		if(!patches_map[indexes[i]].size())
		{
			tmp = indexes[i];
			while(i < f && !patches_map[indexes[f]].size()) f--;
			
			if(patches_map[indexes[f]].size())
			{
				indexes[i] = indexes[f];
				indexes[f] = tmp;
				i++;
				f--;
			}
		}
		else i++;
	}

	return i;
}

void mesh_denoising(che * mesh, vector<index_t> & points, const size_t & freq, size_t & rt, const size_t & m, size_t & M, const distance_t & f, const bool & learn)
{
	dictionary_learning_process(mesh, points, freq, rt, m, M, f, 0, learn);
}

void mesh_inpaiting(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn)
{
	dictionary_learning_process(mesh, points, freq, rt, m, M, f, 2, learn);
}

void mesh_super_resolution(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn)
{
	dictionary_learning_process(mesh, points, freq, rt, m, M, f, 1, learn);
}

void mesh_iterative_inpaiting(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn)
{
	dictionary_learning_process(mesh, points, freq, rt, m, M, f, 3, learn);
}

} // mdict

