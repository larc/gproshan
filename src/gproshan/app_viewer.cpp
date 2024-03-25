#include <gproshan/app_viewer.h>

#include <gproshan/scenes/scene.h>
#include <gproshan/pointcloud/knn.h>

#include <random>
#include <queue>


using namespace gproshan::mdict;


// geometry processing and shape analysis framework
namespace gproshan {


app_viewer::app_viewer(const char * title, const int width, const int height): viewer(title, width, height) {}

app_viewer::~app_viewer()
{
	for(auto & m: meshes)
		delete *m;
}

int app_viewer::main(int nargs, const char ** args)
{
	if(nargs < 2)
	{
		printf("%s [mesh_paths.(off,obj,ply)]\n", args[0]);
		return 0;
	}

	TIC(time)
	for(int i = 1; i < nargs; ++i)
		add_mesh(che::load_mesh(args[i]));
	TOC(time)
	update_status_message("meshes loaded in %.3fs", time);

	init();
	run();

	return 0;
}

void app_viewer::init()
{
	add_menu("Point Cloud",
	{
		add_process("KNN", process_knn),
		add_process("Compute Normals", process_compute_normals)
	});

	add_menu("Scenes",
	{
		add_process("Scan Scene", process_simulate_scanner),
		add_process("RT Scatter Points", process_scatter)
	});

	add_menu("Geometry",
	{
		add_process("Sampling 4points", process_sampling_4points),
		add_process("2D Convex Hull", process_convex_hull, GLFW_KEY_H),
		add_process("Connected Components", process_connected_components, GLFW_KEY_O),
		add_process("Gaussian curvature", process_gaussian_curvature, GLFW_KEY_K),
		add_process("Edge Collapse", process_edge_collapse, GLFW_KEY_Q),
		add_process("Multiplicate", process_multiplicate_vertices, GLFW_KEY_M),
		add_process("Delete vertices", process_delete_vertices, GLFW_KEY_SLASH),
		add_process("Delete non-manifold vertices", process_delete_non_manifold_vertices, GLFW_KEY_MINUS)
	});

	add_menu("Fairing",
	{
		add_process("Fairing Taubin", process_fairing_taubin, GLFW_KEY_F),
		add_process("Fairing Spectral", process_fairing_spectral, GLFW_KEY_E)
	});

	add_menu("Geodesics",
	{
		add_process("Geodesics", process_geodesics, GLFW_KEY_G),
		add_process("Geodesic Farthest Point Sampling", process_farthest_point_sampling, GLFW_KEY_S),
		add_process("Geodesic Voronoi", process_voronoi, GLFW_KEY_V),
		add_process("Toplesets", process_compute_toplesets, GLFW_KEY_T)
	});

	add_menu("Sparse Coding",
	{
		add_process("Mesh Sparse Coding", process_msparse_coding, GLFW_KEY_U),
		add_process("MDICT Patch", process_mdict_patch, GLFW_KEY_J),
		add_process("MDICT Mask", process_mask, GLFW_KEY_D),
		add_process("PC reconstruction", process_pc_reconstruction, GLFW_KEY_L)
	});

	add_menu("Features",
	{
		add_process("Eigenfunctions", process_eigenfuntions, GLFW_KEY_2),
		add_process("Descriptors", process_descriptor_heatmap, GLFW_KEY_3),
		add_process("Key Points", process_key_points, GLFW_KEY_4),
		add_process("Key Components", process_key_components, GLFW_KEY_5)
	});

	add_menu("Hole Filling",
	{
		add_process("Poisson: Membrane surface", process_poisson_laplacian_1, GLFW_KEY_X),
		add_process("Poisson: Thin-plate surface", process_poisson_laplacian_2, GLFW_KEY_Y),
		add_process("Poisson: Minimum variation surface", process_poisson_laplacian_3, GLFW_KEY_Z),
		add_process("Fill hole: planar mesh", process_fill_holes, GLFW_KEY_6),
		add_process("Fill hole: biharmonic splines", process_fill_holes_biharmonic_splines, GLFW_KEY_7)
	});

	add_menu("Others",
	{
		add_process("Select multiple vertices", process_select_multiple, GLFW_KEY_SEMICOLON),
		add_process("Threshold", process_threshold, GLFW_KEY_BACKSLASH),
		add_process("Noise", process_noise, GLFW_KEY_N),
		add_process("Black noise", process_black_noise, GLFW_KEY_P)
	});
}


// Point Cloud

bool app_viewer::process_knn(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static int alg = 0;
	static int k = 9;
	ImGui::Combo("algorithm", &alg, "select\0grid\0k3tree\0\0");
	ImGui::InputInt("k", &k);

	if(ImGui::Button("Run"))
	{
		auto query = mesh.selected;
		if(!size(query)) query.push_back(0);

		mesh.selected.clear();

		switch(alg)
		{
			case 1:
			{
				knn::grid grid(&mesh->point(0), mesh->n_vertices, mesh.model_mat);
				for(const index_t p: query)
				{
					for(const index_t v: grid(mesh.model_mat * (mesh->point(p), 1), k))
						mesh.selected.push_back(v);
				}
			}
			break;

			case 2:
			{
				knn::k3tree k3tree(&mesh->point(0), mesh->n_vertices, k);
				for(const index_t p: query)
				{
					const int * result = k3tree(p);
					for(int i = 0; i < k; ++i)
						mesh.selected.push_back(result[i]);
				}
			}
			break;

			default: break;
		}
	}

	return true;
}

bool app_viewer::process_compute_normals(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	gproshan_log_var(mesh->n_vertices);

	// TODO

	return false;
}


// Scenes

bool app_viewer::process_simulate_scanner(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static size_t n_rows = 5000;
	static size_t n_cols = 2000;
	static const size_t n_min = 100;
	static const size_t n_max = 10000;
	static vertex cam_pos = {0, 0, 0};
	static bool dist_error = true;
	static bool save_jpg = false;

	if(size(view->sphere_points) != 1)
	{
		view->sphere_points.clear();
		view->sphere_points.push_back(cam_pos);
	}

	ImGui::SliderScalar("n_rows", ImGuiDataType_U64, &n_rows, &n_min, &n_max);
	ImGui::SliderScalar("n_cols", ImGuiDataType_U64, &n_cols, &n_min, &n_max);
	if(	ImGui::InputFloat("cam_pos.x", &cam_pos.x(), 0.001, 2, "%.3lf") ||
		ImGui::InputFloat("cam_pos.y", &cam_pos.y(), 0.001, 2, "%.3lf") ||
		ImGui::InputFloat("cam_pos.z", &cam_pos.z(), 0.001, 2, "%.3lf") )
	{
		view->sphere_points[0] = cam_pos;
	}

	ImGui::Checkbox("dist_error", &dist_error);
	ImGui::Checkbox("save_jpg", &save_jpg);

	if(ImGui::Button("Scan"))
	{
		che * ptx_mesh = save_jpg	? scanner_ptx_jpg(mesh.rt_embree, n_rows, n_cols, cam_pos, mesh->filename)
									: scanner_ptx(mesh.rt_embree, n_rows, n_cols, cam_pos, dist_error);

		che_ptx::write_file(ptx_mesh, mesh->filename, n_rows, n_cols);
		ptx_mesh->filename = mesh->filename + ".ptx";
		view->add_mesh(ptx_mesh);
	}

	return true;
}

bool app_viewer::process_scatter(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;

	rt::eval_hit h;
	std::vector<vertex> scatter(100);

	rt::random<float> rnd(0xABCDEF);
	for(vertex & v: scatter)
		h.scatter_diffuse(v, rnd);

	view->add_mesh(new che(scatter.data(), scatter.size(), nullptr, 0));

	return false;
}


// Geometry
bool app_viewer::process_sampling_4points(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static size_t n = 10;
	static std::vector<vertex> points;

	ImGui::InputScalar("resolution", ImGuiDataType_U64, &n);

	if(ImGui::Button("Clean")) points.clear();

	ImGui::SameLine();
	if(size(mesh.selected) > 3)
	{
		if(ImGui::Button("Add Samples"))
		{
			std::vector<vertex> samples = sampling_4points(n, mesh.selected_point(0),
															mesh.selected_point(1),
															mesh.selected_point(2),
															mesh.selected_point(3)
															);

			for(auto & p: samples)
				points.push_back(p);

			mesh.selected.clear();
		}
	}
	else ImGui::Text("Select 4 points.");

	ImGui::Text("Total added points: %lu", size(points));
	if(ImGui::Button("Create Point Cloud"))
		view->add_mesh(new che(points.data(), size(points), nullptr, 0));

	return true;
}

bool app_viewer::process_convex_hull(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	convex_hull ch(&mesh->point(0), mesh->n_vertices);
	mesh.selected = ch;

	return false;
}

bool app_viewer::process_connected_components(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	float * label = &mesh->heatmap(0);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		label[v] = -1;

	int nc = 0;
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(label[v] < 0)
		{
			++nc;

			std::queue<index_t> q;
			q.push(v);
			label[v] = nc;

			while(!q.empty())
			{
				for(const index_t v: mesh->link(q.front()))
					if(label[v] < 0)
					{
						label[v] = nc;
						q.push(v);
					}

				q.pop();
			}
		}
	}

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		label[v] /= nc;

	view->update_status_message("found %d connected components", nc);
	mesh.update_vbo_heatmap();

	return false;
}

bool app_viewer::process_gaussian_curvature(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	float g, g_max = -INFINITY, g_min = INFINITY;
	vertex a, b;

	arma::fvec gv(mesh->n_vertices);

	#pragma omp parallel for private(g, a, b) reduction(max: g_max) reduction(min: g_min)
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		g = 0;
		for(const index_t he: mesh->star(v))
		{
			a = mesh->vertex_he(he_next(he)) - mesh->point(v);
			b = mesh->vertex_he(he_prev(he)) - mesh->point(v);
			g += acos(dot(a, b) / (norm(a) * norm(b)));
		}
		//gv(v) = (2 * M_PI - g) / mesh->area_vertex(v);
		gv(v) = mesh->mean_curvature(v);

		g_max = std::max(g_max, gv(v));
		g_min = std::min(g_min, gv(v));
	}

	g = g_max - g_min;
	gproshan_log_var(g);
	gproshan_log_var(g_min);
	gproshan_log_var(g_max);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		gv(v) = (gv(v) + g_min) / g;

	float gm = mean(gv);
	float gs = var(gv);

	gproshan_debug_var(gm);
	gproshan_debug_var(gs);

	auto f = [&](float x, float a = 4) -> float
	{
		if(x < gm - a * gs) return 0;
		if(x > gm + a * gs) return 1;
		return (x - gm) / (2 * a * gs) + 0.5;
	};

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		gv(v) = f(gv(v));

	mesh.update_vbo_heatmap(gv.memptr());

	return false;
}

bool app_viewer::process_edge_collapse(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	index_t levels;
	std::cin >> levels;

	TIC(view->time) simplification sampling(mesh, levels); TOC(view->time)
	gproshan_debug_var(view->time);

	//if(view->n_meshes < 2)
	//	view->add_mesh(new che(*mesh));

	return false;
}

bool app_viewer::process_multiplicate_vertices(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	mesh->multiplicate_vertices();
	mesh->update_normals();
	mesh.update();

	mesh.log_info();

	return false;
}

bool app_viewer::process_delete_vertices(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	if(!size(mesh.selected)) return true;

	gproshan_debug(removing vertex);
	mesh->remove_vertices(mesh.selected);
	mesh.selected.clear();
	gproshan_debug(removing vertex);

	return false;
}

bool app_viewer::process_delete_non_manifold_vertices(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	gproshan_debug(removing vertex);
	mesh->remove_non_manifold_vertices();
	gproshan_debug(removing vertex);

	return false;
}


// Fairing

bool app_viewer::process_fairing_spectral(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static std::vector<vertex> vertices;
	static size_t min_neigs = 1;
	static size_t max_neigs = std::min(1000lu, mesh->n_vertices);
	static fairing_spectral fair(max_neigs);

	if(ImGui::SliderScalar("n_eigs", ImGuiDataType_U64, &fair.n_eigs, &min_neigs, &max_neigs))
	{
		if(size(vertices) != mesh->n_vertices)
		{
			vertices.resize(mesh->n_vertices);
			memcpy(vertices.data(), &mesh->point(0), mesh->n_vertices * sizeof(vertex));
		}
		else mesh->update_vertices(vertices.data());

		fair.run(mesh);

		mesh->update_vertices(fair.new_vertices());
		mesh->update_normals();

		mesh.update_vbo_geometry();
		mesh.update_vbo_normal();
	}

	return true;
}

bool app_viewer::process_fairing_taubin(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static std::vector<vertex> vertices;
	static fairing_taubin fair(0);

	if(ImGui::InputFloat("step", &fair.step, 0.001))
	{
		if(size(vertices) != mesh->n_vertices)
		{
			vertices.resize(mesh->n_vertices);
			memcpy(vertices.data(), &mesh->point(0), mesh->n_vertices * sizeof(vertex));
		}
		else mesh->update_vertices(vertices.data());

		fair.run(mesh);

		mesh->update_vertices(fair.new_vertices());
		mesh->update_normals();

		mesh.update_vbo_geometry();
		mesh.update_vbo_normal();
	}

	return true;
}


// Geodesics

bool app_viewer::process_geodesics(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static geodesics::params params;

#ifdef GPROSHAN_CUDA
	ImGui::Combo("alg", (int *) &params.alg, "FM\0PTP_CPU\0HEAT_METHOD\0PTP_GPU\0HEAT_METHOD_GPU\0\0");
#else
	ImGui::Combo("alg", (int *) &params.alg, "FM\0PTP_CPU\0HEAT_METHOD\0\0");
#endif // GPROSHAN_CUDA

	if(ImGui::Button("Run"))
	{
		if(!size(mesh.selected))
			mesh.selected.push_back(0);

		params.dist_alloc = &mesh->heatmap(0);

		TIC(view->time)
			geodesics dm(mesh, mesh.selected, params);
		TOC(view->time)
		view->update_status_message("geodesics time: %.3fs", view->time);

		dm.normalize();
		mesh.update_vbo_heatmap();
	}

	return true;
}

bool app_viewer::process_farthest_point_sampling(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static int n = 10;
	static float radio;

	ImGui::SliderInt("samples", &n, 1, mesh->n_vertices / 6);
	ImGui::Text("radio: %.3f", radio);

	if(ImGui::Button("Run"))
	{
		TIC(view->time)
		load_sampling(mesh.selected, radio, mesh, n);
		TOC(view->time)
		gproshan_log_var(view->time);
	}

	return true;
}

bool app_viewer::process_voronoi(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	geodesics::params params;
	params.cluster = true;

#ifdef GPROSHAN_CUDA
	params.alg = geodesics::PTP_GPU;
#endif

	TIC(view->time)
	geodesics ptp(mesh, mesh.selected, params);
	TOC(view->time)

	gproshan_log_var(view->time);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		mesh->heatmap(i) = ptp.clusters[i];
		mesh->heatmap(i) /= size(mesh.selected) + 1;
	}

	mesh.update_vbo_heatmap();

	return false;
}

bool app_viewer::process_compute_toplesets(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	if(!size(mesh.selected))
		mesh.selected.push_back(0);

	index_t * toplesets = new index_t[mesh->n_vertices];
	index_t * sorted = new index_t[mesh->n_vertices];
	std::vector<index_t> limites;
	mesh->compute_toplesets(toplesets, sorted, limites, mesh.selected);

	size_t n_toplesets = size(limites) - 1;

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(toplesets[v] < n_toplesets)
			mesh->heatmap(v) = float(toplesets[v]) / (n_toplesets);
	}

	mesh.update_vbo_heatmap();

	gproshan_debug_var(n_toplesets);

	delete [] toplesets;
	delete [] sorted;

	return false;
}


// Mesh Sparse Coding

bool app_viewer::process_msparse_coding(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static msparse_coding::params params;
	static size_t n = 12;

	assert(sizeof(ImGuiDataType_U64) != sizeof(size_t));

	ImGui::InputFloat("nyquist_factor", &patch::nyquist_factor, 0.01, 0.01, "%.2lf");
	ImGui::InputScalar("basis", ImGuiDataType_U64, &n);
	ImGui::InputScalar("atoms", ImGuiDataType_U64, &params.n_atoms);
	ImGui::InputFloat("delta", &params.delta, 0.001, 0.1, "%.3lf");
	ImGui::InputFloat("proj_thres", &params.sum_thres, 1.001, 0.1, "%.6lf");
	ImGui::InputFloat("area_thres", &params.area_thres, 0.001, 0.1, "%6lf");
	ImGui::Checkbox("learn", &params.learn);

	if(ImGui::Button("Run"))
	{
		basis_dct phi(n);
		msparse_coding msc(mesh, &phi, params);

		float max_error = msc.execute();
		gproshan_log_var(max_error);

		mesh->update_heatmap(msc);
		mesh->update_normals();
	}

	return true;
}

bool app_viewer::process_mdict_patch(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	TIC(view->time)
	index_t * toplevel = new index_t[mesh->n_vertices];
	size_t avg_nvp = 0;

	vertex vdir;
	patch p;
	float mean_edge = mesh->mean_edge();
	for(auto & v: mesh.selected)
	{
		p.init(mesh, v, msparse_coding::T, msparse_coding::T * mean_edge, toplevel);
		for(auto & u: p.vertices)
			mesh->heatmap(u) = 1;

		vdir = {p.T(0, 0), p.T(0, 1), p.T(0, 2)};
		view->vectors.push_back(mesh->point(v));
		view->vectors.push_back(mesh->point(v) + 3 * mean_edge * vdir);

		vdir = {p.T(1, 0), p.T(1, 1), p.T(1, 2)};
		view->vectors.push_back(mesh->point(v));
		view->vectors.push_back(mesh->point(v) + 3 * mean_edge * vdir);

		vdir = {p.T(2, 0), p.T(2, 1), p.T(2, 2)};
		view->vectors.push_back(mesh->point(v));
		view->vectors.push_back(mesh->point(v) + 3 * mean_edge * vdir);

		view->vectors.push_back(mesh->point(v));
		view->vectors.push_back(mesh->point(v) + 3 * mean_edge * mesh->normal(v));

		avg_nvp += size(p.vertices);
	}

	avg_nvp /= size(mesh.selected);
	gproshan_log_var(avg_nvp);

	delete [] toplevel;
	TOC(view->time)
	gproshan_log_var(view->time);

	return false;
}

bool app_viewer::process_mask(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static msparse_coding::params params;
	static size_t n = 12;

	assert(sizeof(ImGuiDataType_U64) != sizeof(size_t));

	ImGui::InputScalar("basis", ImGuiDataType_U64, &n);
	ImGui::InputScalar("atoms", ImGuiDataType_U64, &params.n_atoms);
	ImGui::InputFloat("delta", &params.delta, 0.001, 0.1, "%.3lf");
	ImGui::InputFloat("proj_thres", &params.sum_thres, 1.001, 0.1, "%.6lf");
	ImGui::InputFloat("area_thres", &params.area_thres, 0.001, 0.1, "%6lf");
	ImGui::Checkbox("learn", &params.learn);

	if(ImGui::Button("Run"))
	{
		basis_dct phi(n);
		msparse_coding msc(mesh, &phi, params);

		msc.init_radial_feature_patches();
		//dict.init_voronoi_patches();
		mesh->update_heatmap(msc);
		std::string f_points = tmp_file_path(std::string(msc) + ".rsampl");

		arma::fvec points_out;
		gproshan_debug_var(f_points);
		points_out.load(f_points);
		gproshan_debug_var(size(points_out));

		for(index_t i = 0; i < points_out.size(); ++i)
			mesh.selected.push_back(points_out(i));

		mesh->update_normals();
	}

	return true;
}

bool app_viewer::process_pc_reconstruction(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static msparse_coding::params params;
	static size_t n = 12;
	static float percentage_size = 100;
	static float radio_factor = 1;

	assert(sizeof(ImGuiDataType_U64) != sizeof(size_t));

	ImGui::InputScalar("basis", ImGuiDataType_U64, &n);
	ImGui::InputScalar("atoms", ImGuiDataType_U64, &params.n_atoms);
	ImGui::InputFloat("delta", &params.delta, 0.001, 0.1, "%.3lf");
	ImGui::InputFloat("proj_thres", &params.sum_thres, 1.001, 0.1, "%.6lf");
	ImGui::InputFloat("area_thres", &params.area_thres, 0.001, 0.1, "%6lf");
	ImGui::Checkbox("learn", &params.learn);

	ImGui::InputFloat("percentage_size", &percentage_size, 100, 10, "%.3f");
	ImGui::InputFloat("radio_factor", &radio_factor, 1, 0.1, "%.3f");

	if(ImGui::Button("Run"))
	{
		basis_dct phi(n);
		msparse_coding msc(mesh, &phi, params);

		view->add_mesh(msc.point_cloud_reconstruction(percentage_size, radio_factor));
	}

	return true;
}


// Features

bool app_viewer::process_eigenfuntions(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static size_t n_eigs = 20;
	ImGui::InputScalar("n_eigs", ImGuiDataType_U64, &n_eigs);

	if(ImGui::Button("Run"))
	{
		arma::sp_fmat L, A;
		arma::fvec eigval;
		arma::fmat eigvec;

		TIC(view->time)
			n_eigs = eigs_laplacian(mesh, eigval, eigvec, L, A, n_eigs);
		TOC(view->time)

		view->update_status_message("computing %lu eigs in %.3fs", n_eigs, view->time);

		for(index_t k = 0; k < n_eigs; ++k)
		{
			if(k)
			{
				if(!view->add_mesh(new che(*mesh)))
					break;
			}

			view->idx_selected_mesh = k;
			che_viewer & mesh = view->selected_mesh();

			eigvec.col(k) -= eigvec.col(k).min();
			eigvec.col(k) /= eigvec.col(k).max();

			#pragma omp parallel for
			for(index_t v = 0; v < mesh->n_vertices; ++v)
				mesh->heatmap(v) = eigvec(v, k);

			mesh.update_vbo();
		}

		view->idx_selected_mesh = 0;
	}

	return true;
}

bool app_viewer::process_descriptor_heatmap(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static int n_eigs = 50;
	static bool status = true;
	static descriptor::signature sig = descriptor::GPS;

	ImGui::Combo("descriptor", (int *) &sig, "GPS\0HKS\0WKS\0\0");
	ImGui::InputInt("n_eigs", &n_eigs);
	if(!status) ImGui::TextColored({1, 0, 0, 1}, "Error computing features.");

	if(ImGui::Button("Run"))
	{
		descriptor features(sig, mesh, n_eigs);

		if(features)
		{
			status = true;
			n_eigs = features.n_eigs();

			float max_s = 0;
			#pragma omp parallel for reduction(max: max_s)
			for(index_t v = 0; v < mesh->n_vertices; ++v)
			{
				mesh->heatmap(v) = features(v);
				max_s = std::max(max_s, mesh->heatmap(v));
			}

			#pragma omp parallel for
			for(index_t v = 0; v < mesh->n_vertices; ++v)
				mesh->heatmap(v) /= max_s;

			mesh.update_vbo_heatmap();
		}
		else status = false;
	}

	return true;
}

bool app_viewer::process_key_points(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static float percent = 0.1;
	if(ImGui::InputFloat("percent", &percent, 0.01, 0.1, "%.2f"))
	{
		key_points kps(mesh, percent);
		mesh.selected = kps;
	}

	return true;
}

bool app_viewer::process_key_components(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static float radio = 0.25;
	if(ImGui::InputFloat("radio", &radio, 0.01, 0.1, "%.2f"))
	{
		// use the mesh selected points as key points.
		key_components kcs(mesh, mesh.selected, radio);

		#pragma omp parallel for
		for(index_t v = 0; v < mesh->n_vertices; ++v)
			mesh->heatmap(v) = (float) kcs(v) / kcs;

		mesh.update_vbo_heatmap();
	}

	return true;
}


// Hole Filling

bool app_viewer::process_poisson(viewer * p_view, const index_t k)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	size_t old_n_vertices = mesh->n_vertices;
	delete [] fill_all_holes(mesh);

	TIC(view->time) poisson(mesh, old_n_vertices, k); TOC(view->time)
	gproshan_log_var(view->time);

	mesh.update();

	return false;
}

bool app_viewer::process_poisson_laplacian_1(viewer * p_view)
{
	return process_poisson(p_view, 1);
}

bool app_viewer::process_poisson_laplacian_2(viewer * p_view)
{
	return process_poisson(p_view, 2);
}

bool app_viewer::process_poisson_laplacian_3(viewer * p_view)
{
	return process_poisson(p_view, 3);
}

bool app_viewer::process_fill_holes(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	// TODO
	//	fill_all_holes(mesh);
	/*********************************************************************/
	che * fill_mesh = new che(*mesh);
	for(index_t & v: mesh->bounds())
	{
		const std::vector<index_t> & vbounds = mesh->boundary(v);
		std::vector<vertex> vertices;
		std::vector<index_t> trigs;

		vertex center;
		for(const index_t v: vbounds)
		{
			vertices.push_back(mesh->point(v));
			center += mesh->point(v);
		}
		center /= size(vbounds);

		std::priority_queue<std::pair<float, index_t> > front;
		std::vector<uvec2> neigs(size(vertices));
/*
		auto bprev = [&](const index_t v) -> index_t &
		{
			return neigs[v].x();
		};
		auto bnext = [&](const index_t v) -> index_t &
		{
			return neigs[v].y();
		};
*/
		auto push = [&](const uvec3 & p)
		{
			neigs[p.x()] = {p.y(), p.z()};
			const float angle = 21;
			if(angle <= M_PI)
				front.push({angle, p.x()});
		};

		index_t nv = size(vertices);
		push({0, nv - 1, 1});
		for(index_t i = 1; i < size(vertices) - 1; ++i)
			push({i, i - 1, i + 1});
		push({nv - 1, nv - 2, 0});

		std::vector<bool> border;
		border.assign(true, size(vertices));

//		float angle;
		index_t v0 = NIL;
		index_t v1 = NIL;
		index_t v2 = NIL;
		while(!front.empty())
		{
//			angle = front.top().first;

			if(!(border[v0] && border[v1] && border[v2]))
				continue;


		}

		// vertices.push_back(center);

		che * old = fill_mesh;
		che * hole = new che(vertices.data(), size(vertices), trigs.data(), size(trigs) / 3);
		fill_mesh = old->merge(hole, vbounds);
		delete old;
		delete hole;
	}

	view->add_mesh(fill_mesh);
	/*********************************************************************/

	return false;
}

bool app_viewer::process_fill_holes_biharmonic_splines(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	size_t old_n_vertices, n_vertices = mesh->n_vertices;
	size_t n_holes = 0; // FIX_BOUND mesh->n_borders;

	const auto & [border_vertices, holes] = fill_all_holes_meshes(mesh);
	if(!holes) return true;

	index_t k = 2;

	for(index_t h = 0; h < n_holes; ++h)
		if(holes[h])
		{
			old_n_vertices = n_vertices;
			biharmonic_interp_2(mesh, old_n_vertices, n_vertices += holes[h]->n_vertices - size(border_vertices[h]), border_vertices[h], k);
			delete holes[h];
		}

	delete [] holes;
	delete [] border_vertices;

	mesh.update();

	return false;
}


// Others

bool app_viewer::process_select_multiple(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	static char line[128] = "";

	ImGui::InputText("select", line, sizeof(line));

	if(ImGui::Button("Add"))
	{
		std::stringstream ss(line);
		index_t v;
		while(ss >> v)
			mesh.selected.push_back(v);
	}

	return true;
}

bool app_viewer::process_threshold(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	for(index_t v = 0; v < mesh->n_vertices; ++v)
		mesh->heatmap(v) = mesh->heatmap(v) > 0.5 ? 1 : 0.5;

	return false;
}

bool app_viewer::process_noise(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	std::default_random_engine generator;
	std::uniform_int_distribution<int> d_mod_5(0, 4);
	std::uniform_int_distribution<int> d_mod_1000(0, 999);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		float r = float(d_mod_1000(generator)) / 200000;
		mesh->point(v) += (!d_mod_5(generator)) * r * mesh->normal(v);
	}

	mesh->update_normals();

	return false;
}

bool app_viewer::process_black_noise(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->selected_mesh();

	std::default_random_engine generator;
	std::uniform_int_distribution<int> d_mod_5(0, 4);
	std::uniform_int_distribution<int> d_mod_1000(0, 999);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		float r = float(d_mod_1000(generator)) / 200000;
		mesh->point(v) += (!d_mod_5(generator)) * r * mesh->normal(v);
	}

	mesh->update_normals();

	return false;
}


} // namespace gproshan

