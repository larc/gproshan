#include <gproshan/mesh/che.h>

#include <gproshan/mesh/che_obj.h>
#include <gproshan/mesh/che_off.h>
#include <gproshan/mesh/che_ply.h>
#include <gproshan/mesh/che_ptx.h>
#include <gproshan/mesh/che_xyz.h>
#include <gproshan/mesh/che_pts.h>
#include <gproshan/mesh/che_pcd.h>
#include <gproshan/mesh/che_img.h>

#include <gproshan/scenes/scene.h>

#include <cassert>
#include <cstring>
#include <cmath>
#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan {


size_t & che::rw(const size_t & n)
{
	return const_cast<size_t&>(n);
}

const che::options che::default_opts;


che::che(const che & mesh, const index_t * sorted, const che::options & opts)
{
	filename = mesh.filename;

	if(!alloc(mesh.n_vertices, mesh.n_trigs, opts))
		return;

	rw(n_edges)	= mesh.n_edges;

	if(!sorted)
	{
		memcpy(GT, mesh.GT, n_vertices * sizeof(vertex));
		memcpy(EVT, mesh.EVT, n_vertices * sizeof(index_t));
		memcpy(VT, mesh.VT, n_half_edges * sizeof(index_t));
	}
	else
	{
		index_t * inv = new index_t[n_vertices];

		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
		{
			GT[v] = mesh.GT[sorted[v]];
			inv[sorted[v]] = v;
		}

		#pragma omp parallel for
		for(index_t he = 0; he < n_half_edges; ++he)
		{
			const index_t v = mesh.VT[he];
			VT[he] = inv[v];
			if(mesh.EVT[v] == he)
				EVT[inv[v]] = he;
		}

		delete [] inv;
	}

	memcpy(OT, mesh.OT, n_half_edges * sizeof(index_t));

	if(opts.edges)
	{
		memcpy(ET, mesh.ET, n_edges * sizeof(index_t));
		memcpy(EHT, mesh.EHT, n_half_edges * sizeof(index_t));
	}

	if(opts.normals)
	{
		if(!sorted)
		{
			memcpy(VN, mesh.VN, n_vertices * sizeof(vertex));
		}
		else
		{
			#pragma omp parallel for
			for(index_t v = 0; v < n_vertices; ++v)
				VN[v] = mesh.VN[sorted[v]];
		}
	}

	if(opts.colors)
	{
		if(!sorted)
		{
			memcpy(VC, mesh.VC, n_vertices * sizeof(rgb_t));
			memcpy(VHC, mesh.VHC, n_vertices * sizeof(float));
		}
		else
		{
			#pragma omp parallel for
			for(index_t v = 0; v < n_vertices; ++v)
			{
				VC[v] = mesh.VC[sorted[v]];
				VHC[v] = mesh.VHC[sorted[v]];
			}
		}
	}
}

che::che(const size_t nv, const size_t nf)
{
	alloc(nv, nf);
}

che::che(const vertex * vertices, const index_t nv, const index_t * trigs, const index_t nf)
{
	init(vertices, nv, trigs, nf);
}

che::~che()
{
	free();
}


void che::reload()
{
	free();
	init(filename);
}

mat4 che::normalize_sphere(const float r) const
{
	vertex center;
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
	{
		#pragma omp critical
		center += GT[v];
	}
	center /= n_vertices;

	float mean_dist = 0;
	#pragma omp parallel for reduction(+: mean_dist)
	for(index_t v = 0; v < n_vertices; ++v)
		mean_dist += length(GT[v] - center);
	mean_dist /= n_vertices;

	float sigma_dist = 0;
	#pragma omp parallel for reduction(+: sigma_dist)
	for(index_t v = 0; v < n_vertices; ++v)
	{
		const float diff = mean_dist - length(GT[v] - center);
		sigma_dist += diff * diff;
	}
	sigma_dist = sqrt(sigma_dist / n_vertices);

	mat4 model_mat;

	const float scale = r / (mean_dist + sigma_dist);
	model_mat(0, 0) = model_mat(1, 1) = model_mat(2, 2) = scale;

	center *= -scale;
	model_mat(0, 3) = center.x();
	model_mat(1, 3) = center.y();
	model_mat(2, 3) = center.z();
	model_mat(3, 3) = 1;

	return model_mat;
}

mat4 che::normalize_box(const float side) const
{
	vertex pmin = INFINITY;
	vertex pmax = 0;

	for(index_t v = 0; v < n_vertices; ++v)
	{
		const vertex & p = point(v);

		pmin.x() = std::min(pmin.x(), p.x());
		pmin.y() = std::min(pmin.y(), p.y());
		pmin.z() = std::min(pmin.z(), p.z());

		pmax.x() = std::max(pmax.x(), p.x());
		pmax.y() = std::max(pmax.y(), p.y());
		pmax.z() = std::max(pmax.z(), p.z());
	}

	mat4 model_mat;

	const float scale = side / std::max({pmax.x() - pmin.x(), pmax.y() - pmin.y(), pmax.z() - pmin.z()});
	model_mat(0, 0) = model_mat(1, 1) = model_mat(2, 2) = scale;

	const vertex & translate = - scale * (pmax + pmin) / 2;
	model_mat(0, 3) = translate.x();
	model_mat(1, 3) = translate.y();
	model_mat(2, 3) = translate.z();
	model_mat(3, 3) = 1;

	return model_mat;
}

///< vcommon correspond to the first vertices of the mesh with indices to the main mesh (this)
che * che::merge(const che * mesh, const std::vector<index_t> & vcommon)
{
	const size_t n_vcommon = size(vcommon);
	const size_t n_vnew = mesh->n_vertices - n_vcommon;

	che * new_mesh = new che(n_vertices + n_vnew, n_trigs + mesh->n_trigs);

	memcpy(new_mesh->GT, GT, sizeof(vertex) * n_vertices);
	memcpy(new_mesh->GT + n_vertices, mesh->GT + n_vcommon, sizeof(vertex) * n_vnew);
	memcpy(new_mesh->VN, VN, sizeof(vertex) * n_vertices);
	memcpy(new_mesh->VN + n_vertices, mesh->VN + n_vcommon, sizeof(vertex) * n_vnew);
	memcpy(new_mesh->VC, VC, sizeof(rgb_t) * n_vertices);
	memcpy(new_mesh->VC + n_vertices, mesh->VC + n_vcommon, sizeof(rgb_t) * n_vnew);
	memcpy(new_mesh->VHC, VHC, sizeof(float) * n_vertices);
	memcpy(new_mesh->VHC + n_vertices, mesh->VHC + n_vcommon, sizeof(float) * n_vnew);

	memcpy(new_mesh->VT, VT, sizeof(index_t) * n_half_edges);

	index_t * tVT = new_mesh->VT + n_half_edges;
	for(index_t he = 0; he < mesh->n_half_edges; ++he)
		tVT[he] = mesh->VT[he] < size(vcommon) ? vcommon[mesh->VT[he]] : mesh->VT[he] + n_vertices - n_vcommon;

	new_mesh->update_evt_ot_et();
	new_mesh->update_eht();

	return new_mesh;
}

void che::update_vertices(const vertex * positions, const size_t n, const index_t v_i)
{
	if(!positions) return;
	memcpy(GT + v_i, positions, sizeof(vertex) * (!n ? n_vertices : n));
}

void che::update_heatmap(const float * hm)
{
	if(!hm)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			VHC[v] = 0.45;	// default heatmap value

		return;
	}

	memcpy(VHC, hm, n_vertices * sizeof(float));
	scale_hm = normalize(VHC, n_vertices);
}

void che::update_normals()
{
	if(!n_trigs) return;

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
	{
		vertex & n = VN[v];

		n = 0;
		for(const index_t he: star(v))
			n += area_trig(he_trig(he)) * normal_he(he);

		n /= norm(n);
	}
}

void che::invert_normals()
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		VN[v] = - VN[v];
}

void che::multiplicate_vertices()
{
	vertex * old_GT = GT;
	index_t * old_VT = VT;
	size_t nv = n_vertices;
	size_t nf = n_trigs;

	GT = nullptr;
	VT = nullptr;

	free();
	alloc(nv + nf, 3 * nf);

	memcpy(GT, old_GT, nv * sizeof(vertex));

	#pragma omp parallel for
	for(index_t f = 0; f < nf; ++f)
	{
		const index_t v = nv + f;
		const index_t old_he = f * 3;
		const index_t he = 3 * f * 3;

		GT[v] = (GT[old_VT[old_he]] + GT[old_VT[old_he + 1]] + GT[old_VT[old_he + 2]]) / 3;

		VT[he] = old_VT[old_he];
		VT[he + 1] = old_VT[old_he + 1];
		VT[he + 2] = v;

		VT[he + 3] = old_VT[old_he + 1];
		VT[he + 4] = old_VT[old_he + 2];
		VT[he + 5] = v;

		VT[he + 6] = old_VT[old_he + 2];
		VT[he + 7] = old_VT[old_he];
		VT[he + 8] = v;
	}

	delete [] old_GT;
	delete [] old_VT;

	update_evt_ot_et();
	update_eht();

	for(index_t e = 0; e < n_edges; ++e)
	{
		const index_t he = ET[e];
		if(!(he % 3) && OT[he] != NIL)
			flip(e);
	}
}

void che::remove_vertices(const std::vector<index_t> & vertices)
{
	if(!size(vertices)) return;

	gproshan_debug(removing vertex);
	for(index_t v: vertices)
	{
		for(const index_t he: star(v))
		{
			VT[he] = NIL;
			VT[he_prev(he)] = NIL;
			VT[he_next(he)] = NIL;

			gproshan_debug_var(he);
			gproshan_debug_var(he_next(he));
			gproshan_debug_var(he_prev(he));
		}

		gproshan_debug_var(EVT[v]);
		EVT[v] = NIL;
	}
	/* save in vectors */
	std::vector<vertex> new_vertices;
	std::vector<index_t> removed;
	std::vector<index_t> new_trigs; // each 3

	gproshan_debug(removing vertex);
	for(index_t v = 0; v < n_vertices; ++v)
	{
		if(EVT[v] != NIL)
			new_vertices.push_back(GT[v]);
		else
			removed.push_back(v);
	}
	removed.push_back(n_vertices);

	gproshan_debug_var(size(removed));
	gproshan_debug_var(removed[0]);
	index_t r = 1;
	index_t d = 1;
	for(index_t v = removed[0] + 1; v < n_vertices; ++v)
	{
		if(v < removed[r])
		{
			for(const index_t he: star(v))
				if(VT[he] != NIL) VT[he] = v - d;
		}
		else if(v == removed[r])
		{
			++d;
			++r;
		}
	}

	for(index_t he = 0; he < n_half_edges; ++he)
		if(VT[he] != NIL)
			new_trigs.push_back(VT[he]);
		else gproshan_error_var(he);

	gproshan_debug_var(size(new_vertices));
	gproshan_debug_var(size(new_trigs));
	gproshan_debug(removing vertex);
	free();
	gproshan_debug(removing vertex);
	init(new_vertices.data(), size(new_vertices), new_trigs.data(), size(new_trigs) / 3);
	gproshan_debug(removing vertex);
}

void che::remove_non_manifold_vertices()
{
	for( index_t he = 0; he < n_half_edges; he+=3)
	{
		if(EVT[VT[he]] == NIL || EVT[VT[he+1]] == NIL || EVT[VT[he+2]] == NIL)
		{
			VT[he] = NIL;
			VT[he + 1] = NIL;
			VT[he + 2] = NIL;
		}
	}

	/* save in vectors */
	std::vector<vertex> new_vertices;
	std::vector<index_t> removed;
	std::vector<index_t> new_trigs; // each 3

	gproshan_debug(removing vertex);
	for(index_t v = 0; v < n_vertices; ++v)
	{
		if(EVT[v] != NIL)
			new_vertices.push_back(GT[v]);
		else
			removed.push_back(v);
	}
	removed.push_back(n_vertices);

	gproshan_debug_var(size(removed));
	gproshan_debug_var(removed[0]);
	index_t r = 1;
	index_t d = 1;
	for(index_t v = removed[0] + 1; v < n_vertices; ++v)
	{
		if(v < removed[r])
		{
			for(const index_t he: star(v))
				if(VT[he] != NIL) VT[he] = v - d;
		}
		else if(v == removed[r])
		{
			++d;
			++r;
		}
	}

	for(index_t he = 0; he < n_half_edges; ++he)
		if(VT[he] != NIL)
			new_trigs.push_back(VT[he]);
		else gproshan_error_var(he);

	gproshan_debug_var(size(new_vertices));
	gproshan_debug_var(size(new_trigs));
	gproshan_debug(removing vertex);
	free();
	gproshan_debug(removing vertex);
	init(new_vertices.data(), size(new_vertices), new_trigs.data(), size(new_trigs) / 3);
	gproshan_debug(removing vertex);
}

void che::set_head_vertices(index_t * head, const size_t n)
{
	for(index_t v, i = 0; i < n; ++i)
	{
		v = head[i];

		for(index_t j = i + 1; j < n; ++j)
			if(i == head[j])
			{
				head[j] = v;
				break;
			}

		std::swap(GT[v], GT[i]);

		for(const index_t he: star(v))
			VT[he] = i;

		for(const index_t he: star(i))
			VT[he] = i;

		std::swap(EVT[v], EVT[i]);
	}
}


vertex che::normal_trig(const index_t f) const
{
	return normal_he(f * 3);
}

vertex che::normal_he(const index_t he) const
{
	const vertex & a = GT[VT[he]];
	const vertex & b = GT[VT[he_next(he)]];
	const vertex & c = GT[VT[he_prev(he)]];

	return normalize(cross(b - a, c - a));
}

vertex che::gradient_he(const index_t he, const float * f) const
{
	index_t i = VT[he];
	index_t j = VT[he_next(he)];
	index_t k = VT[he_prev(he)];

	const vertex & xi = GT[i];
	const vertex & xj = GT[j];
	const vertex & xk = GT[k];

	const vertex & n = normal_he(he);

	vertex pij = cross(n, xj - xi);
	vertex pjk = cross(n, xk - xj);
	vertex pki = cross(n, xi - xk);

	return normalize(f[i] * pjk + f[j] * pki + f[k] * pij);
}

dvec3 che::gradient_he(const index_t he, const double * f) const
{
	index_t i = VT[he];
	index_t j = VT[he_next(he)];
	index_t k = VT[he_prev(he)];

	const vertex & gi = GT[i];
	const vertex & gj = GT[j];
	const vertex & gk = GT[k];

	const dvec3 xi = {gi.x(), gi.y(), gi.z()};
	const dvec3 xj = {gj.x(), gj.y(), gj.z()};
	const dvec3 xk = {gk.x(), gk.y(), gk.z()};

	const vertex & n = normal_he(he);
	const dvec3 dn = {n.x(), n.y(), n.z()};

	dvec3 pij = cross(dn, xj - xi);
	dvec3 pjk = cross(dn, xk - xj);
	dvec3 pki = cross(dn, xi - xk);

	return normalize(f[i] * pjk + f[j] * pki + f[k] * pij);
}

vertex che::gradient(const index_t v, const float * f)
{
	vertex g;
	float area, area_star = 0;

	for(const index_t he: star(v))
	{
		area = area_trig(he_trig(he));
		area_star += area;
		g += area * gradient_he(he, f);
	}

	return g / area_star;
}


float che::heatmap_scale() const
{
	return scale_hm;
}

void che::heatmap_scale(const float shm)
{
	scale_hm = shm;
}

float che::heatmap_scale(const index_t v) const
{
	assert(v < n_vertices);
	return scale_hm * VHC[v];
}


index_t che::twin_he(const index_t he) const
{
	assert(he < n_half_edges);
	return OT[he];
}

index_t che::edge_u(const index_t e) const
{
	assert(e < n_edges);
	return VT[ET[e]];
}

index_t che::edge_v(const index_t e) const
{
	assert(e < n_edges);
	return VT[he_next(ET[e])];
}

index_t che::edge_he_0(const index_t e) const
{
	assert(e < n_edges);
	return ET[e];
}

index_t che::edge_he_1(const index_t e) const
{
	assert(e < n_edges);
	return OT[ET[e]];
}

const vertex & che::vertex_he(const index_t he) const
{
	assert(he < n_half_edges);
	return GT[VT[he]];
}

const vertex & che::vertex_edge_u(const index_t e) const
{
	assert(e < n_edges);
	return GT[VT[ET[e]]];
}

const vertex & che::vertex_edge_v(const index_t e) const
{
	assert(e < n_edges);
	return GT[VT[he_next(ET[e])]];
}

index_t che::evt(const index_t v) const
{
	assert(v < n_vertices);
	return EVT[v];
}


std::vector<index_t> che::link(const index_t v) const
{
	assert(v < n_vertices);

	std::vector<index_t> vlink;

	if(is_vertex_bound(v))
		vlink.push_back(VT[he_next(EVT[v])]);

	for(const index_t he: star(v))
		vlink.push_back(VT[he_prev(he)]);

	return vlink;
}

void che::edge_collapse(const std::vector<index_t> & sort_edges)
{
	gproshan_error_var(size(sort_edges));
	// TODO
}

void che::compute_toplesets(index_t * toplesets, index_t * sorted, std::vector<index_t> & limits, const std::vector<index_t> & sources, const index_t k)
{
	if(!size(sources)) return;

	memset(toplesets, -1, sizeof(index_t) * n_vertices);

	index_t level = 0;

	index_t p = 0;
	for(const index_t s: sources)
	{
		sorted[p++] = s;
		toplesets[s] = level;
	}

	limits.push_back(0);
	for(index_t i = 0; i < p; ++i)
	{
		const index_t v = sorted[i];

		if(toplesets[v] > level)
		{
			if(++level > k) break;
			limits.push_back(i);
		}

		for(const index_t u: link(v))
		{
			if(toplesets[u] == NIL)
			{
				toplesets[u] = toplesets[v] + 1;
				sorted[p++] = u;
			}
		}
	}

	assert(p <= n_vertices);
	limits.push_back(p);
}


///< return a vector of indices of one vertex per boundary
std::vector<index_t> che::bounds() const
{
	if(!n_trigs) return {};
	if(!manifold) return {};

	std::vector<index_t> vbounds;

	bool * is_bound = new bool[n_vertices];
	memset(is_bound, 0, sizeof(bool) * n_vertices);

	for(index_t v = 0; v < n_vertices; ++v)
		if(!is_bound[v] && is_vertex_bound(v))
		{
			vbounds.push_back(v);

			for(const index_t b: boundary(v))
				is_bound[b] = true;
		}

	delete [] is_bound;

	return vbounds;
}

///< return a vector of the indices of the boundary where v belongs
std::vector<index_t> che::boundary(const index_t v) const
{
	std::vector<index_t> vbound;

	index_t he_end = EVT[v];
	index_t he = he_end;
	do
	{
		vbound.push_back(VT[he]);
		he = EVT[VT[he_next(he)]];
	}
	while(he != NIL && he != he_end);

	return vbound;
}

bool che::is_vertex_bound(const index_t v) const
{
	assert(v < n_vertices);
	return EVT[v] != NIL && OT[EVT[v]] == NIL;
}

bool che::is_edge_bound(const index_t e) const
{
	assert(e < n_edges);
	return OT[ET[e]] == NIL;
}


const std::string che::name() const
{
	index_t p = filename.find_last_of('/');
	index_t q = filename.find_last_of('.');
	return filename.substr(p + 1, q - p - 1);
}

const std::string che::name_size() const
{
	return name() + "_" + std::to_string(n_vertices);
}

const std::string che::filename_size() const
{
	return filename + "_" + std::to_string(n_vertices);
}


size_t che::genus() const
{
	size_t g = n_vertices - n_edges + n_trigs;
	return (g - 2) / (-2);
}

size_t che::memory() const
{
	return	sizeof(*this) +
			n_vertices * (2 * sizeof(vertex) + sizeof(index_t) + sizeof(float) + sizeof(rgb_t)) +
			sizeof(index_t) * (3 * n_half_edges + n_edges) +
			size(filename);
}

size_t che::max_degree() const
{
	size_t d, md = 0;

	#pragma omp parallel for private(d) reduction(max: md)
	for(index_t v = 0; v < n_vertices; ++v)
	{
		d = 0;
		for([[maybe_unused]] const index_t he: star(v)) ++d;
		d += is_vertex_bound(v);
		md = std::max(md, d);
	}

	return md;
}

float che::quality() const
{
	float q = 0;

	#pragma omp parallel for reduction(+: q)
	for(index_t t = 0; t < n_trigs; ++t)
		q += pdetriq(t) > 0.6; // is confederating good triangle

	return q * 100 / n_trigs;
}

float che::mean_edge() const
{
	float m = 0;

	#pragma omp parallel for reduction(+: m)
	for(index_t e = 0; e < n_edges; ++e)
		m += norm(GT[VT[ET[e]]] - GT[VT[he_next(ET[e])]]);

	return m / n_edges;
}

float che::area_surface() const
{
	float area = 0;

	#pragma omp parallel for reduction(+: area)
	for(index_t i = 0; i < n_trigs; ++i)
		area += area_trig(i);

	return area;
}

bool che::is_manifold() const
{
	return manifold;
}

bool che::is_scene() const
{
	return false;
}

bool che::is_pointcloud() const
{
	return n_trigs == 0;
}


void che::flip(const index_t e)
{
	index_t ha = ET[e];
	index_t hb = OT[ha];

	if(hb == NIL)
		return;

	index_t va = VT[ha];
	index_t vb = VT[hb];
	index_t vc = VT[he_prev(ha)];
	index_t vd = VT[he_prev(hb)];

	index_t et_pa = EHT[he_prev(ha)];
	index_t et_na = EHT[he_next(ha)];
	index_t et_pb = EHT[he_prev(hb)];
	index_t et_nb = EHT[he_next(hb)];

	index_t ot_pa = OT[he_prev(ha)];
	index_t ot_na = OT[he_next(ha)];
	index_t ot_pb = OT[he_prev(hb)];
	index_t ot_nb = OT[he_next(hb)];

	VT[he_prev(ha)] = vb;
	VT[ha] = vc;
	VT[he_next(ha)] = vd;
	VT[he_prev(hb)] = va;
	VT[hb] = vd;
	VT[he_next(hb)] = vc;

	if(ot_pa != NIL) OT[ot_pa] = he_next(hb);
	if(ot_na != NIL) OT[ot_na] = he_prev(ha);
	if(ot_pb != NIL) OT[ot_pb] = he_next(ha);
	if(ot_nb != NIL) OT[ot_nb] = he_prev(hb);

	OT[he_prev(ha)] = ot_na;
	OT[he_next(ha)] = ot_pb;
	OT[he_prev(hb)] = ot_nb;
	OT[he_next(hb)] = ot_pa;

	ET[et_pa] = he_prev(ha);
	ET[et_na] = he_next(ha);
	ET[et_pb] = he_prev(hb);
	ET[et_nb] = he_next(hb);

	EHT[he_prev(ha)] = EHT[OT[he_prev(ha)]] = et_pa;
	EHT[he_next(ha)] = EHT[OT[he_next(ha)]] = et_na;
	EHT[he_prev(hb)] = EHT[OT[he_prev(hb)]] = et_pb;
	EHT[he_next(hb)] = EHT[OT[he_next(hb)]] = et_nb;

	if(EVT[va] == he_next(hb) || EVT[va] == ha) EVT[va] = he_prev(hb);
	if(EVT[vb] == he_next(ha) || EVT[vb] == hb) EVT[vb] = he_prev(ha);
	if(EVT[vc] == he_prev(ha)) EVT[vc] = he_next(hb);
	if(EVT[vd] == he_prev(hb)) EVT[vd] = he_next(ha);
}

float che::cotan(const index_t he) const
{
	if(he == NIL) return 0;

	vertex a = GT[VT[he]] - GT[VT[he_prev(he)]];
	vertex b = GT[VT[he_next(he)]] - GT[VT[he_prev(he)]];

	return dot(a, b) / norm(cross(a, b));
}

// https://www.mathworks.com/help/pde/ug/pdetriq.html
// 4*sqrt(3)*a
// q = ----------------
// h1^2+h2^2+h3^2
float che::pdetriq(const index_t t) const
{
	index_t he = t * 3;
	float h[3] = {
					norm(GT[VT[he_next(he)]] - GT[VT[he]]),
					norm(GT[VT[he_prev(he)]] - GT[VT[he_next(he)]]),
					norm(GT[VT[he]] - GT[VT[he_prev(he)]])
					};
	return (4 * sqrt(3) * area_trig(t)) / (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
}

float che::area_trig(const index_t t) const
{
	index_t he = t * 3;
	vertex a = GT[VT[he_next(he)]] - GT[VT[he]];
	vertex b = GT[VT[he_prev(he)]] - GT[VT[he]];

	return norm(cross(a, b)) / 2;
}

float che::area_vertex(const index_t v) const
{
	float area_star = 0;
	for(const index_t he: star(v))
		area_star += area_trig(he_trig(he));

	return area_star / 3;
}

// The Gauss-Bonnet Scheme
float che::mean_curvature(const index_t v) const
{
	float h = 0;
	float a = 0;

	for(const index_t he: star(v))
	{
		a += area_trig(he_trig(he));
		h += norm(GT[VT[he_next(he)]] - GT[v]) * dot(normal(v), normal_he(he));
	}

	return 0.75 * h / a;
}


void che::init(const vertex * vertices, const index_t n_v, const index_t * trigs, const index_t n_f)
{
	alloc(n_v, n_f);

	memcpy(GT, vertices, n_vertices * sizeof(vertex));
	memcpy(VT, trigs, n_half_edges * sizeof(index_t));

	update_evt_ot_et();
	update_eht();
}

void che::init(const std::string & file)
{
	filename = file;
	read_file(filename);

	update_evt_ot_et();
	update_eht();
}

bool che::alloc(const size_t n_v, const size_t n_f, const che::options & opt)
{
	if(!n_v) return false;

	rw(n_vertices)		= n_v;
	rw(n_trigs)			= n_f;
	rw(n_half_edges)	= 3 * n_trigs;
	rw(n_edges)			= n_half_edges;				// max number of edges

	GT	= new vertex[n_vertices];
	EVT	= new index_t[n_vertices];

	if(n_half_edges)
	{
		VT	= new index_t[n_half_edges];
		OT	= new index_t[n_half_edges];

		if(opt.edges)
		{
			ET	= new index_t[n_half_edges];
			EHT	= new index_t[n_half_edges];
		}
	}

	if(opt.normals)
		VN	= new vertex[n_vertices];

	if(opt.colors)
	{
		VC	= new rgb_t[n_vertices];
		VHC	= new float[n_vertices];
		update_heatmap();
	}

	return true;
}

void che::free()
{
	delete [] GT;	GT = nullptr;
	delete [] EVT;	EVT = nullptr;
	delete [] VT;	VT = nullptr;
	delete [] OT;	OT = nullptr;
	delete [] ET;	ET = nullptr;
	delete [] EHT;	EHT = nullptr;
	delete [] VN;	VN = nullptr;
	delete [] VC;	VC = nullptr;
	delete [] VHC;	VHC = nullptr;
}


void che::read_file(const std::string & ) {}		/* virtual */


void che::update_evt_ot_et()
{
	if(!n_trigs) return;

	memset(EVT, -1, sizeof(index_t) * n_vertices);
	memset(OT, -1, sizeof(index_t) * n_half_edges);

	std::vector<index_t> vnhe;
	vnhe.assign(n_vertices, 0);

	for(index_t he = 0; he < n_half_edges; ++he)
		++vnhe[VT[he]];

	std::vector<index_t *> vhe(n_vertices);
	vhe[0] = new index_t[n_half_edges];
	for(index_t v = 1; v < n_vertices; ++v)
		vhe[v] = vhe[v - 1] + vnhe[v - 1];

	vnhe.assign(n_vertices, 0);
	for(index_t he = 0; he < n_half_edges; ++he)
		vhe[VT[he]][vnhe[VT[he]]++] = he;

	size_t ne = 0;
	for(index_t he = 0; he < n_half_edges; ++he)
	{
		const index_t u = VT[he];
		const index_t v = VT[he_next(he)];

		EVT[u] = he;

		if(OT[he] == NIL)
		{
			if(ET) ET[ne++] = he;

			index_t ohe = NIL;
			for(index_t j = 0; j < vnhe[v]; ++j)
				if(VT[he_next(vhe[v][j])] == u)
				{
					ohe = vhe[v][j];
					break;
				}

			if(ohe != NIL && OT[ohe] == NIL)
			{
				OT[he] = ohe;
				OT[ohe] = he;
			}
		}
	}

	delete [] vhe[0];

	rw(n_edges) = ne;


	for(index_t he = 0; he < n_half_edges; ++he)
		if(OT[he] == NIL && EVT[VT[he]] != NIL)
		{
			if(OT[EVT[VT[he]]] == NIL && EVT[VT[he]] != he)
			{
				manifold = false;
				EVT[VT[he]] = NIL;
			}
			else EVT[VT[he]] = he;
		}
}

void che::update_eht()
{
	#pragma omp parallel for
	for(index_t e = 0; e < n_edges; ++e)
	{
		EHT[ET[e]] = e;
		if(OT[ET[e]] != NIL)
			EHT[OT[ET[e]]] = e;
	}
}


std::vector<index_t> che::trig_convex_polygon(const index_t * P, const size_t n)
{
	std::vector<index_t> trigs;
	trigs.reserve(3 * (n - 2));

	index_t a = n - 1;
	index_t b = 0;
	index_t c = 1;
	while(a > c)
	{
		trigs.push_back(P[a]);
		trigs.push_back(P[b]);
		trigs.push_back(P[c]);

		b = (b + 1 == c) ? a-- : c++;
	}

	return trigs;
}

che * che::load_mesh(const std::string & file_path)
{
	size_t pos = file_path.rfind('.');
	assert(pos != std::string::npos);

	std::string extension = file_path.substr(pos + 1);

	if(extension == "obj")
	{
		scene * sc = new scene(file_path);
		if(!sc->is_scene())
		{
			delete sc;
			return new che_obj(file_path);
		}
		return sc;
	}
	if(extension == "off") return new che_off(file_path);
	if(extension == "ply") return new che_ply(file_path);
	if(extension == "ptx") return new che_ptx(file_path);
	if(extension == "xyz") return new che_xyz(file_path);
	if(extension == "txt") return new che_xyz(file_path);
	if(extension == "pts") return new che_pts(file_path);
	if(extension == "pcd") return new che_pcd(file_path);

	return new che_img(file_path);
}


} // namespace gproshan

