#include "mesh/che.h"

#include "mesh/kdtree.h"
#include "include_arma.h"

#include <cassert>
#include <cstring>
#include <cmath>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


size_t & rw(const size_t & n)
{
	return const_cast<size_t&>(n);
}

index_t trig(const index_t & he)
{
	if(he == NIL) return NIL;
	return he / che::mtrig;
}

index_t next(const index_t & he)
{
	if(he == NIL) return NIL;
	return che::mtrig * trig(he) + (he + 1) % che::mtrig;
}

index_t prev(const index_t & he)
{
	if(he == NIL) return NIL;
	return che::mtrig * trig(he) + (he + che::mtrig - 1) % che::mtrig;
}

CHE::CHE(const che * mesh)
{
	n_vertices = mesh->n_vertices;
	n_faces = mesh->n_faces;
	n_half_edges = mesh->n_half_edges;

	GT = (vertex_cu *) mesh->GT;
	VN = (vertex_cu *) mesh->VN;
	VC = mesh->VC;
	VT = mesh->VT;
	OT = mesh->OT;
	EVT = mesh->EVT;
}

che::che(const che & mesh)
{
	filename = mesh.filename;

	alloc(mesh.n_vertices, mesh.n_faces);
	rw(n_edges)	= mesh.n_edges;

	memcpy(GT, mesh.GT, n_vertices * sizeof(vertex));
	memcpy(VT, mesh.VT, n_half_edges * sizeof(index_t));
	memcpy(OT, mesh.OT, n_half_edges * sizeof(index_t));
	memcpy(EVT, mesh.EVT, n_vertices * sizeof(index_t));
	memcpy(ET, mesh.ET, n_edges * sizeof(index_t));
	memcpy(EHT, mesh.EHT, n_half_edges * sizeof(index_t));
	memcpy(VN, mesh.VN, n_vertices * sizeof(vertex));
	memcpy(VC, mesh.VC, n_vertices * sizeof(rgb_t));
	memcpy(VHC, mesh.VHC, n_vertices * sizeof(real_t));
}

che::che(const size_t & n_v, const size_t & n_f)
{
	alloc(n_v, n_f);
}

che::che(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f)
{
	init(vertices, n_v, faces, n_f);
}

che::~che()
{
	free();
}

vector<index_t> che::star(const index_t & v) const
{
	assert(v >= n_vertices);

	vector<index_t> vstar;
	for_star(he, this, v)
		vstar.push_back(he);

	return vstar;
}

vector<index_t> che::link(const index_t & v) const
{
	assert(v >= n_vertices);

	vector<index_t> vlink;

	if(is_vertex_bound(v))
		vlink.push_back(VT[next(EVT[v])]);

	for_star(he, this, v)
		vlink.push_back(VT[prev(he)]);

	return vlink;
}

///< return a vector of indices of one vertex per boundary
vector<index_t> che::bounds() const
{
	if(!n_faces) return {};
	if(!manifold) return {};

	vector<index_t> vbounds;

	bool * is_bound = new bool[n_vertices];
	memset(is_bound, 0, sizeof(bool) * n_vertices);

	for(index_t v = 0; v < n_vertices; ++v)
		if(!is_bound[v] && is_vertex_bound(v))
		{
			vbounds.push_back(v);

			for_boundary(he, this, v)
				is_bound[VT[he]] = true;
		}

	delete [] is_bound;

	return vbounds;
}

///< return a vector of the indices of the boundary where v belongs
vector<index_t> che::boundary(const index_t & v) const
{
	vector<index_t> vbound;

	for_boundary(he, this, v)
		vbound.push_back(VT[he]);

	return vbound;
}

bool che::is_vertex_bound(const index_t & v) const
{
	assert(v < n_vertices && EVT[v] < n_half_edges);
	return EVT[v] != NIL && OT[EVT[v]] == NIL;
}

bool che::is_edge_bound(const index_t & e) const
{
	return OT[ET[e]] == NIL;
}

void che::flip(const index_t & e)
{
	index_t ha = ET[e];
	index_t hb = OT[ha];

	if(hb == NIL)
		return;

	index_t va = VT[ha];
	index_t vb = VT[hb];
	index_t vc = VT[prev(ha)];
	index_t vd = VT[prev(hb)];

	index_t et_pa = EHT[prev(ha)];
	index_t et_na = EHT[next(ha)];
	index_t et_pb = EHT[prev(hb)];
	index_t et_nb = EHT[next(hb)];

	index_t ot_pa = OT[prev(ha)];
	index_t ot_na = OT[next(ha)];
	index_t ot_pb = OT[prev(hb)];
	index_t ot_nb = OT[next(hb)];

	VT[prev(ha)] = vb;
	VT[ha] = vc;
	VT[next(ha)] = vd;
	VT[prev(hb)] = va;
	VT[hb] = vd;
	VT[next(hb)] = vc;

	if(ot_pa != NIL) OT[ot_pa] = next(hb);
	if(ot_na != NIL) OT[ot_na] = prev(ha);
	if(ot_pb != NIL) OT[ot_pb] = next(ha);
	if(ot_nb != NIL) OT[ot_nb] = prev(hb);

	OT[prev(ha)] = ot_na;
	OT[next(ha)] = ot_pb;
	OT[prev(hb)] = ot_nb;
	OT[next(hb)] = ot_pa;

	ET[et_pa] = prev(ha);
	ET[et_na] = next(ha);
	ET[et_pb] = prev(hb);
	ET[et_nb] = next(hb);

	EHT[prev(ha)] = EHT[OT[prev(ha)]] = et_pa;
	EHT[next(ha)] = EHT[OT[next(ha)]] = et_na;
	EHT[prev(hb)] = EHT[OT[prev(hb)]] = et_pb;
	EHT[next(hb)] = EHT[OT[next(hb)]] = et_nb;

	if(EVT[va] == next(hb) || EVT[va] == ha) EVT[va] = prev(hb);
	if(EVT[vb] == next(ha) || EVT[vb] == hb) EVT[vb] = prev(ha);
	if(EVT[vc] == prev(ha)) EVT[vc] = next(hb);
	if(EVT[vd] == prev(hb)) EVT[vd] = next(ha);
}

// https://www.mathworks.com/help/pde/ug/pdetriq.html
// 4*sqrt(3)*a
// q = ----------------
// h1^2+h2^2+h3^2
real_t che::pdetriq(const index_t & t) const
{
	index_t he = t * che::mtrig;
	real_t h[3] = {
						*(GT[VT[next(he)]] - GT[VT[he]]),
						*(GT[VT[prev(he)]] - GT[VT[next(he)]]),
						*(GT[VT[he]] - GT[VT[prev(he)]])
					};
	return (4 * sqrt(3) * area_trig(t)) / (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
}

real_t che::quality()
{
	real_t q = 0;

	#pragma omp parallel for reduction(+: q)
	for(index_t t = 0; t < n_faces; ++t)
		q += pdetriq(t) > 0.6; // is confederating good triangle

	return q * 100 / n_faces;
}

real_t che::area_trig(const index_t & t) const
{
	index_t he = t * che::mtrig;
	vertex a = GT[VT[next(he)]] - GT[VT[he]];
	vertex b = GT[VT[prev(he)]] - GT[VT[he]];

	return *(a * b) / 2;
}

real_t che::area_vertex(const index_t & v) const
{
	real_t area_star = 0;
	for_star(he, this, v)
		area_star += area_trig(trig(he));

	return area_star / 3;
}

real_t che::area_surface() const
{
	real_t area = 0;

	#pragma omp parallel for reduction(+: area)
	for(index_t i = 0; i < n_faces; ++i)
		area += area_trig(i);

	return area;
}

void che::update_heatmap(const real_t * hm, real_t max_color)
{
	if(!hm)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; ++v)
			VHC[v] = 0.45;

		return;
	}

	if(max_color < numeric_limits<real_t>::epsilon())
	{
		#pragma omp parallel for reduction(max: max_color)
		for(index_t v = 0; v < n_vertices; ++v)
			if(hm[v] < INFINITY)
				max_color = max(hm[v], max_color);
	}

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		VHC[v] = hm[v] / max_color;
}

const che::rgb_t & che::rgb(const index_t & v) const
{
	assert(VC && v < n_vertices);
	return VC[v];
}

che::rgb_t & che::rgb(const index_t & v)
{
	assert(VC && v < n_vertices);
	return VC[v];
}

vertex che::color(const index_t & v) const
{
	assert(VC && v < n_vertices);
	return VC[v];
}

vertex che::shading_color(const index_t & f, const float & u, const float & v, const float & w) const
{
	index_t he = f * che::mtrig;
	return VC ? u * color(VT[he]) + v * color(VT[he + 1]) + w * color(VT[he + 2]) : rgb_t();
}

const real_t & che::heatmap(const index_t & v) const
{
	assert(VHC && v < n_vertices);
	return VHC[v];
}

real_t & che::heatmap(const index_t & v)
{
	assert(VHC && v < n_vertices);
	return VHC[v];
}

void che::update_normals()
{
	if(!n_faces) return;

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
	{
		vertex & n = VN[v];

		n = 0;
		for_star(he, this, v)
			n += area_trig(trig(he)) * normal_he(he);

		n /= *n;
	}
}

const vertex & che::normal(const index_t & v) const
{
	assert(VN && v < n_vertices);
	return VN[v];
}

vertex & che::normal(const index_t & v)
{
	assert(VN && v < n_vertices);
	return VN[v];
}

vertex che::shading_normal(const index_t & f, const float & u, const float & v, const float & w) const
{
	index_t he = f * che::mtrig;

	return {u * VN[VT[he]] + v * VN[VT[he + 1]] + w * VN[VT[he + 2]]};
}

vertex che::normal_trig(const index_t & f) const
{
	return normal_he(f * che::mtrig);
}

vertex che::normal_he(const index_t & he) const
{
	vertex n = (GT[VT[next(he)]] - GT[VT[he]]) * (GT[VT[prev(he)]] - GT[VT[he]]);
	return n / *n;
}

vertex che::gradient_he(const index_t & he, const real_t *const & f) const
{
	index_t i = VT[he];
	index_t j = VT[next(he)];
	index_t k = VT[prev(he)];

	vertex xi = GT[i];
	vertex xj = GT[j];
	vertex xk = GT[k];

	vertex n = normal_he(he);

	real_t A2 = area_trig(trig(he)) * 2;

	vertex pij = n * (xj - xi);
	vertex pjk = n * (xk - xj);
	vertex pki = n * (xi - xk);

	vertex g = ( f[i] * pjk + f[j] * pki + f[k] * pij ) / A2;
	return g / *g;
}

vertex che::gradient(const index_t & v, const real_t *const & f)
{
	vertex g;
	real_t area, area_star = 0;

	for_star(he, this, v)
	{
		area = area_trig(trig(he));
		area_star += area;
		g += area * gradient_he(he, f);
	}

	return g / area_star;
}

vertex che::barycenter(const index_t & t) const
{
	vertex bc;
	index_t tmp, he = t * che::mtrig;
	tmp = he;

	do
	{
		bc += GT[VT[he]];
		he = next(he);
	}
	while(he != tmp);

	return bc / che::mtrig;
}

real_t che::cotan(const index_t & he) const
{
	if(he == NIL) return 0;

	vertex a = GT[VT[he]] - GT[VT[prev(he)]];
	vertex b = GT[VT[next(he)]] - GT[VT[prev(he)]];

	return (a, b) / *(a * b);
}

real_t che::mean_edge() const
{
	real_t m = 0;

	#pragma omp parallel for reduction(+: m)
	for(index_t e = 0; e < n_edges; ++e)
		m += *(GT[VT[ET[e]]] - GT[VT[next(ET[e])]]);

	return m / n_edges;
}

size_t che::memory() const
{
	return sizeof(*this) + n_vertices * (sizeof(vertex) + sizeof(index_t)) + filename.size()
						+ sizeof(index_t) * (3 * n_half_edges + n_edges);
}

size_t che::genus() const
{
	size_t g = n_vertices - n_edges + n_faces;
	return (g - 2) / (-2);
}

// The Gauss-Bonnet Scheme
real_t che::mean_curvature(const index_t & v)
{
	real_t h = 0;
	real_t a = 0;

	for_star(he, this, v)
	{
		a += area_trig(trig(he));
		h += *(GT[VT[next(he)]] - GT[v]) * (normal(v), normal_he(he));
	}

	return 0.75 * h / a;
}

void che::normalize()
{
	vertex center;

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
	{
		#pragma omp critical
		center += GT[v];
	}

	center /= n_vertices;

	real_t max_norm = 0;

	#pragma omp parallel for reduction(max : max_norm)
	for(index_t v = 0; v < n_vertices; ++v)
	{
		GT[v] -= center;
		max_norm = std::max(max_norm, *GT[v]);
	}

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		GT[v] /= max_norm;
}

bool che::is_pointcloud() const
{
	return n_faces == 0;
}

bool che::is_manifold() const
{
	return manifold;
}

const index_t & che::vt(const index_t & he) const
{
	assert(he < n_half_edges);
	return VT[he];
}

const vertex & che::gt(const index_t & v) const
{
	assert(v < n_vertices);
	return GT[v];
}

const vertex & che::gt_vt(const index_t & he) const
{
	assert(he < n_half_edges);
	return GT[VT[he]];
}

const vertex & che::gt_vt_next_evt(const index_t & v) const
{
	assert(v < n_vertices);
	return GT[VT[next(EVT[v])]];
}

const vertex & che::gt_e(const index_t & e, const bool & op)
{
	assert(e < n_edges);
	return op ? GT[VT[next(ET[e])]] : GT[VT[ET[e]]];
}

const index_t & che::vt_e(const index_t & e, const bool & op)
{
	assert(e < n_edges);
	return op ? VT[next(ET[e])] : VT[ET[e]];
}

const index_t & che::et(const index_t & e) const
{
	assert(e < n_edges);
	return ET[e];
}

const index_t & che::ot_et(const index_t & e) const
{
	assert(e < n_edges);
	return OT[ET[e]];
}

const index_t & che::ot(const index_t & he) const
{
	assert(he < n_half_edges);
	return OT[he];
}

const index_t & che::ot_evt(const index_t & v) const
{
	assert(v < n_vertices);
	return OT[EVT[v]];
}

const index_t & che::evt(const index_t & v) const
{
	assert(v < n_vertices);
	return EVT[v];
}

size_t che::max_degree() const
{
	size_t d, md = 0;

	#pragma omp parallel for private(d) reduction(max: md)
	for(index_t v = 0; v < n_vertices; ++v)
	{
		d = 0;
		for_star(he, this, v) ++d;
		d += is_vertex_bound(v);
		md = max(md, d);
	}

	return md;
}

vertex & che::get_vertex(index_t v)
{
	return GT[v];
}

void che::set_vertices(const vertex *const& positions, size_t n, const index_t & v_i)
{
	if(!positions) return;
	if(!n) n = n_vertices;
	memcpy(GT + v_i, positions, sizeof(vertex) * n);
}

const string che::filename_size() const
{
	return filename + "_" + to_string(n_vertices);
}

const string che::name() const
{
	index_t p = filename.find_last_of('/');
	index_t q = filename.find_last_of('.');
	return filename.substr(p + 1, q - p - 1);
}

const string che::name_size() const
{
	return name() + "_" + to_string(n_vertices);
}

void che::reload()
{
	free();
	init(filename);
}

void che::compute_toplesets(index_t *& toplesets, index_t *& sorted, vector<index_t> & limits, const vector<index_t> & sources, const index_t & k)
{
	if(!sources.size()) return;

	memset(toplesets, -1, sizeof(index_t) * n_vertices);

	index_t level = 0;

	index_t p = 0;
	for(const index_t & s: sources)
	{
		sorted[p++] = s;

		if(toplesets[s] == NIL)
			toplesets[s] = level;
	}

	limits.push_back(0);
	for(index_t i = 0; i < p; ++i)
	{
		const index_t & v = sorted[i];

		if(toplesets[v] > level)
		{
			if(++level > k) break;

			limits.push_back(i);
		}

		for(const index_t & u: link(v))
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

void che::multiplicate_vertices()
{
	vertex * old_GT = GT;
	index_t * old_VT = VT;
	size_t nv = n_vertices;
	size_t nf = n_faces;

	GT = nullptr;
	VT = nullptr;

	free();
	alloc(nv + nf, 3 * nf);

	memcpy(GT, old_GT, nv * sizeof(vertex));

	#pragma omp parallel for
	for(index_t f = 0; f < nf; ++f)
	{
		const index_t & v = nv + f;
		const index_t & old_he = f * che::mtrig;
		const index_t & he = 3 * f * che::mtrig;

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
		const index_t & he = ET[e];
		if(!(he % 3) && OT[he] != NIL)
			flip(e);
	}
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
	vector<vertex> new_vertices;
	vector<index_t> removed;
	vector<index_t> new_faces; // each 3

	gproshan_debug(removing vertex);
	for(index_t v = 0; v < n_vertices; ++v)
	{
		if(EVT[v] != NIL)
			new_vertices.push_back(GT[v]);
		else
			removed.push_back(v);
	}
	removed.push_back(n_vertices);

	gproshan_debug_var(removed.size());
	gproshan_debug_var(removed[0]);
	index_t r = 1;
	index_t d = 1;
	for(index_t v = removed[0] + 1; v < n_vertices; ++v)
	{
		if(v < removed[r])
		{
			for_star(he, this, v)
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
			new_faces.push_back(VT[he]);
		else gproshan_error_var(he);

	gproshan_debug_var(new_vertices.size());
	gproshan_debug_var(new_faces.size());
	gproshan_debug(removing vertex);
	free();
	gproshan_debug(removing vertex);
	init(new_vertices.data(), new_vertices.size(), new_faces.data(), new_faces.size() / che::mtrig);
	gproshan_debug(removing vertex);
}

void che::remove_vertices(const vector<index_t> & vertices)
{
	if(!vertices.size()) return;

	gproshan_debug(removing vertex);
	for(index_t v: vertices)
	{
		for_star(he, this, v)
		{
			VT[he] = NIL;
			VT[prev(he)] = NIL;
			VT[next(he)] = NIL;

			gproshan_debug_var(he);
			gproshan_debug_var(next(he));
			gproshan_debug_var(prev(he));
		}

		gproshan_debug_var(EVT[v]);
		EVT[v] = NIL;
	}
	/* save in vectors */
	vector<vertex> new_vertices;
	vector<index_t> removed;
	vector<index_t> new_faces; // each 3

	gproshan_debug(removing vertex);
	for(index_t v = 0; v < n_vertices; ++v)
	{
		if(EVT[v] != NIL)
			new_vertices.push_back(GT[v]);
		else
			removed.push_back(v);
	}
	removed.push_back(n_vertices);

	gproshan_debug_var(removed.size());
	gproshan_debug_var(removed[0]);
	index_t r = 1;
	index_t d = 1;
	for(index_t v = removed[0] + 1; v < n_vertices; ++v)
	{
		if(v < removed[r])
		{
			for_star(he, this, v)
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
			new_faces.push_back(VT[he]);
		else gproshan_error_var(he);

	gproshan_debug_var(new_vertices.size());
	gproshan_debug_var(new_faces.size());
	gproshan_debug(removing vertex);
	free();
	gproshan_debug(removing vertex);
	init(new_vertices.data(), new_vertices.size(), new_faces.data(), new_faces.size() / che::mtrig);
	gproshan_debug(removing vertex);
}

void che::merge(const che * mesh, const vector<index_t> & com_vertices)
{
//	write_file("big.off");
//	mesh->write_file("small.off");
gproshan_debug(fill_holes);
	size_t ncv = com_vertices.size();
	bool is_open = mesh->VT[next(mesh->EVT[0])] >= ncv;
	gproshan_debug_var(is_open);

	size_t nv = n_vertices + mesh->n_vertices - ncv;
	size_t nf = n_faces + mesh->n_faces;
	size_t nh = n_half_edges + mesh->n_half_edges;
	size_t ne = n_edges + mesh->n_edges - (ncv - is_open);

	vertex * aGT = new vertex[nv];
	index_t * aVT = new index_t[nh];
	index_t * aOT = new index_t[nh];
	index_t * aEVT = new index_t[nv];
	index_t * aET = new index_t[ne];
	index_t * aEHT = new index_t[nh];

gproshan_debug(fill_holes);
	memcpy(aGT, GT, sizeof(vertex) * n_vertices);
	memcpy(aGT + n_vertices, mesh->GT + ncv, sizeof(vertex) * (nv - n_vertices));

	memcpy(aVT, VT, sizeof(index_t) * n_half_edges);

	index_t * t_aVT = aVT + n_half_edges;
	for(index_t he = 0; he < mesh->n_half_edges; ++he)
		t_aVT[he] = mesh->VT[he] < ncv ? com_vertices[mesh->VT[he]] : mesh->VT[he] + n_vertices - ncv;
gproshan_debug(fill_holes);

	memcpy(aOT, OT, sizeof(index_t) * n_half_edges);
gproshan_debug(fill_holes);

	index_t * t_aOT = aOT + n_half_edges;
	for(index_t he = 0; he < mesh->n_half_edges; ++he)
		t_aOT[he] = mesh->OT[he] != NIL ? mesh->OT[he] + n_half_edges : NIL;
gproshan_debug(fill_holes);

	for(index_t v, he_v, he_i, i = 0; i < ncv; ++i)
	{
		he_i = mesh->EVT[i];
		if(he_i != NIL && mesh->VT[next(he_i)] < ncv)
		{
			v = com_vertices[mesh->VT[next(he_i)]];
			he_v = EVT[v];
			aOT[he_v] = he_i + n_half_edges;
			aOT[aOT[he_v]] = he_v;
		}
	}

gproshan_debug(fill_holes);
	memcpy(aEVT, EVT, sizeof(index_t) * n_vertices);
gproshan_debug(fill_holes);
	if(is_open)
		aEVT[com_vertices[0]] = mesh->EVT[0] != NIL ? mesh->EVT[0] + n_half_edges : NIL;

gproshan_debug(fill_holes);
	index_t * t_aEVT = aEVT + n_vertices;
	for(index_t v = ncv; v < mesh->n_vertices; ++v)
		t_aEVT[v - ncv] = mesh->EVT[v] != NIL ? mesh->EVT[v] + n_half_edges : NIL;
gproshan_debug(fill_holes);

	memcpy(aET, ET, sizeof(index_t) * n_edges);
gproshan_debug(fill_holes);

	bool * common_edge = new bool[mesh->n_edges];
	memset(common_edge, 0, sizeof(bool) * mesh->n_edges);
	for(index_t he_i, i = 0; i < ncv; ++i)
	{
		he_i = mesh->EVT[i];
		if(he_i != NIL && mesh->VT[next(he_i)] < ncv)
			common_edge[mesh->EHT[he_i]] = true;
	}

	index_t ae = n_edges;
	for(index_t e = 0; e < mesh->n_edges; ++e)
		if(!common_edge[e])
			aET[ae++] = mesh->ET[e] + n_half_edges;

	gproshan_debug_var(ae == ne);
	gproshan_debug_var(ae);
	gproshan_debug_var(ne);
	assert(ae == ne);
	delete [] common_edge;
gproshan_debug(fill_holes);
	free();

gproshan_debug(fill_holes);
	GT	= aGT;
	VT	= aVT;
	OT	= aOT;
	EVT	= aEVT;
	ET	= aET;
	EHT	= aEHT;

	rw(n_vertices)		= nv;
	rw(n_faces)			= nf;
	rw(n_half_edges)	= nh;
	rw(n_edges)			= ne;

	update_eht();
}

void che::set_head_vertices(index_t * head, const size_t & n)
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

		swap(GT[v], GT[i]);

		for_star(he, this, v)
			VT[he] = i;
		for_star(he, this, i)
			VT[he] = v;

		swap(EVT[v], EVT[i]);
	}
}

index_t che::link_intersect(const index_t & v_a, const index_t & v_b)
{
	index_t intersect = 0;

	for(index_t & a: link(v_a))
	for(index_t & b: link(v_b))
		if(a == b) ++intersect;

	return intersect;
}

void che::edge_collapse(const index_t *const & sort_edges)
{
}

void che::init(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f)
{
	alloc(n_v, n_f);

	memcpy(GT, vertices, n_vertices * sizeof(vertex));
	memcpy(VT, faces, n_half_edges * sizeof(index_t));

	update_evt_ot_et();
	update_eht();
}

void che::init(const string & file)
{
	filename = file;
	read_file(filename);

	update_evt_ot_et();
	update_eht();
}

void che::alloc(const size_t & n_v, const size_t & n_f)
{
	rw(n_vertices)		= n_v;
	rw(n_faces)			= n_f;
	rw(n_half_edges)	= che::mtrig * n_faces;
	rw(n_edges)			= n_half_edges;				// max number of edges

	if(n_vertices)		GT	= new vertex[n_vertices];
	if(n_half_edges)	VT	= new index_t[n_half_edges];
	if(n_half_edges)	OT	= new index_t[n_half_edges];
	if(n_vertices)		EVT	= new index_t[n_vertices];
	if(n_half_edges)	ET	= new index_t[n_half_edges];
	if(n_half_edges)	EHT	= new index_t[n_half_edges];

	if(n_vertices)		VN	= new vertex[n_vertices];
	if(n_vertices)		VC	= new rgb_t[n_vertices];
	if(n_vertices)		VHC	= new real_t[n_vertices];

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		VC[v] = rgb_t();

	update_heatmap();
}

void che::free()
{
	delete [] GT;	GT = nullptr;
	delete [] VT;	VT = nullptr;
	delete [] OT;	OT = nullptr;
	delete [] EVT;	EVT = nullptr;
	delete [] ET;	ET = nullptr;
	delete [] EHT;	EHT = nullptr;
	delete [] VN;	VN = nullptr;
	delete [] VC;	VC = nullptr;
	delete [] VHC;	VHC = nullptr;
}

void che::read_file(const string & ) {}		/* virtual */

void che::update_evt_ot_et()
{
	if(!n_faces) return;

	memset(EVT, -1, sizeof(index_t) * n_vertices);
	memset(OT, -1, sizeof(index_t) * n_half_edges);

	vector<index_t> vnhe;
	vnhe.assign(n_vertices, 0);

	for(index_t he = 0; he < n_half_edges; ++he)
		++vnhe[VT[he]];

	vector<index_t *> vhe(n_vertices);
	vhe[0] = new index_t[n_half_edges];
	for(index_t v = 1; v < n_vertices; ++v)
		vhe[v] = vhe[v - 1] + vnhe[v - 1];

	vnhe.assign(n_vertices, 0);
	for(index_t he = 0; he < n_half_edges; ++he)
		vhe[VT[he]][vnhe[VT[he]]++] = he;

	size_t ne = 0;
	for(index_t ohe, he = 0; he < n_half_edges; ++he)
	{
		const index_t & u = VT[he];
		const index_t & v = VT[next(he)];

		EVT[u] = he;

		if(OT[he] == NIL)
		{
			ohe = NIL;
			for(index_t j = 0; j < vnhe[v]; ++j)
			{
				index_t & h = vhe[v][j];
				if(VT[next(h)] == u)
				{
					ohe = h;
					break;
				}
			}

			if(ohe != NIL && OT[ohe] == NIL)
			{
				ET[ne++] = he;
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


// static

vector<index_t> che::trig_convex_polygon(const index_t * P, const size_t & n)
{
	vector<index_t> trigs;

	trigs.reserve(che::mtrig * (n - 2));
	for(index_t i = 2; i < n; ++i)
	{
		trigs.push_back(P[0]);
		trigs.push_back(P[i - 1]);
		trigs.push_back(P[i]);
	}

	return trigs;
}


} // namespace gproshan

