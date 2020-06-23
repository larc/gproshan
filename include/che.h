#ifndef CHE_H
#define CHE_H

#include "include.h"
#include "vertex.h"

#include <vector>
#include <string>

#define for_star(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->ot(prev(he))) != stop ? he : NIL)
#define for_border(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->evt(mesh->vt(next(he)))) != stop ? he : NIL)


// geometry processing and shape analysis framework
namespace gproshan {


typedef std::vector<index_t> star_t;		// star (vector of he)
typedef std::vector<index_t> link_t;		// link (vector of he)

index_t trig(const index_t & he);
index_t next(const index_t & he);
index_t prev(const index_t & he);

struct corr_t;

class che
{
	public:
		static const size_t P = 3; ///< default polygon size 3, triangular meshes

	protected:
		std::string filename_;

		size_t n_vertices_		= 0;
		size_t n_faces_			= 0;
		size_t n_half_edges_	= 0;
		size_t n_edges_			= 0;
		size_t n_borders_		= 0;

		vertex * GT		= nullptr;	///< geometry table			: v		-> vertex
		index_t * VT	= nullptr;	///< vertex table (faces)	: he	-> v
		index_t * OT	= nullptr;	///< opposite table			: he	-> he
		index_t * EVT	= nullptr;	///< extra vertex table		: v		-> he
		index_t * ET	= nullptr;	///< edge table				: e		-> he
		index_t * EHT	= nullptr;	///< extra half edge table	: he	-> e
		index_t * BT	= nullptr;	///< boundary table			: b 	-> v
		
		vertex * VN		= nullptr;	///< vertex normals			: v		-> normal(v)
		real_t * VC		= nullptr;	///< vertex color			: v		-> color(v)

		bool manifold = true;

	public:
		che(const size_t & n_v = 0, const size_t & n_f = 0);
		che(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		che(const che & mesh);		
		virtual ~che();
		
		void star(star_t & s, const index_t & v);
		void link(link_t & l, const index_t & v);
		void border(std::vector<index_t> & border, const index_t & b);
		bool is_border_v(const index_t & v) const;
		bool is_border_e(const index_t & e) const;
		void flip(const index_t & e);
		real_t pdetriq(const index_t & t) const;
		real_t quality();
		real_t real_trig(const index_t & t) const;
		real_t area_vertex(const index_t & v);
		real_t area_surface() const;
		void update_colors(const real_t * vcolor = nullptr);
		const real_t & color(const index_t & v) const;
		real_t & color(const index_t & v);
		void update_normals();
		const vertex & normal(const index_t & v) const;
		vertex & normal(const index_t & v);
		vertex shading_normal(const index_t & f, const float & u, const float & v, const float & w) const;
		vertex normal_trig(const index_t & f) const;
		vertex normal_he(const index_t & he) const;
		vertex gradient_he(const index_t & he, const real_t *const & f) const;
		vertex gradient(const index_t & v, const real_t *const & f);
		vertex barycenter(const index_t & t) const;
		vertex corr_vertex(corr_t & corr) const;
		real_t cotan(const index_t & he) const;
		real_t mean_edge() const;
		size_t memory() const;
		size_t genus() const;

		void normalize();
		bool is_pointcloud() const;
		bool is_manifold() const;
		bool is_pointcloud() const;
		const index_t & vt(const index_t & he) const;
		const vertex & gt(const index_t & v) const;
		const vertex & gt_vt(const index_t & he) const;
		const vertex & gt_vt_next_evt(const index_t & v) const;
		const vertex & gt_e(const index_t & e, const bool & op = false);
		const index_t & vt_e(const index_t & e, const bool & op = false);
		const index_t & et(const index_t & e) const;
		const index_t & ot_et(const index_t & e) const;
		const index_t & ot(const index_t & he) const;
		const index_t & ot_evt(const index_t & v) const;
		const index_t & evt(const index_t & v) const;
		const index_t & bt(const index_t & b) const;
		const size_t & n_vertices() const;
		const size_t & n_faces() const;
		const size_t & n_half_edges() const;
		const size_t & n_edges() const;
		const size_t & n_borders() const;
		size_t max_degree() const;
		vertex & get_vertex(index_t v);
		void set_vertices(const vertex *const& positions, size_t n = 0, const index_t & v_i = 0);
		void set_filename(const std::string & f);
		const std::string & filename() const;
		const std::string filename_size() const;
		const std::string name() const;
		const std::string name_size() const;
		void reload();
		void compute_toplesets(index_t *& rings, index_t *& sorted, std::vector<index_t> & limites, const std::vector<index_t> & sources, const index_t & k = NIL);
		void multiplicate_vertices();
		void remove_non_manifold_vertices();
		void remove_vertices(const std::vector<index_t> & vertices);
		void merge(const che * mesh, const std::vector<index_t> & com_vertices);
		void set_head_vertices(index_t * head, const size_t & n);
		index_t link_intersect(const index_t & v_a, const index_t & v_b);
		corr_t * edge_collapse(const index_t *const & sort_edges, const vertex *const & normals);
		corr_t find_corr(const vertex & v, const vertex & n, const std::vector<index_t> & triangles);

	protected:
		void delete_me();
		void init(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		void init(const std::string & file);
		void init(const size_t & n_v, const size_t & n_f);
		virtual void read_file(const std::string & file);

	private:
		void update_evt_ot_et();
		void update_eht();
		void update_bt();

	friend struct CHE;
};

struct vertex_cu;

struct CHE
{
	size_t n_vertices;
	size_t n_faces;
	size_t n_half_edges;

	vertex_cu * GT;
	index_t * VT;
	index_t * OT;
	index_t * EVT;

	CHE(che * mesh);
};

struct corr_t
{
	index_t t;
	vertex alpha;

	corr_t()
	{
		t = NIL;
	}

	void init(const index_t & he)
	{
		t = trig(he) * che::P;
		alpha[0] = he == t;
		alpha[1] = he == next(t);
		alpha[2] = he == prev(t);
		t /= che::P;
	}
};


} // namespace gproshan

#endif // CHE_H

