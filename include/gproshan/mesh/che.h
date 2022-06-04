#ifndef CHE_H
#define CHE_H

#include <gproshan/include.h>
#include <gproshan/mesh/vertex.h>

#include <vector>
#include <string>

#define for_star(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->ot(prev(he))) != stop ? he : NIL)
#define for_boundary(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->evt(mesh->vt(next(he)))) != stop ? he : NIL)


// geometry processing and shape analysis framework
namespace gproshan {


size_t & rw(const size_t & n);
index_t trig(const index_t & he);
index_t next(const index_t & he);
index_t prev(const index_t & he);


class che
{
	public:
		enum mesh_type : unsigned char { mtrig = 3, mquad = 4 };	///< meshes_types
		struct rgb_t
		{
			unsigned char r = 230;
			unsigned char g = 240;
			unsigned char b = 250;

			unsigned char & operator [] (const index_t & i)
			{
				return (&r)[i];
			}

			operator vertex () const
			{
				return {float(r) / 255, float(g) / 255, float(b) / 255};
			}
		};

		const size_t n_vertices		= 0;
		const size_t n_faces		= 0;
		const size_t n_half_edges	= 0;
		const size_t n_edges		= 0;

		std::string filename;		///< get and set data member

	protected:
		vertex * GT		= nullptr;	///< geometry table			: v		-> vertex
		index_t * VT	= nullptr;	///< vertex table (faces)	: he	-> v
		index_t * OT	= nullptr;	///< opposite table			: he	-> he
		index_t * EVT	= nullptr;	///< extra vertex table		: v		-> he
		index_t * ET	= nullptr;	///< edge table				: e		-> he
		index_t * EHT	= nullptr;	///< extra half edge table	: he	-> e

		vertex * VN		= nullptr;	///< vertex normals			: v		-> normal(v)
		rgb_t * VC		= nullptr;	///< vertex color			: v		-> color(v)
		real_t * VHC	= nullptr;	///< vertex color heat map	: v		-> heatmap(v)

		bool manifold = true;

	public:
		che(const che & mesh);
		che(const size_t & n_v = 0, const size_t & n_f = 0);
		che(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		virtual ~che();

		std::vector<index_t> star(const index_t & v) const;
		std::vector<index_t> link(const index_t & v) const;
		std::vector<index_t> bounds() const;
		std::vector<index_t> boundary(const index_t & v) const;
		bool is_vertex_bound(const index_t & v) const;
		bool is_edge_bound(const index_t & e) const;
		void flip(const index_t & e);
		real_t pdetriq(const index_t & t) const;
		real_t quality() const;
		real_t area_trig(const index_t & t) const;
		real_t area_vertex(const index_t & v) const;
		real_t area_surface() const;
		void update_heatmap(const real_t * hm = nullptr);
		const rgb_t & rgb(const index_t & v) const;
		rgb_t & rgb(const index_t & v);
		vertex color(const index_t & v) const;
		vertex shading_color(const index_t & f, const float & u, const float & v, const float & w) const;
		const real_t & heatmap(const index_t & v) const;
		real_t & heatmap(const index_t & v);
		void update_normals();
		void invert_normals();
		const vertex & normal(const index_t & v) const;
		vertex & normal(const index_t & v);
		vertex shading_normal(const index_t & f, const float & u, const float & v, const float & w) const;
		vertex normal_trig(const index_t & f) const;
		vertex normal_he(const index_t & he) const;
		vertex gradient_he(const index_t & he, const real_t *const & f) const;
		vertex gradient(const index_t & v, const real_t *const & f);
		vertex barycenter(const index_t & t) const;
		real_t cotan(const index_t & he) const;
		real_t mean_edge() const;
		size_t memory() const;
		size_t genus() const;
		real_t mean_curvature(const index_t & v);

		void normalize();
		bool is_pointcloud() const;
		bool is_manifold() const;
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
		size_t max_degree() const;
		vertex & point(index_t v);
		void set_vertices(const vertex *const& positions, size_t n = 0, const index_t & v_i = 0);
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
		void edge_collapse(const index_t *const & sort_edges);

	protected:
		void init(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		void init(const std::string & file);
		void free();
		void alloc(const size_t & n_v, const size_t & n_f);
		virtual void read_file(const std::string & file);

	private:
		void update_evt_ot_et();
		void update_eht();

	public:
		static std::vector<index_t> trig_convex_polygon(const index_t * P, const size_t & n);

	friend struct CHE;
};

struct vertex_cu;

struct CHE
{
	size_t n_vertices = 0;
	size_t n_faces = 0;
	size_t n_half_edges = 0;

	vertex_cu * GT	= nullptr;
	vertex_cu * VN	= nullptr;
	che::rgb_t * VC	= nullptr;
	index_t * VT	= nullptr;
	index_t * OT	= nullptr;
	index_t * EVT	= nullptr;

	CHE() = default;
	CHE(const che * mesh);
};


} // namespace gproshan

#endif // CHE_H

