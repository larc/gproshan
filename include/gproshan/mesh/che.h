#ifndef CHE_H
#define CHE_H

#include <gproshan/include.h>
#include <gproshan/geometry/vec.h>
#include <gproshan/geometry/mat.h>

#include <vector>
#include <string>


// geometry processing and shape analysis framework
namespace gproshan {


using vertex = vec3;

class che
{
	protected:
		static size_t & rw(const size_t & n);

	private:
		class star_he;

	public:
		enum mesh_type : unsigned char { mtrig = 3, mquad = 4 };	///< meshes_types

		struct rgb_t
		{
			unsigned char r = 230;
			unsigned char g = 240;
			unsigned char b = 250;

			rgb_t() = default;
			rgb_t(const float & fr, const float & fg, const float & fb);
			unsigned char & operator [] (const index_t & i);
			operator vertex () const;
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
		real_t * VHC	= nullptr;	///< vertex color heatmap	: v		-> heatmap(v)

		bool manifold = true;

	public:
		che(const che & mesh);
		che(const size_t & n_v = 0, const size_t & n_f = 0);
		che(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		virtual ~che();

		// vertex access geometry methods to xyz point values, normals, and gradient
		const vertex & point(const index_t & v) const;
		vertex & point(const index_t & v);
		const vertex & normal(const index_t & v) const;
		vertex & normal(const index_t & v);
		vertex shading_normal(const index_t & f, const float & u, const float & v) const;
		vertex normal_trig(const index_t & f) const;
		vertex normal_he(const index_t & he) const;
		vertex gradient_he(const index_t & he, const real_t * f) const;
		vertex gradient(const index_t & v, const real_t * f);

		// vertex color methods
		const real_t & heatmap(const index_t & v) const;
		real_t & heatmap(const index_t & v);
		const rgb_t & rgb(const index_t & v) const;
		rgb_t & rgb(const index_t & v);
		vertex color(const index_t & v) const;
		vertex shading_color(const index_t & f, const float & u, const float & v) const;

		// update methods
		void reload();
		mat4 normalize_sphere(const real_t & r = 1) const;
		mat4 normalize_box(const real_t & side = 2) const;
		che * merge(const che * mesh, const std::vector<index_t> & com_vertices);
		void update_vertices(const vertex * positions, const size_t & n = 0, const index_t & v_i = 0);
		void update_heatmap(const real_t * hm = nullptr);
		void update_normals();
		void invert_normals();
		void multiplicate_vertices();
		void remove_vertices(const std::vector<index_t> & vertices);
		void remove_non_manifold_vertices();
		void set_head_vertices(index_t * head, const size_t & n);

		// half edge access methods triangular faces and navigation
		const index_t & halfedge(const index_t & he) const;
		const index_t & twin_he(const index_t & he) const;
		const index_t &	edge_u(const index_t & e) const;
		const index_t &	edge_v(const index_t & e) const;
		const index_t &	edge_he_0(const index_t & e) const;
		const index_t &	edge_he_1(const index_t & e) const;
		const vertex & vertex_he(const index_t & he) const;
		const vertex & vertex_edge_u(const index_t & e) const;
		const vertex & vertex_edge_v(const index_t & e) const;
		const index_t & evt(const index_t & v) const;

		// topology methods
		che::star_he star(const index_t & v) const;
		std::vector<index_t> link(const index_t & v) const;
		void edge_collapse(const std::vector<index_t> & sort_edges);
		void compute_toplesets(index_t *& rings, index_t *& sorted, std::vector<index_t> & limites, const std::vector<index_t> & sources, const index_t & k = NIL);

		// boundary methods
		std::vector<index_t> bounds() const;
		std::vector<index_t> boundary(const index_t & v) const;
		bool is_vertex_bound(const index_t & v) const;
		bool is_edge_bound(const index_t & e) const;

		// file, name, and system methods
		const std::string name() const;
		const std::string name_size() const;
		const std::string filename_size() const;

		// mesh information methods
		size_t genus() const;
		size_t memory() const;
		size_t max_degree() const;
		real_t quality() const;
		real_t mean_edge() const;
		real_t area_surface() const;
		bool is_manifold() const;
		bool is_pointcloud() const;

		// operation methods
		void flip(const index_t & e);
		real_t cotan(const index_t & he) const;
		real_t pdetriq(const index_t & t) const;
		real_t area_trig(const index_t & t) const;
		real_t area_vertex(const index_t & v) const;
		real_t mean_curvature(const index_t & v) const;

	protected:
		void init(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		void init(const std::string & file);
		void alloc(const size_t & n_v, const size_t & n_f);
		void free();

		virtual void read_file(const std::string & file);

	private:
		void update_evt_ot_et();
		void update_eht();

	public:
		static std::vector<index_t> trig_convex_polygon(const index_t * P, const size_t & n);

	friend struct CHE;
};

class che::star_he
{
	class iterator;

	const che * mesh;
	const index_t & v;

	public:
		star_he(const che * p_mesh, const index_t & p_v);
		iterator begin() const;
		iterator end() const;
};

class che::star_he::iterator
{
	const che * mesh;
	index_t he;
	const index_t & he_end;

	public:
		iterator(const che * p_mesh, const index_t & p_he, const index_t & p_he_end);
		iterator & operator ++ ();
		bool operator != (const iterator & it) const;
		const index_t & operator * ();
};


// che halfedge functions

inline index_t trig(const index_t & he)
{
	if(he == NIL) return NIL;
	return he / che::mtrig;
}

inline index_t next(const index_t & he)
{
	if(he == NIL) return NIL;
	return che::mtrig * trig(he) + (he + 1) % che::mtrig;
}

inline index_t prev(const index_t & he)
{
	if(he == NIL) return NIL;
	return che::mtrig * trig(he) + (he + che::mtrig - 1) % che::mtrig;
}


// simple che data structure
struct CHE
{
	size_t n_vertices = 0;
	size_t n_faces = 0;
	size_t n_half_edges = 0;

	vertex * GT	= nullptr;
	vertex * VN	= nullptr;
	che::rgb_t * VC	= nullptr;
	index_t * VT	= nullptr;
	index_t * OT	= nullptr;
	index_t * EVT	= nullptr;

	CHE() = default;
	CHE(const che * mesh);
};


} // namespace gproshan

#endif // CHE_H

