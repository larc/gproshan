#ifndef CHE_H
#define CHE_H

#include <vector>
#include <string>

#include "include.h"
#include "vertex.h"

#define for_star(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->ot(prev(he))) != stop ? he : NIL)
#define for_border(he, mesh, v) for(index_t stop = mesh->evt(v), he = mesh->evt(v); he != NIL; he = (he = mesh->evt(mesh->vt(next(he)))) != stop ? he : NIL)


const size_t P = 3;

using namespace std;

typedef vector<index_t> star_t;		//star (vector of he)
typedef vector<index_t> link_t;		//link (vector of he)

index_t trig(const index_t & he);
index_t next(const index_t & he);
index_t prev(const index_t & he);

struct corr_t
{
	index_t t;
	vertex alpha;
	
	void init(const index_t & he)
	{
		t = trig(he) * P;
		alpha[0] = he == t;
		alpha[1] = he == next(t);
		alpha[2] = he == prev(t);
		t /= P;
	}
};

class che
{
	protected:
		string filename_;

		size_t n_vertices_;
		size_t n_faces_;
		size_t n_half_edges_;
		size_t n_edges_;
		size_t n_borders_;

		vertex * GT;	//geometry table		v	-> vertex
		index_t * VT;	//vertex table (faces)	he	-> v
		index_t * OT;	//opposite table		he	-> he
		index_t * EVT;	//extra vertex table	v	-> he
		index_t * ET;	//edge table			e	-> he
		index_t * EHT;	//extra half edge table	he	-> e
		index_t * BT;	//boundary table		b 	-> v
		
		bool manifold;

	public:
		virtual ~che();
		void star(star_t & s, const index_t & v);
		void link(link_t & l, const index_t & v);
		void border(vector<index_t> & border, const index_t & b);
		bool is_border_v(const index_t & v) const;
		bool is_border_e(const index_t & e) const;
		void flip(const index_t & e);
		vertex_t pdetriq(const index_t & t) const;
		percent_t quality();
		area_t area_trig(const index_t & t) const;
		area_t area_vertex(const index_t & v);
		area_t area_surface() const;
		vertex normal_he(const index_t & he) const;
		vertex normal(const index_t & v);
		vertex gradient_he(const index_t & he, const distance_t *const & f) const;
		vertex gradient(const index_t & v, const distance_t *const & f);
		vertex barycenter(const index_t & t) const;
		area_t cotan(const index_t & he) const;
		vertex_t mean_edge() const;
		size_t memory() const;
		size_t genus() const;

		void normalize();
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
		const size_t & n_vertices() const;
		const size_t & n_faces() const;
		const size_t & n_half_edges() const;
		const size_t & n_edges() const;
		const size_t & n_borders() const;
		vertex & get_vertex(index_t v);
		void set_vertices(const vertex *const& positions, size_t n = 0, const index_t & v_i = 0);
		const string & filename() const;
		void set_filename(const string & f);
		string name() const;
		void reload();
		void sort_by_rings(index_t *& rings, index_t *& sorted, vector<index_t> & limites, const vector<index_t> & sources, const index_t & k = NIL);
		void multiplicate_vertices();
		void remove_non_manifold_vertices();
		void remove_vertices(const vector<index_t> & vertices);
		void merge(const che * mesh, const vector<index_t> & com_vertices);
		void set_head_vertices(index_t * head, const size_t & n);
		index_t link_intersect(const index_t & v_a, const index_t & v_b);
		corr_t * edge_collapse(const index_t *const & sort_edges);
		corr_t find_corr(const vertex & v, const vector<index_t> & triangles);

	protected:
		void delete_me();
		void init(const vertex * vertices, const index_t & n_v, const index_t * faces, const index_t & n_f);
		void init(const string & file);
		void init(const size_t & n_v, const size_t & n_f);
		virtual void read_file(const string & file) = 0;

	public:
		virtual void write_file(const string & file) const = 0;

	private:
		void update_evt_ot_et();
		void update_eht();
		void update_bt();

	friend struct CHE;
};

struct vertex_cu;

struct CHE
{
	length_t n_vertices;
	length_t n_faces;
	length_t n_half_edges;

	vertex_cu * GT;
	index_t * VT;
	index_t * OT;
	index_t * EVT;

	CHE(che * mesh);
};

#endif // CHE_H

