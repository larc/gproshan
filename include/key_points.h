#ifndef KEY_POINTS_H
#define KEY_POINTS_H

#include "che.h"

#include <utility>


typedef std::pair<real_t, index_t> real_idx_t;

class key_points
{
	private:
		size_t n_faces;
		size_t n_vertices;
		real_idx_t * face_areas;
		index_t * kps;
		bool * is_kp;
		size_t n_kps;

	public:
		key_points(che * mesh, const real_t & percent = 0.10);
		~key_points();
		const index_t & operator[](const index_t & i) const;
		const bool & operator()(const index_t & i) const;
		const size_t & size() const;

	private:
		void compute_kps(che * mesh);
};

#endif // KEY_POINTS_H

