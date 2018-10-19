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

	public:
		key_points(che * mesh);
		~key_points();
		const index_t & operator[](const index_t & i);

	private:
		void compute_kps(che * mesh);
};

#endif // KEY_POINTS_H

