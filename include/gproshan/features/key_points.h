#ifndef KEY_POINTS_H
#define KEY_POINTS_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class key_points
{
	private:
		std::vector<index_t> kps;
		std::vector<bool> is_kp;

	public:
		key_points(che * mesh, const float percent = 0.10);
		operator const std::vector<index_t> & () const;

	private:
		void compute_kps_areas(che * mesh, const float percent);
};


} // namespace gproshan

#endif // KEY_POINTS_H

