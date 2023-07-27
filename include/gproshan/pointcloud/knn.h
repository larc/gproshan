#ifndef KNN_H
#define KNN_H

#include <gproshan/geometry/vec.h>
#include <gproshan/geometry/mat.h>

#include <vector>
#include <unordered_map>


inline gproshan::uvec3 hash(gproshan::vec3 p, const float & res = 1000)
{
	p = (res - 1) * (0.5f * p + 1);
	return {(unsigned int) p.x(), (unsigned int) p.y(), (unsigned int) p.z()};
}


template<>
struct std::hash<gproshan::uvec3>
{
	std::size_t operator () (const gproshan::uvec3 & p) const noexcept
	{
		return p.x() * 1e6 + p.y() * 1e3 + p.z();
	}
};


// geometry processing and shape analysis framework
namespace gproshan {


using point = vec3;


class knn	// grid
{
	private:
		std::vector<point> points;
		std::unordered_map<uvec3, std::vector<index_t> > grid;

	public:
		knn(const point * pc, const size_t & n_points, const mat4 & transform);
		virtual ~knn() = default;

	private:
};


} // namespace gproshan

#endif // KNN_H

