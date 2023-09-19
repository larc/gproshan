#ifndef KNN_H
#define KNN_H

#include <gproshan/geometry/vec.h>
#include <gproshan/geometry/mat.h>

#include <vector>
#include <unordered_map>

#include <flann/flann.hpp>


inline gproshan::uvec3 hash(gproshan::vec3 p, const float & res)
{
	p = (res - 1) * (0.5f * p + 1);
	return {(unsigned int) p.x(), (unsigned int) p.y(), (unsigned int) p.z()};
}


template<>
struct std::hash<gproshan::uvec3>
{
	std::size_t operator () (const gproshan::uvec3 & p) const noexcept
	{
		return p.x() * 1e12 + p.y() * 1e6 + p.z();
	}
};


// geometry processing and shape analysis framework
namespace gproshan::knn {


using point = vec3;


class grid
{
	private:
		float res = 1000;
		std::vector<point> points;
		std::unordered_map<uvec3, std::vector<index_t> > voxels;

	public:
		grid(const point * pc, const size_t & n_points, const mat4 & transform);
		~grid() = default;

		std::vector<index_t> operator () (const point & p, int k);
};


class k3tree
{
	private:
		flann::Matrix<int> indices;
	
	public:
		k3tree(const point * pc, const size_t & n_points, const size_t & k = 8, const std::vector<point> & query = {});
		~k3tree();

		int * operator () (const index_t & i);
};


} // namespace gproshan

#endif // KNN_H

