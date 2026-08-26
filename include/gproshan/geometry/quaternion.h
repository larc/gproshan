#ifndef QUATERNION_H
#define QUATERNION_H

#include <gproshan/geometry/vec.h>

#include <ostream>


// geometry processing and shape analysis framework
namespace gproshan {


class quaternion
{
	private:
		float s = 0;
		vec3 v;

	public:
		quaternion(const vec3 & v = {});
		quaternion(float s, const vec3 & v = {});

		operator const vec3 & () const;
		float & operator [] (int index);
		float operator [] (int index) const;

		quaternion operator + (const quaternion & q) const;
		quaternion operator - (const quaternion & q) const;
		quaternion operator - () const;
		quaternion operator * (float c) const;
		quaternion operator / (float c) const;
		quaternion operator * (const quaternion & q) const;

		quaternion conj() const;
		quaternion inv() const;
		float norm() const;
		float norm2() const;

	friend std::ostream & operator << (std::ostream & os, const quaternion & q);
	friend std::istream & operator >> (std::istream & is, quaternion & q);
};


float norm(const quaternion & q);
quaternion normalize(const quaternion & q);
quaternion operator * (float c, const quaternion & q);


} // namespace gproshan

#endif // QUATERNION_H

