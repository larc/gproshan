#ifndef CAMERA_H
#define CAMERA_H

#include <gproshan/mesh/quaternion.h>
#include <gproshan/geometry/mat.h>


// geometry processing and shape analysis framework
namespace gproshan {


class camera
{
	private:
		quaternion p_click	= 1;
		quaternion p_drag	= 1;
		quaternion p_last	= 1;
		quaternion r_last	= 1;

	public:
		quaternion eye;
		quaternion pos		= vertex{0, 0, -3.14};
		quaternion front	= vertex{0, 0, 1};
		quaternion up		= vertex{0, 1, 0};

	public:
		static mat4 perspective(const real_t & fovy, const real_t & aspect, const real_t & near, const real_t & far);
		mat4 look_at(const quaternion & r);
		quaternion current_rotation() const;
		void mouse(const bool & press, const double & x, const double & y, const int & w, const int & h);
		void motion(const double & x, const double & y, const int & w, const int & h);
		void zoom_in();
		void zoom_out();
		real_t zoom() const;

	private:
		quaternion click_to_sphere(const double & x, const double & y, const int & w, const int & h);

	friend std::ostream & operator << (std::ostream & os, const camera & cam);
	friend std::istream & operator >> (std::istream & is, camera & cam);
};


} // namespace gproshan

#endif // CAMERA_H

