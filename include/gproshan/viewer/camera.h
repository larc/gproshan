#ifndef CAMERA_H
#define CAMERA_H

#include "mesh/quaternion.h"

#include <glm/glm.hpp>


// geometry processing and shape analysis framework
namespace gproshan {


class camera
{
	private:
		quaternion p_click;
		quaternion p_drag;
		quaternion p_last;
		quaternion r_last;

	public:
		quaternion eye;
		quaternion center;
		quaternion up;
		double zoom;

	public:
		camera();
		glm::mat4 look_at(const quaternion & r);
		quaternion current_rotation() const;
		void mouse(const bool & press, const double & x, const double & y, const int & w, const int & h);
		void motion(const double & x, const double & y, const int & w, const int & h);
		void zoom_in();
		void zoom_out();

	private:
		quaternion click_to_sphere(const double & x, const double & y, const int & w, const int & h);

	friend std::ostream & operator << (std::ostream & os, const camera & cam);
	friend std::istream & operator >> (std::istream & is, camera & cam);
};


} // namespace gproshan

#endif // CAMERA_H

