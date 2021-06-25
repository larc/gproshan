#include "viewer/camera.h"

#include <cmath>

#include <glm/gtc/matrix_transform.hpp>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


camera::camera(): p_click(1), p_drag(1), p_last(1), r_last(1), zoom(2) {}

glm::mat4 camera::look_at(const quaternion & r)
{
	eye		= vertex(0, 0, -zoom);
	center	= vertex(0, 0, 0);
	up		= vertex(0, 1, 0);

	eye = r.conj() * eye * r;

	return glm::lookAt( glm_vec3(eye),
						glm_vec3(r.conj() * center * r),
						glm_vec3(r.conj() * up * r));
}

quaternion camera::click_to_sphere(const double & x, const double & y, const int & w, const int & h)
{
	quaternion p(0, 2 * x / w - 1, 2 * y / h - 1, 0);

	if(p.norm2() > 1)
	{
		p.normalize();
		p.im().z = 0;
	}
	else
	{
		p.im().z = sqrt(1 - p.norm2());
	}

	return p;
}

quaternion camera::current_rotation() const
{
	return (p_drag * p_click.conj()) * r_last;
}

void camera::mouse(const bool & press, const double & x, const double & y, const int & w, const int & h)
{
	if(press)
	{
		p_click = p_drag = p_last = click_to_sphere(x, y, w, h);
	}
	else
	{
		r_last = p_drag * p_click.conj() * r_last;
		p_click = p_drag = 1.;
	}
}

void camera::motion(const double & x, const double & y, const int & w, const int & h)
{
	p_last = p_drag;
	p_drag = click_to_sphere(x, y, w, h);
}

void camera::zoom_in()
{
	zoom -= 0.02;
}

void camera::zoom_out()
{
	zoom += 0.02;
}

ostream & operator << (ostream & os, const camera & cam)
{
	return os << cam.p_click << "\n"
			<< cam.p_drag << "\n"
			<< cam.p_last << "\n"
			<< cam.r_last << "\n"
			<< cam.zoom << "\n";
}

istream & operator >> (istream & is, camera & cam)
{
	return is >> cam.p_click
			>> cam.p_drag
			>> cam.p_last
			>> cam.r_last
			>> cam.zoom;
}


} // namespace gproshan

