#include "viewer/camera.h"

#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


camera::camera(): p_click(1), p_drag(1), p_last(1), r_last(1), zoom(1) {}

quaternion camera::click_to_sphere(int x, int y, int w, int h)
{
	quaternion p(	0.,
					2. * (double) x / (double) w - 1.,
					2. * (double) y / (double) h - 1.,
					0.);

	if(p.norm2() > 1.)
	{
		p.normalize();
		p.im().z = 0.;
	}
	else
	{
		p.im().z = sqrt(1. - p.norm2());
	}

	return p;
}

quaternion camera::current_rotation() const
{
	return (p_drag * p_click.conj()) * r_last;
}

void camera::mouse(int, int state, int x, int y, int w, int h)
{
	quaternion momentum = 1;

	if(state == GLFW_PRESS)
		p_click = p_drag = p_last = click_to_sphere(x, y, w, h);

	if(state == GLFW_RELEASE)
	{
		double timeSinceDrag = (clock() - t_last) / (double) CLOCKS_PER_SEC;

		if(timeSinceDrag < .1)
		{
			momentum = p_drag * p_last.conj();
			momentum = (.03 * momentum + .97).unit();
		}
		else
		{
			momentum = 1.;
		}

		r_last = p_drag * p_click.conj() * r_last;
		p_click = p_drag = 1.;
	}
}

void camera::motion(int x, int y, int w, int h)
{
	t_last = clock();
	p_last = p_drag;
	p_drag = click_to_sphere(x, y, w, h);
}

void camera::zoom_in()
{
	zoom -= 0.01;
}

void camera::zoom_out()
{
	zoom += 0.01;
}

ostream & operator << (ostream & os, const camera & cam)
{
	return os << cam.p_click << "\n"
			<< cam.p_drag << "\n"
			<< cam.p_last << "\n"
			<< cam.r_last << "\n"
			<< cam.t_last << "\n"
			<< cam.zoom << "\n";
}

istream & operator >> (istream & is, camera & cam)
{
	return is >> cam.p_click
			>> cam.p_drag
			>> cam.p_last
			>> cam.r_last
			>> cam.t_last
			>> cam.zoom;
}


} // namespace gproshan

