#include "camera.h"

#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


camera::camera()
: pClick(1.),
	pDrag(1.),
	pLast(1.),
	rLast(1.),
	momentum(1.),
	zoom(1.)
{}

quaternion camera::clickToSphere(int x, int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	int w = viewport[2];
	int h = viewport[3];

	x %= w;
	y %= h;

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
	return (pDrag * pClick.conj()) * rLast;
}

void camera::mouse(int , int state, int x, int y)
{
	if(state == GLFW_PRESS)
	{
		pClick = pDrag = pLast = clickToSphere(x, y);
		momentum = 1.;
	}
	if(state == GLFW_RELEASE)
	{
		double timeSinceDrag = (clock() - tLast) / (double) CLOCKS_PER_SEC;

		if(timeSinceDrag < .1)
		{
			momentum = pDrag * pLast.conj();
			momentum = (.03 * momentum + .97).unit();
		}
		else
		{
			momentum = 1.;
		}

		rLast = pDrag * pClick.conj() * rLast;
		pClick = pDrag = 1.;
	}
}

void camera::motion(int x, int y)
{
	tLast = clock();
	pLast = pDrag;
	pDrag = clickToSphere(x, y);
}

void camera::idle()
{
	// get time since last idle event
	static int t0 = clock();
	int t1 = clock();
	double dt = (t1-t0) / (double) CLOCKS_PER_SEC;

	rLast = momentum * rLast;
	momentum = ((1.-.5*dt) * momentum + .5*dt).unit();

	zoom += vZoom*dt;
	vZoom *= max(0., 1.-5.*dt);

	t0 = t1;
}

void camera::zoom_in()
{
	zoom -= 0.01;
}

void camera::zoom_out()
{
	zoom += 0.01;
}


} // namespace gproshan

