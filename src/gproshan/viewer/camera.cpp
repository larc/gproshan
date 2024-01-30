#include <gproshan/viewer/camera.h>

#include <cmath>


// geometry processing and shape analysis framework
namespace gproshan {


mat4 camera::look_at(const quaternion & r)
{
	eye = r.conj() * pos * r;

	vec3 Z = r.conj() * front * r;
	Z = normalize(Z);
	vec3 Y = r.conj() * up * r;
	Y = normalize(Y);
	vec3 X = cross(Z, Y);

	mat4 view;
	view[0] = {X, -dot(X, eye.v)};
	view[1] = {Y, -dot(Y, eye.v)};
	view[2] = {-Z, dot(Z, eye.v)};
	view[3] = {0, 0, 0, 1};

	return view;
}

mat4 camera::perspective()
{
	return perspective(fovy, aspect, near, far);
}

mat4 camera::perspective(const real_t fovy, const real_t aspect, const real_t near, const real_t far)
{
	const real_t tan_fovy_2 = std::tan(fovy * M_PI / 360);

	mat4 P;
	P(0, 0) = 1 / (aspect * tan_fovy_2);
	P(1, 1) = 1 / tan_fovy_2;
	P(2, 2) = - (far + near) / (far - near);
	P(2, 3) = -2 * far * near / (far - near);
	P(3, 2) = -1;

	return P;
}

quaternion camera::click_to_sphere(const double x, const double y, const int w, const int h)
{
	quaternion p(0, 2 * x / w - 1, 2 * y / h - 1, 0);

	if(p.norm2() > 1)
	{
		p.normalize();
		p.im().z() = 0;
	}
	else
	{
		p.im().z() = sqrt(1 - p.norm2());
	}

	return p;
}

quaternion camera::current_rotation() const
{
	return p_drag * p_click.conj() * r_last;
}

void camera::mouse(const bool press, const double x, const double y, const int w, const int h)
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

void camera::motion(const double x, const double y, const int w, const int h)
{
	p_last = p_drag;
	p_drag = click_to_sphere(x, y, w, h);
}

void camera::zoom_in()
{
	pos.v.z() += 0.02;
}

void camera::zoom_out()
{
	pos.v.z() -= 0.02;
}

real_t camera::zoom() const
{
	return -pos.v.z();
}

std::ostream & operator << (std::ostream & os, const camera & cam)
{
	return os << cam.p_click << "\n"
			<< cam.p_drag << "\n"
			<< cam.p_last << "\n"
			<< cam.r_last << "\n"
			<< cam.pos << "\n";
}

std::istream & operator >> (std::istream & is, camera & cam)
{
	return is >> cam.p_click
			>> cam.p_drag
			>> cam.p_last
			>> cam.r_last
			>> cam.pos;
}


} // namespace gproshan

