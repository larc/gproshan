#ifndef VIEWER_CAMERA_H
#define VIEWER_CAMERA_H

#include "mesh/quaternion.h"

#include "viewer/include_opengl.h"


// geometry processing and shape analysis framework
namespace gproshan {


class camera
{
	private:
		quaternion p_click;
		quaternion p_drag;
		quaternion p_last;
		quaternion r_last;
		int t_last;

	public:
		double zoom;

	public:
		camera();
		void mouse(int button, int state, int x, int y);
		void motion(int x, int y);
		void idle();
		void zoom_in();
		void zoom_out();
		quaternion current_rotation() const;
	
	private:
		quaternion click_to_sphere(int x, int y);
	
	friend std::ostream & operator << (std::ostream & os, const camera & cam);
	friend std::istream & operator >> (std::istream & is, camera & cam);
};


} // namespace gproshan

#endif

