#ifndef VIEWER_CAMERA_H
#define VIEWER_CAMERA_H

#include "quaternion.h"

#include "include_opengl.h"


// geometry processing and shape analysis framework
namespace gproshan {


class camera
{
	private:
		quaternion pClick;
		quaternion pDrag;
		quaternion pLast;
		quaternion rLast;
		quaternion momentum;
		int tLast;
		double vZoom;

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
};


} // namespace gproshan

#endif

