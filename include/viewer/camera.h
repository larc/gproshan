#ifndef VIEWER_CAMERA_H
#define VIEWER_CAMERA_H

#include "quaternion.h"

#include "include_opengl.h"


// geometry processing and shape analysis framework
namespace gproshan {


class camera
{
	public:
		quaternion pClick;
		quaternion pDrag;
		quaternion pLast;
		quaternion rLast;
		quaternion momentum;
		int tLast;
		double zoom, vZoom;

	public:
		camera();
		quaternion clickToSphere(int x, int y);
		void setView() const;
		void mouse(int button, int state, int x, int y);
		void motion(int x, int y);
		void idle();
		void zoomIn();
		void zoomOut();
		quaternion currentRotation() const;
};


} // namespace gproshan

#endif

