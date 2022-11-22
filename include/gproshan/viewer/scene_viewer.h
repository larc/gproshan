#ifndef SCENE_VIEWER_H
#define SCENE_VIEWER_H

#include <gproshan/viewer/che_viewer.h>


// geometry processing and shape analysis framework
namespace gproshan {


class scene_viewer: public che_viewer
{
	public:
		virtual void draw(shader & program);
};


} // namespace gproshan

#endif // SCENE_VIEWER_H

