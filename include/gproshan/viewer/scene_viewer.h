#ifndef SCENE_VIEWER_H
#define SCENE_VIEWER_H

#include <gproshan/viewer/che_viewer.h>
#include <gproshan/scenes/scene.h>


// geometry processing and shape analysis framework
namespace gproshan {


class scene_viewer: public che_viewer
{
	private:
		scene * sc = nullptr;

	public:
		scene_viewer(scene * p_sc);
		void draw(shader & program);
};


} // namespace gproshan

#endif // SCENE_VIEWER_H

