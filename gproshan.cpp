#include "app_viewer.h"
#include "geometry/convex_hull.h"

using namespace std;

int main(int nargs, const char ** args)
{
	std::vector<gproshan::vertex> vv{	{0, 2, 0},
										{1, 0, 0},
										{2, 5, 0},
										{3, 3, 0},
										{4, -1, 0},
										{2, 2, 0},
										{5, 3, 0},
										{6, 4, 0}	};
			
	gproshan::convex_hull ch(vv);
	const std::vector<gproshan::vertex> & p = ch;
	for(auto v: p)
		gproshan_log_var(v);

	gproshan::app_viewer app;
	return app.main(nargs, args);
}

