#include "testkeypoints.h"

testkeypoints::testkeypoints(percent_t _pkps)
{
	pkps = _pkps;
}

void testkeypoints::run(string database, string sub)
{
	string file;
	ifstream is(PATH_DATA + database);
	while(is>>file)
	{
		off shape(PATH_DATA + file);
		size_t nkps = shape.get_nvertices()*pkps;

		keypoints kps(shape);
		ofstream oskp(PATH_KPS + sub + "/"+ getFileName(file) + "kp");
		kps.print(oskp, nkps);
		oskp.close();
	}
}
