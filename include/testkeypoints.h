#ifndef TESTKEYPOINTS_H
#define TESTKEYPOINTS_H

#include "test.h"
#include "keypoints.h"

class testkeypoints : public test
{
	private:
		percent_t pkps;

	private:
		void run(string database, string sub);

	public:
		testkeypoints(percent_t _pkps = 0);
};

#endif // TESTKEYPOINTS_H
