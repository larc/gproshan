#include "test.h"

test::test()
{
	time_total = 0;
}

test::~test()
{
	time_total = 0;
}

void test::execute(string database, string sub)
{
	run(database, sub);
}

string getFileName(string file)
{
	return file.substr(0, file.find_last_of(".") + 1);
}
