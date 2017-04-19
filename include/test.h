#ifndef TEST_H
#define TEST_H

#include <string>

using namespace std;

class test
{
	private:
		time_t time_total;

	private:
		virtual void run(string database, string sub) = 0;

	public:
		test();
		virtual ~test();
		void execute(string database, string sub);
};

string getFileName(string file);

#endif // TEST_H
