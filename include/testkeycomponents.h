#ifndef TESTKEYCOMPONENTS_H
#define TESTKEYCOMPONENTS_H

#include "test.h"
#include "keycomponents.h"

class testkeycomponents : public test
{
	private:
		size_t k;
		percent_t percent_kps;

	private:
		void run(string database, string sub);

	public:
		testkeycomponents(size_t k_ = 45, percent_t pkps = 0.10);
		void configure(size_t k_, percent_t pkps);
		void one_test(string file, string sub);
		void one_test_fm(string file, string sub);
		~testkeycomponents();

};

#endif // TESTKEYCOMPONENTS_H

