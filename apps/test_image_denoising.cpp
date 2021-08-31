#include "mdict/image_denoising.h"

int main(int nargs, const char ** args)
{
	for(int i = 1; i < nargs; ++i)
		gproshan::mdict::test_image_denoising(args[i]);

	return 0;
}

