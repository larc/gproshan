#include "mdict/denoising.h"
#include "che_off.h"

int main(int nargs, const char ** args)
{
 

    string file = args[1];
    gproshan::mdict::test_image_denoising(args[1]);


	return 0;
}
