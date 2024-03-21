#ifndef CHE_CUDA_H
#define CHE_CUDA_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class che_cuda : public che
{
	che * device = nullptr;

	public:
		che_cuda(const che * mesh = nullptr, const che::options & opts = che::default_opts);
		~che_cuda();

		operator const che * () const;
};


} // namespace gproshan

#endif // CHE_CUDA_H

