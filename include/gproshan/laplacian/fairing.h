#ifndef FAIRING_H
#define FAIRING_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class fairing
{
	protected:
		vertex * vertices = nullptr;

	public:
		fairing() = default;
		virtual ~fairing();
		void run(che * mesh);
		const vertex * new_vertices();

	protected:
		virtual void compute(che * mesh) = 0;
};


} // namespace gproshan

#endif // FAIRING_H

