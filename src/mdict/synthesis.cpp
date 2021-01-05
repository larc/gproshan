#include "mdict/synthesis.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


synthesis::synthesis(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const real_t & _f, const bool & _learn, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _learn, _plot)
{
}

real_t synthesis::execute()
{
	TIC(d_time) init_sampling(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) init_patches(); TOC(d_time)
	gproshan_debug_var(d_time);

	//TIC(d_time) learning(); TOC(d_time)
	/*
	string name;
	d_message(Dictionary name:)
	cin>>name;
	string f_dict = "tmp/" + name + ".dict";
	debug(f_dict)
	d_message(loading dictionary)
	if(!A.load(f_dict))
	{ 
		d_message(This dictionary does not exist Bye) return;
	}
	gproshan_debug_var(d_time);*/
	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	gproshan_debug_var(d_time);

	return 0;
}


} // namespace gproshan::mdict

