#include "d_basis.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void basis::plot_basis()
{
	string file = tmp_file_path("basis.gpi");
	ofstream os(file);
	os << "set term qt size 1000,1000;" << endl;
	os << "set isosamples 50,50;" << endl;
	os << "set parametric;" << endl;
	os << "set vrange [-"<< 0 << ":" << radio <<"];" << endl;
	os << "set urange [-pi:pi];" << endl;
	os << "unset key;" << endl;
	os << "set pm3d at b;" << endl;
	os << "unset colorbox;" << endl;

	plot_basis(os);

	os << "unset multiplot;" << endl;
	os << "pause -1;" << endl;

	os.close();

	file = "gnuplot -persist " + file + " &";

	debug(system(file.c_str()));
}

void basis::plot_atoms(const a_mat & A)
{
	size_t K = A.n_rows;
	size_t m = A.n_cols;
	size_t s = sqrt(m);
	s += !(s * s == K);

	string file = tmp_file_path("atoms.gpi");
	ofstream os(file);

	os << "set term qt size 1000,1000;" << endl;
	os << "set multiplot layout " << s << "," << s << " rowsfirst scale 1.2;" << endl;
	os << "set isosamples 25,25;" << endl;
	os << "set parametric;" << endl;
	os << "set vrange [-"<< 0 << ":" << radio <<"];" << endl;
	os << "set urange [-pi:pi];" << endl;
	os << "unset key;" << endl;
	os << "set pm3d at b;" << endl;
	os << "unset colorbox;" << endl;

	for(index_t i = 0; i < m; i++)
	{
		os << "splot v * cos(u), v * sin(u), 0 ";
		plot_atoms(os, A.col(i));
		os << ";" << endl;
	}
	os << "unset multiplot;" << endl;
	os << "pause -1;" << endl;

	os.close();

	file = "gnuplot -persist " + file + " &";

	debug(system(file.c_str()));
}


} // namespace gproshan::mdict

