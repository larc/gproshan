#include <gproshan/mdict/basis.h>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


basis::basis(const real_t & r, const size_t & d): _radio(r), _dim(d) {}

real_t & basis::radio()
{
	return _radio;
}

const size_t & basis::dim() const
{
	return _dim;
}

void basis::plot_basis()
{
	string file = tmp_file_path("basis.gpi");
	ofstream os(file);

	os << "set term qt size 1000,1000;" << endl;
	os << "set isosamples 50,50;" << endl;
	os << "set parametric;" << endl;
	os << "set vrange [-"<< 0 << ":" << _radio <<"];" << endl;
	os << "set urange [-pi:pi];" << endl;
	os << "unset key;" << endl;
	os << "set pm3d at b;" << endl;
	os << "unset colorbox;" << endl;

	plot_basis(os);

	os << "unset multiplot;" << endl;
	os << "pause -1;" << endl;

	os.close();

	file = "gnuplot -persist " + file + " &";

	system(file.c_str());
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
	os << "set vrange [-"<< 0 << ":" << _radio <<"];" << endl;
	os << "set urange [-pi:pi];" << endl;
	os << "unset key;" << endl;
	os << "set pm3d at b;" << endl;
	os << "unset colorbox;" << endl;

	for(index_t i = 0; i < m; ++i)
	{
		os << "splot v * cos(u), v * sin(u), 0 ";
		plot_atoms(os, A.col(i));
		os << ";" << endl;
	}

	os << "unset multiplot;" << endl;
	os << "pause -1;" << endl;

	os.close();

	file = "gnuplot -persist " + file + " &";

	system(file.c_str());
}

void basis::plot_patch(const a_mat & A, const a_mat & xyz, const index_t & p)
{
	a_mat tmp = xyz.t();
	string data = tmp_file_path("xyz_" + to_string(p) + ".dat");
	tmp.save(data.c_str(), arma::arma_ascii);

	size_t K = A.n_rows;
	size_t m = A.n_cols;
	size_t s = sqrt(m);
	s += !(s * s == K);

	string file = tmp_file_path("atoms_patch_"+ to_string(p) + ".gpi");
	ofstream os(file);

	os << "set term qt size 1000,1000;" << endl;
	os << "set multiplot layout " << s << "," << s << " rowsfirst scale 1.2;" << endl;
	os << "set isosamples 25,25;" << endl;
	os << "set parametric;" << endl;
	os << "set vrange [-"<< 0 << ":" << _radio <<"];" << endl;
	os << "set urange [-pi:pi];" << endl;
	os << "unset key;" << endl;
	os << "set pm3d at b;" << endl;
	os << "unset colorbox;" << endl;
	os << "splot \"xyz_" << to_string(p) << ".dat\" u 1:2:3 with points palette pointsize 2 pointtype 7,";

	for(index_t i = 0; i < m; ++i)
	{
		os << " v * cos(u), v * sin(u), 0 ";
		plot_atoms(os, A.col(i));
		os << ";" << endl;
	}

	os << "unset multiplot;" << endl;
	os << "pause -1;" << endl;

	os.close();

	//file = "gnuplot -persist " + file + " &";

	//system(file.c_str());
}


} // namespace gproshan::mdict
