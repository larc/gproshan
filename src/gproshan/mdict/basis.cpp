#include <gproshan/mdict/basis.h>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


basis::basis(const real_t r, const size_t & d): _radio(r), _dim(d) {}

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
	std::string file = tmp_file_path("basis.gpi");
	std::ofstream os(file);

	os << "set term qt size 1000,1000;" << std::endl;
	os << "set isosamples 50,50;" << std::endl;
	os << "set parametric;" << std::endl;
	os << "set vrange [-"<< 0 << ":" << _radio <<"];" << std::endl;
	os << "set urange [-pi:pi];" << std::endl;
	os << "unset key;" << std::endl;
	os << "set pm3d at b;" << std::endl;
	os << "unset colorbox;" << std::endl;

	plot_basis(os);

	os << "unset multiplot;" << std::endl;
	os << "pause -1;" << std::endl;

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

	std::string file = tmp_file_path("atoms.gpi");
	std::ofstream os(file);

	os << "set term qt size 1000,1000;" << std::endl;
	os << "set multiplot layout " << s << "," << s << " rowsfirst scale 1.2;" << std::endl;
	os << "set isosamples 25,25;" << std::endl;
	os << "set parametric;" << std::endl;
	os << "set vrange [-"<< 0 << ":" << _radio <<"];" << std::endl;
	os << "set urange [-pi:pi];" << std::endl;
	os << "unset key;" << std::endl;
	os << "set pm3d at b;" << std::endl;
	os << "unset colorbox;" << std::endl;

	for(index_t i = 0; i < m; ++i)
	{
		os << "splot v * cos(u), v * sin(u), 0 ";
		plot_atoms(os, A.col(i));
		os << ";" << std::endl;
	}

	os << "unset multiplot;" << std::endl;
	os << "pause -1;" << std::endl;

	os.close();

	file = "gnuplot -persist " + file + " &";

	system(file.c_str());
}

void basis::plot_patch(const a_mat & A, const a_mat & xyz, const index_t p)
{
	a_mat tmp = xyz.t();
	std::string data = tmp_file_path("xyz_" + std::to_string(p) + ".dat");
	tmp.save(data.c_str(), arma::arma_ascii);

	size_t K = A.n_rows;
	size_t m = A.n_cols;
	size_t s = sqrt(m);
	s += !(s * s == K);

	std::string file = tmp_file_path("atoms_patch_"+ std::to_string(p) + ".gpi");
	std::ofstream os(file);

	os << "set term qt size 1000,1000;" << std::endl;
	os << "set multiplot layout " << s << "," << s << " rowsfirst scale 1.2;" << std::endl;
	os << "set isosamples 25,25;" << std::endl;
	os << "set parametric;" << std::endl;
	os << "set vrange [-"<< 0 << ":" << _radio <<"];" << std::endl;
	os << "set urange [-pi:pi];" << std::endl;
	os << "unset key;" << std::endl;
	os << "set pm3d at b;" << std::endl;
	os << "unset colorbox;" << std::endl;
	os << "splot \"xyz_" << std::to_string(p) << ".dat\" u 1:2:3 with points palette pointsize 2 pointtype 7,";

	for(index_t i = 0; i < m; ++i)
	{
		os << " v * cos(u), v * sin(u), 0 ";
		plot_atoms(os, A.col(i));
		os << ";" << std::endl;
	}

	os << "unset multiplot;" << std::endl;
	os << "pause -1;" << std::endl;

	os.close();

	//file = "gnuplot -persist " + file + " &";

	//system(file.c_str());
}


} // namespace gproshan::mdict

