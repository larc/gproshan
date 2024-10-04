#ifndef GRAPH_H
#define GRAPH_H


// geometry processing and shape analysis framework
namespace gproshan {


class union_find
{
	private:
		unsigned * sets = nullptr;
		unsigned n = 0;

	public:
		union_find(unsigned n);
		~union_find();

		unsigned find(const unsigned x);
		bool merge(unsigned x, unsigned y);
};


} // namespace gproshan


#endif // GRAPH_H

