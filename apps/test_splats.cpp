#include <gproshan/include.h>

#include <vector>
#include <unordered_set>
#include <unordered_map>


using namespace gproshan;


struct row
{
	void * view;
	size_t mnv;
	void * splat;
	size_t pnv;
	size_t pnt;
	size_t ns;
	double tknn;
	double tseg;
	double tvor;
	double tisp;
	double time;
};

int main()
{
	std::vector<std::string> pcs;
	std::unordered_set<std::string> set_pcs;

	char str[32];
	while(scanf("%s", str) != EOF)
	{
		pcs.emplace_back(str);
		set_pcs.insert(str);
	}


	std::unordered_map<std::string, row> table[2];

	FILE * fp = fopen(tmp_file_path("rt_build_times").c_str(), "r");

	row r;
	index_t rt;
	while(fscanf(fp, "%s", str) != EOF)
	{
		fscanf(fp, "%p %s %lu %*d %u", &r.view, str, &r.mnv, &rt);
		if(rt < 3)
		{
			fscanf(fp, "%*f");
			continue;
		}

		fscanf(fp, "%p %lu %lu %lu", &r.splat, &r.pnv, &r.pnt, &r.ns);
		fscanf(fp, "%lf %lf %lf %lf %lf", &r.tknn, &r.tseg, &r.tvor, &r.tisp, &r.time);

		table[rt - 4][str] = r;
	}

	fclose(fp);

	for(auto & s: pcs)
	{
		const row & e = table[0][s];
		const row & o = table[1][s];

		printf("%20s & embree & %16lu & %16lu & %16lu & %16lu \\\\\n", s.c_str(), e.mnv, e.pnv, e.pnt, e.ns);
		printf("%20s & optix & %16lu & %16lu & %16lu & %16lu \\\\\\hline\n", s.c_str(), o.mnv, o.pnv, o.pnt, o.ns);
	}

	return 0;
}

