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
		if(rt < 4)
		{
			fscanf(fp, "%*f");
			continue;
		}

		fscanf(fp, "%p %lu %lu %lu", &r.splat, &r.pnv, &r.pnt, &r.ns);
		fscanf(fp, "%lf %lf %lf %lf %lf", &r.tknn, &r.tseg, &r.tvor, &r.tisp, &r.time);

		if(set_pcs.find(str) != set_pcs.end())
			table[rt - 4][str] = r;
	}

	fclose(fp);


	fp = fopen("results/table_splats.tex", "w");
	for(const auto & s: pcs)
	{
		const row & e = table[0][s];
		const row & o = table[1][s];

		fprintf(fp, "\\filename{%20s} & embree & %16lu & %16lu & %16lu & %16lu \\\\        %% %12f %12f %12f %12f %12f\n",
					s.c_str(), e.mnv, e.pnv, e.pnt, e.ns, e.tknn, e.tseg, e.tvor, e.tisp, e.time);
		fprintf(fp, "\\filename{%20s} &  optix & %16lu & %16lu & %16lu & %16lu \\\\\\hline  %% %12f %12f %12f %12f %12f\n\n\n", 
					s.c_str(), o.mnv, o.pnv, o.pnt, o.ns, o.tknn, o.tseg, o.tvor, o.tisp, o.time);

		printf("cp %s_%p results/frametime_%s_embree\n", tmp_file_path("frametime").c_str(), e.view, s.c_str());
		printf("cp %s_%p results/frametime_%s_optix\n", tmp_file_path("frametime").c_str(), o.view, s.c_str());

		printf("cp %s_%p results/histogram_%s_embree\n", tmp_file_path("histogram").c_str(), e.splat, s.c_str());
		printf("cp %s_%p results/histogram_%s_optix\n", tmp_file_path("histogram").c_str(), o.splat, s.c_str());
	}
	fclose(fp);

	return 0;
}

