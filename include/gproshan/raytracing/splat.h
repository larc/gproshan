#ifndef SPLAT_H
#define SPLAT_H


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat : public raytracing
{
	private:
		che * pc = nullptr;
		std::vector<int> morton_codes;
		std::vector<index_t> start_splats;

	public:
		splat(const che * mesh);
		virtual ~splat();

		virtual void render(vec4 * img, const render_params & params, const bool & flat);

	protected:
		void init_splats_mesh(const che * mesh, std::vector<index_t> & vertices);
		void build_splats_ch(const che * mesh, const std::vector<index_t> & vertices);
};


} // namespace gproshan

#endif // SPLAT_H

