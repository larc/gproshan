#ifndef SHADER_H
#define SHADER_H

#include <string>
#include <map>

#include <gproshan/viewer/include_opengl.h>

// geometry processing and shape analysis framework
namespace gproshan {


class shader
{
	protected:
		GLuint program = 0;
		std::map<std::string, GLint> uniform;
		bool linked = false;

	public:
		shader() = default;
		~shader();
		const GLint & operator () (const std::string & name);
		operator GLuint () const;
		void load_vertex(const std::string & filename);
		void load_fragment(const std::string & filename);
		void load_geometry(const std::string & filename);
		void enable();
		void disable() const;

	protected:
		bool load(GLenum shader_type, const std::string & filename);
		bool read_source(const std::string & filename, std::string & source);
};


} // namespace gproshan

#endif // SHADER_H

