#ifndef SHADER_H
#define SHADER_H

#include <gproshan/viewer/include_opengl.h>
#include <gproshan/geometry/vec.h>
#include <gproshan/geometry/mat.h>

#include <string>


// geometry processing and shape analysis framework
namespace gproshan {


class shader
{
	protected:
		GLuint program = 0;
		bool linked = false;

	public:
		shader() = default;
		~shader();
		operator GLuint () const;
		void load_vertex(const std::string & filename);
		void load_fragment(const std::string & filename);
		void load_geometry(const std::string & filename);
		void enable();
		void disable() const;

		void uniform(const std::string & name, bool value);
		void uniform(const std::string & name, int value);
		void uniform(const std::string & name, unsigned int value);
		void uniform(const std::string & name, float value);
		void uniform(const std::string & name, const vec3 & value);
		void uniform(const std::string & name, const vec4 & value);
		void uniform(const std::string & name, const mat3 & value);
		void uniform(const std::string & name, const mat4 & value);

	protected:
		bool load(GLenum shader_type, const std::string & filename);
		bool read_source(const std::string & filename, std::string & source);
};


} // namespace gproshan

#endif // SHADER_H

