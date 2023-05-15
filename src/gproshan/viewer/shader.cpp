#include <gproshan/viewer/shader.h>

#include <gproshan/include.h>

#include <cassert>
#include <fstream>
#include <sstream>


// geometry processing and shape analysis framework
namespace gproshan {


shader::~shader()
{
	glDeleteProgram(program);
}

shader::operator GLuint() const
{
	return program;
}

void shader::load_vertex(const std::string & filename)
{
	load(GL_VERTEX_SHADER, filename);
}

void shader::load_fragment(const std::string & filename)
{
	load(GL_FRAGMENT_SHADER, filename);
}

void shader::load_geometry(const std::string & filename)
{
	load(GL_GEOMETRY_SHADER_EXT, filename);
}

void shader::enable()
{
	if(!linked)
	{
		glLinkProgram(program);
		linked = true;
	}

	glUseProgram(program);
}

void shader::disable() const
{
	glUseProgram(0);
}

void shader::uniform(const std::string & name, bool value)
{
	glProgramUniform1i(program, glGetUniformLocation(program, name.c_str()), value);
}

void shader::uniform(const std::string & name, int value)
{
	glProgramUniform1i(program, glGetUniformLocation(program, name.c_str()), value);
}

void shader::uniform(const std::string & name, unsigned int value)
{
	glProgramUniform1ui(program, glGetUniformLocation(program, name.c_str()), value);
}

void shader::uniform(const std::string & name, float value)
{
	glProgramUniform1f(program, glGetUniformLocation(program, name.c_str()), value);
}

void shader::uniform(const std::string & name, const vec3 & value)
{
	glProgramUniform3fv(program, glGetUniformLocation(program, name.c_str()), 1, &value[0]);
}

void shader::uniform(const std::string & name, const vec4 & value)
{
	glProgramUniform4fv(program, glGetUniformLocation(program, name.c_str()), 1, &value[0]);
}

void shader::uniform(const std::string & name, const mat3 & value)
{
	glProgramUniformMatrix3fv(program, glGetUniformLocation(program, name.c_str()), 1, true, &value[0][0]);
}

void shader::uniform(const std::string & name, const mat4 & value)
{
	glProgramUniformMatrix4fv(program, glGetUniformLocation(program, name.c_str()), 1, true, &value[0][0]);
}

bool shader::load(GLenum shader_type, const std::string & filename)
{
	std::string source;

	if(!read_source(filename, source))
	{
		std::cerr << "Not load shader file: " << filename << std::endl;
		return false;
	}

	if(!program)
		program = glCreateProgram();

	// glDetachShader(program, shader); when/where can we need this?

	GLuint shader = glCreateShader(shader_type);

	const char * source_c_str = source.c_str();
	int size = source.size();

	glShaderSource(shader, 1, &(source_c_str), &size);

	glCompileShader(shader);

	GLint compile_status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);

	if(compile_status == GL_TRUE)
	{
		glAttachShader(program, shader);
		linked = false;
	}
	else
	{
		GLsizei maxLength = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

		if(maxLength > 0)
		{
			GLchar* infoLog = new char[maxLength];
			GLsizei length;

			glGetShaderInfoLog(shader, maxLength, &length, infoLog);

			std::cerr << filename << " GLSL Error: " << infoLog << std::endl;

			delete[] infoLog;
		}

		glDeleteShader(shader);
		return false;
	}

	glDeleteShader(shader);
	return true;
}

bool shader::read_source(const std::string & filename, std::string & source)
{
	std::ifstream is(filename);

	if(!is.is_open())
		return false;

	source = "";

	std::string line, include;
	while(getline(is, line))
	{
		std::stringstream ss(line);

		ss >> include;
		if(include == "#include")
		{
			ss >> include;
			if(read_source(shaders_path(include), include))
				source += include + '\n';
		}
		else
			source += line + '\n';
	}

	is.close();

	return true;
}


} // namespace gproshan

