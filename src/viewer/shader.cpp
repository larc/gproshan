#include "viewer/shader.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


shader::~shader()
{
	if(program) glDeleteProgram(program);
}

const GLint & shader::operator () (const string & name)
{
	if(uniform.find(name) != uniform.end())
		uniform[name] = glGetUniformLocation(program, name.c_str());
	
	return uniform[name];
}

shader::operator GLuint() const
{
	return program;
}

void shader::load_vertex(const char * filename)
{
	load(GL_VERTEX_SHADER, filename);
}

void shader::load_fragment(const char * filename)
{
	load(GL_FRAGMENT_SHADER, filename);
}

void shader::load_geometry(const char * filename)
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

bool shader::load(GLenum shader_type, const char * filename)
{
	string source;

	if(!read_source(filename, source))
	{
		cerr << "Not load shader file: " << filename << endl;
		return false;
	}

	if(program == 0)
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

			cerr << filename << " GLSL Error: " << infoLog << endl;

			delete[] infoLog;
		}
		
		glDeleteShader(shader);
		return false;
	}

	glDeleteShader(shader);
	return true;
}

bool shader::read_source(const char * filename, std::string & source)
{
	ifstream is(filename);

	if(!is.is_open())
		return false;

	source = "";
	
	string line, include;
	while(getline(is, line))
	{
		stringstream ss(line);
		
		ss >> include;
		if(include == "#include")
		{
			ss >> include;
			if(read_source(include.c_str(), include))
				source += include + '\n';
		}
		else
			source += line + '\n';
	}
	
	is.close();

	return true;
}


} // namespace gproshan

