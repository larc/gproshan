#include "shader.h"

#include <fstream>
#include <iostream>

using namespace std;

shader::shader(): vertexshader(0),
					fragmentshader(0),
					geometryshader(0),
					program(0),
					linked(false)
{
}

shader::~shader()
{
	if(program) glDeleteProgram(program);

	if(vertexshader) glDeleteShader(vertexshader);
	if(fragmentshader) glDeleteShader(fragmentshader);
	if(geometryshader) glDeleteShader(geometryshader);
}

void shader::loadVertex(const char * filename)
{
	load(GL_VERTEX_SHADER, filename, vertexshader);
}

void shader::loadFragment(const char * filename)
{
	load(GL_FRAGMENT_SHADER, filename, fragmentshader);
}

void shader::loadGeometry(const char * filename)
{
	#ifdef GL_GEOMETRY_SHADER_EXT
		load(GL_GEOMETRY_SHADER_EXT, filename, geometryshader);
	#else
		cerr << "Error: geometry shaders not supported!" << endl;
	#endif
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

shader::operator GLuint() const
{
	return program;
}

void shader::load(GLenum shaderType, const char * filename, GLuint & _shader)
{
	string source;
	
	if(!readSource(filename, source))
	{
		return;
	}

	if(program == 0)
	{
		program = glCreateProgram();
	}

	if(_shader != 0)
	{
		glDetachShader(program, _shader);
	}

	_shader = glCreateShader(shaderType);
	const char * source_c_str = source.c_str();
	int size = source.size();
	glShaderSource(_shader, 1, &(source_c_str), &size);

	glCompileShader(_shader);
	GLint compileStatus;
	glGetShaderiv(_shader, GL_COMPILE_STATUS, &compileStatus);

	if(compileStatus == GL_TRUE)
	{
		glAttachShader(program, _shader);
		linked = false;
	}
	else
	{
		GLsizei maxLength = 0;
		glGetShaderiv(_shader, GL_INFO_LOG_LENGTH, &maxLength);

		if(maxLength > 0)
		{
			GLchar* infoLog = new char[ maxLength ];
			GLsizei length;

			glGetShaderInfoLog(_shader, maxLength, &length, infoLog);

			cerr << filename << " GLSL Error: " << infoLog << endl;

			delete[] infoLog;
		}
	}
}

bool shader::readSource(const char * filename, std::string & source)
{
	source = "";

	ifstream in(filename);
	if(!in.is_open())
	{
		cerr << "Error: could not open shader file ";
		cerr << filename;
		cerr << " for input!" << endl;
		return false;
	}

	string line;
	while(getline(in, line))
	{
		source += line + '\n';
	}

	return true;
}

