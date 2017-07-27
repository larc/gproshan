#ifndef VIEWER_SHADER_H
#define VIEWER_SHADER_H

#include <GLES3/gl3.h>
#include <GL/glut.h>
#include <string>

class shader
{
	protected:
		GLuint vertexshader;
		GLuint fragmentshader;
		GLuint geometryshader;
		GLuint program;
		bool linked;

	public:
		shader();
		~shader();
		void loadVertex( const char* filename );
		void loadFragment( const char* filename );
		void loadGeometry( const char* filename );
		void enable();
		void disable() const;
		operator GLuint() const;

	protected:
		void load(GLenum shaderType, const char * filename, GLuint & shader);
		bool readSource( const char * filename, std::string & source);
};

#endif
