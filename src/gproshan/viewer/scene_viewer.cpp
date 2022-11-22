#include <gproshan/viewer/scene_viewer.h>


// geometry processing and shape analysis framework
namespace gproshan {


scene_viewer::scene_viewer(scene * p_sc): che_viewer(p_sc), sc(p_sc)
{
}

void scene_viewer::draw(shader & program)
{
	glProgramUniformMatrix4fv(program, program("model_mat"), 1, true, &model_mat[0][0]);
	glProgramUniform1ui(program, program("idx_colormap"), idx_colormap);
	glProgramUniform1i(program, program("render_flat"), render_flat);
	glProgramUniform1i(program, program("render_lines"), render_lines);
	glProgramUniform1i(program, program("render_wireframe"), render_triangles);

	glPolygonMode(GL_FRONT_AND_BACK, render_wireframe ? GL_LINE : GL_FILL);

	program.enable();

	glBindVertexArray(vao);
	glDrawArrays(GL_TRIANGLES, 0, mesh->n_vertices);
	glBindVertexArray(0);

	program.disable();

}


} // namespace gproshan

