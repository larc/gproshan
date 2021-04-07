#version 410 core

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec3 vs_position[];
in vec3 vs_normal[];
in vec3 vs_mesh_color[];
in float vs_color[];

out vec3 gs_position;
out vec3 gs_normal;
out vec3 gs_mesh_color;
out float gs_color;

noperspective out vec3 edge_dist;

void main()
{
	float a = length(vs_position[1] - vs_position[2]);
	float b = length(vs_position[2] - vs_position[0]);
	float c = length(vs_position[1] - vs_position[0]);

	float alpha = acos((b * b + c * c - a * a) / (2 * b * c));
	float beta = acos((a * a + c * c - b * b) / (2 * a * c));

	float ha = abs(c * sin(beta));
	float hb = abs(c * sin(alpha));
	float hc = abs(b * sin(alpha));


	gs_position = vs_position[0];
	gs_normal = vs_normal[0];
	gs_mesh_color = vs_mesh_color[0];
	gs_color = vs_color[0];
	edge_dist = vec3(ha, 0, 0);
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();

	gs_position = vs_position[1];
	gs_normal = vs_normal[1];
	gs_mesh_color = vs_mesh_color[1];
	gs_color = vs_color[1];
	edge_dist = vec3(0, hb, 0);
	gl_Position = gl_in[1].gl_Position;
	EmitVertex();

	gs_position = vs_position[2];
	gs_normal = vs_normal[2];
	gs_mesh_color = vs_mesh_color[2];
	gs_color = vs_color[2];
	edge_dist = vec3(0, 0, hc);
	gl_Position = gl_in[2].gl_Position;
	EmitVertex();

	EndPrimitive();
}

