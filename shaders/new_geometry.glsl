#version 450 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;
uniform mat4  projection;  //from application program.
out vec3      normal;      //to fragment shader.
 
void main (void) {
  vec3        vector1;
  vec3        vector2;
 
  gl_Position = gl_in[0].gl_Position;
  vector1 = gl_in[1].gl_Position.xyz - gl_Position.xyz;
  vector2 = gl_in[2].gl_Position.xyz - gl_Position.xyz;
  normal = normalize (cross (vector1, vector2));
  gl_Position = projection * gl_Position;
  EmitVertex ();
 
  gl_Position = projection * gl_in[1].gl_Position;
  EmitVertex ();
 
  gl_Position = projection * gl_in[2].gl_Position;
  EmitVertex ();
  EndPrimitive ();
}
