#version 410 core

in vec2 UV;

out vec4 color;

uniform sampler2D render_tex;

void main()
{
	color = texture(render_tex, UV).rgba;
}

