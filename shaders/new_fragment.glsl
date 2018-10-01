#version 450 core

uniform vec3 eye;
uniform vec3 light;
uniform bool is_flat;
uniform bool lines;

in vec3 position;
in vec3 normal;
in float color;

layout(location = 0) out vec4 FragColor;

float diffuse( vec3 N, vec3 L )
{
	return max( 0., dot( N, L ));
}

float specular( vec3 N, vec3 L, vec3 E )
{
	const float shininess = 8.;
	vec3 R = 2.*dot(L,N)*N - L;
	return pow( max( 0., dot( R, E )), shininess );
}

float fresnel( vec3 N, vec3 E )
{
	const float sharpness = 10.;
	float NE = max( 0., dot( N, E ));
	return pow( sqrt( 1. - NE*NE ), sharpness );
}

void main()
{
	// color
	float d = 1. - color;
	float r = (1. - d*d) * .8;
	float g = (1. - (2. * (d - .5)) * (2. * (d - .5))) * .7;
	float b = (1. - (1. - d) * (1. - d));
	vec3 vcolor = vec3(r, g, b);

	// lines
	if(lines)
	{
		float h = color;
		h = h * 30.;
		h = h - floor( h );
		h = (1. / (1. + exp(-100.*(h - .55)))) + (1. / (1. + exp(-100.*(-h + .45))));
		h = 1. - h;
		vcolor.xyz = vec3(h, h, h) + (1. - h) * vcolor.xyz;
	}

 	vec3 N;

 	if(is_flat)
	{
		vec3 X = dFdx(position);
		vec3 Y = dFdy(position);
		vec3 normal_flat = cross(X,Y);
		N = normalize( normal_flat );
	}
	else
	{
		N = normalize( normal );
	}
	
	vec3 L = normalize( light - position );
	vec3 E = normalize( eye - position );
	vec3 R = 2.*dot(L,N)*N - L;
	vec3 one = vec3( 1., 1., 1. );

	FragColor.rgb = diffuse(N,L) * vcolor + .5 * specular(N,L,E) * one + .5 * fresnel(N,E) * one;
	FragColor.a = 1.;
}

