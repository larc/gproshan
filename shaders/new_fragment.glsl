#version 460 core

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
	const float shininess = 4.;
	vec3 R = 2.*dot(L,N)*N - L;
	return pow( max( 0., dot( R, E )), shininess );
}

float fresnel( vec3 N, vec3 E )
{
	const float sharpness = 10.;
	float NE = max( 0., dot( N, E ));
	return pow( sqrt( 1. - NE*NE ), sharpness );
}

//https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_CB-PuBu.frag
float colormap_red(float x)
{
	if (x < 0.7520372909206926)
		return (((9.68615208861418E+02 * x - 1.16097242960380E+03) * x + 1.06173672031378E+02) * x - 1.68616613530379E+02) * x + 2.56073136099945E+02;
	else 
		return -1.20830453148990E+01 * x + 1.44337397593436E+01;
}

float colormap_green(float x)
{
	if (x < 0.7485333535031721)
		return (((-4.58537247030064E+02 * x + 5.67323181593790E+02) * x - 2.56714665792882E+02) * x - 1.14205365680507E+02) * x + 2.47073841488433E+02;
	else
		return ((-2.99774273328017E+02 * x + 4.12147041403012E+02) * x - 2.49880079288168E+02) * x + 1.93578601034431E+02;
}

float colormap_blue(float x)
{
	if (x < 0.7628468501376879)
		return ((-5.44257972228224E+01 * x + 2.70890554876532E+01) * x - 9.12766750739247E+01) * x + 2.52166182860177E+02;
	else
		return (((4.55621137729287E+04 * x - 1.59960900638524E+05) * x + 2.09530452721547E+05) * x - 1.21704642900945E+05) * x + 2.66644674068694E+04;
}

vec3 colormap(float x)
{
	float r = clamp(colormap_red(x) / 255.0, .0, .9);
	float g = clamp(colormap_green(x) / 255.0, .0, .9);
	float b = clamp(colormap_blue(x) / 255.0, .0, .9);
	return vec3(r, g, b);
}
//winter
vec3 colormap2(float x) {
    return vec3(0.0, clamp(x, 0.0, 1.0), clamp(-0.5 * x + 1.0, 0.0, 1.0));
}

void main()
{
	// color
	/*
	float d = 1. - color;
	float r = (1. - d*d) * .8;
	float g = (1. - (2. * (d - .5)) * (2. * (d - .5))) * .7;
	float b = (1. - (1. - d) * (1. - d));
	vec3 vcolor = vec3(r, g, b);
	*/
	vec3 vcolor = colormap(color);


	// lines
	if(lines)
	{
		float h = color;
		h = h * 40.;
		h = h - floor( h );
		h = (1. / (1. + exp(-100.*(h - .55)))) + (1. / (1. + exp(-100.*(-h + .45))));
		h = 1. - h;
		vcolor.xyz = vec3(0, 0, 0) + (1. - h) * vcolor.xyz;
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

	FragColor.rgb = diffuse(N,L) * vcolor + .1 * specular(N,L,E) * one + .5 * fresnel(N,E) * one;
	FragColor.a = 1.;
}

