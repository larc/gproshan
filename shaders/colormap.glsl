// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_CB-PuBu.frag

float colormap_red_0(float x)
{
	if (x < 0.7520372909206926)
		return (((9.68615208861418E+02 * x - 1.16097242960380E+03) * x + 1.06173672031378E+02) * x - 1.68616613530379E+02) * x + 2.56073136099945E+02;
	else 
		return -1.20830453148990E+01 * x + 1.44337397593436E+01;
}

float colormap_green_0(float x)
{
	if (x < 0.7485333535031721)
		return (((-4.58537247030064E+02 * x + 5.67323181593790E+02) * x - 2.56714665792882E+02) * x - 1.14205365680507E+02) * x + 2.47073841488433E+02;
	else
		return ((-2.99774273328017E+02 * x + 4.12147041403012E+02) * x - 2.49880079288168E+02) * x + 1.93578601034431E+02;
}

float colormap_blue_0(float x)
{
	if (x < 0.7628468501376879)
		return ((-5.44257972228224E+01 * x + 2.70890554876532E+01) * x - 9.12766750739247E+01) * x + 2.52166182860177E+02;
	else
		return (((4.55621137729287E+04 * x - 1.59960900638524E+05) * x + 2.09530452721547E+05) * x - 1.21704642900945E+05) * x + 2.66644674068694E+04;
}

vec3 colormap_0(float x)
{
	float r = clamp(colormap_red_0(x) / 255.0, .0, .9);
	float g = clamp(colormap_green_0(x) / 255.0, .0, .9);
	float b = clamp(colormap_blue_0(x) / 255.0, .0, .9);
	return vec3(r, g, b);
}


// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_CB-YIOrBr.frag

float colormap_red_1(float x)
{
	return ((((1.30858855846896E+03 * x - 2.84649723684787E+03) * x + 1.76048857883363E+03) * x - 3.99775093706324E+02) * x + 2.69759225316811E+01) * x + 2.54587325383574E+02;
}

float colormap_green_1(float x)
{
	return ((((-8.85605750526301E+02 * x + 2.20590941129997E+03) * x - 1.50123293069936E+03) * x + 2.38490009587258E+01) * x - 6.03460495073813E+01) * x + 2.54768707485247E+02;
}

float colormap_blue_1(float x)
{
	if (x < 0.2363454401493073)
		return (-3.68734834041388E+01 * x - 3.28163398692792E+02) * x + 2.27342862588147E+02;
	else if (x < 0.7571054399013519)
		return ((((1.60988309475108E+04 * x - 4.18782706486673E+04) * x + 4.14508040221340E+04) * x - 1.88926043556059E+04) * x + 3.50108270140290E+03) * x - 5.28541997751406E+01;
	else
		return 1.68513761929930E+01 * x - 1.06424668227935E+01;
}

vec3 colormap_1(float x)
{
	float r = clamp(colormap_red_1(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green_1(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue_1(x) / 255.0, 0.0, 1.0);
	return vec3(r, g, b);
}


// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_Blue-Red_2.frag

float colormap_red_2(float x)
{
	if (x < 0.75) 
		return 1012.0 * x - 389.0;
	else
		return -1.11322769567548E+03 * x + 1.24461193212872E+03;
}

float colormap_green_2(float x)
{
	if (x < 0.5)
		return 1012.0 * x - 129.0;
	else
		return -1012.0 * x + 899.0;
}

float colormap_blue_2(float x)
{
	if (x < 0.25)
		return 1012.0 * x + 131.0;
	else
		return -1012.0 * x + 643.0;
}

vec3 colormap_2(float x)
{
	float r = clamp(colormap_red_2(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green_2(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue_2(x) / 255.0, 0.0, 1.0);
	return vec3(r, g, b);
}

vec3 colormap(uint i, float x)
{
	if(i == 0) return colormap_0(x);
	if(i == 1) return colormap_1(x);
	if(i == 2) return colormap_2(x);
}

