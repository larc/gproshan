// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_CB-PuBu.frag

float colormap_red_1(float x)
{
	if(x < 0.7520372909206926)
		return(((9.68615208861418E+02 * x - 1.16097242960380E+03) * x + 1.06173672031378E+02) * x - 1.68616613530379E+02) * x + 2.56073136099945E+02;
	else
		return -1.20830453148990E+01 * x + 1.44337397593436E+01;
}

float colormap_green_1(float x)
{
	if(x < 0.7485333535031721)
		return(((-4.58537247030064E+02 * x + 5.67323181593790E+02) * x - 2.56714665792882E+02) * x - 1.14205365680507E+02) * x + 2.47073841488433E+02;
	else
		return((-2.99774273328017E+02 * x + 4.12147041403012E+02) * x - 2.49880079288168E+02) * x + 1.93578601034431E+02;
}

float colormap_blue_1(float x)
{
	if(x < 0.7628468501376879)
		return((-5.44257972228224E+01 * x + 2.70890554876532E+01) * x - 9.12766750739247E+01) * x + 2.52166182860177E+02;
	else
		return(((4.55621137729287E+04 * x - 1.59960900638524E+05) * x + 2.09530452721547E+05) * x - 1.21704642900945E+05) * x + 2.66644674068694E+04;
}

vec3 colormap_1(float x)
{
	float r = clamp(colormap_red_1(x) / 255.0, .0, .9);
	float g = clamp(colormap_green_1(x) / 255.0, .0, .9);
	float b = clamp(colormap_blue_1(x) / 255.0, .0, .9);
	return vec3(r, g, b);
}


// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_CB-YIOrBr.frag

float colormap_red_2(float x)
{
	return((((1.30858855846896E+03 * x - 2.84649723684787E+03) * x + 1.76048857883363E+03) * x - 3.99775093706324E+02) * x + 2.69759225316811E+01) * x + 2.54587325383574E+02;
}

float colormap_green_2(float x)
{
	return((((-8.85605750526301E+02 * x + 2.20590941129997E+03) * x - 1.50123293069936E+03) * x + 2.38490009587258E+01) * x - 6.03460495073813E+01) * x + 2.54768707485247E+02;
}

float colormap_blue_2(float x)
{
	if(x < 0.2363454401493073)
		return(-3.68734834041388E+01 * x - 3.28163398692792E+02) * x + 2.27342862588147E+02;
	else if(x < 0.7571054399013519)
		return((((1.60988309475108E+04 * x - 4.18782706486673E+04) * x + 4.14508040221340E+04) * x - 1.88926043556059E+04) * x + 3.50108270140290E+03) * x - 5.28541997751406E+01;
	else
		return 1.68513761929930E+01 * x - 1.06424668227935E+01;
}

vec3 colormap_2(float x)
{
	float r = clamp(colormap_red_2(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green_2(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue_2(x) / 255.0, 0.0, 1.0);
	return vec3(r, g, b);
}


// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_Blue-Red_3.frag

float colormap_red_3(float x)
{
	if(x < 0.75)
		return 1012.0 * x - 389.0;
	else
		return -1.11322769567548E+03 * x + 1.24461193212872E+03;
}

float colormap_green_3(float x)
{
	if(x < 0.5)
		return 1012.0 * x - 129.0;
	else
		return -1012.0 * x + 899.0;
}

float colormap_blue_3(float x)
{
	if(x < 0.25)
		return 1012.0 * x + 131.0;
	else
		return -1012.0 * x + 643.0;
}

vec3 colormap_3(float x)
{
	float r = clamp(colormap_red_3(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green_3(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue_3(x) / 255.0, 0.0, 1.0);
	return vec3(r, g, b);
}


// https://github.com/kbinani/colormap-shaders/blob/master/shaders/glsl/IDL_CB-Set3.frag

float colormap_red_4(float x)
{
	if(x < 0.09082479229584027)
		return 1.24879652173913E+03 * x + 1.41460000000000E+02;
	else if(x < 0.1809653122266933)
		return -7.21339920948626E+02 * x + 3.20397233201581E+02;
	else if(x < 0.2715720097177793)
		return 6.77416996047422E+02 * x + 6.72707509881444E+01;
	else if(x < 0.3619607687891861)
		return -1.36850782608711E+03 * x + 6.22886666666710E+02;
	else if(x < 0.4527609316115322)
		return 1.38118774703557E+03 * x - 3.72395256916997E+02;
	else if(x < 0.5472860687991931)
		return -7.81436521739194E+02 * x + 6.06756521739174E+02;
	else if(x < 0.6360981817705944)
		return 8.06836521739242E+02 * x - 2.62483188405869E+02;
	else if(x < 0.8158623444475089)
		return -3.49616157878512E+02 * x + 4.73134258402717E+02;
	else if(x < 0.9098023786863947)
		return 1.72428853754953E+02 * x + 4.72173913043111E+01;
	else
		return 5.44142292490101E+02 * x - 2.90968379446626E+02;
}

float colormap_green_4(float x)
{
	if(x < 0.08778161310534617)
		return 4.88563478260870E+02 * x + 2.10796666666667E+02;
	else if(x < 0.2697669137324175)
		return -6.96835646006769E+02 * x + 3.14852913968545E+02;
	else if(x < 0.3622079895714037)
		return 5.40799130434797E+02 * x - 1.90200000000068E+01;
	else if(x < 0.4519795462045253)
		return 3.23774703557373E+01 * x + 1.65134387351785E+02;
	else if(x < 0.5466820192751115)
		return 4.43064347826088E+02 * x - 2.04876811594176E+01;
	else if(x < 0.6368889369442862)
		return -1.83472332015826E+02 * x + 3.22028656126484E+02;
	else if(x < 0.728402572416003)
		return 1.27250988142231E+02 * x + 1.24132411067220E+02;
	else if(x < 0.8187333479165154)
		return -9.82116600790428E+02 * x + 9.32198616600708E+02;
	else if(x < 0.9094607880855196)
		return 1.17713438735149E+03 * x - 8.35652173912769E+02;
	else
		return 2.13339920948864E+01 * x + 2.15502964426857E+02;
}

float colormap_blue_4(float x)
{
	if(x < 0.09081516507716858)
		return -2.27937391304345E+02 * x + 1.99486666666666E+02;
	else if(x < 0.1809300436999751)
		return 4.33958498023703E+02 * x + 1.39376482213440E+02;
	else if(x < 0.2720053156712806)
		return -1.14300000000004E+03 * x + 4.24695652173923E+02;
	else if(x < 0.3616296568054424)
		return 1.08175889328072E+03 * x - 1.80450592885399E+02;
	else if(x < 0.4537067088757783)
		return -1.22681999999994E+03 * x + 6.54399999999974E+02;
	else if(x < 0.5472726179445029)
		return 8.30770750988243E+01 * x + 6.00909090909056E+01;
	else if(x < 0.6374811920489858)
		return 1.36487351778676E+03 * x - 6.41401185770872E+02;
	else if(x < 0.7237636846906381)
		return -1.27390769230737E+02 * x + 3.09889230769173E+02;
	else if(x < 0.8178226469606309)
		return -3.01831168831021E+02 * x + 4.36142857142782E+02;
	else if(x < 0.9094505664375214)
		return 8.47622811970801E+01 * x + 1.19977978543158E+02;
	else
		return -9.06117391304296E+02 * x + 1.02113405797096E+03;
}

vec3 colormap_4(float x)
{
	float r = clamp(colormap_red_4(x) / 255.0, 0.0, 1.0);
	float g = clamp(colormap_green_4(x) / 255.0, 0.0, 1.0);
	float b = clamp(colormap_blue_4(x) / 255.0, 0.0, 1.0);
	return vec3(r, g, b);
}


vec3 colormap(uint i, float x)
{
	if(i == 1) return colormap_1(x);
	if(i == 2) return colormap_2(x);
	if(i == 3) return colormap_3(x);
	if(i == 4) return colormap_4(x);
	return vec3(0, 0, 0);
}

