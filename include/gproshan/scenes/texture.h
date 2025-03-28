#ifndef TEXTURE_H
#define TEXTURE_H

#include <gproshan/geometry/vec.h>

#include <string>


// geometry processing and shape analysis framework
namespace gproshan {


struct texture
{
	unsigned char * data = nullptr;
	unsigned int width = 0;
	unsigned int height = 0;
	unsigned int spectrum = 0;

	__host_device__
	texture() = default;
	texture(const std::string & file);

	__host_device__
	operator bool () const
	{
		return data != nullptr;
	}

	__host_device__
	vec4 operator () (const vec2 & coord) const
	{
		const int i = (width + int(coord.x() * (width - 1))) % width;
		const int j = (height + int(coord.y() * (height - 1))) % height;
		const int k = (j * width + i) * spectrum;

		unsigned char * tex = data + k;

		vec4 v;
		for(unsigned int i = 0; i < spectrum; ++i)
			v[i] = float(tex[i]) / 255;

		if(spectrum == 1)
			v[3] = v[2] = v[1] = v[0];

		if(spectrum == 3)
			v[3] = 1;

		return v;
	}
};


} // namespace gproshan


#endif // TEXTURE_H

