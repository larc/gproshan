#include "scenes/scanner.h"
#include <CImg.h>
#include <thread>

using namespace cimg_library;


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


che * scanner_ptx(const raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam)
{
	std::vector<vertex> ptx;

	return new che(ptx.data(), ptx.size(), nullptr, 0);
}

che * scanner_ptx(const che * mesh, raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam)
{
	std::vector<vertex> vertices;
	std::vector<che::rgb_t> vertices_color;
	vertices.reserve(n_cols * n_rows);
	vertices_color.reserve(n_cols * n_rows);

	glm::vec3 cam_pos = glm_vec3(cam);
	glm::vec3 p, n_v;

	const real_t  r = 1;
	index_t v_idx;
	float distance;

	
	gproshan_log("init");

	const real_t delta_phi = (2 * M_PI) / n_rows;
	const real_t delta_theta = M_PI / n_cols;


	for(real_t phi = 0; phi < 2 * M_PI - 0.5 * delta_phi; phi += delta_phi)
	for(real_t theta = delta_theta; theta < M_PI - 0.5 * delta_theta + delta_theta; theta += delta_theta)
	//for(real_t phi = 0; phi < 2 * M_PI - 0.5 * delta_phi; phi += delta_phi)
	{
		// p is the direction of the ray
		p = glm::vec3( r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta) ) - cam_pos;
		// quering the closest point in the mesh to the hit point and the distance 
		auto [v_idx, distance] = rt->cast_ray_intersect_depth(cam_pos, glm::normalize(p - cam_pos));
		
		if(v_idx == NIL)
		{
			vertices.push_back( {0, 0, 0} );
			vertices_color.push_back({0,0,0});
		}
		else
		{
			n_v = cam_pos + p * distance;
			vertices.push_back( vertex(n_v.x, n_v.y, n_v.z) );
			vertices_color.push_back( mesh->rgb(v_idx) );
		}
	}

	gproshan_log_var(n_cols * n_rows);
	gproshan_log_var(vertices.size());

	gproshan_log_var(n_cols * n_rows == vertices.size());

	che * mesh_ptx = new che(vertices.data(), vertices.size(), nullptr, 0);
	memcpy(&mesh_ptx->rgb(0), vertices_color.data(), vertices_color.size() * sizeof(che::rgb_t));


	CImg<unsigned char> img((unsigned char *) vertices_color.data(), 3, n_cols, n_rows);
	img.permute_axes("zycx");
	img.save((mesh->name().substr(0, mesh->name().size()-4) + ".jpg").c_str());


	std::thread([](CImg<real_t> img) { img.display(); }, img).detach();

	
	return mesh_ptx;
	//return new che(ptx.data(), ptx.size(), nullptr, 0);
}


} // namespace gproshan

