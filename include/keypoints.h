#ifndef KEYPOINTS_H
#define KEYPOINTS_H

#include <map>

#include "off.h"

using namespace std;

class keypoints
{
	private:
		size_t n_faces;
		areapair_t * areas;

		size_t n_vertices;
		areapair_t * stds;
		areapair_t * means;
		size_t * kps;
		bool * is_kp;

	public:
		keypoints(off & shape);
		~keypoints();
		size_t * get_keypoints();
		bool * get_iskeypoint();
		size_t get_nro_faces(size_t i);
		area_t get_mean_area(size_t i);
		void print(ostream & os, size_t n);

	private:
		void calcular_areas(off & shape);
		void normalizar_areas(off & shape);
		void calcular_stds(off & shape);
		void keypoints_faces(off & shape);
		void keypoints_vertices(off & shape);
};

#endif // KEYPOINTS_H
