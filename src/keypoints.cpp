#include "keypoints.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>

keypoints::keypoints(off & shape)
{
	n_faces = shape.get_nfaces();
	areas = new areapair_t[n_faces];

	n_vertices = shape.get_nvertices();
	stds = new areapair_t[n_vertices];
	means = new areapair_t[n_vertices];
	kps = new size_t[n_vertices];
	is_kp = new bool[n_vertices];

	//keypoints_vertices(shape);
	keypoints_faces(shape);
}

keypoints::~keypoints()
{
	delete [] areas;
	delete [] stds;
	delete [] means;
	delete [] kps;
	delete [] is_kp;
}

size_t * keypoints::get_keypoints()
{
	return kps;
}

bool * keypoints::get_iskeypoint()
{
	return is_kp;
}

size_t keypoints::get_nro_faces(size_t i)
{
	return means[i].second;
}

void keypoints::print(ostream & os, size_t n)
{
	if(n >= n_vertices) return;
	for(size_t i = 0; i < n; i++)
		os<<kps[i]<<" ";
	os<<endl;
}

area_t keypoints::get_mean_area(size_t i)
{
	return means[i].first;
}

void keypoints::calcular_areas(off & shape)
{
	vertex a, b, c;

	for(size_t i = 0; i < n_vertices; i++)
	{
		means[i].first = 0;
		means[i].second = 0;
	}

	for(size_t i = 0; i < n_faces; i++)
	{
		a = shape(shape[i][0]);
		b = shape(shape[i][1]);
		c = shape(shape[i][2]);

		areas[i].first = *((a - b) * (a - c)) / 2;
		areas[i].second = i;

		means[shape[i][0]].second++;
		means[shape[i][1]].second++;
		means[shape[i][2]].second++;
	}
}

void keypoints::normalizar_areas(off & shape)
{
	vertex a, b, c;
	area_t R;

	for(size_t i = 0; i < n_faces; i++)
	{
		a = shape(shape[i][0]);
		b = shape(shape[i][1]);
		c = shape(shape[i][2]);

		R = *(a-b)**(b-c)**(c-a)/(4*areas[i].first);
		areas[i].first /= M_PI*R*R;

		means[shape[i][0]].first += areas[i].first;
		means[shape[i][1]].first += areas[i].first;
		means[shape[i][2]].first += areas[i].first;
	}
}

void keypoints::calcular_stds(off & shape)
{
	calcular_areas(shape);
	normalizar_areas(shape);

	for(size_t i = 0; i < n_vertices; i++)
	{
		means[i].first /= means[i].second;
		stds[i].first = 0;
		stds[i].second = i;
	}

	area_t tmp;
	for(size_t i = 0; i < n_faces; i++)
	{
		tmp = means[shape[i][0]].first - areas[i].first;
		stds[shape[i][0]].first += tmp*tmp;
		tmp = means[shape[i][1]].first - areas[i].first;
		stds[shape[i][1]].first += tmp*tmp;
		tmp = means[shape[i][2]].first - areas[i].first;
		stds[shape[i][2]].first += tmp*tmp;
	}

	for(size_t i = 0; i < n_vertices; i++)
		stds[i].first = sqrt(stds[i].first)/(means[i].second - 1);
}

void keypoints::keypoints_faces(off & shape)
{
	calcular_areas(shape);
	//normalizar_areas(shape);
	//calcular_stds(shape);
	sort(areas, areas + n_faces);

	memset(is_kp, 0, sizeof(bool)*n_vertices);

	size_t p = 0;
	for(size_t i = 0; i < n_faces; i++)
	{
		if(!is_kp[shape[areas[i].second][0]])
		{
			kps[p++] = shape[areas[i].second][0];
			is_kp[shape[areas[i].second][0]] = 1;
		}
		if(!is_kp[shape[areas[i].second][1]])
		{
			kps[p++] = shape[areas[i].second][1];
			is_kp[shape[areas[i].second][1]] = 1;
		}
		if(!is_kp[shape[areas[i].second][2]])
		{
			kps[p++] = shape[areas[i].second][2];
			is_kp[shape[areas[i].second][2]] = 1;
		}
	}
}

void keypoints::keypoints_vertices(off & shape)
{
	calcular_stds(shape);
	sort(stds, stds + n_vertices);

	for(size_t i = 0; i < n_vertices; i++)
		kps[i] = stds[n_vertices - i - 1].second;
}

