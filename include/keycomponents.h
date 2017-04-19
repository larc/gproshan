#ifndef KEYCOMPONENTS_H
#define KEYCOMPONENTS_H

#include <queue>

#include "off.h"
#include "keypoints.h"
#include "fastmarching.h"

#define NSTD 1.5

class keycomponents
{
	private:
		size_t n_vertices;
		size_t ** components;
		size_t * components_index;
		size_t * marcados;
		size_t n_keypoints;
		size_t * v_keypoints;
		bool * is_keypoint;
		size_t k;
		distance_t radio;
		size_t n_components;
		queue<size_t*> * m_index;

	public:
		keycomponents(keypoints & kps, off & shape, size_t k_ = 45, percent_t percent = 0.10);
		keycomponents(keypoints & kps, off & shape, fastmarching & fm, distance_t radio_ = 0.15, percent_t percent = 0.10);
		~keycomponents();
		void print(ostream & os);
		void print_fm(ostream & os);
		size_t get_ncomponents();

	private:
		void calcular_krings(off & shape);
		void calcular_krings(off & shape, size_t p, size_t * index, size_t kr);
		void calcular_rings_fm(off & shape, fastmarching & fm);
		void agrupar_kpoints(off & shape);
		void agrupar_kpoints(off & shape, size_t p, size_t * index, size_t kr, bool * keypoints);
		void contar_components();

	public:
		void new_keypoints(off & shape, keypoints & kps, ostream & os, percent_t p);
};

void join(size_t * n, size_t a, size_t b);
size_t find(size_t * n, size_t a);

#endif // KEYCOMPONENTS_H
