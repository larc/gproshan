/*! \file */

#ifndef APP_VIEWER_H
#define APP_VIEWER_H

#include "viewer/viewer.h"
#include "include.h"
#include "che_off.h"
#include "fastmarching.h"
#include "laplacian.h"
#include "che_off.h"
#include "dijkstra.h"
#include "geodesics.h"
#include "fairing_taubin.h"
#include "fairing_spectral.h"
#include "sampling.h"
#include "che_fill_hole.h"
#include "che_poisson.h"
#include "d_mesh_apps.h"
#include "decimation.h"
#include "denoising.h"
#include "d_basis_dct.h"
#include "d_basis_cosine.h"

void viewer_process_fairing_taubin();
void viewer_process_fairing_spectral();

void viewer_process_fastmarching();
void viewer_process_geodesics();
void viewer_process_farthest_point_sampling();
void viewer_process_farthest_point_sampling_radio();
void viewer_process_fastmarching_cpu();
void viewer_process_fastmarching_gpu();
void viewer_sort_by_rings();
void viewer_process_voronoi();

void viewer_process_denoising();
void viewer_process_super_resolution();
void viewer_process_inpaiting();
void viewer_process_iterative_inpaiting();

void viewer_process_gps();
void viewer_process_hks();
void viewer_process_wks();

void viewer_process_poisson(const index_t & lk);
void viewer_process_poisson_laplacian_1();
void viewer_process_poisson_laplacian_2();
void viewer_process_poisson_laplacian_3();

void viewer_process_thresold();
void viewer_process_noise();
void viewer_process_multiplicate_vertices();
void viewer_process_fill_holes();
void viewer_process_delete_vertices();
void viewer_process_fill_holes_test();
void viewer_process_delete_non_manifold_vertices();
void viewer_process_fill_holes_biharmonic_splines();
void viewer_process_gaussian_curvature();
void viewer_process_edge_collapse();
void viewer_select_multiple();

int viewer_main(int nargs, char ** args);

#endif //APP_VIEWER_H

