#ifndef APP_VIEWER_H
#define APP_VIEWER_H

#include "viewer/viewer.h"
#include "include.h"
#include "che_off.h"
#include "laplacian.h"
#include "che_off.h"
#include "dijkstra.h"
#include "geodesics.h"
#include "fairing_taubin.h"
#include "fairing_spectral.h"
#include "sampling.h"
#include "che_fill_hole.h"
#include "che_poisson.h"
#include "decimation.h"
#include "mdict/denoising.h"
#include "mdict/super_resolution.h"
#include "mdict/inpainting.h"
#include "mdict/d_basis_dct.h"
#include "mdict/patch.h"
#include "key_points.h"
#include "key_components.h"

void viewer_process_fairing_taubin();
void viewer_process_fairing_spectral();

void viewer_process_fastmarching();
void viewer_process_geodesics_fm();
void viewer_process_geodesics_ptp_gpu();
void viewer_process_geodesics_ptp_cpu();
void viewer_process_geodesics_heat_flow();
void viewer_process_geodesics_heat_flow_gpu();
void viewer_process_farthest_point_sampling();
void viewer_process_farthest_point_sampling_radio();
void viewer_compute_toplesets();
void viewer_process_voronoi();

void viewer_process_mdict_patch();
void viewer_process_denoising();
void viewer_process_super_resolution();
void viewer_process_inpaiting();
void viewer_process_iterative_inpaiting();

void viewer_process_functional_maps();
void viewer_process_gps();
void viewer_process_hks();
void viewer_process_wks();
void viewer_process_key_points();
void viewer_process_key_components();

void viewer_process_poisson(const index_t & lk);
void viewer_process_poisson_laplacian_1();
void viewer_process_poisson_laplacian_2();
void viewer_process_poisson_laplacian_3();

void viewer_process_thresold();
void viewer_process_noise();
void viewer_process_black_noise();
void viewer_process_multiplicate_vertices();
void viewer_process_fill_holes();
void viewer_process_delete_vertices();
void viewer_process_fill_holes_test();
void viewer_process_delete_non_manifold_vertices();
void viewer_process_fill_holes_biharmonic_splines();
void viewer_process_gaussian_curvature();
void viewer_process_edge_collapse();
void viewer_select_multiple();

int viewer_main(int nargs, const char ** args);

#endif //APP_VIEWER_H

