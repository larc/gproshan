/*
 * Base on the code in https://github.com/larc/dgpdec-course
 *
 * Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow
 * Keenan Crane, Clarisse Weischedel, Max Wardetzky
 * To appear at ACM Transactions on Graphics
 */

#ifndef HEAT_FLOW_H
#define HEAT_FLOW_H

#include "che.h"
#include "include_arma.h"

#include <cholmod.h>

distance_t * heat_flow(che * mesh, const vector<index_t> & sources);

void compute_divergence(che * mesh, const a_mat & u, a_mat & div);

#endif // HEAT_FLOW_H

