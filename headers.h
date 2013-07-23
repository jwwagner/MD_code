//headers.h
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include <time.h>

#include "prototypes.h"

//PBC
#define NO_PERIODICAL_BD	0
#define PERIODICAL_BD		1

// renew of neighbour list
#define NO_RENEW			0
#define RENEW				1

//control
#define NOT_YET_CALCULATION_MODE_FOR_SYSTEM	    -1
#define GOING_ON_CALCULATION_MODE_FOR_SYSTEM	0
#define FINISHED_CALCULATION_MODE_FOR_SYSTEM	1

#define NO_NEIGHBOUR_IN_PARTICLE	0
#define NEIGHBOUR_IN_PARTICLE		1

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

typedef struct {
    int atoms;
    int copies;
} Topology;

typedef struct {
    int restart_flag, ensemble_flag;
    int n_steps;
    int molecule_types;
    int num_particles;
    int total_cells, cellx, celly, cellz;
    int counter;
    int integrator;
    int parallel_flag;

    double nose_velocity, pressure_s, pressure_v;
    double kb, nose_hoover, nose_Q, npt_pressure, npt_mass, npt_nose;
    double temp, pressure, actual_pressure, delta_t;
    double cutoff, cutoff2;
    double box_x, box_y, box_z;
    double angle1, angle2, angle3;
    double energy, volume;
    double neighbor_factor;
    double kinetic_energy, potential_energy;

    Topology* types;
} Control;


typedef struct {
    char type[3];
    double sigma;
    double epsilon;
    double x_rel;
    double y_rel;
    double z_rel;
    double charge;
    double mass;

    double cutoff_potential, cutoff_derivative;

} Prototypes;

typedef struct {
    char type[3];
    double sigma;
    double epsilon;
    double x, y ,z;
    double oldx, oldy, oldz;
    double vx, vy, vz;
    double fx, fy, fz;
    double oldfx, oldfy, oldfz;//velocity verlet only
    double potential;
    double charge;
    double mass;
    int cell_id;
    int style;

    double cutoff_potential, cutoff_derivative;
} Particles;

typedef struct {
    int* cell;
    int length;
} NeighborList;
