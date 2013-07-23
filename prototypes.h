//prototypes.h

void    FreeMemory(void);   //closing.c

double    cpu_time(void);     //cpu_time.c

//compute forces
void    calc_forces(void);              //forces.c
void    compute_IMF(void);              //forces.c
void    compute_IMF_with_neighbors(void);//force.c
void    compute_IMF_same_cell(int);     //force.c
void    compute_IMF_in_cells(int, int); //force.c
void    total_energy(void);             //force.c
void    zero_forces_and_potential(void);// forces.c

//neighborlist
void    build_neighbor_list(void); //forces.c
void    initialize_neighbor_list(void); //forces.c

//initializiation and file reads
void    create_particles(void); //initialize.c
void    read_controls(void); //initialize.c
void    initialize_system(void); //initialize.c
void    initialize_internally(void);//intitialize.c
void    initialize_from_file(void);//initialize.c
void    initialize_coordinates_from_file(void); //initialize.c
void    initialize_topology_from_file(void); //initialize.c

//integration
void    time_step(void);          //integrate.c
void    do_verlet(void);          //integrate.c
void    do_velocity_verlet(void); //integrate.c
void    do_NVT_velocity_verlet(void); //integrate.c
void    do_bendersen_thermostat(void); //integrate.c
void    do_nose_hoover(void);     //integrate.c
void    do_NPT_nose_hoover(void); //integrate.c
void    do_NPT_bendersen(void);   //integrate.c
void    do_leapfrog_verlet(void); //integrate.c
void    wrap_PBC(void);           //integrate.c
void do_NPT_AT_method(void);      //integrate.c
void do_NPT_AT_method_Temp(void); //integrate.c
void do_NPT_AT_method_Pressure(void); //integrate.c
void do_NPT_berendsen_method_Temp(void);
void do_NPT_berendsen_method_Pressure(void);

//free memory
void    ErrExit (int); //err.c

//outputs
void interval_output(void); //output.c
void intialize_VMD_output(void); //output.c
void VMD_output(void);  //output.c
void initialize_energy_output(void); //output.c
void energy_output(void); //output.c
