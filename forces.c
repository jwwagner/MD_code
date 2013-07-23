//forces.c

void calc_forces(void)
{
//	initialization_for_computing_forces();
//	compute_bond_information();

//	if(nebrNow == RENEW){

//		BuildNebrList();

//		if(flagEntaglementForce == ENTANGLEMENT_FORCE){
//			BuildNebrList_for_entanglement();
//		}

//		saving_position_for_neighbouring_list();
//		nebrNow = NO_RENEW;
//	}

//	compute_intra_particle_forces();

	//compute_IMF();
	printf("calculating forces %d\n", System.counter);
	zero_forces_and_potential();
	//printf("IMF\n");
	compute_IMF_with_neighbors();
	//printf("total energy\n");
    total_energy();
//	nebrNow = check_neighbouring_list();
}

void zero_forces_and_potential(void)
{
    int i;

    for(i = 0; i < System.num_particles; i++)
    {
        Particle[i].fx = 0;
        Particle[i].fy = 0;
        Particle[i].fz = 0;
        Particle[i].potential = 0;
    }

    System.actual_pressure = 0.0;
}

void compute_IMF(void)
{
    int num = System.num_particles;
    int i, j;
    double r2, rx, ry, rz, rij;
    double force, potential;

    //normally need neighbor list _HERE_ we use the full particle list

    //loop through particle list
    for(i = 0; i < num; i++)
    {
        for(j = i + 1; j < num; j++)
        {
            rx = Particle[i].x - Particle[j].x;
            ry = Particle[i].y - Particle[j].y;
            rz = Particle[i].z - Particle[j].z;

            //printf("%lf\n", (rx*rx + ry*ry + rz*rz));
            //check each for nearest image
            if(abs(rx) > (System.box_x / 2.0) )
            {
                if(rx > 0) rx -= System.box_x / 2.0;
                else rx += System.box_x / 2.0;
            }

            if( abs(ry) > (System.box_y/ 2.0) )
            {
                if(ry > 0) ry -= System.box_y / 2.0;
                else ry += System.box_y / 2.0;
            }

            if( abs(rz) > (System.box_z/ 2.0) )
            {
                if(rz > 0) rz -= System.box_z / 2.0;
                else rz += System.box_z / 2.0;
            }

            r2 = rx*rx + ry*ry + rz*rz;
            //printf("%lf for %d %d interactions\n", r2, i, j);
            if(r2 < System.cutoff2)
            {
                rij = sqrt(r2);

                force = - 48.0 * Particle[i].epsilon / r2 * ( pow(Particle[i].sigma,12) / pow(r2,6) - pow(Particle[i].sigma,6) / pow(r2,3)) - Particle[i].cutoff_derivative * rij;
                potential = 4 * Particle[i].epsilon * ( pow(Particle[i].sigma,12) / pow(r2,6) - pow(Particle[i].sigma,6) / pow(r2,3)) - Particle[i].cutoff_potential - Particle[i].cutoff_derivative * (rij - System.cutoff);

                Particle[i].fx += force * rx;
                Particle[i].fy += force * ry;
                Particle[i].fz += force * rz;
                Particle[i].potential += potential;

                Particle[j].fx += -force * rx;
                Particle[j].fy += -force * ry;
                Particle[j].fz += -force * rz;
                Particle[j].potential += potential;

                System.actual_pressure += rx * Particle[i].fx + ry * Particle[i].fy + rz * Particle[i].fz;

               //printf("force %d %d %lf \t %lf \t %lf \t %lf\n", i, j, Particle[i].fx, Particle[i].fy, Particle[i].fz, r2);
            }
            //else no new forces added
        }
    }
    System.actual_pressure += 2.0 * System.kinetic_energy;
    System.actual_pressure /= (3.0 * System.volume);
}

void compute_IMF_with_neighbors(void)
{
    int i,j;
    int w, x, y, z;
    int plus_x, plus_y, plus_z;
    int minus_y, minus_z;
    //for each cell
    //printf("total cells = %d\n", System.total_cells);
    for(i = 0; i < System.total_cells; i++)
    {
        if(Neighbor[i].length != 0)
        {

        //printf("computing same cell forces for cell %d\n", i);
        //calculate forces with_in the same cell
        compute_IMF_same_cell(i);
        //printf("computed same cell forces for cell %d\n", i);
        //calculate neighboring cell interactions via template

        //determine cell of i
        x = i % System.cellx;

        w = i - x;
        if(w > 0) y = (w / System.cellx) % System.celly;
        else y = 0;

        w -= y*System.cellx;
        if(w > 0) z = w / (System.cellx * System.celly);
        else z = 0;

        //go through neighbor list
        //increase x
        if(x == System.cellx - 1) plus_x = 1 - System.cellx;
        else plus_x = 1;

        //if(x == 0) minus_x = System.cellx - 1;
        //else minus_x = -1;

        if(y == System.celly - 1) plus_y = System.cellx *(1 - System.celly);
        else plus_y = System.cellx;

        if(y == 0) minus_y = System.cellx * (System.celly - 1);
        else minus_y = - System.cellx;

        if(z == System.cellz - 1) plus_z = System.cellx * System.celly * (1 - System.cellz);
        else plus_z = System.cellx * System.celly;

        if(z == 0) minus_z = System.cellx * System.celly * (System.cellz - 1);
        else minus_z = - System.cellx * System.celly;

        //go through the 13 combinations of neighbor cells

/*
        if(i == 4)
        {
            printf("1, %d\n", i + plus_x);
            printf("2, %d\n", i + plus_x + plus_y);
            printf("3, %d\n", i + plus_x + plus_z);
            printf("4, %d\n", i + plus_x + minus_y);
            printf("5, %d\n", i + plus_x + minus_z);
            printf("6, %d\n", i + plus_x + plus_y + minus_z);
            printf("7, %d\n", i + plus_x + minus_y + plus_z);
            printf("8, %d\n", i + plus_x + minus_y + minus_z);
            printf("9, %d\n", i + plus_x + plus_y + plus_z);
            printf("10, %d\n", i + plus_y + plus_z);
            printf("11, %d\n", i + plus_y);
            printf("12, %d\n", i + plus_y + minus_z);
            printf("13, %d\n", i + plus_z);

            printf("+x = %d, +y = %d, -y  = %d, +z = %d, -z = %d from x = %d, y = %d, z = %d\n", plus_x, plus_y, minus_y, plus_z, minus_z, x, y ,z);
        }
*/

        compute_IMF_in_cells(i, i + plus_x);                    //x+1   y   z   1
        compute_IMF_in_cells(i, i + plus_x + plus_y);           //x+1   y+1 z   2
        compute_IMF_in_cells(i, i + plus_x + plus_z);           //x+1   y   z+1 3
        compute_IMF_in_cells(i, i + plus_x + minus_y);          //x+1   y-1 z   4
        compute_IMF_in_cells(i, i + plus_x + minus_z);          //x+1   y   z-1 5
        compute_IMF_in_cells(i, i + plus_x + plus_y + minus_z); //x+1   y+1 z-1 6
        compute_IMF_in_cells(i, i + plus_x + minus_y + plus_z); //x+1   y-1 z+1 7
        compute_IMF_in_cells(i, i + plus_x + minus_y + minus_z);//x+1   y-1 z-1 8
        compute_IMF_in_cells(i, i + plus_x + plus_y + plus_z);  //x+1   y+1 z+1 9
        compute_IMF_in_cells(i, i + plus_y + plus_z);           //x     y+1 z+1 10
        compute_IMF_in_cells(i, i + plus_y);                    //x     y+1 z   11d
        compute_IMF_in_cells(i, i + plus_y + minus_z);          //x     y+1 z-1 12
        compute_IMF_in_cells(i, i + plus_z);                    //x     y   z+1 13
        }
    }
    System.actual_pressure += 2.0 * System.kinetic_energy;
    System.actual_pressure /= (3.0 * System.volume);
}

void compute_IMF_same_cell(int host)
{
    int ii, jj, i, j;
    double r2, rx, ry, rz, rij;
    double force, potential;

    for(ii = 0; ii < Neighbor[host].length - 1; ii++)
    {
        i = Neighbor[host].cell[ii];

        for(jj = ii + 1; jj < Neighbor[host].length; jj++)
        {
            j = Neighbor[host].cell[jj];
            //printf("cell %d: interaction between %d and %d\n", host, i, j);
            rx = Particle[i].x - Particle[j].x;
            ry = Particle[i].y - Particle[j].y;
            rz = Particle[i].z - Particle[j].z;

            //printf("%lf\n", (rx*rx + ry*ry + rz*rz));
            //check each for nearest image
            if(abs(rx) > (System.box_x / 2.0) )
            {
                if(rx > 0) rx -= System.box_x;
                else rx += System.box_x;
            }

            if( abs(ry) > (System.box_y/ 2.0) )
            {
                if(ry > 0) ry -= System.box_y;
                else ry += System.box_y;
            }

            if( abs(rz) > (System.box_z/ 2.0) )
            {
                if(rz > 0) rz -= System.box_z;
                else rz += System.box_z;
            }

            r2 = rx*rx + ry*ry + rz*rz;

            if(r2 < System.cutoff2)
            {
                rij = sqrt(r2);

                force = - 48.0 * Particle[i].epsilon / r2 * ( pow(Particle[i].sigma,12) / pow(r2,6) - pow(Particle[i].sigma,6) / pow(r2,3)) - Particle[i].cutoff_derivative * rij;
                potential = 4 * Particle[i].epsilon * ( pow(Particle[i].sigma,12) / pow(r2,6) - pow(Particle[i].sigma,6) / pow(r2,3)) - Particle[i].cutoff_potential - Particle[i].cutoff_derivative * (rij - System.cutoff);

                Particle[i].fx += force * rx;
                Particle[i].fy += force * ry;
                Particle[i].fz += force * rz;
                Particle[i].potential += potential;

                Particle[j].fx += -force * rx;
                Particle[j].fy += -force * ry;
                Particle[j].fz += -force * rz;
                Particle[j].potential += potential;

                System.actual_pressure += rx * Particle[i].fx + ry * Particle[i].fy + rz * Particle[i].fz;
            }
            //otherwise force is zero

            //printf("force i=%d j=%d fx=%lf fy=%lf fz=%lf r2=%lf\n", i, j, Particle[i].fx*pow(10,5), Particle[i].fy*pow(10,5), Particle[i].fz*pow(10,5), r2);
        }
    }
}

void compute_IMF_in_cells(int host, int other)
{
    if(Neighbor[other].length != 0)
    {

    //printf("compute pair cells %d %d \n", host, other);
    int ii, jj, i, j;
    double r2, rx, ry, rz, rij;
    double force, potential;

    for(ii = 0; ii < Neighbor[host].length; ii++)
    {
        i = Neighbor[host].cell[ii];
        for(jj = 0; jj < Neighbor[other].length; jj++)
        {
            j = Neighbor[other].cell[jj];

            rx = Particle[i].x - Particle[j].x;
            ry = Particle[i].y - Particle[j].y;
            rz = Particle[i].z - Particle[j].z;

            //printf("%lf\n", (rx*rx + ry*ry + rz*rz));
            //check each for nearest image
            if(abs(rx) > (System.box_x / 2.0) )
            {
                if(rx > 0) rx -= System.box_x;
                else rx += System.box_x;
            }

            if( abs(ry) > (System.box_y/ 2.0) )
            {
                if(ry > 0) ry -= System.box_y;
                else ry += System.box_y;
            }

            if( abs(rz) > (System.box_z/ 2.0) )
            {
                if(rz > 0) rz -= System.box_z;
                else rz += System.box_z;
            }

            r2 = rx*rx + ry*ry + rz*rz;

            if(r2 < System.cutoff2)
            {
                rij = sqrt(r2);

                force = - 48.0 * Particle[i].epsilon / r2 * ( pow(Particle[i].sigma,12) / pow(r2,6) - pow(Particle[i].sigma,6) / pow(r2,3)) - Particle[i].cutoff_derivative * rij;
                potential = 4 * Particle[i].epsilon * ( pow(Particle[i].sigma,12) / pow(r2,6) - pow(Particle[i].sigma,6) / pow(r2,3)) - Particle[i].cutoff_potential - Particle[i].cutoff_derivative * (rij - System.cutoff);

                Particle[i].fx += force * rx;
                Particle[i].fy += force * ry;
                Particle[i].fz += force * rz;
                Particle[i].potential += potential;

                Particle[j].fx += -force * rx;
                Particle[j].fy += -force * ry;
                Particle[j].fz += -force * rz;
                Particle[j].potential += potential;

                System.actual_pressure += rx * Particle[i].fx + ry * Particle[i].fy + rz * Particle[i].fz;
            }
            //otherwise force is zero

            //printf("force %d %d %lf \t %lf \t %lf \t %lf\n", i, j, Particle[i].fx, Particle[i].fy, Particle[i].fz, r2);
        }
    }
    }
}

//neighborlist initialization and repopulation
void initialize_neighbor_list(void)
{
    double inv_cell_size = 1 / (System.neighbor_factor * System.cutoff);

    System.cellx = ceil(System.box_x * inv_cell_size);
    System.celly = ceil(System.box_y * inv_cell_size);
    System.cellz = ceil(System.box_z * inv_cell_size);
    System.total_cells = System.cellx * System.celly * System.cellz;
    //printf("totalcells:%d from x:%d of %lf y:%d of %lf z:%d of %lf\n", System.total_cells, System.cellx, System.box_x, System.celly, System.box_y, System.cellz, System.box_z);

    //printf("CELLS: x %d y %d z %d", System.cellx, System.celly, System.cellz);
    Neighbor = malloc( System.total_cells * sizeof(NeighborList));
    int i;
    for(i=0; i < System.total_cells; i++)
    {
        Neighbor[i].length = 0;
    }
}

void build_neighbor_list(void)
{
    int i, x, y, z;
    int j;
    double inv_cell_size = 1/ (System.neighbor_factor * System.cutoff);

    //printf("\nBuild neighbor list with %d cells \n\n", System.total_cells);
    //reset Neighbor counts


    for(i = 0; i < System.total_cells; i++)
    {
        Neighbor[i].length = 0;

        if(System.counter != 0)
        {
            free(Neighbor[i].cell);
        }
    }

    if(System.ensemble_flag == 2)
    {

            //determine num number of cells
            j = ceil(System.box_x * inv_cell_size) + ceil(System.box_y * inv_cell_size) + ceil(System.box_z * inv_cell_size);

            if(j != System.total_cells)
            {
                free(Neighbor);
                initialize_neighbor_list();
                printf("re_intialized neighbor list with %d cells\n", System.total_cells);
            }
    }

    //label each particle by cell
    for(i = 0; i < System.num_particles; i++)
    {
        x = (int) (Particle[i].x * inv_cell_size);
        y = (int) (Particle[i].y * inv_cell_size);
        z = (int) (Particle[i].z * inv_cell_size);

        Particle[i].cell_id = x + y * System.cellx + z * System.cellx * System.celly;
        Neighbor[Particle[i].cell_id].length++;
        //printf("particle %d is in cell %d\n", i, Particle[i].cell_id);
        //printf("cells: x = %d, y = %d, z = %d with total cells = %d\n", x, y, z, System.total_cells);
    }

    //printf("calculated cells\n");
    //Allocate neighbor_list array space
    for(i = 0; i < System.total_cells; i++)
    {
        //allocate memory for cell index list
        //printf("Neighbor[%d].cell is size %d ints \n", i, Neighbor[i].length);
        Neighbor[i].cell = malloc(Neighbor[i].length * sizeof(int));

        //reset length for space counter
        Neighbor[i].length = 0;
    }
    //printf("allocated cells\n");

    //Populate neighbor_list indecies
    for(i = 0; i < System.num_particles; i++)
    {
        //printf("add %d to Neighbor[%d].cell[%d]\n", i, Particle[i].cell_id, Neighbor[Particle[i].cell_id].length);
        Neighbor[ Particle[i].cell_id ].cell[ Neighbor[ Particle[i].cell_id ].length ] = i;
        Neighbor[ Particle[i].cell_id ].length++;
    }
    //printf("assigned cells\n");
    //printf("\n");
    for(i = 0; i < System.total_cells; i++)
    {
        for(j = 0; j < Neighbor[i].length; j++ )
        {
            //printf("Neighbor %d has %d atoms\n", i, Neighbor[i].length);
            //printf("particle id is %d\n", Neighbor[i].cell[j]);
        }
    }
    //printf("\n");
}

void total_energy(void)
{
    int i;
    System.kinetic_energy = 0.0;
    System.potential_energy = 0.0;

    for(i = 0; i < System.num_particles; i++)
    {
        System.kinetic_energy += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
        System.potential_energy += Particle[i].potential;
    }
   // printf("\nKinetic energy = %lf and Potential Energy = %lf\n\n", System.kinetic_energy, System.potential_energy);
   // printf("Pressure is %lf\n", System.actual_pressure);

    //printf("output at counter %d\n", System.counter);
    if(System.counter == 0)
    {
        initialize_energy_output();
        printf("output initialized\n");
    }
    else energy_output();
}
