//initialize.c

void    initialize_system(void)
{
    int file_info = 1;
    System.kb = 8.31 * pow(10, -3);
    System.pressure = 0.0;

    //read control variables from file
    read_controls();


    //get coordinates and forces/velocities
    initialize_from_file();
    initialize_neighbor_list();
    intialize_VMD_output();
}

void read_controls(void)
{
   FILE* fr;
    //Control System;

    char line[80];
    char junk;
    System.integrator = 1;
    System.nose_hoover = 0.0;
    System.npt_pressure = 0.0;

    fr = fopen("001.in", "rt"); //need to set intput file to something

    fgets(line,80,fr);
    sscanf(line, "%d", &System.restart_flag);
    fgets(line,80,fr);
    sscanf(line, "%d", &System.ensemble_flag);
    fgets(line,80,fr);
    sscanf(line, "%d", &System.integrator);

    fgets(line,80,fr);
    sscanf(line, "%s", &junk);

    fgets(line,80,fr);
    sscanf(line, "%lf", &System.temp);
    fgets(line,80,fr);
    sscanf(line, "%lf", &System.pressure);

    fgets(line,80,fr);
    sscanf(line, "%s", &junk);

    fgets(line,80,fr);
    sscanf(line, "%d", &System.n_steps);
    fgets(line,80,fr);
    sscanf(line, "%lf", &System.delta_t);

    fgets(line,80,fr);
    sscanf(line, "%s", &junk);

    fgets(line,80,fr);
    sscanf(line, "%lf", &System.cutoff);
    fgets(line,80,fr);
    sscanf(line, "%lf", &System.neighbor_factor);

    fclose(fr);

    System.cutoff2 = System.cutoff * System.cutoff;

    printf("from control file *.in \n");
    printf("%d %d \n", System.restart_flag, System.ensemble_flag);
    printf("%lf %lf \n", System.temp, System.pressure);
    printf("%d %lf %lf \n", System.n_steps, System.delta_t, System.cutoff);
    printf("finished *.in \n");
}

void initialize_from_file(void)
{
   //read in topology file
   initialize_topology_from_file();

    //create types based on topology
    create_particles();

   //read in coordinate file
   initialize_coordinates_from_file();
}

void initialize_internally(void)
{
//ignore
}

void initialize_topology_from_file(void)
{
    FILE* fr;

    char line[80];
    char junk;
    int i, j;
    //double other, other2, other3, other4, other5, other6;

    printf("\n");
    printf("reading from topology file \n");


    fr = fopen("001.top", "rt"); //need to set intput file to something
    fgets(line,80,fr);
    sscanf(line,"%d", &System.molecule_types);   // need to define array for next part
    printf("%d molecule types\n", System.molecule_types);

    //allocate arrays for types
    System.types = malloc(System.molecule_types * sizeof(int));

    fgets(line,80,fr);
    sscanf(line, "%s", &junk);

    //allocate array for Prototype
    Prototype = malloc(System.molecule_types * sizeof(Prototypes));
    //printf("sizeof(Prototypes) = %d", sizeof(Prototypes));

    for(i=0; i < System.molecule_types; i++)
    {
        fgets(line,80,fr);
        sscanf(line, "%d %d", &System.types[i].atoms, &System.types[i].copies);
    }


    fgets(line,80,fr);
    sscanf(line, "%s", &junk);

    for (i=0; i < System.molecule_types; i++)
    {
        printf("System.types[%d].atoms = %d \n", i, System.types[i].atoms);
        for (j=0; j < System.types[i].atoms; j++)
        {
            fgets(line,80,fr);
            //printf("%s \n", line);
            sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf", Prototype[i].type, &Prototype[i].sigma, &Prototype[i].epsilon, &Prototype[i].x_rel, &Prototype[i].y_rel, &Prototype[i].z_rel, &Prototype[i].charge, &Prototype[i].mass);
            printf("read %s %lf %lf %lf %lf %lf %lf %lf \n", Prototype[i].type, Prototype[i].sigma, Prototype[i].epsilon, Prototype[i].x_rel, Prototype[i].y_rel, Prototype[i].z_rel, Prototype[i].charge, Prototype[i].mass);

            //scale epsilon over kbt
            Prototype[i].epsilon *= System.kb;
            Prototype[i].cutoff_potential = 4 * Prototype[i].epsilon * ( pow(Prototype[i].sigma,12) / pow(System.cutoff, 12) - pow(Prototype[i].sigma,6) / pow(System.cutoff, 6) );
            Prototype[i].cutoff_derivative =  48 * Prototype[i].epsilon * (  pow(Prototype[i].sigma,12) / pow(System.cutoff, 13) - pow(Prototype[i].sigma,6) / (2 * pow(System.cutoff, 7) ) );
        }

        fgets(line,80,fr);
        sscanf(line, "%s", &junk);
    }

    printf("closing file\n");
    fclose(fr);
    printf("file closed\n");

    printf("%d molecule types \n", System.molecule_types);

    for(i=0; i < System.molecule_types; i++)
    {
        printf("%d atoms and %d molecules \n", System.types[i].atoms,  System.types[i].copies);
    }



    for (i=0; i < System.molecule_types; i++)
    {
        for (j=0; j < System.types[i].atoms; j++)
        {
            printf("%s \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", Prototype[i].type, Prototype[i].sigma, Prototype[i].epsilon, Prototype[i].x_rel, Prototype[i].y_rel, Prototype[i].z_rel, Prototype[i].charge, Prototype[i].mass);
        }
    }

    printf("finished topology file *.top \n");

}

void create_particles(void)
{
    int n = 0;
    int i, j, m;
    double ke = 0.0;
    double v_range = 0.0;
    double net_x = 0.0;
    double net_y = 0.0;
    double net_z = 0.0;

    //determine total number of particles
    for(i = 0; i < System.molecule_types; i++)
    {
        n += System.types[i].copies;
    }

    //allocate particle array
    Particle = malloc(n * sizeof(Particles));

    //initialize defaults for each particle by type
    m = 0;
    for(i = 0; i < System.molecule_types; i++)
    {
        for(j = 0; j < System.types[i].copies; j++)
        {
            *Particle[m].type = *Prototype[i].type;
            Particle[m].style = i;
            Particle[m].sigma = Prototype[i].sigma;
            Particle[m].epsilon = Prototype[i].epsilon;
            Particle[m].x = Prototype[i].x_rel;
            Particle[m].y = Prototype[i].y_rel;
            Particle[m].z = Prototype[i].z_rel;
            Particle[m].charge = Prototype[i].charge;
            Particle[m].mass = Prototype[i].mass;
            Particle[m].cutoff_potential = Prototype[i].cutoff_potential;
            Particle[m].cutoff_derivative = Prototype[i].cutoff_derivative;
            m++;
        }
    }

    //initialize velocities randomly
    for(i = 0; i < n; i++)
    {
        v_range = sqrt( 3.0 * System.kb * System.temp  / Particle[i].mass);
        Particle[i].vx = v_range * (Random() * 2.0 - 1.0); //these velocities should be initialized randomly ...
        Particle[i].vy = v_range * (Random() * 2.0 - 1.0);
        Particle[i].vz = v_range * (Random() * 2.0 - 1.0);

        net_x += Particle[i].vx;
        net_y += Particle[i].vy;
        net_z += Particle[i].vz;

        printf("vx = %lf, vy = %lf, vz = %lf with i = %d and n = %d\n", Particle[i].vx, Particle[i].vy, Particle[i].vz, i, n);
        ke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz) / n;
    }
    printf("\n1st KE = %lf with excess X: %lf Y: %lf Z: %lf\n", ke, net_x, net_y, net_z);
    ke = 0.0;
    net_x /= n;

    for(i = 0; i < n; i++)
    {
        Particle[i].vx -= net_x;
        Particle[i].vy -= net_y;
        Particle[i].vz -= net_z;
        ke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz) / n;
    }
    printf("next KE = %lf with %lf\n", ke, (1.50 * System.kb * System.temp));
    net_x = sqrt( 1.50 * System.kb * System.temp / ke);
    ke = 0.0;

    for(i = 0; i < n; i++)
    {
        Particle[i].vx *= net_x;
        Particle[i].vy *= net_x;
        Particle[i].vz *= net_x;
        ke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz) / n;
    }
    printf("final KE = %lf with target %lf\n", ke, (1.5 * System.kb * System.temp));
    System.kinetic_energy = ke;
}


void initialize_coordinates_from_file(void)
{
    FILE* fr;

    char line[80];
    char junk;
    int i;
    char type[3];

    printf("\nstart reading coordinate file 001.xyz\n");

    fr = fopen("001.xyz", "rt"); //need to set intput file to something

    printf("open file \n");
    fgets(line,80,fr);
    sscanf(line, "%d", &System.num_particles);
    printf("read particles %s = %d\n", line, System.num_particles);

    fgets(line,80,fr);
    sscanf(line, "%s", &junk);
    fgets(line,80,fr);
    sscanf(line,"%lf %lf %lf %lf %lf %lf", &System.box_x, &System.box_y, &System.box_z, &System.angle1, &System.angle2, &System.angle3);
    printf("%lf %lf %lf %lf %lf %lf\n",System.box_x, System.box_y, System.box_z, System.angle1, System.angle2, System.angle3);

    System.volume = System.box_x * System.box_y * System.box_z;

    fgets(line,80,fr);
    sscanf(line, "%s", &junk);

    printf("System.num_particles = %d \n", System.num_particles);
    for(i = 0; i < System.num_particles; i++)
    {
        fgets(line,80,fr);
        sscanf(line, "%s %lf %lf %lf", type, &Particle[i].x, &Particle[i].y, &Particle[i].z);

        printf("LINE%d: %s \t %lf \t %lf \t %lf\n", i, type, Particle[i].x, Particle[i].y, Particle[i].z);
    }

    fclose(fr);

    printf("finished coordinate files\n");

}
