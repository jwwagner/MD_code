//integrate.c

//NVE: Verlet, Velocity Verlet
//NVT: Velocity Rescaling, Berendsen, Nose-Hoover
//NPT: Nose-Poincare-Anderson (NPA) with different NPT subsets

void time_step(void)
{
    if(System.ensemble_flag == 1) // NVT
    {
        if(System.integrator == 0) do_nose_hoover();                // Nose-Hoover thermostat
        else {
             if(System.integrator == 1) do_bendersen_thermostat();  // Bendersen Thermostat
             else do_NVT_velocity_verlet();                         // Simple velocity rescaling
             }
    }

    else
    {
        if(System.ensemble_flag == 2) // NPT
        {
            //if(System.integrator == 0) do_NPT_nose_hoover();
            //else do_NPT_bendersen();
            if(System.integrator == 0)
                {
                    if(System.counter % 2 == 0) do_NPT_berendsen_method_Temp();
                    else do_NPT_berendsen_method_Pressure();
                }
            if(System.integrator == 1)
                {
                    if(System.counter % 2 == 0) do_NPT_AT_method_Temp();
                    else do_NPT_AT_method_Pressure();
                }
            if(System.integrator == 2) do_NPT_AT_method();
        }
        else {                        // NVE
            if(System.integrator == 0)
            {
                do_verlet();
            }
            else
            {
                do_velocity_verlet();
            }
        }
   }


    //do_leapfrog_verlet();
    //do_predictor_corrector();
    //printf("finished integrating\n"); /**/

    wrap_PBC();

    //printf("PBC wrapped\n");/**/

}


void    do_verlet(void)
{
    int i;
    double tempx, tempy, tempz;

    if(System.counter == 0)
    {
        for(i=0; i < System.num_particles; i++)
        {
            Particle[i].oldx = Particle[i].x;
            Particle[i].oldy = Particle[i].y;
            Particle[i].oldz = Particle[i].z;

            //temp is "old" coordinates
            tempx = Particle[i].x - Particle[i].vx * System.delta_t;
            tempy = Particle[i].y - Particle[i].vy * System.delta_t;
            tempz = Particle[i].z - Particle[i].vz * System.delta_t;

            Particle[i].x = 2.0 * Particle[i].x - tempx + System.delta_t * System.delta_t * Particle[i].fx / Particle[i].mass;
            Particle[i].y = 2.0 * Particle[i].y - tempy + System.delta_t * System.delta_t * Particle[i].fy / Particle[i].mass;
            Particle[i].z = 2.0 * Particle[i].z - tempz + System.delta_t * System.delta_t * Particle[i].fz / Particle[i].mass;

            Particle[i].vx = (Particle[i].oldx - tempx) / (2.0 * System.delta_t);
            Particle[i].vy = (Particle[i].oldy - tempy) / (2.0 * System.delta_t);
            Particle[i].vz = (Particle[i].oldz - tempz) / (2.0 * System.delta_t);
        }
    }
    else
    {
        for(i=0; i < System.num_particles; i++)
        {
            //printf("x = %lf, vx = %lf, ax = %lf\n", Particle[i].x, Particle[i].vx, Particle[i].fx);/**/
            //temp is "new" coordinates
            tempx = 2.0 * Particle[i].x - Particle[i].oldx + System.delta_t * System.delta_t * Particle[i].fx / Particle[i].mass;
            tempy = 2.0 * Particle[i].y - Particle[i].oldy + System.delta_t * System.delta_t * Particle[i].fy / Particle[i].mass;
            tempz = 2.0 * Particle[i].z - Particle[i].oldz + System.delta_t * System.delta_t * Particle[i].fz / Particle[i].mass;

            //make sure old coordinates are wraped
            if(Particle[i].oldx > System.box_x) Particle[i].oldx -= System.box_x;
            if(Particle[i].oldy > System.box_y) Particle[i].oldy -= System.box_x;
            if(Particle[i].oldz > System.box_z) Particle[i].oldz -= System.box_x;

            if(Particle[i].oldx < 0.0) Particle[i].oldx += System.box_x;
            if(Particle[i].oldy < 0.0) Particle[i].oldy += System.box_x;
            if(Particle[i].oldz < 0.0) Particle[i].oldz += System.box_x;


            Particle[i].vx = (Particle[i].oldx - tempx) / (2.0 * System.delta_t);
            Particle[i].vy = (Particle[i].oldy - tempy) / (2.0 * System.delta_t);
            Particle[i].vz = (Particle[i].oldz - tempz) / (2.0 * System.delta_t);

            Particle[i].oldx = Particle[i].x;
            Particle[i].oldy = Particle[i].y;
            Particle[i].oldz = Particle[i].z;

            Particle[i].x = tempx;
            Particle[i].y = tempy;
            Particle[i].z = tempz;
        }
    }
}

void do_velocity_verlet(void)
{
    int i;
    double tempx, tempy, tempz;

    if(System.counter == 0) //start-up
    {
        for(i=0; i < System.num_particles; i++)
        {
            //temp is "old" coordinates
            Particle[i].oldx = Particle[i].x - Particle[i].vx * System.delta_t;
            Particle[i].oldy = Particle[i].y - Particle[i].vy * System.delta_t;
            Particle[i].oldz = Particle[i].z - Particle[i].vz * System.delta_t;

            Particle[i].x += Particle[i].vx * System.delta_t + 0.5 * System.delta_t * System.delta_t * Particle[i].fx / Particle[i].mass;
            Particle[i].y += Particle[i].vy * System.delta_t + 0.5 * System.delta_t * System.delta_t * Particle[i].fy / Particle[i].mass;
            Particle[i].z += Particle[i].vz * System.delta_t + 0.5 * System.delta_t * System.delta_t * Particle[i].fz / Particle[i].mass;

            //use normal Verlet as proxy for start-up velocity since its hard to get previous accelleration/force
            Particle[i].vx = (Particle[i].oldx - Particle[i].x) / (2.0 * System.delta_t);
            Particle[i].vy = (Particle[i].oldy - Particle[i].y) / (2.0 * System.delta_t);
            Particle[i].vz = (Particle[i].oldz - Particle[i].z) / (2.0 * System.delta_t);

            Particle[i].oldfx = Particle[i].fx;
            Particle[i].oldfy = Particle[i].fy;
            Particle[i].oldfz = Particle[i].fz;
        }
    }
    else
    {
        for(i=0; i < System.num_particles; i++)
        {
            //temp is "new" coordinates
            //printf("x = %lf, vx = %lf, ax = %lf\n", Particle[i].x, Particle[i].vx, Particle[i].fx);
            Particle[i].x += Particle[i].vx * System.delta_t + System.delta_t * System.delta_t * Particle[i].fx / (2.0 * Particle[i].mass);
            Particle[i].y += Particle[i].vy * System.delta_t + System.delta_t * System.delta_t * Particle[i].fy / (2.0 * Particle[i].mass);
            Particle[i].z += Particle[i].vz * System.delta_t + System.delta_t * System.delta_t * Particle[i].fz / (2.0 * Particle[i].mass);

            Particle[i].vx += (Particle[i].fx + Particle[i].oldfx) * System.delta_t / (2.0 * Particle[i].mass);
            Particle[i].vy += (Particle[i].fy + Particle[i].oldfy) * System.delta_t / (2.0 * Particle[i].mass);
            Particle[i].vz += (Particle[i].fz + Particle[i].oldfz) * System.delta_t / (2.0 * Particle[i].mass);

            Particle[i].oldfx = Particle[i].fx;
            Particle[i].oldfy = Particle[i].fy;
            Particle[i].oldfz = Particle[i].fz;
        }
    }
}

void    do_leapfrog_verlet(void)
{

}

void do_NVT_velocity_verlet(void)
{
    int i;
    double tempx, tempy, tempz;

    System.kinetic_energy = 0.0;

    for(i = 0; i < System.num_particles; i++)
    {
        System.kinetic_energy += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    double factor = sqrt(System.temp / (2.0 * System.kinetic_energy /(3.0 * System.num_particles * System.kb) ) );

    //printf("factor %lf with temp %lf and relative %lf\n", factor, System.temp, (2.0 * System.kinetic_energy / ( 3.0 * System.num_particles * System.kb)) );/**/
    if(System.counter == 0) //start-up
    {
        for(i=0; i < System.num_particles; i++)
        {
            //temp is "old" coordinates
            Particle[i].oldx = Particle[i].x - Particle[i].vx * System.delta_t;
            Particle[i].oldy = Particle[i].y - Particle[i].vy * System.delta_t;
            Particle[i].oldz = Particle[i].z - Particle[i].vz * System.delta_t;

            Particle[i].x += Particle[i].vx * System.delta_t + 0.5 * System.delta_t * System.delta_t * Particle[i].fx / Particle[i].mass;
            Particle[i].y += Particle[i].vy * System.delta_t + 0.5 * System.delta_t * System.delta_t * Particle[i].fy / Particle[i].mass;
            Particle[i].z += Particle[i].vz * System.delta_t + 0.5 * System.delta_t * System.delta_t * Particle[i].fz / Particle[i].mass;

            //use normal Verlet as proxy for start-up velocity since its hard to get previous accelleration/force
            //temp rescaling not needed for first step because velocities already satisfy this
            Particle[i].vx = (Particle[i].oldx - Particle[i].x) / (2.0 * System.delta_t);
            Particle[i].vy = (Particle[i].oldy - Particle[i].y) / (2.0 * System.delta_t);
            Particle[i].vz = (Particle[i].oldz - Particle[i].z) / (2.0 * System.delta_t);

            Particle[i].oldfx = Particle[i].fx;
            Particle[i].oldfy = Particle[i].fy;
            Particle[i].oldfz = Particle[i].fz;
        }
    }
    else
    {
        for(i=0; i < System.num_particles; i++)
        {
            //temp is "new" coordinates
            //printf("x = %lf, vx = %lf, ax = %lf\n", Particle[i].x, Particle[i].vx, Particle[i].fx);
            Particle[i].x += Particle[i].vx * System.delta_t + System.delta_t * System.delta_t * Particle[i].fx / (2.0 * Particle[i].mass);
            Particle[i].y += Particle[i].vy * System.delta_t + System.delta_t * System.delta_t * Particle[i].fy / (2.0 * Particle[i].mass);
            Particle[i].z += Particle[i].vz * System.delta_t + System.delta_t * System.delta_t * Particle[i].fz / (2.0 * Particle[i].mass);

            Particle[i].vx = Particle[i].vx * factor + (Particle[i].fx + Particle[i].oldfx) * System.delta_t / (2.0 * Particle[i].mass);
            Particle[i].vy = Particle[i].vy * factor + (Particle[i].fy + Particle[i].oldfy) * System.delta_t / (2.0 * Particle[i].mass);
            Particle[i].vz = Particle[i].vz * factor + (Particle[i].fz + Particle[i].oldfz) * System.delta_t / (2.0 * Particle[i].mass);

            Particle[i].oldfx = Particle[i].fx;
            Particle[i].oldfy = Particle[i].fy;
            Particle[i].oldfz = Particle[i].fz;
        }
    }
}

void    do_bendersen_thermostat(void)    //CONSIDER CORRECTING Q VALUE IN TERMS OF TAU
{
    int i,j;
    double newke = 0.0;
    System.kinetic_energy = 0.0;
    System.potential_energy = 0.0;
    for(i = 0; i < System.num_particles; i++)
    {
        System.kinetic_energy += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
        System.potential_energy += Particle[i].potential;
    }

    if(System.counter == 0)
    {
        System.nose_Q = 3.0 * ((double)System.num_particles) * System.kb * System.temp * (1.0) * Particle[1].mass * Particle[1].sigma * Particle[1].sigma / Particle[1].epsilon;
        //printf("B:Initialize Nose-Hoover as %lf \n", System.nose_hoover);
        //printf("Initialize Nose-Hoover as %lf with N:%d T:%lf M:%lf S:%lf E:%lf \n", System.nose_hoover, System.num_particles, System.temp, Particle[1].mass, Particle[1].sigma, Particle[1].epsilon);
    }

    //integrate nose_hoover term
    System.nose_hoover = System.delta_t * (System.kinetic_energy - 1.5 * System.kb * System.num_particles * System.temp) / System.nose_Q;
    //printf("Nose Hoover is %lf\n", System.nose_hoover);

    //integrate position and momentum
    for(i = 0; i < System.num_particles;i++)
    {
        //position
        Particle[i].x += Particle[i].vx * System.delta_t;
        Particle[i].y += Particle[i].vy * System.delta_t;
        Particle[i].z += Particle[i].vz * System.delta_t;

        //velocity
        Particle[i].vx += Particle[i].fx / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vx * System.delta_t;
        Particle[i].vy += Particle[i].fy / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vy * System.delta_t;
        Particle[i].vz += Particle[i].fz / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vz * System.delta_t;
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    //printf("Target: %lf OldKE: %lf NewKE: %lf B-NH: %lf\n",(1.50 * System.num_particles * System.kb * System.temp), System.kinetic_energy, newke, System.nose_hoover );
    //printf("used KE %lf and Target %lf\n", System.kinetic_energy, (1.5 * System.kb * System.num_particles * System.temp) );

}

void do_nose_hoover(void) // IMPLEMENT NOSE HOOVER ALGORITHM
{
    int i,j;
    double newke = 0.0;
    System.kinetic_energy = 0.0;
    System.potential_energy = 0.0;
    for(i = 0; i < System.num_particles; i++)
    {
        System.kinetic_energy += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
        System.potential_energy += Particle[i].potential;
    }
    //integrate nose_hoover term
    if(System.counter == 0)
    {
        System.nose_Q = 3.0 * ((double)System.num_particles) * System.kb * System.temp * (0.09622 * 0.09622) * Particle[1].mass * Particle[1].sigma * Particle[1].sigma / Particle[1].epsilon;
        //printf("NH:Initialize Nose-Hoover as %lf \n", System.nose_hoover);
        //printf("Initialize Nose-Hoover as %lf with N:%d T:%lf M:%lf S:%lf E:%lf \n", System.nose_hoover, System.num_particles, System.temp, Particle[1].mass, Particle[1].sigma, Particle[1].epsilon);
    }

    System.nose_hoover += System.delta_t * (2.0 * System.kinetic_energy - 3.0 * System.kb * System.num_particles * System.temp) / System.nose_Q;
    //printf("Nose Hoover is %lf\n", System.nose_hoover);

    //integrate position and momentum
    for(i = 0; i < System.num_particles;i++)
    {
        //position
        Particle[i].x += Particle[i].vx * System.delta_t;
        Particle[i].y += Particle[i].vy * System.delta_t;
        Particle[i].z += Particle[i].vz * System.delta_t;

        //velocity
        Particle[i].vx += Particle[i].fx / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vx * System.delta_t;
        Particle[i].vy += Particle[i].fy / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vy * System.delta_t;
        Particle[i].vz += Particle[i].fz / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vz * System.delta_t;
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    //printf("Target: %lf OldKE: %lf NewKE: %lf NH: %lf\n",(1.50 * System.num_particles * System.kb * System.temp), System.kinetic_energy, newke, System.nose_hoover );
    //printf("used KE %lf and Target %lf\n", System.kinetic_energy, (1.5 * System.kb * System.num_particles * System.temp) );

}

void do_NPT_nose_hoover(void)
{
    int i = 0;
    double dVolume = 0.0;
    double dPressure = 0.0;
    double dLength = 0.0;
    double newke = 0.0;
    System.kinetic_energy = 0.0;
    System.npt_mass = System.pressure;

    for(i = 0; i < System.num_particles; i++)
    {
        System.kinetic_energy += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    dPressure = System.actual_pressure - System.pressure - System.npt_pressure;
    dVolume = dPressure / System.npt_mass;
    dLength = 1.0 + pow( dVolume, (1.0 / 3.0) );
    printf("dVolume %lf, dPressure %lf, dLength %lf \n at pressure %lf with target %lf and npt %lf\n\n", dVolume, dPressure, dLength, System.actual_pressure, System.pressure, System.npt_pressure);
    //printf("target pressure is %lf\n", System.pressure  );
    if(System.counter == 0)
    {
        System.nose_Q = 3.0 * ((double)System.num_particles) * System.kb * System.temp * (0.09622 * 0.09622) * Particle[1].mass * Particle[1].sigma * Particle[1].sigma / Particle[1].epsilon;
        System.nose_hoover = 0.0;
        System.npt_pressure = 0.0;
        System.npt_nose = 0.0;

        printf("B:Initialize NPT-Nose-Hoover as %lf, dP: %lf, dV: %lf\n", System.nose_hoover, dPressure, dVolume);
        //printf("Initialize Nose-Hoover as %lf with N:%d T:%lf M:%lf S:%lf E:%lf \n", System.nose_hoover, System.num_particles, System.temp, Particle[1].mass, Particle[1].sigma, Particle[1].epsilon);
    }

    //integrate position and momentum
    for(i = 0; i < System.num_particles;i++)
    {
        //position
        Particle[i].x += Particle[i].vx * System.delta_t + Particle[i].x * System.delta_t * dVolume / (3.0 * System.volume);   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y += Particle[i].vy * System.delta_t + Particle[i].y * System.delta_t * dVolume / (3.0 * System.volume);   //
        Particle[i].z += Particle[i].vz * System.delta_t + Particle[i].z * System.delta_t * dVolume / (3.0 * System.volume);   //

        //velocity
        Particle[i].vx += Particle[i].fx / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vx * System.delta_t - Particle[i].vx * dVolume / (3.0 * System.volume);
        Particle[i].vy += Particle[i].fy / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vy * System.delta_t - Particle[i].vy * dVolume / (3.0 * System.volume);
        Particle[i].vz += Particle[i].fz / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vz * System.delta_t - Particle[i].vz * dVolume / (3.0 * System.volume);
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }
    printf("Target: %lf OldKE: %lf NewKE: %lf NH: %lf\n",(1.50 * System.num_particles * System.kb * System.temp), System.kinetic_energy, newke, System.nose_hoover );
    printf("used KE %lf and Target %lf\n", System.kinetic_energy, (1.5 * System.kb * System.num_particles * System.temp) );

    //integrate nose_hoover term & integrate temperature
    System.npt_nose += System.delta_t * (2.0 * System.kinetic_energy + System.npt_nose * System.npt_nose / System.nose_Q - 3.0 * System.num_particles * System.kb * System.temp);
    System.nose_hoover += System.delta_t * System.npt_nose / System.npt_nose;

    System.nose_hoover = System.delta_t * (System.kinetic_energy - 1.5 * System.kb * System.num_particles * System.temp) / System.nose_Q;
    //printf("Nose Hoover is %lf\n", System.nose_hoover);


    //integrate volume
    System.volume += dVolume * System.delta_t;
    System.box_x *= dLength;
    System.box_y *= dLength;
    System.box_z *= dLength;
    printf("new box V:%lf, x:%lf, y:%lf, z:%lf\n", System.volume, System.box_x, System.box_y, System.box_z);

    //integrate pressure
    System.npt_pressure += dPressure * System.delta_t;
}

void do_NPT_bendersen(void)
{
    int i = 0;
    double dVolume = 0.0;
    double dPressure = 0.0;
    double dLength = 0.0;
    double newke = 0.0;
    System.kinetic_energy = 0.0;

    for(i = 0; i < System.num_particles; i++)
    {
        System.kinetic_energy += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }


    dPressure = System.actual_pressure - System.pressure - System.npt_pressure;
    dVolume = dPressure / System.npt_mass;
    dLength = 1.0 + pow( System.volume, (1.0 / 3.0) );

    if(System.counter == 0)
    {
        System.nose_Q = 3.0 * ((double)System.num_particles) * System.kb * System.temp * (0.09622 * 0.09622) * Particle[1].mass * Particle[1].sigma * Particle[1].sigma / Particle[1].epsilon;
        System.nose_hoover = 0.0;
        System.npt_pressure = 0.0;
        System.npt_nose = 0.0;

        printf("B:Initialize NPT-Nose-Hoover as %lf, dP: %lf, dV: %lf\n", System.nose_hoover, dPressure, dVolume);
        //printf("Initialize Nose-Hoover as %lf with N:%d T:%lf M:%lf S:%lf E:%lf \n", System.nose_hoover, System.num_particles, System.temp, Particle[1].mass, Particle[1].sigma, Particle[1].epsilon);
    }

    //integrate nose_hoover term & integrate temperature
    System.npt_nose += System.delta_t * (2.0 * System.kinetic_energy + System.npt_nose * System.npt_nose / System.nose_Q - 3.0 * System.num_particles * System.kb * System.temp);
    System.nose_hoover = System.delta_t * System.npt_nose / System.npt_nose;

    System.nose_hoover = System.delta_t * (System.kinetic_energy - 1.5 * System.kb * System.num_particles * System.temp) / System.nose_Q;
    //printf("Nose Hoover is %lf\n", System.nose_hoover);

    //integrate position and momentum
    for(i = 0; i < System.num_particles;i++)
    {
        //position
        Particle[i].x += Particle[i].vx * System.delta_t + Particle[i].x * System.delta_t * dVolume / (3.0 * System.volume);   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y += Particle[i].vy * System.delta_t + Particle[i].y * System.delta_t * dVolume / (3.0 * System.volume);   //
        Particle[i].z += Particle[i].vz * System.delta_t + Particle[i].z * System.delta_t * dVolume / (3.0 * System.volume);   //

        //velocity
        Particle[i].vx += Particle[i].fx / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vx * System.delta_t - Particle[i].vx * dVolume / (3.0 * System.volume);
        Particle[i].vy += Particle[i].fy / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vy * System.delta_t - Particle[i].vy * dVolume / (3.0 * System.volume);
        Particle[i].vz += Particle[i].fz / Particle[i].mass * System.delta_t - System.nose_hoover * Particle[i].vz * System.delta_t - Particle[i].vz * dVolume / (3.0 * System.volume);
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }
    printf("Target: %lf OldKE: %lf NewKE: %lf NH: %lf\n",(1.50 * System.num_particles * System.kb * System.temp), System.kinetic_energy, newke, System.nose_hoover );
    printf("used KE %lf and Target %lf\n", System.kinetic_energy, (1.5 * System.kb * System.num_particles * System.temp) );

    //integrate volume
    System.volume += dVolume * System.delta_t;
    System.box_x *= dLength;
    System.box_y *= dLength;
    System.box_z *= dLength;

    //integrate pressure
    System.npt_pressure += dPressure * System.delta_t;
}

void    wrap_PBC(void)
{
    int i;

    for(i = 0; i < System.num_particles; i++)
    {
        //printf("%lf \t %lf \t %lf\n", Particle[i].x, Particle[i].y, Particle[i].z);

        if(Particle[i].x < 0.0)
            {
                while(Particle[i].x < 0.0)
                {
                    Particle[i].x += System.box_x;
                    //printf("x+: box is %lf now coord %lf\n", System.box_x, Particle[i].x);
                }
            }
        if(Particle[i].y < 0.0)
            {
                while(Particle[i].y < 0.0)
                {
                    Particle[i].y += System.box_y;
                    //printf("y+: box is %lf now coord %lf\n", System.box_y, Particle[i].y);
                }
            }
        if(Particle[i].z < 0.0)
            {
                while(Particle[i].z < 0.0)
                {
                    Particle[i].z += System.box_z;
                    //printf("z+: box is %lf now coord %lf\n", System.box_z, Particle[i].z);
                }
            }

        if(Particle[i].x > System.box_x)
            {
                while(Particle[i].x > System.box_x)
                {
                    Particle[i].x -= System.box_x;
                    //printf("x-: box is %lf now coord %lf\n", System.box_z, Particle[i].x);
                }
            }
        if(Particle[i].y > System.box_y)
            {
                while(Particle[i].y > System.box_y)
                {
                    Particle[i].y -= System.box_y;
                    //printf("y-: box is %lf now coord %lf\n", System.box_y, Particle[i].y);
                }
            }
        if(Particle[i].z > System.box_z)
            {
                while(Particle[i].z > System.box_z)
                {
                    Particle[i].z -= System.box_z;
                    //printf("z-: box is %lf now coord %lf\n", System.box_z, Particle[i].z);
                }
            }
    }
}

void do_NPT_AT_method(void)
{
    int i;
    double newke = 0.0;
    if(System.counter == 0)
    {
        System.nose_hoover = 3.0 * System.num_particles * System.kb * System.temp;
        System.nose_velocity = 0.0;
        System.pressure_v = 0.0;
        System.npt_mass = System.volume;
        System.nose_Q = System.nose_hoover;
        System.pressure_s = 0.0;
    }


    //determine derivatives
    //integrate pressure values (v and s)
    System.pressure_v += - (System.actual_pressure - System.pressure - System.pressure_v * System.nose_velocity / System.nose_hoover) * System.delta_t;
    printf("actual: %lf, pressure: %lf, pv: %lf, nose_v: %lf, nose_hoover: %lf\n", System.actual_pressure, System.pressure, System.pressure_v, System.nose_velocity, System.nose_hoover);
    printf("pressure_s: %lf, nose_Q: %lf, volume: %lf\n", System.pressure_s, System.nose_Q, System.volume);
    System.pressure_s += (2.0 * System.kinetic_energy  + System.pressure_v * System.pressure_v / System.nose_Q - 3.0 * System.num_particles * System.kb * System.temp) * System.delta_t;


    double dVolume = System.pressure_v / System.npt_mass;
    System.nose_velocity = System.nose_hoover * System.pressure_s / System.nose_Q;

    //integrate nose_hoover
    System.nose_hoover *= (1.0 + System.nose_velocity * System.delta_t);

    //integrate volume and change box sizes
    System.volume = System.volume * (1.0 + dVolume);
    double dx = pow((1.0 + dVolume), (1.0/3.0) );
    System.box_x *= dx;
    System.box_y *= dx;
    System.box_z *= dx;


    //integrate position and velocity
    for(i=0; i < System.num_particles; i++)
    {
        //position
        Particle[i].x = Particle[i].vx * System.delta_t + Particle[i].x * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y = Particle[i].vy * System.delta_t + Particle[i].y * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   //
        Particle[i].z = Particle[i].vz * System.delta_t + Particle[i].z * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   //

        //velocity
        Particle[i].vx = Particle[i].fx / Particle[i].mass * System.delta_t + Particle[i].vx * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover -  dVolume  * System.delta_t / (3.0 * System.volume) );
        Particle[i].vy = Particle[i].fy / Particle[i].mass * System.delta_t + Particle[i].vy * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover -  dVolume  * System.delta_t / (3.0 * System.volume) );
        Particle[i].vz = Particle[i].fz / Particle[i].mass * System.delta_t + Particle[i].vz * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover -  dVolume  * System.delta_t / (3.0 * System.volume) );
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

}

void do_NPT_AT_method_Temp(void)
{

int i;
    double newke = 0.0;
    if(System.counter == 0)
    {
        System.nose_hoover = 3.0 * System.num_particles * System.kb * System.temp;
        System.nose_velocity = 0.0;
        System.pressure_v = 0.0;
        System.npt_mass = System.volume;
        System.nose_Q = System.nose_hoover;
        System.pressure_s = 0.0;
    }


    //determine derivatives
    //integrate pressure values (v and s)
    System.pressure_s += (2.0 * System.kinetic_energy  + System.pressure_v * System.pressure_v / System.nose_Q - 3.0 * System.num_particles * System.kb * System.temp) * System.delta_t;
    System.nose_velocity = System.nose_hoover * System.pressure_s / System.nose_Q;

    //integrate position and velocity
    for(i=0; i < System.num_particles; i++)
    {
        //position
        Particle[i].x = Particle[i].vx * System.delta_t + Particle[i].x ;   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y = Particle[i].vy * System.delta_t + Particle[i].y;   //
        Particle[i].z = Particle[i].vz * System.delta_t + Particle[i].z;   //

        //velocity
        Particle[i].vx = Particle[i].fx / Particle[i].mass * System.delta_t + Particle[i].vx * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover);
        Particle[i].vy = Particle[i].fy / Particle[i].mass * System.delta_t + Particle[i].vy * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover);
        Particle[i].vz = Particle[i].fz / Particle[i].mass * System.delta_t + Particle[i].vz * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover);
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    //integrate nose_hoover
    System.nose_hoover *= (1.0 + System.nose_velocity * System.delta_t);
}


void do_NPT_AT_method_Pressure(void)
{
    int i;
    double newke;

    System.pressure_v += - (System.actual_pressure - System.pressure - System.pressure_v * System.nose_velocity / System.nose_hoover) * System.delta_t;
    double dVolume = System.pressure_v / System.npt_mass;

    //integrate volume and change box sizes
    System.volume = System.volume * (1.0 + dVolume);
    double dx = pow((1.0 + dVolume), (1.0/3.0) );
    System.box_x *= dx;
    System.box_y *= dx;
    System.box_z *= dx;

    printf("actual: %lf, pressure: %lf, pv: %lf, nose_v: %lf, nose_hoover: %lf\n", System.actual_pressure, System.pressure, System.pressure_v, System.nose_velocity, System.nose_hoover);
    printf("pressure_s: %lf, nose_Q: %lf, volume: %lf\n", System.pressure_s, System.nose_Q, System.volume);

    //integrate position and velocity
    for(i=0; i < System.num_particles; i++)
    {
        //position
        Particle[i].x = Particle[i].vx * System.delta_t + Particle[i].x * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y = Particle[i].vy * System.delta_t + Particle[i].y * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   //
        Particle[i].z = Particle[i].vz * System.delta_t + Particle[i].z * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   //

        //velocity
        Particle[i].vx = Particle[i].fx / Particle[i].mass * System.delta_t + Particle[i].vx * (1.0 -  dVolume  * System.delta_t / (3.0 * System.volume) );
        Particle[i].vy = Particle[i].fy / Particle[i].mass * System.delta_t + Particle[i].vy * (1.0 -  dVolume  * System.delta_t / (3.0 * System.volume) );
        Particle[i].vz = Particle[i].fz / Particle[i].mass * System.delta_t + Particle[i].vz * (1.0 -  dVolume  * System.delta_t / (3.0 * System.volume) );
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    printf("at %d, dV: %lf \n", System.counter, dVolume);
    printf("pv: %lf, npt_mass: %lf \n", System.pressure_v, System.npt_mass);


}

void do_NPT_berendsen_method_Temp(void)
{

int i;
    double newke = 0.0;
    if(System.counter == 0)
    {
        System.nose_hoover = 3.0 * System.num_particles * System.kb * System.temp * 0.95;
        System.nose_velocity = 0.0;
        System.pressure_v = 0.0;
        System.npt_mass = System.pressure;
        System.nose_Q = System.pressure * System.pressure * 0.95;
        System.pressure_s = 0.0;
    }


    //determine derivatives
    //integrate pressure values (v and s)
    System.pressure_s = (2.0 * System.kinetic_energy  + System.pressure_v * System.pressure_v / System.nose_Q - 3.0 * System.num_particles * System.kb * System.temp) * 2.0 * System.delta_t;
    System.nose_velocity = System.nose_hoover * System.pressure_s / System.nose_Q;

    //integrate nose_hoover
    System.nose_hoover += System.nose_velocity * 2.0 * System.delta_t;

    //integrate position and velocity
    for(i=0; i < System.num_particles; i++)
    {
        //position
        Particle[i].x = Particle[i].vx * System.delta_t + Particle[i].x ;   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y = Particle[i].vy * System.delta_t + Particle[i].y;   //
        Particle[i].z = Particle[i].vz * System.delta_t + Particle[i].z;   //

        //velocity
        Particle[i].vx = Particle[i].fx / Particle[i].mass * System.delta_t + Particle[i].vx * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover);
        Particle[i].vy = Particle[i].fy / Particle[i].mass * System.delta_t + Particle[i].vy * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover);
        Particle[i].vz = Particle[i].fz / Particle[i].mass * System.delta_t + Particle[i].vz * (1.0 - System.nose_velocity * System.delta_t / System.nose_hoover);
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

}

void do_NPT_berendsen_method_Pressure(void)
{
    int i;
    double newke;

    System.pressure_v = (System.actual_pressure - System.pressure) * 2.0 * System.delta_t;
    double dVolume = System.pressure_v / System.npt_mass;

    //integrate volume and change box sizes
    System.volume = System.volume * (1.0 + dVolume);
    double dx = pow((1.0 + dVolume), (1.0/3.0) );
    System.box_x *= dx;
    System.box_y *= dx;
    System.box_z *= dx;

    printf("actual: %lf, pressure: %lf, pv: %lf, nose_v: %lf, nose_hoover: %lf\n", System.actual_pressure, System.pressure, System.pressure_v, System.nose_velocity, System.nose_hoover);
    printf("pressure_s: %lf, nose_Q: %lf, volume: %lf\n", System.pressure_s, System.nose_Q, System.volume);

    //integrate position and velocity
    for(i=0; i < System.num_particles; i++)
    {
        //position
        Particle[i].x = Particle[i].vx * System.delta_t + Particle[i].x * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   // or Particle[i].x * delta_t * dLength  / 3.0 ???
        Particle[i].y = Particle[i].vy * System.delta_t + Particle[i].y * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   //
        Particle[i].z = Particle[i].vz * System.delta_t + Particle[i].z * (1.0 + System.delta_t * dVolume / (3.0 * System.volume) );   //

        //velocity
        Particle[i].vx = Particle[i].fx / Particle[i].mass * System.delta_t + Particle[i].vx * (1.0 -  dVolume  * System.delta_t / (3.0 * System.volume) );
        Particle[i].vy = Particle[i].fy / Particle[i].mass * System.delta_t + Particle[i].vy * (1.0 -  dVolume  * System.delta_t / (3.0 * System.volume) );
        Particle[i].vz = Particle[i].fz / Particle[i].mass * System.delta_t + Particle[i].vz * (1.0 -  dVolume  * System.delta_t / (3.0 * System.volume) );
        newke += 0.5 * Particle[i].mass * (Particle[i].vx * Particle[i].vx + Particle[i].vy * Particle[i].vy + Particle[i].vz * Particle[i].vz);
    }

    //printf("at %d, dV: %lf \n", System.counter, dVolume);
    //printf("pv: %lf, npt_mass: %lf \n", System.pressure_v, System.npt_mass);


}
