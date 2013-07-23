//output.c

void interval_output(void)
{
    VMD_output();
}

void intialize_VMD_output(void)
{
    FILE* fr;

    fr = fopen("coord.xyz", "w+");

    fclose(fr);

}

void VMD_output(void)
{
    FILE* fr;
    int i;

    fr = fopen("coord.xyz", "a");
    fprintf(fr, "%d\n\n", System.num_particles);

    for(i = 0; i < System.num_particles; i++)
    {
        fprintf(fr, "%s %lf %lf %lf \n", Prototype[ Particle[i].style ].type, Particle[i].x, Particle[i].y, Particle[i].z);
    }
    //fprintf(fr, "\n");

    fclose(fr);
}

void initialize_energy_output(void)
{
    FILE* fr;

    fr = fopen("energy.txt", "w+");
    fprintf(fr, "#\t Total\t Kinetic\t Potential\t Enthalpy\n");
    fprintf(fr, "%d \t %lf \t %lf \t %lf\t %lf\n", System.counter, (System.kinetic_energy + System.potential_energy), System.kinetic_energy, System.potential_energy, (System.kinetic_energy + System.potential_energy + System.actual_pressure * System.volume));

    fclose(fr);

    fr = fopen("pressure.txt", "w+");
    fprintf(fr, "#\t pressure\t volume\n");
    fclose(fr);


    fr = fopen("temperature.txt", "w+");
    fprintf(fr, "#\t temp(K)\n");
    fclose(fr);

}

void energy_output(void)
{
    FILE* fr;

    fr = fopen("energy.txt", "a");
    fprintf(fr, "%d \t %lf \t %lf \t %lf \t %lf\n", System.counter, (System.kinetic_energy + System.potential_energy), System.kinetic_energy, System.potential_energy,  (System.kinetic_energy + System.potential_energy + System.actual_pressure * System.volume));

    fclose(fr);

    fr = fopen("pressure.txt", "a");
    fprintf(fr, "%d \t %lf \t %lf\n", System.counter, System.actual_pressure, System.volume);

    fclose(fr);



    fr = fopen("temperature.txt", "a");
    fprintf(fr, "%d \t %lf \n", System.counter, System.kinetic_energy * 2.0 / (3.0 * System.num_particles * System.kb));

    fclose(fr);
}
