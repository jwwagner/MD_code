//closing.c

void FreeMemory(void)
{
    long i;

    for(i = 0; i < System.total_cells; i++)
    {
        Neighbor[i].length = 0;

        if(System.counter != 0)
        {
            free(Neighbor[i].cell);
        }
    }

    free(Neighbor);
    //printf("free Neighbor\n");
    free(Prototype);
    //printf("free Prototype\n");
    free(Particle);
    //printf("free Particle\n");
    free(System.types);
    printf("free System.types\n");

    printf("end\n");
}
