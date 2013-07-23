//POPULATE 1.6
//By Jacob Wagner
//Translated from FORTRAN 90 to C by Jacob Wagner
//
//Creates a .XYZ file for input to MD programs.
//Places a user defined number of molecules at given spacing with random orientations
//into a box, while aboiding collision/overlap complying with hard-body exclusion.

//Random number generator: Unifrom random number generator from http://www.cs.wm.edu/~va/software/park/ see rngs.c/rngs.h
//
// Date April 15, 2013

//includes
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include <time.h>

#include "rngs.h"
#include "rngs.c"


typedef struct
{
    char label[2];
    double sigma;
    double epsilon;
    double mass;
    double charge;
    double relx;
    double rely;
    double relz;
    double sigma2;

} Particles;

typedef struct
{
    Particles* atom;
    int num_atoms;
    int copies;
} Prototypes;

Prototypes* Prototype;

double absolute(double);

int main(int argc, char **argv)
{
    //declare variables
    char outname[] = "000";
    char output_file[7];
    int ntot, xmax, ymax, zmax, seed;
    int i, j, k, l, spot, prev, species, Nspecies, atom_loc;
    double xsize, ysize, zsize;
    double randx, randy, randz, rand_theta, rand_phi, rand_psi, zero;
    int overlap;
    double rotation[3][3];
    //double* atoms;
    double* locations;
    int* species_info;
    //char* species_label;
    int* locations_label;
    FILE* fr;
    char line[100];

    //initialize variable
    xsize = 0.0;
    ysize = 0.0;
    zsize = 0.0;
    ntot = 0;
    xmax = 0;
    ymax = 0;
    zmax = 0;
    seed = 0;
    i = 0;
    j = 0;
    k = 0;
    l = 0;
    spot = 0;
    prev = 0;
    species = 0;
    Nspecies = 0;
    atom_loc = 0;
    randx = 0.0;
    randy = 0.0;
    randz = 0.0;
    rand_theta = 0.0;
    rand_phi = 0.0;
    zero = 0.0;
    overlap = 0;


    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            rotation[i][j] = 0.0;
        }
    }

    //read input file
    fr = fopen("create.dat", "rt");

    //read label
    fgets(line,100,fr);//8
    sscanf(line, "%s", outname);
    //fscanf(fr, "%s", outname);
    printf("outname is %s\n", outname);

    fgets(line,100,fr);
    sscanf(line, "%lf %lf %lf", &xsize, &ysize, &zsize);
    printf("box size is %lf %lf %lf\n", xsize, ysize, zsize);
    fgets(line,100,fr);
    //molecule information
    fgets(line,100,fr);
    sscanf(line, "%d", &species);
    printf("species is %d \n", species);
    fgets(line,100,fr);
    //allocate space for prototype
    Prototype = malloc(species * sizeof(Prototypes));
    //species_info = malloc( species * 2 * sizeof(int));

    for(i = 0; i < species; i++)
    {
        fgets(line,100,fr);
        sscanf(line, "%d %d", &Prototype[i].num_atoms, &Prototype[i].copies);
        printf("Prototype %d: %d and %d\n", i, Prototype[i].num_atoms, Prototype[i].copies);
        ntot += Prototype[i].num_atoms * Prototype[i].copies;
        Nspecies += Prototype[i].num_atoms;
    }

    fgets(line,100,fr);

    //allocate prototype space
    for(i = 0; i < species; i++)
    {
            Prototype[i].atom = malloc(Prototype[i].num_atoms * sizeof(Particles));
    }

    //allocate species_specifications

    //atoms = malloc(Nspecies * 7 * sizeof(int));             //replace with Particle[i].atom[j]
    //species_label = malloc(Nspecies * 2 * sizeof(char));    //replace with Particle[i].atom[j]
    printf("ntot is %d\n", ntot);
    locations = malloc(ntot * 3 * sizeof(double));
    locations_label = malloc(ntot * 2 * sizeof(int)); //does this work

    k = 0;

    for( i = 0; i < species; i++)
    {
        for( j = 0; j < Prototype[i].num_atoms; j++) //change sites
        {
            fgets(line,100,fr);
            //printf("%s\n", line);
            sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf", Prototype[i].atom[j].label, &Prototype[i].atom[j].sigma, &Prototype[i].atom[j].epsilon, &Prototype[i].atom[j].mass, &Prototype[i].atom[j].charge, &Prototype[i].atom[j].relx, &Prototype[i].atom[j].rely, &Prototype[i].atom[j].relz);
            printf("ATOM %d: %s %lf %lf %lf %lf %lf %lf %lf\n", j, Prototype[i].atom[j].label, Prototype[i].atom[j].sigma, Prototype[i].atom[j].epsilon, Prototype[i].atom[j].mass, Prototype[i].atom[j].charge, Prototype[i].atom[j].relx, Prototype[i].atom[j].rely, Prototype[i].atom[j].relz);
            Prototype[i].atom[j].sigma2 = Prototype[i].atom[j].sigma / 2.0;
            //printf("%lf sigma2\n", Prototype[i].atom[j].sigma2);
            k++;
        }
        fgets(line,100,fr);
    }
    fclose(fr);

    //generate random loccation and orientation
    xmax = (int) (floor( xsize / 2.0));
    ymax = (int) (floor( ysize / 2.0));
    zmax = (int) (floor( zsize / 2.0));

    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            rotation[i][j] = 0.0;
        }
    }

    printf("start loop for species\n");
    for(i = 0; i < species; i++)
    {
        spot = Prototype[i].num_atoms;
        j = 0;
        while( j < Prototype[i].copies)
        {
            randx = (Random() - 0.50) * xsize;
            randy = (Random() - 0.50) * ysize;
            randz = (Random() - 0.50) * zsize;
            //printf("coor: %lf %lf %lf\n", randx, randy, randz);

            rand_theta = Random() * 2.0 * 3.14159;
            rand_phi = Random() * 1.0 * 3.14159;
            rand_psi = Random() * 1.0 * 3.14159;

            if(Prototype[i].num_atoms == 1) //no rotation
            {
                //place atoms on wall
                k=0;
                /*
                //call random numbers
                randx = (Random() - 0.50) * xsize;
                randy = (Random() - 0.50) * ysize;
                randz = (Random() - 0.50) * zsize;
                */
                locations[(prev + k)*3] = randx;
                locations[(prev + k)*3 + 1] = randy;
                locations[(prev + k)*3 + 2] = randz;
                locations_label[(prev + k)*2] = i;
                locations_label[(prev + k)*2 + 1] = spot - 1;

                printf("%lf %lf %lf for i=%d j=%d at k=%d and prev=%d\n", locations[(prev + k)*3], locations[(prev + k)*3 +1], locations[(prev + k)*3 + 2], i, j, k,prev);

                if(prev >= 1) //check for overlaps
                {
                    for(l = 0; l < prev; l++)
                    {
                        /*
                        printf("check %d vs %d interaction at %lf %lf %lf\n", l, k+prev, locations[l*3], locations[l*3 + 1], locations[l*3 + 2]);
                        printf("%lf vs\n", Prototype[i].atom[spot].sigma2);
                        printf("%lf x\n", absolute(locations[l*3] - locations[prev*3]));
                        printf("%lf y\n", absolute(locations[l*3 + 1] - locations[prev*3 + 1]));
                        printf("%lf z\n", absolute(locations[l*3 + 2] - locations[prev*3 + 2]));
                        */
                       if( (absolute(locations[l*3] - locations[prev*3]) <= Prototype[i].atom[spot].sigma2) && (absolute(locations[l*3 + 1] - locations[prev*3 + 1]) <= Prototype[i].atom[j].sigma2) && (absolute(locations[l*3 + 2] - locations[prev*3 + 2]) <= Prototype[i].atom[j].sigma2) )
                       {
                           overlap = 1;
                           printf("REJECT overlap\n");
                           break;
                       }
                       //printf("next\n");
                    }
                }
                //check if the particle is in the box
                //printf("in the box?\n");
                if( (absolute(locations[(prev + k)*3] > xmax)) || (absolute(locations[(prev + k)*3 + 1] > ymax)) || (absolute(locations[(prev + k)*3 + 2] > zmax )) )
                {
                    overlap = 1;
                    printf("REJECT box\n");
                }

                if(overlap == 1) //reset overlap and redo iteration
                {
                    overlap = 0;
                    printf("reject\n");
                }
                else
                {
                    //self-increment while loop
                    //printf("\naccepted: prev=%d, j=%d, and k=%d\n", prev, j, k);
                    j++;
                    prev++;
                }
                //printf("moving on\n");
            }
            else //molecular system
            {
                printf("multi-system %d\n", prev+k);
                /////

                k=0;
                while(k <= Prototype[i].copies)
                {
                    //calculate rotation matrix
                    rotation[0][0] =   cos(rand_theta) * cos(rand_phi);
                    rotation[0][1] = - sin(rand_theta) * cos(rand_psi) + sin(rand_theta) * sin(rand_theta) * sin(rand_psi);
                    rotation[0][2] =   sin(rand_theta) * cos(rand_phi) * sin(rand_phi) + cos(rand_psi);

                    rotation[1][0] =   cos(rand_theta) * sin(rand_psi);
                    rotation[1][1] =   cos(rand_phi)   * cos(rand_psi) + sin(rand_theta) * sin(rand_phi) * sin(rand_psi);
                    rotation[1][2] = - cos(rand_psi)   * sin(rand_phi) + sin(rand_theta) * cos(rand_phi) * sin(rand_psi);

                    rotation[2][0] = - sin(rand_theta);
                    rotation[2][1] =   cos(rand_theta) * sin(rand_phi);
                    rotation[2][2] =   cos(rand_theta) * cos(rand_phi);

                    //call random numbers
                    randx = (Random() - 0.50) * xsize;
                    randy = (Random() - 0.50) * ysize;
                    randz = (Random() - 0.50) * zsize;

                    //shift molecule coordinates
                    for(j = 0; j < Prototype[i].num_atoms; j++)
                    {
                        locations[(prev + k*Prototype[i].num_atoms + j)*3] = randx + Prototype[i].atom[j].relx * rotation[0][0] + Prototype[i].atom[j].rely * rotation[0][1] + Prototype[i].atom[j].relz * rotation[0][2];
                        locations[(prev + k*Prototype[i].num_atoms + j)*3 + 1] = randy + Prototype[i].atom[j].relx * rotation[1][0] + Prototype[i].atom[j].rely * rotation[1][1] + Prototype[i].atom[j].relz * rotation[1][2];
                        locations[(prev + k*Prototype[i].num_atoms + j)*3 + 2] = randz + Prototype[i].atom[j].relx * rotation[2][0] + Prototype[i].atom[j].rely * rotation[2][1] + Prototype[i].atom[j].relz * rotation[2][2];
                        //strcpy(&locations_label[prev + k*Prototype[i].num_atoms + j], Prototype[i].atom[j].label);
                    }

                    printf("%lf %lf %lf for i=%d j=%d at k=%d and prev=%d\n", locations[(prev + k)*3], locations[(prev + k)*3 +1], locations[(prev + k)*3 + 2], i, j, k,prev);

                    if(prev >= 1) //check for overlaps
                    {
                      for(j = 0; j < Prototype[i].num_atoms; j++)
                      {
                            for(l = 1; l < prev; l++)
                            {
                                if( ((abs(locations[l*3] - locations[(k+prev + j)*3])) <= Prototype[i].atom[j].sigma / 2.0) && ((abs(locations[l*3 + 1] - locations[(k+prev + j)*3 + 1])) <= Prototype[i].atom[j].sigma / 2.0) && ((abs(locations[l*3 + 2] - locations[(k+prev + j)*3 + 2])) <= Prototype[i].atom[j].sigma / 2.0) )
                                {
                                    overlap = 1;
                                }
                            }
                      }
                    }
                    //check if ne particle is in the box
                    if( (abs(locations[(prev + k)*3] > xmax)) || (abs(locations[(prev + k)*3 + 1] > ymax)) || (abs(locations[(prev + k)*3 + 2] > zmax )) )
                    {
                        overlap = 1;
                    }

                    if(overlap == 1) //reset overlap and redo iteration
                    {
                        overlap = 0;
                    }
                    else
                    {
                        //self-increment while loop
                        printf("accepted %d\n", k);
                        j++;
                        k++;
                    }
                }
                prev += k * Prototype[i].num_atoms;

                ////
                j++;
                k++;
                prev += spot;
            }
        }
        atom_loc += spot;
    }
    printf("output\n");
    //create output file -- XYZ file
    strcpy(output_file, outname);
    //strcat(output_file, outname);
    strcat(output_file, ".xyz");
    printf("file is %s and outname = %s\n", output_file, outname);

    fr = fopen(output_file, "w+");
    fprintf(fr,"%d\n",ntot);
    fprintf(fr,"\n");
    fprintf(fr,"%lf\t%lf\t%lf\n", xsize, ysize, zsize);
    fprintf(fr, "\n");
    printf("done with header; onto %d outputs\n", ntot);
    //printf("locations_label = %s", locations_label[i]);
    for(i = 0; i < ntot; i++)
    {
        fprintf(fr,"%s\t%lf\t%lf\t%lf\n", Prototype[locations_label[i*2]].atom[locations_label[i*2 + 1]].label, locations[i*3], locations[i*3 + 1], locations[i*3 + 2]);
    }
    fclose(fr);

    //create output file -- TOP file
    //strcpy(output_file, "");
    //strcpy(output_file, outname);
    output_file[4] = 't';
    output_file[5] = 'o';
    output_file[6] = 'p';
    //strcat(output_file, ".top");
    printf("output_file is %s \n", output_file);

    fr = fopen(output_file, "w+");

    fprintf(fr, "%d\n\n", species);
    for(i = 0; i < species; i++)
    {
        fprintf(fr, "%d\t%d\n", Prototype[i].num_atoms, Prototype[i].copies);
    }

    fprintf(fr,"\n");

    for(i = 0; i < species; i++)
    {
        for(j = 0; j < Prototype[i].num_atoms; j++)
        {
            fprintf(fr, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Prototype[i].atom[j].label, Prototype[i].atom[j].sigma, Prototype[i].atom[j].epsilon, Prototype[i].atom[j].relx, Prototype[i].atom[j].rely, Prototype[i].atom[j].relz, Prototype[i].atom[j].mass);
        }
    }

    fclose(fr);

    //free allocated memory
    for(i = 0; i < species; i++)
    {
        free(Prototype[i].atom);
    }
    free(Prototype);

    //free(species_info);
    //free(atoms);
    //free(species_label);
    free(locations);
    free(locations_label);

    printf("done");
    return 0;
}

double absolute(double arg)
{
  if(arg < 0)
  {
      return (-arg);
  }
  else
  {
    return arg;
  }
}
