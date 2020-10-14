#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_SECURE_NO_WARNINGS  1
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void coleso_add_function(char* FUNCNAME, int* ID);
void coleso_set_parameter(char* PARAM, char* VALUE);
void coleso_read_set(int ID);
void coleso_read_file(int ID, char* FILENAME);
void coleso_init(int ID);
void coleso_pointvalue(int ID, double T, double* C, double* V);

void solution_using_prescribing_parameters() {
    int ID; // handle for the solution
    double t; // time
    double c[3]; // coordinates
    double V[1]; // solution (for this case the solution is an 1-component function)
    FILE* F;

    // Creating an object for the solution
    coleso_add_function("4peak", &ID);

    // Setting solution parameters
    // Xmin and Xmax prescribes a linear coordinate mapping
    coleso_set_parameter("Xmin", "0.0");
    coleso_set_parameter("Xmax", "1.0");
    // Set additional space of size 0.5 between pulses
    coleso_set_parameter("Period", "1.5");

    // Applying parameters to the solution object
    coleso_read_set(ID);

    // Initializing of the solution
    // (for this particular case this call has no effect)
    coleso_init(ID);

    // Calculating the solution and printing it to test1.dat
    F = fopen("test1.dat", "wt");
    if(F==NULL) { printf("Error writing to test1.dat\n"); return; }
    t = 0.4;
    c[0] = c[1] = c[2] = 0.0;
    for(c[0]=0.0; c[0] <= 3.0; c[0] += 0.0078125) {
        coleso_pointvalue(ID, t, c, V);
        fprintf(F, "%e %e\n", c[0], V[0]);
    }
    fclose(F);
}

void solution_using_parameters_in_file() {
    int ID; // handle for the solution
    double t; // time
    double c[3]; // coordinates
    double V[1]; // solution (for this case the solution is an 1-component function)
    FILE* F; // file descriptor

    // First write the parameters to a file
    F = fopen("test.txt", "wt");
    if(F==NULL) { printf("Error writing to test.txt\n"); return; }
    // Xmin and Xmax prescribes a linear coordinate mapping
    fprintf(F, "Xmin 0.0\n");
    fprintf(F, "Xmax 1.0\n");
    // Set additional space of size 0.5 between pulses
    fprintf(F, "Period 1.5\n");
    fclose(F);

    // Creating an object for the solution
    coleso_add_function("4peak", &ID);

    // Reading solution parameters from a file
    coleso_read_file(ID, "test.txt");

    // Initializing of the solution
    // (for this particular case this call has no effect)
    coleso_init(ID);

    // Calculating the solution and printing it to test2.dat
    F = fopen("test2.dat", "wt");
    if(F==NULL) { printf("Error writing to test2.dat\n"); return; }
    t = 0.4;
    c[0] = c[1] = c[2] = 0.0;
    for(c[0]=0.0; c[0] <= 3.0; c[0] += 0.0078125) {
        coleso_pointvalue(ID, t, c, V);
        fprintf(F, "%e %e\n", c[0], V[0]);
    }
    fclose(F);
}

int main(int argc, char** argv) {
    solution_using_prescribing_parameters();
    solution_using_parameters_in_file();
}
