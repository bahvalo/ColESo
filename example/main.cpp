#include "coleso.h"

void solution_using_prescribing_parameters() {
    // Creating an object for the solution
    s_4peak S;

    // Setting solution parameters
    // Xmin and Xmax prescribes a linear coordinate mapping
    S.Xmin = 0.0;
    S.Xmax = 1.0;
    // Set additional space of size 0.5 between pulses
    S.Period = 1.5;

    // Initializing of the solution
    // (for this particular case this call has no effect)
    S.Init();

    // Calculating the solution and printing it to test1.dat
    FILE* F = fopen("test1.dat", "wt");
    if(F==NULL) { printf("Error writing to test1.dat\n"); return; }
    double t = 0.4;
    double c[3] = {0.0, 0.0, 0.0};
    for(c[0]=0.0; c[0] <= 3.0; c[0] += 0.0078125) {
        double V[1]; // solution (for this case the solution is an 1-component function)
        S.PointValue(t, c, V);
        fprintf(F, "%e %e\n", c[0], V[0]);
    }
    fclose(F);
}

void solution_using_parameters_in_file() {
    // First write the parameters to a file
    FILE* F = fopen("test.txt", "wt");
    if(F==NULL) { printf("Error writing to test.txt\n"); return; }
    // Xmin and Xmax prescribes a linear coordinate mapping
    fprintf(F, "Xmin 0.0\n");
    fprintf(F, "Xmax 1.0\n");
    // Set additional space of size 0.5 between pulses
    fprintf(F, "Period 1.5\n");
    fclose(F);

    // Creating an object for the solution
    s_4peak S;

    // Reading solution parameters from the file
    S.ReadParamsFromFile("test.txt");

    // Initializing of the solution
    // (for this particular case this call has no effect)
    S.Init();

    // Calculating the solution and printing it to test1.dat
    F = fopen("test2.dat", "wt");
    if(F==NULL) { printf("Error writing to test1.dat\n"); return; }
    double t = 0.4;
    double c[3] = {0.0, 0.0, 0.0};
    for(c[0]=0.0; c[0] <= 3.0; c[0] += 0.0078125) {
        double V[1]; // solution (for this case the solution is an 1-component function)
        S.PointValue(t, c, V);
        fprintf(F, "%e %e\n", c[0], V[0]);
    }
    fclose(F);
}

int main(int, char**) {
    solution_using_prescribing_parameters();
    solution_using_parameters_in_file();
}
