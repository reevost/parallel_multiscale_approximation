#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <cblas.h>

unsigned int prod (unsigned int end, unsigned int begin){
    if (begin == end){
        return 1;
    }
    // here I made the splitting to avoid problems if mistakenly the number are exchanged.
    if (begin < end){
        unsigned int cumulative_prod = 1;
        for (unsigned int prod_index = end; prod_index > begin; prod_index--){
            cumulative_prod *= prod_index;
        }
        return cumulative_prod;
    }
    else {
        unsigned int cumulative_prod = 1;
        for (unsigned int prod_index = begin; prod_index > end; prod_index--) {
            cumulative_prod *= prod_index;
        }
        return cumulative_prod;
    }
}

double square_bracket_operator(int alpha, int grade) {
    if (grade == -1){
        return 1.0/(alpha+1);
    }
    if (grade == 0){
        return 1;
    }
    else{
        return prod(alpha, alpha-grade);
    }
}

double curly_bracket_operator(int nu, int grade){
    if (grade == 0){
        return 1;
    }
    if (nu == 0){
        //error handling
        return 0;
    }
    else{
        return prod(nu+grade-1, nu-1);
    }

}

#pragma clang diagnostic push // to avoid errors notifications
#pragma ide diagnostic ignored "misc-no-recursion"
double beta_for_wendland_function(int j, int k_plus_one, int nu){
    if (j == 0 && k_plus_one == 0){
        return 1;
    }
    else{
        double coefficient_sum = 0;
        int k = k_plus_one -1;
        if (j-1 < 0){
            for (int n = 0; n < k_plus_one ; n++) {
                coefficient_sum += beta_for_wendland_function(n, k, nu) * square_bracket_operator(n+1, n-j+1)/
                                   curly_bracket_operator(nu+2*k-n+1, n-j+2);
            }
            return coefficient_sum;
        }
        else{
            for (int n = j-1; n < k_plus_one ; n++) {
                coefficient_sum += beta_for_wendland_function(n, k, nu) * square_bracket_operator(n+1, n-j+1)/
                                   curly_bracket_operator(nu+2*k-n+1, n-j+2);
            }
            return coefficient_sum;
        }
    }
}
#pragma clang diagnostic pop

double wendland_function(double r, double k, int d) {
    // r "points" where the function is evaluated. r is supposed to be ||x-y|| i.e. the norm of some point difference, so r in R^+.
    // k is degree of the function, who is C^2k.
    // d is the dimension of the embedding space.
    if (k == (int) k) {
        int nu = d/2 + k + 1; // optimal value for the purpose of the class function. NOLINT(cppcoreguidelines-narrowing-conversions,bugprone-integer-division)
        double progress_evaluation = 0;
        for (int n = 0; n < k + 1; n++) {
            progress_evaluation += beta_for_wendland_function(n, k, nu) * pow(r, n) * pow(fmax(1 - r, 0), (nu + 2 * k + 1)); // NOLINT(cppcoreguidelines-narrowing-conversions,bugprone-integer-division)
        }
        return progress_evaluation;
    }
    else {
        // we expect even dimension here
        double estimate_integration_value;
        int nu = d/2 + k + 1.0/2; //NOLINT(cppcoreguidelines-narrowing-conversions,bugprone-integer-division)
        if (1 > r && r >= 0) {
            double result = 0; // integrate.quad(lambda t: t * (1 - t) * *nu * (t * *2 - r * *2) * *(k - 1) / (special.gamma(k) * 2 * *(k - 1)), r, 1)
            estimate_integration_value = result;
        }
        else {
            estimate_integration_value = 0;
        }
        return estimate_integration_value;
    }
}

double separation_distance(double ** points, unsigned int dim_, unsigned int number_of_points_){
    // printf("separation distance routine with:\n points dim: %d, points numbers: %d\n", dim_, number_of_points_);
    double minimum = 1;
    for (int i=0; i<number_of_points_-1; i++){
        for (int j=i+1; j<number_of_points_; j++) {
            // printf("%d,%d\n",i,j);
            double temp[dim_];
            for (int d=0; d<dim_; d++){
                temp[d] = points[i][d] - points[j][d];
            }
            // printf("[%f, %f]\n", temp[0], temp[1]);
            double p_dist = cblas_dnrm2((int) dim_, temp, 1);
            // printf("%f\n",p_dist);
            if (p_dist<minimum){
                minimum = p_dist;
            }
        }
    }
    return minimum;
}

double fill_distance(double ** points, unsigned int dim_, unsigned int number_of_points_, double threshold){
    double min_vector[dim_], max_vector[dim_];
    int partial_temp_points[dim_+1];
    int total_temp_points = 1;
    partial_temp_points[0] = 1;
    for (int d=0; d<dim_; d++) {
        min_vector[d] = 0;
        max_vector[d] = 1;
        /* if not in the standard settings [0,1]^d
        min_vector[d] = points[0][d];
        max_vector[d] = points[0][d];
        for (int i=0; i<number_of_points_; i++){
            if (points[i][d] < min_vector[d]){
                min_vector[d]= points[i][d];
            }
            else if (points[i][d] > max_vector[d]){
                max_vector[d] = points[i][d];
            }
        }*/
        total_temp_points *= floor((max_vector[d]-min_vector[d])/threshold+1);
        partial_temp_points[d+1] = total_temp_points;
    }
    double temp_grid_point[dim_];
    double maximum = 0;
    for (int temp_i=0; temp_i < total_temp_points; temp_i++){
        // test the grid point generator: printf("[ "); for (int d=0; d<dim_; d++) {printf("%f ", temp_grid_point[d]);} printf("]\n");
        double minimum = 1; // should be changed to the maximum possible value
        for (int p_ind = 0; p_ind < number_of_points_; p_ind++){
            double temp_array[dim_];
            for (int p_dim = 0; p_dim < dim_; p_dim++){
                temp_array[p_dim] = temp_grid_point[p_dim] - points[p_ind][p_dim];
            }
            double p_dist = cblas_dnrm2((int) dim_, temp_array, 1);
            if (p_dist < minimum){
                minimum = p_dist;
            }
        }
        // increase the grid by one step
        temp_grid_point[0] += threshold;
        for (int d=1; d<dim_; d++) {
            if ((temp_i+1)%partial_temp_points[d]==0){
                temp_grid_point[d] += threshold;
                if (d != 0){
                    temp_grid_point[d-1] = min_vector[d-1];
                }
            }
        } // end increase

        if (minimum > maximum){
            maximum = minimum;
        }
    }
    return maximum;
}

int main(){
    clock_t start_time = clock();
    FILE *file_pointer;
    file_pointer = fopen("/app/home/lotf/CLionProjects/parallel_multiscale_approximation/data_2D.csv", "r");
    char data_line[100]; // enough to have all read digits in memory

    // find out the dimension
    fgets(data_line, sizeof data_line, file_pointer);
    char * token;
    int temp_dim = 0;
    token = strtok(data_line, ",");
    while (token != NULL)
    {
        temp_dim ++;
        token = strtok(NULL, ",");
    }
    const unsigned int dim = temp_dim;
    printf("dim: %d\n", dim);

    int temp_number_of_points =1;
    while (fgets(data_line, sizeof data_line, file_pointer)){
        temp_number_of_points++;
    }
    const unsigned int total_number_of_points = temp_number_of_points;
    printf("number of points: %d\n", total_number_of_points);

    // back flow to avoid data loss is needed!
    fseek(file_pointer, 0, SEEK_SET);

    // allocate the memory for the matrix of points.
    int index = 0, temp_index;
    double **data = (double **) malloc(total_number_of_points * sizeof (double *));
    for (temp_index = 0; temp_index < total_number_of_points; temp_index ++)
    {
        data[temp_index] = (double *) malloc(dim * sizeof (double));
    }

    // read values from file and store them in the matrix defined above.
    while (fgets(data_line, sizeof data_line, file_pointer))
    {
        token = strtok(data_line, ",");
        temp_dim = 0;
        while (token != NULL)
        {
            data[index][temp_dim] = strtod(token, NULL);
            // printf("dd %f, dim %d, index %d\n", data[index][temp_dim], temp_dim, index);
            token = strtok(NULL, ",");
            temp_dim ++;
        }
        index ++;

    }
    fclose(file_pointer);
    clock_t read_time = clock()-start_time;
    printf("Data read in %f sec\n", (double) read_time/CLOCKS_PER_SEC);


    // here I have successfully imported the values (double format) from the file into the 2-pointer data!
    double sep = separation_distance(data, dim, total_number_of_points);
    clock_t separation_time = clock()- read_time;
    printf("Separation distance: %f, evaluated in %f sec.\n", sep, (double) separation_time/CLOCKS_PER_SEC);

    double fill = fill_distance(data, dim, total_number_of_points, 0.01);
    clock_t fill_time = clock()- separation_time;
    printf("Fill distance: %f, evaluated in %f sec.\n", fill, (double) fill_time/CLOCKS_PER_SEC);




    free(data);
    return 0;
}