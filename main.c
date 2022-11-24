# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>
# include <math.h>
# include "multiscale_structure.h"
# include "wendland_functions.h"

int main(){
    clock_t start_time = clock();
    char problem_setting[40] = "2D";
    FILE *file_pointer;
    char f_path[100];
    snprintf(f_path, 100, "/app/home/lotf/CLionProjects/parallel_multiscale_approximation/evaluated_data_%s.csv", problem_setting);
    file_pointer = fopen(f_path, "r");
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
    const int dim = temp_dim-1;
    printf("dim: %d\n", dim); // here we assume function of with real values, i.e. 1-dimensional.
    // find out the number of points
    int temp_number_of_points =1;
    while (fgets(data_line, sizeof data_line, file_pointer)){
        temp_number_of_points++;
    }
    const unsigned int total_number_of_points = temp_number_of_points;
    printf("number of points: %d\n", total_number_of_points);
    // size of the point matrix found

    // back flow to avoid data loss is needed!
    fseek(file_pointer, 0, SEEK_SET);

    // allocate the memory for the matrix of points.
    int index = 0, temp_index;
    double **data = (double **) malloc(total_number_of_points * sizeof (double *));
    for (temp_index = 0; temp_index < total_number_of_points; temp_index ++)
    {
        data[temp_index] = (double *) malloc((dim+1) * sizeof (double));
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
    double sep_all = separation_distance(data, dim, total_number_of_points);
    clock_t separation_time = clock()- read_time;
    printf("Separation distance: %f, evaluated in %f sec.\n", sep_all, (double) separation_time/CLOCKS_PER_SEC);

    double fill_all = fill_distance(data, dim, total_number_of_points, 0.001);
    clock_t fill_time = clock()- separation_time;
    printf("Fill distance: %f, evaluated in %f sec.\n", fill_all, (double) fill_time/CLOCKS_PER_SEC);

    const int number_of_levels = 4; int points_on_level[number_of_levels]; double mu = 0.5;
    double *** data_nested_structure;
    data_nested_structure = data_multilevel_structure(data, total_number_of_points, dim, number_of_levels, points_on_level, mu, 0.25,
                                                      true);

    for (int i=0;i<number_of_levels;i++){
        printf("\nlevel %d.\n", i+1);
        double sep = separation_distance(data_nested_structure[i], dim, points_on_level[i]);
        separation_time = clock()- fill_time;
        printf("Separation distance: %f, evaluated in %f sec.\n", sep, (double) separation_time/CLOCKS_PER_SEC);

        double fill = fill_distance(data_nested_structure[i], dim, points_on_level[i], 0.001);
        fill_time = clock()- separation_time;
        printf("Fill distance: %f, evaluated in %f sec.\n", fill, (double) fill_time/CLOCKS_PER_SEC);
    }
    /// Function that makes the interpolation using the multiscale approximation technique (see holger paper) ==> returns a 2-pointer, pointer[i] points to the array of alphas for level i.
    // PARAMETERS ========================================================================================================================
    // solving technique: "cg", "gmres" or "single_iteration". The first one solve with conjugate gradient iteratively the linear systems Ax=b level by level,
    // the second solve A x = b using gmres, while in both A = the big triangular matrix.
    double eps = 0.0001; // convergence threshold;
    int wendland_coefficients[2] = {1,3}; //A couple of integers (n, k) associated with the wanted Wendland compactly supported radial basis function phi_(n, k).
    double nu=1; // The coefficient that define the radius (delta) of the compactly supported rbf, delta_j = nu * mu_j, where mu_j is the mesh norm at the level j.
    int matrix_cumulative_size[number_of_levels+1]; matrix_cumulative_size[0]=0;
    for(int j=0; j<number_of_levels;j++){matrix_cumulative_size[j+1]=matrix_cumulative_size[j]+points_on_level[j];} // dimension of the big matrix T_n
    // At a certain point the controls for the parameters types and shapes should be included.
    double solution_vector[matrix_cumulative_size[number_of_levels]];

    char * solving_technique = "cg";
    // char * solving_technique = "gmres"
    if (solving_technique[0] == 'c') {
        printf("\nStarting solving process with conjugate gradient\n");
        for (int level = 0; level < number_of_levels; level++) {
            clock_t cg_time = clock();
            printf("\nStep %d out of %d\n", level+1, number_of_levels);
            const int points_on_this_level = points_on_level[level]; // since I used it a lot, I decide to "waste" some space do have a quicker access.
            /// initialize the parameters
            double delta_j = nu * fill_distance(data_nested_structure[level], dim, points_on_this_level, 0.001);
            unsigned int cg_j = 0; // cg iterations
            // initializing the list of solution points.
            double alpha_j_list[points_on_this_level + 1][points_on_this_level];
            for (int i=0; i < points_on_this_level; i++) { alpha_j_list[0][i] = 0.5; } // chose starting point
            // evaluate and save A as symmetric matrix
            double A_upper_packed[points_on_this_level * (points_on_this_level + 1) / 2];
            // I have to create a vector to store the difference of x_i-x_j, to use it later to compute r:= ||(x_i-x_j)/delta_j||_2
            double temp_vector[dim]; // I will also use it later to update the rhs, i.e. rhs = f_j-sum_{i=0}^{cg_j-1} B_ij@alpha_i
            for (int row = 0; row < points_on_this_level; row++) {
                for (int col = row; col < points_on_this_level; col++) {
                    // copy in the temp vector the point x_i
                    cblas_dcopy(dim, data_nested_structure[level][col], 1, temp_vector, 1);
                    // subtract the point x_j to temp vector
                    cblas_daxpy(dim, -1, data_nested_structure[level][row], 1, temp_vector, 1);
                    // store values on the matrix in packed form (ord=RowMajor)
                    A_upper_packed[row * (2*points_on_this_level-1-row)/2 + col] =
                            wendland_function(cblas_dnrm2(dim, temp_vector, 1) / delta_j, wendland_coefficients[0],
                                              wendland_coefficients[1]) / pow(delta_j, dim);
                }
            } //end storing packed matrix
            double residual_list[points_on_this_level+1][points_on_this_level]; // Residual. Since I do know a priori the maximum number of iteration just initialize 2-array.
            double direction_j[points_on_this_level]; // The direction along with I will apply cg. I do not want to take track of it during the algorithm, so I initialize it as a single array.
            // set the residual_0 = the rhs. See below.
            for (int i = 0; i < points_on_this_level; i++) {residual_list[0][i] = data_nested_structure[level][i][dim];}
            // Above we set the residual_0 = the rhs. In this way, the below call of dspmv for the creation of the residual (r_0 = b-Ax) will store the value of directly in r_0, avoiding the creation of a vector f_level just for few steps.
            printf("r_0 = "); for (int pp = 0; pp < points_on_this_level; ++pp) {printf("%f, ", residual_list[0][pp]);} printf("\n");

            /// update the rhs
            for (int k = 0; k < level; k++) {
                printf("sol_k = "); for (int pp = 0; pp < points_on_level[k]; ++pp) {printf("%f, ", solution_vector[matrix_cumulative_size[k]+pp]);} printf("\n");
                // generate B_kj
                double *B_jk = (double *) malloc(points_on_this_level * points_on_level[k] * sizeof(double));
                double delta_k = delta_j * pow(mu, k - level);
                for (int row = 0; row < points_on_this_level; row++) {
                    for (int col = row; col < points_on_level[k]; col++) {
                        // copy in the temp vector the point x_i^(level)
                        cblas_dcopy(dim, data_nested_structure[level][row], 1, temp_vector, 1);
                        // subtract the point x_j^(k) to temp vector
                        cblas_daxpy(dim, -1, data_nested_structure[k][col], 1, temp_vector, 1);
                        // store values on the matrix in packed form (ord=RowMajor)
                        B_jk[row * points_on_this_level + col] =
                                wendland_function(cblas_dnrm2(dim, temp_vector, 1) / delta_k, wendland_coefficients[0],
                                                  wendland_coefficients[1]) / pow(delta_k, dim);
                    }
                }
                printf("sol_k = "); for (int pp = 0; pp < points_on_level[k]; ++pp) {printf("%f, ", solution_vector[matrix_cumulative_size[k]+pp]);} printf("\n");

                // update the rhs: subtract to the rhs (now stored in r_0) B_jk@alpha_k
                cblas_dgemv(CblasRowMajor, CblasConjNoTrans, points_on_this_level, points_on_level[k], -1, B_jk,
                            points_on_this_level, &solution_vector[matrix_cumulative_size[k]], 1, 1,
                            residual_list[0], 1);
                printf("sol_k = "); for (int pp = 0; pp < points_on_level[k]; ++pp) {printf("%f, ", solution_vector[matrix_cumulative_size[k]+pp]);} printf("\n");

                free(B_jk);
            }
            // r_0 is the residual at j = 0, therefore: r_0 = rhs-Ax_0. since we stored the rhs in r_0, we just need to update it with dspmv.
            cblas_dspmv(CblasRowMajor, CblasUpper, points_on_this_level, -1, A_upper_packed, alpha_j_list[cg_j], 1, 1,
                        residual_list[0], 1); // b-Ax_0
            // the direction at j = 0 is defined as the starting residual, therefore we just copy the values.
            cblas_dcopy(points_on_this_level, residual_list[0], 1, direction_j, 1);
            printf("r_0 = dir_0 = "); for (int pp = 0; pp < points_on_this_level; ++pp) {printf("%f, ", direction_j[pp]);} printf("\n");
            /// Conjugate Gradient routine. Reference NLA Wendland 2017, algorithm 22, pag 195. all comments on the steps of the routine follows this algorithm notation.
            double t_j[points_on_this_level], a_j, b_j; // parameters of the cg routine, since at every cycle is overwritten better initialize it outside to avoid memory loss
            while (cblas_dnrm2(points_on_this_level, residual_list[cg_j], 1) > eps) {
                // \t_j = A@\p_j
                cblas_dspmv(CblasRowMajor, CblasUpper, points_on_this_level, 1, A_upper_packed, direction_j, 1, 0, t_j, 1);
                //a_j = \frac{||\r_j||_2^2}{<\t_j,\p_j>_2}
                a_j = pow(cblas_dnrm2(points_on_this_level, residual_list[cg_j], 1), 2) /
                      cblas_ddot(points_on_this_level, t_j, 1, direction_j, 1);
                // \x_{cg_j+1} = \x_{cg_j} + a_j \p_j;
                cblas_dcopy(points_on_this_level, alpha_j_list[cg_j], 1, alpha_j_list[cg_j + 1], 1);
                cblas_daxpy(points_on_this_level, a_j, direction_j, 1, alpha_j_list[cg_j + 1], 1);
                // \r_{cg_j+1} = \r_{cg_j} + a_j \p_j;
                cblas_dcopy(points_on_this_level, residual_list[cg_j], 1, residual_list[cg_j + 1], 1);
                cblas_daxpy(points_on_this_level, -a_j, t_j, 1, residual_list[cg_j + 1], 1);
                // b_j = \frac{||\r_{cg_j+1}||_2^2}{||\r_j||_2^2}
                b_j = pow(cblas_dnrm2(points_on_this_level, residual_list[cg_j + 1], 1) /
                          cblas_dnrm2(points_on_this_level, residual_list[cg_j], 1), 2);
                // \p_{cg_j+1} = \r_{cg_j+1} + b_j \p_j
                cblas_dscal(points_on_this_level, b_j, direction_j, 1);
                cblas_daxpy(points_on_this_level, 1, residual_list[cg_j + 1], 1, direction_j, 1);
                cg_j++;
                printf("cg iteration %d: \n2-norm of the residual: %f\n", cg_j,
                       cblas_dnrm2(points_on_this_level, residual_list[cg_j], 1));
                // PRINTS AND DEBUG OF CG
                /*printf("t_j = "); for (int pp = 0; pp < points_on_this_level; ++pp) {printf("%f, ", t_j[pp]);} printf("\n");
                printf("a_j = %f\n", a_j);
                printf("x_j+1 = "); for (int pp = 0; pp < points_on_this_level; ++pp) {printf("%f, ", alpha_j_list[cg_j][pp]);} printf("\n");
                printf("r_j+1 = "); for (int pp = 0; pp < points_on_this_level; ++pp) {printf("%f, ", residual_list[cg_j][pp]);} printf("\n");
                printf("b_j = %f\n", b_j);
                printf("dir_j+1 = "); for (int pp = 0; pp < points_on_this_level; ++pp) {printf("%f, ", direction_j[pp]);} printf("\n");*/
            }// endwhile
            // Now I have a list of solutions and residuals for every step. However, I simply store the result, i.e. alpha[cg_j] on the solution vector
            for (int i = 0; i < points_on_this_level; i++) {
                solution_vector[matrix_cumulative_size[level] +i] = alpha_j_list[cg_j][i];
            }
            clock_t cg_time_delta = clock()- cg_time;
            printf("Conjugate gradient terminated in %f sec.\n", (double) cg_time_delta/CLOCKS_PER_SEC);
        }
    }
    else if (solving_technique[0] == 'g'){
            printf("I will implement gmres non-sequential routine");
        }
    else if (solving_technique[0] == 'j') {
        printf("I will implement jacobi routine");
    }

    printf("\nNow is time to test the accuracy of the approximation and save the interpolant for every level\n");
    // s_j = sum_{i=1}^{N_j} alpha_i^(j) Phi_j(dot, \x_i^(j))

    // domain separation distance
    double threshold = 0.01;
    // If we assume that our domain is within [0,1]^d, then
    double min_vector[dim], max_vector[dim];
    int partial_temp_points[dim+1]; int total_temp_points = 1; partial_temp_points[0] = 1;
    for (int d=0; d<dim; d++) {
        min_vector[d] = 0; max_vector[d] = 1;
        total_temp_points *= floor((max_vector[d]-min_vector[d])/threshold+1);
        partial_temp_points[d+1] = total_temp_points;
    }
    double temp_value; double temp_vector[dim];

    for (int level = 0; level < number_of_levels; ++level) {
        clock_t save_time = clock();
        printf("Level %d", level+1);
        // create the file for data storage
        FILE *evaluated_file_pointer;
        char path[200];
        snprintf(path, 200, "/app/home/lotf/CLionProjects/parallel_multiscale_approximation/%s_domain/approximation_%d_evaluated_in_domain_2D_few.csv", problem_setting, level+1);
        evaluated_file_pointer = fopen(path, "w");

        double temp_grid_point[dim];
        for (int i = 0; i < dim; ++i) { temp_grid_point[i] = min_vector[i]; }
        for (int i = 0; i < pow(1.0 / threshold, dim); i++) {
            temp_value = 0;
            for (int d = 0; d < dim; ++d) {
                fprintf(evaluated_file_pointer, "%f,", temp_grid_point[d]);
            }
            double delta_j = nu * fill_distance(data_nested_structure[level], dim, points_on_level[level], 0.001);
            // evaluate the sum_{i=1}^{N_j} alpha_i^(j) Phi_j(domain_point, \x_i^(j)) for fixed j=level.
            for (int point_index = matrix_cumulative_size[level];
                 point_index < matrix_cumulative_size[level + 1]; ++point_index) {
                // copy in the temp vector the point x_i
                cblas_dcopy(dim, data_nested_structure[level][point_index - matrix_cumulative_size[level]], 1,
                            temp_vector, 1);
                // subtract the point x_j to temp vector
                cblas_daxpy(dim, -1, temp_grid_point, 1, temp_vector, 1);
                // store values on the matrix in packed form (ord=RowMajor)
                temp_value += solution_vector[point_index] *
                              wendland_function(cblas_dnrm2(dim, temp_vector, 1) / delta_j, wendland_coefficients[0],wendland_coefficients[1]) / pow(delta_j, dim);
            }
            fprintf(evaluated_file_pointer, "%f\n", temp_value);

            // increase the grid by one step
            temp_grid_point[0] += threshold;
            for (int d = 1; d < dim; d++) {
                if ((i + 1) % partial_temp_points[d] == 0) {
                    temp_grid_point[d] += threshold;
                    if (d != 0) {
                        temp_grid_point[d - 1] = min_vector[d - 1];
                    }
                }
            } // end increase
        }
        clock_t save_time_delta = clock()- save_time;
        printf(" saved in %f sec.\n", (double) save_time_delta/CLOCKS_PER_SEC);
    }
    free(data);
    return 0;
}