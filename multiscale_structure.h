#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <cblas.h>

#ifndef PARALLEL_MULTISCALE_APPROXIMATION_MULTISCALE_STRUCTURE_H
#define PARALLEL_MULTISCALE_APPROXIMATION_MULTISCALE_STRUCTURE_H

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
    return minimum/2;
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
    double temp_grid_point[dim_]; for (int i = 0; i < dim_; ++i) {temp_grid_point[i]=min_vector[i];}
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
        if (minimum > maximum){
            maximum = minimum;
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


    }
    return maximum;
}

double *** data_multilevel_structure(double ** data,const unsigned int number_of_points,const unsigned int points_dim, const int number_of_levels, int * points_on_level, const double mu,const double starting_mesh_norm, bool nest_flag) {
    /// data is the set of data points from which we want to build our nested sequence of data sets.
    /// number_of_levels is the number of level on which we want to split our data.
    /// mu is a value in (0,1) which describe the relation between the mesh norm of two nested data sets: h(X_j+1) = mu * h(X_j)
    /// starting_mesh_norm is the mesh norm used to find the first set of data.
    /// nest_flag: for generating a nested sequence or not
    double mesh_norm = 2 * starting_mesh_norm;  // for the code purpose it is better to consider twice of the mesh norm
    double *** data_nest_list = (double ***) malloc(number_of_levels * sizeof (double **)); // initialize the empty array where we will put the nested sets
    /// the flag vector help us to take track of the selected point on each level and save them when we have the number of points on that level // IS NECESSARY?
    int data_flag[number_of_points]; for (int i = 1; i < number_of_points; i++) { data_flag[i] = 1; } // create a flag array to keep track of the added points (positive: available, zero: not available, negative:already chosen)
    // initial point.
    double *selected_point = data[0]; data_flag[0] = -(number_of_levels + 1);  // in this way they stay negative in the nest setting

    //double ** new_data = data_nest_list[0]; // future set with the actual mesh norm
    points_on_level[0] = 1; // number of points in the first set, i.e. at the beginning one, the starting point.

    // find the new set based on the starting mesh norm
    while (selected_point) {
        // filter the list of points removing the point that are too close with the selected point
        for (int i = 0; i < number_of_points; i++) {
            if (data_flag[i] > 0) {
                double sum = 0;
                //evaluate the distance between the selected point and the data[i]
                for (int d = 0; d < points_dim; d++) {
                    sum += pow(selected_point[d] - data[i][d], 2);
                }
                // compare the distance with the desired norm, to see if the point is too close from the original
                // instead of sqrt(sum(square)) > value i consider sum(square) > value**2
                if (sum < pow(mesh_norm, 2)) { data_flag[i] = 0;}
            }//endif
        } //endfor

        // now I have set to 0 all the flags of not-eligible points
        // now I have to find a new point between the eligible ones and repeat
        int tmp_index = 0;
        do {
            tmp_index++;
        } while (tmp_index < number_of_points && data_flag[tmp_index] < 1);
        if (tmp_index < number_of_points) { // an eligible point was found
            selected_point = data[tmp_index]; // assign it as new selected point
            data_flag[tmp_index] = -(number_of_levels +1); // in this way they stay negative in the nest setting i.e. they won't be selected anymore
            points_on_level[0]++;
        } else { break; } // no eligible point => end of the run for this level
    }//end first while => we successfully created the first layer
    printf("number of points in the set 1: %d\n", points_on_level[0]);

    // save this list of points on the heap. Specifically, save the address of this selected points.
    data_nest_list[0]  = (double **) malloc(points_on_level[0] * sizeof (double *)); // allocate memory for the number of points
    int temp_counter = 0;
    for (int i=0; i<number_of_points; i++){ // iterate over the points and store point whose flag is negative
        if(data_flag[i]<0){
            data_nest_list[0][temp_counter] = data[i]; // store point reference
            temp_counter++;
        }
    }

    for (int j = 1; j < number_of_levels; j++) {// for every other step
        mesh_norm *= mu; // update the mesh norm
        if (nest_flag) { // nested setting
            // reset the flags: 1 for the non-selected points in the previous round, negative all the others
            for (int i = 0; i < number_of_points; i++) {data_flag[i]++;}

            // filter the list of points zeroing the points that are too close with the point brought from the last set
            for (int i = 0; i < points_on_level[j - 1]; i++) {
                for (int k = 1; k < number_of_points; k++) {
                    if (data_flag[k] > 0) {
                        double sum = 0;
                        //evaluate the distance between the selected point and the data[k]
                        for (int d = 0; d < points_dim; d++) {
                            sum += pow(data_nest_list[j - 1][i][d] - data[k][d], 2);
                        }
                        // compare the distance with the desired norm, to see if the point is too close from the original
                        // instead of sqrt(sum(square)) > value k consider sum(square) > value**2
                        if (sum < pow(mesh_norm, 2)) {data_flag[k] = 0;}
                    }//endif
                } //endfor
            }//endfor
            points_on_level[j] = points_on_level[j-1]; // set the actual number of points on the new level equal to the previous level
        }
        else { // not nested setting
            for (int i = 0; i < number_of_points; i++) {data_flag[i] = 1;}
            points_on_level[j] = 0;
            // set all values to 1, i.e. all points can be selected again, and set the starting number of points to 0.
        }
        // endif-else nest_flag

        // select a new point to begin the new round - OPTION: add more randomness.
        int tmp_index = -1;
        do {
            tmp_index++;
        } while (tmp_index < number_of_points && data_flag[tmp_index] < 1);
        if (tmp_index < number_of_points) { // an eligible point was found
            selected_point = data[tmp_index]; // assign it as new selected point
            data_flag[tmp_index] = -(number_of_levels+1); // in this way they stay negative in the nest setting i.e. they won't be selected anymore
            points_on_level[j]++;
        } else {
            printf("full at level %d\n", j-1);
            break;
        } // if no eligible point is found as first point, then we are necessary in a nested setting, where all points are already assigned. Therefore, we notify this, and we stop the generation of the sets

        // find the new set based on the j-th mesh norm
        while (selected_point) {
            // filter the list of points removing the point that are too close with the selected point
            for (int i = 0; i < number_of_points; i++) {
                if (data_flag[i] > 0) {
                    double sum = 0;
                    //evaluate the distance between the selected point and the data[i]
                    for (int d = 0; d < points_dim; d++) {
                        sum += pow(selected_point[d] - data[i][d], 2);
                    }
                    // compare the distance with the desired norm, to see if the point is too close from the original
                    // instead of sqrt(sum(square)) > value i consider sum(square) > value**2
                    if (sum < pow(mesh_norm, 2)) { data_flag[i] = 0; }
                }//endif
            } //endfor

            // now I have set to 0 all the flags of not-eligible points
            // now I have to find a new point between the eligible ones and repeat
            tmp_index = 0;
            do {
                tmp_index++;
            } while (tmp_index < number_of_points && data_flag[tmp_index] < 1);
            if (tmp_index < number_of_points) { // an eligible point was found
                selected_point = data[tmp_index]; // assign it as new selected point
                data_flag[tmp_index] = -(number_of_levels +1); // in this way they stay negative in the nest setting i.e. they won't be selected anymore
                points_on_level[j]++;
            } else { break; } // no eligible point => end of the run for this level
        }//end first while => we successfully created the first layer
        printf("number of points in the set %d: %d\n",j+1, points_on_level[j]);

        // save this list of points on the heap. Specifically, save the address of this selected points.
        data_nest_list[j]  = (double **) malloc(points_on_level[j] * sizeof (double *)); // allocate memory for the number of points
        temp_counter = 0;
        for (int i=0; i<number_of_points; i++){ // iterate over the points and store point whose flag is negative
            if(data_flag[i]<0){
                data_nest_list[j][temp_counter] = data[i]; // store point reference
                temp_counter++;
            }
        }
    }
    return data_nest_list;
}

#endif //PARALLEL_MULTISCALE_APPROXIMATION_MULTISCALE_STRUCTURE_H
