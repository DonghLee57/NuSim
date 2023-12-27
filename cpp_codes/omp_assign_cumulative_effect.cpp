#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;

int main(){
    int nrow = 4;
    vector<int> vec = {1,2,3,4,5,6,7,8,9,10};
    vector<vector<int>> mat(nrow, vector<int>(vec.size(), 0));

    // Simple assignment
    cout << "Simple - assign 2D matrix" << endl;
    for (int i=0; i<nrow; i++){
        for (int& v : vec){ v -= nrow; }
        mat[i] = vec;
    }
    for (int i=0; i<nrow; i++){
        for (int j=0; j<vec.size(); j++){
            cout << mat[i][j] << " ";
	}
        cout << endl;
    }

    // OpenMP assignment
    cout << "OpenMP - assign 2D matrix" << endl;
    vec = {1,2,3,4,5,6,7,8,9,10};
    #pragma omp parallel for
    for (int i=0; i<nrow; i++){
        for (int& v : vec){ v -= nrow; }
        mat[i] = vec;
    }
    for (int i=0; i<nrow; i++){
        for (int j=0; j<vec.size(); j++){
            cout << mat[i][j] << " ";
	}
        cout << endl;
    }

    /*
    To achieve the same result in the OpenMP block as in the simple loop block,
    while avoiding race conditions, it's important to understand how the data is being modified and accessed by multiple threads.
    */
    cout << "OpenMP - assign 2D matrix changing algorithm to avoid the cumulative effect" << endl;
    vec = {1,2,3,4,5,6,7,8,9,10};
    #pragma omp parallel for
    for (int i=0; i<nrow; i++){
        for (int j=0; j< vec.size(); j++)
        mat[i][j] = vec[j] - nrow*(i+1);
    }
    for (int i=0; i<nrow; i++){
        for (int j=0; j<vec.size(); j++){
            cout << mat[i][j] << " ";
	}
        cout << endl;
    }

    return 0;
}
