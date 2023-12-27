#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;

int main(){
    int nrow = 4;
    int count;
    vector<int> vec = {1,2,3,4,5,6,7,8,9,10};
    vector<vector<int>> mat(nrow, vector<int>(vec.size(), 0));

    // Simple assignment
    cout << "Simple - assign 2D matrix" << endl;
    count = 0;
    for (int i=0; i<nrow; i++){
        for (int j=0; j<vec.size(); j++){
			      count = i*j;
            mat[i][j] = i + j;
			      mat[i][j] = count;
        }
    }
    for (int i=0; i<nrow; i++){
        for (int j=0; j<vec.size(); j++){
            cout << mat[i][j] << " ";
	      }
        cout << endl;
    }
    cout << count << endl;

    // OpenMP assignment
    // !!! Without private(scalar), sometimes wrong values are assigned !!!
    cout << "OpenMP - assign 2D matrix without private" << endl;
    vec = {1,2,3,4,5,6,7,8,9,10};
    #pragma omp parallel for
    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < vec.size(); j++){
	          count = i*j;
            mat[i][j] = i + j;
            mat[i][j] = count;
        }
	      #pragma omp critical
	      cout << i << " " << count << endl;
    }
    cout << "Array: " << endl;
    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < vec.size(); j++){
            cout << mat[i][j] << " ";
	      }
        cout << endl;
    }
    cout << "Scalar: " << count << endl;

    cout << "OpenMP - assign 2D matrix with private" << endl;
    vec = {1,2,3,4,5,6,7,8,9,10};
    #pragma omp parallel for private(count)
    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < vec.size(); j++){
	          count = i*j;
            mat[i][j] = i + j;
            mat[i][j] = count;
        }
	      #pragma omp critical
	      cout << i << " " << count << endl;
    }
    cout << "Array: " << endl;
    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < vec.size(); j++){
            cout << mat[i][j] << " ";
	      }
        cout << endl;
    }
    cout << "Scalar: " << count << endl;
    return 0;
}

/* - Example Output -
Simple - assign 2D matrix
0 1 2 3 4 5 6 7 8 9 
1 2 3 4 5 6 7 8 9 10 
2 3 4 5 6 7 8 9 10 11 
3 4 5 6 7 8 9 10 11 12 
27
OpenMP - assign 2D matrix without private
0 27
1 27
3 27
2 27
Array: !!! Without private(scalar), sometimes wrong values are assigned !!!
1 0 0 0 0 0 0 0 0 0 
12 1 2 3 4 5 6 7 8 9 
0 2 4 6 8 10 0 14 16 0 
0 3 6 9 12 15 18 21 24 27 
Scalar: 27
OpenMP - assign 2D matrix with private
1 9
0 0
3 27
2 18
Array: 
0 0 0 0 0 0 0 0 0 0 
0 1 2 3 4 5 6 7 8 9 
0 2 4 6 8 10 12 14 16 18 
0 3 6 9 12 15 18 21 24 27 
Scalar: 27
*/
