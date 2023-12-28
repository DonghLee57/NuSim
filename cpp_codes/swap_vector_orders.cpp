#include <iostream>
#include <vector>
using namespace std;
void my_swap(vector<vector<vector<int>>>& MAT);

int main(){
    vector<vector<vector<int>>> MAT3D;
    MAT3D = {{{1,2,3},{4,5,6},{7,8,9}},
             {{0,0,0},{0,0,0},{0,0,0}}};

    // Print initial state
    cout << "Initial state" << endl;
    for (int i=0; i < MAT3D.size(); i++){
        for (int j=0; j < MAT3D[i].size(); j++){
            for (int k=0; k < MAT3D[i][j].size(); k++){
                cout << MAT3D[i][j][k] << " ";
            }   
            cout << endl;
        }
        cout << endl;
    }

    // Using 'swap' function in std
    swap(MAT3D[0],MAT3D[1]);
    
  
    // Print final state after swap
    cout << "Final state" << endl;
    for (int i=0; i < MAT3D.size(); i++){
        for (int j=0; j < MAT3D[i].size(); j++){
            for (int k=0; k < MAT3D[i][j].size(); k++){
                cout << MAT3D[i][j][k] << " ";
            }   
            cout << endl;
        }
        cout << endl;
    }

    my_swap(MAT3D);
    // Print final state after swap
    cout << "Final state" << endl;
    for (int i=0; i < MAT3D.size(); i++){
        for (int j=0; j < MAT3D[i].size(); j++){
            for (int k=0; k < MAT3D[i][j].size(); k++){
                cout << MAT3D[i][j][k] << " ";
            }   
            cout << endl;
        }
        cout << endl;
    }
    return 0;
}

void my_swap(vector<vector<vector<int>>>& MAT){
    swap(MAT[0],MAT[1]);
}
