#include "seir.hpp"
#include "node.hpp"
#include "factor.hpp"
#include <iostream>

using namespace std;



vector<double> prob_vector( int N) {
    auto v = vector<double>(N);
    double sum = 0;
    for( int i=0; i<N; i++) {
        sum += v[i] = (double)rand() / RAND_MAX;
    }
    for( int i=0; i<N; v[i++] /= sum);
    return v;
}



void loopy_test() {

    int N = 3;
    Node A(N), B(N), C(N), D(N);

    auto tA = prob_vector(N);
    vector<Node*> nA = {&A};
    TableFactor fA( nA, tA );  

    auto tB = prob_vector(N);
    vector<Node*> nB = {&B};
    TableFactor fB( nB, tB );  

    vector<double> tABC;
    for(int i=0; i<N*N; i++) {
        auto v = prob_vector(N);
        tABC.insert(tABC.end(),v.begin(),v.end());
    }
    vector<Node*> nABC = {&A, &B, &C};
    TableFactor fABC( nABC, tABC );  

    vector<double> tBCD;
    for(int i=0; i<N*N; i++) {
        auto v = prob_vector(N);
        tBCD.insert(tBCD.end(),v.begin(),v.end());
    }
    vector<Node*> nBCD = {&B, &C, &D};
    TableFactor fBCD( nBCD, tBCD );  

    cout << std::fixed; 
    cout.precision(5);

    for( int i=0; i<5; i++) {
        D.update();
        C.update();
        B.update();
        A.update();
        cout << "A:" << normalize(A.message_to()) << endl;
        cout << "B:" << normalize(B.message_to()) << endl;
        cout << "C:" << normalize(C.message_to()) << endl;
        cout << "D:" << normalize(D.message_to()) << endl;
        cout << endl;
    }

    vector<double> PA = tA;
    vector<double> PB = tB;
    vector<double> PC(N, 0.0);
    auto it = tABC.begin();
    for(int i=0; i<N; i++) 
        for(int j=0; j<N; j++)
            for(int k=0; k<N; k++) 
                PC[k] += PA[i] * PB[j] *  (*it++);
    vector<double> PD(N, 0.0);
    it = tBCD.begin();
    for(int i=0; i<N; i++) 
        for(int j=0; j<N; j++)
            for(int k=0; k<N; k++) 
                PD[k] += PB[i] * PC[j] *  (*it++);

    

    cout << "P(A) = " << PA << endl;
    cout << "P(B) = " << PB << endl;
    cout << "P(C) = " << PC << endl;
    cout << "P(D) = " << PD << endl;
}


void seir_state_test() {
    for( auto s: SEIRState::all_states( 10, 15))
        cout << s <<" ";
    cout << endl;




    SEIRNode node1( SEIRState::all_states( 10, 15)); 
    SEIRNode node2( SEIRState::all_states( 10, 15)); 
    VirusLoadNode load(4);


    vector<double> qE = {0.0000000000, 0.05908981283, 0.1656874653, 0.1819578343, 0.154807057,
                0.1198776096, 0.08938884645, 0.06572939883, 0.04819654533,
                0.03543733758, 0.02620080839, 0.01950646727, 0.01463254844,
                0.0110616426, 0.008426626119};

    vector<double> qI = {0.000000000000, 0.000000000000, 0.00000000000, 0.000000000000, 0.000000000000,
                0.0001178655952, 0.0006658439543, 0.002319264193, 0.005825713197, 0.01160465163,
                0.01949056696, 0.02877007836, 0.03842711373, 0.04743309657, 0.05496446107,
                0.06050719418, 0.06386313651, 0.065094874, 0.06444537162, 0.06225794729,
                0.0589104177, 0.05476817903, 0.05015542853, 0.0453410888, 0.04053528452,
                0.03589255717, 0.03151878504, 0.02747963753, 0.02380914891, 0.02051758911,
                0.01759822872, 0.01503287457, 0.0127962154, 0.01085910889, 0.009190974483,
                0.007761463001, 0.006541562648, 0.005504277076};
    double p0 = 0.01;
    double p1 = 0.5;

    SEIRFactor f(qE, qI, p0, p1, node1, node2, load);


    f.message_to( &node1);
}


int main() {
    loopy_test();
    return 0;
}
