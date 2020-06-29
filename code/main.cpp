#include "seir.hpp"
#include "node.hpp"
#include "factor.hpp"
#include <iostream>
#include <list>
#include <memory>

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

double sum( const Message &m) {
    double res = 0.0;
    for( auto i=m.begin(); i!=m.end();++i) res += *i;
    return res;
}
 


void loopy_test() {

/*
small example of loopy BP, node D ends up with incorrect marginals
                  .─.            
                 ( B )           
                 ╱`─'╲           
       .─.      ╱     ╲      .─. 
██────( A )────██     ██────( D )
       `─'      ╲     ╱      `─' 
                 ╲.─.╱           
                 ( C )           
                  `─'            

*/
    int N = 3;
    Node A(N), B(N), C(N), D(N);

    auto tA = prob_vector(N);
    vector<Node*> nA = {&A};
    TableFactor fA( nA, tA );  

    vector<double> tABC;
    for(int i=0; i<N; i++) {
        auto v = prob_vector(N*N);
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
        B.update();
        D.update();
        A.update();
        C.update();
        cout << "iteration " << i+1 << endl;
        cout << "A:" << normalize(*A.message_to()) << endl;
        cout << "B:" << normalize(*B.message_to()) << endl;
        cout << "C:" << normalize(*C.message_to()) << endl;
        cout << "D:" << normalize(*D.message_to()) << endl;
        cout << endl;
    }

    vector<double> PA = tA;
    vector<double> PB(N,0.0);
    vector<double> PC(N, 0.0);
    vector< vector<double>> PBC(N, vector<double>(N,0.0));
    auto it = tABC.begin();
    for(int i=0; i<N; i++) 
        for(int j=0; j<N; j++)
            for(int k=0; k<N; k++) {
                PB[j] += PA[i] *(*it);
                PC[k] += PA[i] *(*it);
                PBC[j][k] += PA[i] *(*it);
                it++;
            }

    vector<double> PD(N, 0.0);
    vector<double> PD_(N, 0.0);
    it = tBCD.begin();
    for(int i=0; i<N; i++) 
        for(int j=0; j<N; j++)
            for(int k=0; k<N; k++) {
                PD[k] += PBC[i][j] *(*it);
                PD_[k] += PB[i]*PC[j] *(*it);
                it++;
            }

    cout << "P(A) = " << PA << endl;
    cout << "P(B) = " << PB << endl;
    cout << "P(C) = " << PC << endl;
    cout << "P(D) = " << PD << "(true values, \\sum_{B,C} P(D|B,C) P(B,C) )" << endl;
    cout << "P(D) = " << PD_ << "(factorizing, \\sum_{B,C} P(D|B,C) P(B) P(C) )" << endl;
}



void seir_state_test() {


    // !!! t=0 based !!!
    vector<double> qE = {0.0, 0.05908981283, 0.1656874653, 0.1819578343, 0.154807057,
                0.1198776096, 0.08938884645, 0.06572939883, 0.04819654533,
                0.03543733758, 0.02620080839, 0.01950646727, 0.01463254844,
                0.0110616426, 0.008426626119};

    vector<double> qI = {0.0, 0.000000000000, 0.00000000000, 0.000000000000, 0.000000000000,
                0.0001178655952, 0.0006658439543, 0.002319264193, 0.005825713197, 0.01160465163,
                0.01949056696, 0.02877007836, 0.03842711373, 0.04743309657, 0.05496446107,
                0.06050719418, 0.06386313651, 0.065094874, 0.06444537162, 0.06225794729,
                0.0589104177, 0.05476817903, 0.05015542853, 0.0453410888, 0.04053528452,
                0.03589255717, 0.03151878504, 0.02747963753, 0.02380914891, 0.02051758911,
                0.01759822872, 0.01503287457, 0.0127962154, 0.01085910889, 0.009190974483,
                0.007761463001, 0.006541562648, 0.005504277076};
    
    double p0 = 0.02;
    double p1 = 0.99;

    double alpha = 0.001;
    double beta = 0.01;

    SEIRStateSpace states( qE.size()-1, qI.size()-1);
    for( auto s: states)
        cout << s <<" ";
    cout << endl;

    SEIRNode node1( states, p1); 
    SEIRNode node2( states, p1); 


    unsigned int S = 3;
    unsigned int T = 30;
    vector<unique_ptr<SEIRNode>> nodes;//, SEIRNode(states));
    for( int u=0; u<S; u++)
        for( int t=0; t<T; t++)
            nodes.emplace_back(new SEIRNode(states, p1));


    vector<unique_ptr<Factor>> factors;
    for( unsigned int s=0; s<S; s++)
        factors.emplace_back(new SEIRInitFactor(*nodes[s]));

    for( unsigned int t=0; t<T-1; t++) {
        for( unsigned int s=0; s<S; s++) {
            if      (t==17 && s==0) factors.emplace_back( new SEIRFactor(qE, qI, p0, p1, *nodes[S*t+s], *nodes[S*(t+1)+s], vector<SEIRNode*>({nodes[S*t+1].get()}) ));
            else if (t==17 && s==1) factors.emplace_back( new SEIRFactor(qE, qI, p0, p1, *nodes[S*t+s], *nodes[S*(t+1)+s], vector<SEIRNode*>({nodes[S*t+0].get(),nodes[S*t+2].get()}) ));
            else if (t==17 && s==2) factors.emplace_back( new SEIRFactor(qE, qI, p0, p1, *nodes[S*t+s], *nodes[S*(t+1)+s], vector<SEIRNode*>({nodes[S*t+1].get()}) ));
            
            else
                factors.emplace_back( new SEIRFactor(qE, qI, p0, p1, *nodes[S*t+s], *nodes[S*(t+1)+s] ));
        }
    }   


    factors.emplace_back( new SEIRTestFactor( *nodes[S*9+0], false, alpha, beta) );
    factors.emplace_back( new SEIRTestFactor( *nodes[S*15+0], true, alpha, beta) );
    factors.emplace_back( new SEIRTestFactor( *nodes[S*19+0], true, alpha, beta) );
    factors.emplace_back( new SEIRTestFactor( *nodes[S*29+2], true, alpha, beta) );

    // for( auto f = factors.begin(); f!=factors.end(); f++) {
    //     cerr << (*f)->_nodes.size() <<"(" << (*f) << "): ";
    //     for( auto n: (*f)->_nodes) cerr << (void*)(n) << " "; 
    //     cerr << endl;
    // }
    cerr << factors.size() << "factors" << endl;

    for( int i=0; i<25; i++) {
        for(auto &node: nodes) 
            node->update();
        for( auto node = nodes.rbegin(); node !=nodes.rend(); ++node)
            (*node)->update();
    }

    cout << std::fixed; 
    cout.precision(3);

    // auto fac = factors[0]->message_to(&nodes[1]);
    // cout << *fac << endl;
    // cout << "sum="<< sum(*fac) << endl;
    for(auto &node: nodes) {
        cout << normalize(basic_states(*node->message_to(), node->states())) <<endl;
    }
    //f.message_to( &node1);

}



void virus_load_test() {

    VirusLoad Px;
    cout << Px << endl;

    // Px.add_source( 0.8, 1.0);
    // Px.add_source( 0.8, 1.0);
    // Px.add_source( 0.8, 1.0);
    Px.add_source( 0.1, 1.0);
    Px.add_source( 0.1, 1.0);
    Px.add_source( 0.1, 1.1);
    Px.add_source( 0.1, 0.95);
    

    cout << Px << endl;

}

int main() {
    loopy_test();
    virus_load_test();
    seir_state_test();
    return 0;
}




/* ToDo:
- strip inner loops of SEIRFactor::message_horizontally, and SEIRFactor::message_vertical 
- make VirusLoad nicer, collecting the load distribution may need explicit scaling of the incoming factor
- cleaup _states, make it global for all SEIRNodes
- implement short_messages (aggregated over the four SEIR states )

*/