#ifndef _NODE_HPP
#define _NODE_HPP

#include "seir.hpp"

#include <vector>
#include <any>

using namespace std;

class Factor;

class Node {
    friend class Factor;

    vector<Factor*> _factors;   // list of associated factors

protected:
    int _N;                     // length of message as discrete distribution
    vector<double> _message;          // local message used for loopy belief propagation
    vector<any> _state;

public:
    Node( int N);

    void update();
};



class SEIRNode: public Node {

public:
    SEIRNode( const vector<SEIRState> &all_states);

};

class VirusLoadNode: public Node {
public:
    VirusLoadNode( int N);
};


#endif