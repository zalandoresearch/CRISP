#ifndef _NODE_HPP
#define _NODE_HPP


#include "seir.hpp"

#include <iostream>
#include <any>

using namespace std;

template <class T>
ostream& operator<<(ostream& os, const vector<T> &v) {
    os << "(";
    for(unsigned int i=0; i!=v.size(); i++) os << v[i] << (i==v.size()-1 ? "": ", ") ;
    os << ")";
    return os;
}

vector<double> normalize(const vector<double> &v_);


class Factor;

class Node {
    friend class Factor;

    vector<Factor*> _factors;            // list of associated factors

protected:
    int _N;                              // length of message as discrete distribution
    vector<vector<double>> _messages;    // local message used for loopy belief propagation
    vector<any> _states;

public:
    Node( int N);
    int size() const {return _N;};
    void update();

    const vector<double> message_to(const Factor *f = 0) const;
};



class SEIRNode: public Node {

public:
    SEIRNode( const vector<SEIRState> &all_states);

};

class VirusLoadNode: public Node {
public:
    VirusLoadNode( unsigned long N);
};


#endif