#ifndef _NODE_HPP
#define _NODE_HPP

#include <any>
#include <memory>
#include <iostream>

using namespace std;

typedef vector<double> Message;
typedef shared_ptr<Message> MessagePtr;

#include "seir.hpp"





template <class T>
ostream& operator<<(ostream& os, const vector<T> &v) {
    os << "[";
    for(unsigned int i=0; i!=v.size(); i++) os << v[i] << (i==v.size()-1 ? "": ", ") ;
    os << "],";
    return os;
}

vector<double> normalize(const vector<double> &v_);


class Factor;

class Node {
    friend class Factor;
public:
    vector<Factor*> _factors;            // list of associated factors

public:
    int _N;                              // length of message as discrete distribution
    vector<MessagePtr> _messages;        // local message used for loopy belief propagation
    vector<any> _states;

public:
    Node( int N);
    int size() const {return _N;};
    void update();

    const vector<any> & states() {return _states;};
    MessagePtr message_to(const Factor *f = 0) const;
};



class SEIRNode: public Node {

public:
    SEIRNode( const vector<SEIRState> &all_states);

};


#endif