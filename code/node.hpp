#ifndef _NODE_HPP
#define _NODE_HPP

#include <memory>
#include <iostream>
#include <vector>

#include "seir.hpp"

using namespace std;


typedef vector<double> Message;
typedef shared_ptr<Message> MessagePtr;



Message basic_states( const Message &message, const SEIRStateSpace &states);



template <class T>
ostream& operator<<(ostream& os, const vector<T> &v) {
    os << "[";
    for(unsigned int i=0; i!=v.size(); i++) os << v[i] << (i==v.size()-1 ? "": ", ") ;
    os << "],";
    return os;
}

//vector<double> normalize(const vector<double> &v_);
Message & normalize( Message &v_);
Message normalize( const Message &v_);


class Factor;

class Node {
    Node( const Node &other) = delete;
    
public:
    vector<Factor*> _factors;            // list of associated factors

public:
    int _N;                              // length of message as discrete distribution
    vector<MessagePtr> _messages;        // local message used for loopy belief propagation
    Message _marginal;

public:
    Node( int N);
    virtual void add_factor( Factor* f);

    int size() const {return _N;};
    virtual void update();

    MessagePtr message_to(const Factor *f = 0) const;

    static int message_counter;
};



class SEIRNode: public Node {
    const SEIRStateSpace &_states;

    double _p1;

    vector<MessagePtr> _z_messages;       
    Message _z_marginal;

public:
    SEIRNode( const SEIRStateSpace &all_states, double p1);
    virtual void add_factor( Factor* f);
    virtual void update();

    const SEIRStateSpace& states() {return _states;};
    MessagePtr infection_message_to( const Factor *f =0) const;
};


#endif