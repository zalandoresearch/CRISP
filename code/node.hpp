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

vector<double> normalize(const vector<double> &v_);


class Factor;

class Node {
    friend class Factor;

    Node( const Node &other) = delete;
    
public:
    vector<Factor*> _factors;            // list of associated factors

public:
    int _N;                              // length of message as discrete distribution
    vector<MessagePtr> _messages;        // local message used for loopy belief propagation
   
public:
    Node( int N);
    int size() const {return _N;};
    void update();

    MessagePtr message_to(const Factor *f = 0) const;

    static int message_counter;
};



class SEIRNode: public Node {
    const SEIRStateSpace &_states;
public:
    SEIRNode( const SEIRStateSpace &all_states);
    const SEIRStateSpace& states() {return _states;};

};


#endif