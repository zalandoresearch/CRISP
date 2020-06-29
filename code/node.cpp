#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"
#include <cassert>




void operator *= (vector<double>& a, const vector<double>& b) {

    assert(a.size() == b.size());
    for( unsigned long i=0; i<a.size(); i++) {
        a[i] *= b[i];
    } 
}

Message normalize(const Message &v_) {
    double sum=0;
    Message v(v_);
    for( auto x: v) sum+=x;
    for( auto &x: v) x/=sum;
    return v;
}

Node::Node( int N): _N(N),  _messages() {}


void Node::update() {

    Message m_update( _N, 1.0);

    for( size_t i=0; i<_factors.size(); i++) {
        _messages[i] = _factors[i]->message_to( this);
        Node::message_counter++;
    }
}


MessagePtr Node::message_to( const Factor *f) const {
    MessagePtr m(new Message(_N, 1.0));
    for( size_t i=0; i<_factors.size(); i++) 
        if( _factors[i] != f)
            (*m) *= (*_messages[i]);
    return m;
}

int Node::message_counter = 0;


SEIRNode::SEIRNode( const SEIRStateSpace &all_states) :
    Node(all_states.size()), _states( all_states) 
{
}

