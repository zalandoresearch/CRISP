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

Node::Node( int N): _N(N),  _messages(), _states(N) {
    for( int i=0; i<N; i++)   
        _states[i] = i;
}


void Node::update() {

    Message m_update( _N, 1.0);

    for( size_t i=0; i<_factors.size(); i++) {
        _messages[i] = _factors[i]->message_to( this);
    }
}


MessagePtr Node::message_to( const Factor *f) const {
    MessagePtr m(new Message(_N, 1.0));
    for( size_t i=0; i<_factors.size(); i++) 
        if( _factors[i] != f)
            (*m) *= (*_messages[i]);
    return m;
}



SEIRNode::SEIRNode( const vector<SEIRState> &all_states) :
    Node(all_states.size()) 
{
    for(size_t i=0; i<all_states.size(); i++)
        _states[i] = all_states[i];
}

