#include "node.hpp"
#include "factor.hpp"
#include "seir.hpp"
#include <cassert>



void operator *= (vector<double>& a, const vector<double>& b) {

    assert(a.size() == b.size());
    for( unsigned long i=0; i<a.size(); i++) {
        a[i] = b[i];
    } 


}



Node::Node( int N): _N(N), _state(N), _message(N, 1.0) {}


void Node::update() {

    vector<double> m_update( _N, 1.0);

    for( auto f: _factors) {
        auto m = f->message_to( this);
        m_update *= m ;
    }

    _message = m_update;
}



SEIRNode::SEIRNode( const vector<SEIRState> &all_states) :
    Node(all_states.size()) 
{
    for(int i=0; i<all_states.size(); i++)
        _state[i] = all_states[i];
}


VirusLoadNode::VirusLoadNode( int N) :
    Node(N) 
{
    for( int i=0; i<N; i++)
        _state[i] = (double) i;
}