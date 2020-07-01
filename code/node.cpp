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

void operator /= (vector<double>& a, const vector<double>& b) {

    assert(a.size() == b.size());
    for( unsigned long i=0; i<a.size(); i++) {
        if( abs(b[i])>1e-10) //numeric_limits<double>::epsilon())    
            a[i] /= b[i];
        // else  { 
        //     cerr << "small a=" << a[i] << endl;
        //     assert(abs(a[i])<=10*numeric_limits<double>::epsilon());
        // }    
    } 
}

Message & normalize( Message &v) {
    double sum=0;
    for( auto x: v) sum+=x;
    for( auto &x: v) x/=sum;
    return v;
}

Message normalize( const Message &v_) {
    Message v(v_);
    double sum=0;
    for( auto x: v) sum+=x;
    for( auto &x: v) x/=sum;
    return v;
}

Node::Node( int N): _N(N),  _messages() , _marginal(N,1.0) {}


void Node::add_factor( Factor* f) {
    _factors.push_back(f);
    _messages.emplace_back( new Message(_N, 1.0) );
}


void Node::update() {

    for( size_t i=0; i<_factors.size(); i++) {
        _marginal /= *(_messages[i]);
        _factors[i]->message_to( this, _messages[i]);
        _marginal *= *(_messages[i]);
        Node::message_counter++;
    }
}


MessagePtr Node::message_to( const Factor *f) const {
    MessagePtr m(new Message(_marginal));
    for( size_t i=0; i<_factors.size(); i++) {
        if( _factors[i] == f)
           (*m) /= (*_messages[i]);
    }
    normalize(*m);
    return m;
}

int Node::message_counter = 0;


SEIRNode::SEIRNode( const SEIRStateSpace &all_states, double p1) :
    Node(all_states.size()), 
    _states( all_states), 
    _z_messages(), 
    _z_marginal(2),
    _p1(p1)
{
}


void SEIRNode::add_factor( Factor* f) {
    Node::add_factor(f);
    _z_messages.emplace_back( new Message(2, 1.0) );
}


void SEIRNode::update() {

    for( size_t i=0; i<_factors.size(); i++) {
        _marginal /= *(_messages[i]);
        _factors[i]->message_to( this, _messages[i]);
        _marginal *= *(_messages[i]);

        _z_marginal /= *(_z_messages[i]);
        fill( _z_messages[i]->begin(), _z_messages[i]->end(), 0.0);
        for( auto s: _states) {
            int j = _states[s];
            if( s.phase()==SEIRState::I)
                (*_z_messages[i])[1] += (*_messages[i])[j];  
            else
                (*_z_messages[i])[0] += (*_messages[i])[j];  
        }
        _z_marginal *= *(_z_messages[i]);


        Node::message_counter++;
    }
}

double min_P = 1.0;
MessagePtr SEIRNode::infection_message_to( const Factor *f) const {
    
    MessagePtr m = message_to(f);
    double p_I = 0.0;
    for( SEIRState s=SEIRState(SEIRState::I,1); +_states.can_continue(s); s = s.next(/* change = */false)) {
        p_I += (*m)[_states[s]];
    }   
    return MessagePtr(new Message({1.0-p_I, p_I}));
}