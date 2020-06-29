
#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>

using namespace std;


Factor::Factor( const vector<Node*> &nodes, const vector<Node*> &child_nodes) : 
    _nodes(nodes) 
{
    if( child_nodes.size()>0) {
        for( auto n: child_nodes) {
            n->_factors.push_back(this);
            n->_messages.emplace_back( new Message(n->_N, 1.0) );
        }
    }
    else {
        for( auto n: _nodes) {
            n->_factors.push_back(this);
            n->_messages.emplace_back( new Message(n->_N, 1.0) );
        }
    }
}


void Factor::message_to( Node* n, MessagePtr to) {

    fill( to->begin(), to->end(), 0.0);

    size_t i_node = 0; 
    for( ; i_node<_nodes.size() && _nodes[i_node]!=n; i_node++){
    }
    assert( i_node<_nodes.size()); // n needs to be in _nodes


    vector< MessagePtr>  incoming_messages;
    vector< Message::const_iterator>   msg_its;

    vector< unsigned int> state_its(_nodes.size());

    for( unsigned int i=0; i!=_nodes.size(); i++ ) {
        incoming_messages.push_back(_nodes[i]->message_to(this));
        msg_its.push_back(incoming_messages[i]->begin());
    }

    MessagePtr outgoing_message = to;
    Message::iterator out_msg_it = outgoing_message->begin();

    while( true) {

        double p = 1.0;
        for(  size_t i=0; i<_nodes.size(); i++) {
            if( i!=i_node) {
                p *= *(msg_its[i]);
            }
        }
        double upd = p * potential( state_its);
        // cout << *(out_msg_it) << " " << upd << " "<< p << " " << endl;
        *(out_msg_it) += upd;

        // generate the next state configuration
        for( int i=_nodes.size()-1; i>=0; i--) {
            msg_its[i]++;
            state_its[i]++;
            if(_nodes[i] == n) out_msg_it++;
            
            if( msg_its[i] == incoming_messages[i]->end()) {
                // if they overrun, wrap them over to the beginning 
                msg_its[i]   = incoming_messages[i]->begin();
                state_its[i] = 0;  
                if(_nodes[i] == n) 
                    out_msg_it = outgoing_message->begin(); 

                // if the first iterator wrapped around we are done
                if( i==0) goto done;           
            }          
            else break;
        }
    } 
    done:
    return;
}

TableFactor::TableFactor( vector<Node*> &nodes, Message tab) :
    Factor(nodes), _tab(tab)
{
}

double TableFactor::potential( const vector<unsigned int> & state_its) {
    assert( state_its.size()== _nodes.size());

    int idx = 0;
    for( unsigned int i=0; i<_nodes.size(); i++) {
        if( i>0) idx *= _nodes[i-1]->size();
        idx += state_its[i];
        //cout << _nodes[i]->size() << " " <<  std::any_cast<int>(*(state_its[i])) << endl;
    }
    //cout << "ïdx: " << idx << endl;
    return _tab[idx];
}


SEIRFactor::SEIRFactor( const Distribution &qE, const Distribution &qI, 
                        double p0, double p1, 
                        SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts,
                        const vector<Node*> &child_nodes
                        ) : 
    Factor(init_helper(in, out, contacts), child_nodes), _piE(init_pi(qE)), _piI(init_pi(qI)), _p0(p0), _p1(p1), _states(in.states()) 
{
}


vector<Node *> SEIRFactor::init_helper(SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts ) {
    vector<Node *> nodes({ &in, &out});
    nodes.insert( nodes.end(), contacts.begin(), contacts.end());
    return nodes;
}


const vector<double> SEIRFactor::init_pi( const Distribution &q) {

    auto res = Message(q.getMaxOutcomeValue()+1);
    for( int i=0; i<res.size(); i++)
        res[i] = exp( q.getLogP(i)-q.getLogPTail(i)); 
    return res;
}



void SEIRFactor::message_forward( MessagePtr to) {


    MessagePtr output_message = to;
    MessagePtr input_message = _nodes[0]->message_to(this);
    

    double p_keep = (1.0-_p0);
    for( size_t i=2; i<_nodes.size(); i++) {
        p_keep *= ((SEIRNode*)_nodes[i])->infection_message_to(this);
    }

    auto it_input_message = input_message->cbegin();
    for( auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) 
    {
        auto in = *it_input_states;
        const double& in_value = *it_input_message;

        switch( in.phase() ) {
            case SEIRState::S: {
                (*output_message)[_states[in.next(/*change = */ true)]] += in_value * (1.0 - p_keep);
                (*output_message)[_states[in.next(/*change = */ false)]] += in_value * p_keep;
            } break;
            case SEIRState::E: {
                (*output_message)[_states[in.next(/*change = */ true)]]  += in_value * _piE[in.days()];
                if( _states.can_continue(in))   
                    (*output_message)[_states[in.next(/*change = */ false)]] += in_value * (1.0-_piE[in.days()]);

            } break;
            case SEIRState::I: {
                (*output_message)[_states[in.next(/*change = */ true)]]  += in_value * _piI[in.days()];
                if( _states.can_continue(in))   
                    (*output_message)[_states[in.next(/*change = */ false)]] += in_value * (1.0-_piI[in.days()]);
            } break;
            case SEIRState::R: {
                (*output_message)[_states[in]] += in_value;
            } break;
        }
        it_input_message++;
    }
}

void SEIRFactor::message_backward( MessagePtr to) {


    MessagePtr input_message = to;
    MessagePtr output_message = _nodes[1]->message_to(this);


    double p_keep = (1.0-_p0);
    for( size_t i=2; i<_nodes.size(); i++) {
        p_keep *= ((SEIRNode*)_nodes[i])->infection_message_to(this);
    }

    for( auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) 
    {
        auto in = *it_input_states;

        switch( in.phase() ) {
            case SEIRState::S: {
                (*input_message)[_states[in]] += (*output_message)[_states[in.next(/*change = */ true)]] * (1.0 - p_keep);
                (*input_message)[_states[in]] += (*output_message)[_states[in.next(/*change = */ false)]] * p_keep;
            } break;
            case SEIRState::E: {
                (*input_message)[_states[in]] += (*output_message)[_states[in.next(/*change = */ true)]] * _piE[in.days()];
                if( _states.can_continue(in))   
                    (*input_message)[_states[in]] += (*output_message)[_states[in.next(/*change = */ false)]] * (1.0-_piE[in.days()]);

            } break;
            case SEIRState::I: {
                (*input_message)[_states[in]] += (*output_message)[_states[in.next(/*change = */ true)]] * _piI[in.days()];
                if( _states.can_continue(in))   
                    (*input_message)[_states[in]] += (*output_message)[_states[in.next(/*change = */ false)]] * (1.0-_piI[in.days()]);
            } break;
            case SEIRState::R: {
                (*input_message)[_states[in]] += (*output_message)[_states[in]];
            } break;
        }
    }
}
        

void SEIRFactor::message_vertical( Node *n, MessagePtr to) {


    auto input_message = _nodes[0]->message_to(this);
    auto output_message = _nodes[1]->message_to(this);

    double p_keep = (1.0-_p0);
    cerr << _nodes.size() << " nodes" << endl;
    for( size_t i=2; i<_nodes.size(); i++) {
        p_keep *= ((SEIRNode*)_nodes[i])->infection_message_to(this);
    }

    Message p_outgoing(2, 0.0);
    for( auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) {
        auto in = *it_input_states;

        switch( in.phase() ) {
            case SEIRState::S: {
                for(int j=0; j<2; j++) {
                   double p_keep_ = p_keep * (j>0 ? (1.0-_p1) : 1);
                    p_outgoing[j] += (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ true)]] * (1.0 - p_keep_);
                    p_outgoing[j] += (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ false)]] * p_keep_;
                }
            } break;
            case SEIRState::E: {
                double p;
                p =  (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ true)]] * _piE[in.days()];
                p += (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ false)]] * (1.0 - _piE[in.days()]);
                p_outgoing[0] += p;
                p_outgoing[1] += p;
            } break;                    
            case SEIRState::I: {
                double p;
                p =  (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ true)]] * _piI[in.days()];
                p += (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ false)]] * (1.0 - _piI[in.days()]);
                p_outgoing[0] += p;
                p_outgoing[1] += p;
            } break;                    
            case SEIRState::R: {
                double p = (*input_message)[_states[in]] * (*output_message)[_states[in.next(/*change = */ false)]];
                p_outgoing[0] += p;
                p_outgoing[1] += p;
            } break;
        }
    }

    p_outgoing = normalize(p_outgoing);
    MessagePtr outgoing_message = to;
    for(int i=0; i<n->_N; i++) {
        (*outgoing_message)[i] = (_states[i].phase()==SEIRState::I ? p_outgoing[1] : p_outgoing[0]);
    }
}


void SEIRFactor::message_to( Node *n, MessagePtr to) {

    fill( to->begin(), to->end(), 0.0);
    if( n==_nodes[0]) return message_backward( to);
    if( n==_nodes[1]) return message_forward( to);
    return message_vertical( n, to);
}


double SEIRFactor::potential( const vector<unsigned int> & ) {

    assert(false); // SEIRFactor implements its own message_to() method
    return 0.0;
}


SEIRInitFactor::SEIRInitFactor( SEIRNode &out, bool patient_zero) :
    Factor({&out}), _patient_zero(patient_zero), _states(out.states())
{
}

double SEIRInitFactor::potential( const vector<unsigned int> &state_its) {
     assert(false); // SEIRInitFactor implements its own message_to() method
     return 0.0;
}

void SEIRInitFactor::message_to( Node *n, MessagePtr to) {

    fill( to->begin(), to->end(), 0.0);
    if( _patient_zero)
        (*to)[_states[SEIRState(SEIRState::E,1)]] = 1.0;
    else
        (*to)[_states[SEIRState(SEIRState::S)]] = 1.0;
}

SEIRTestFactor::SEIRTestFactor( SEIRNode &out, bool positive, double alpha, double beta) :
    Factor({&out}), _states(out.states())
{
    _positive = positive;
    _alpha = alpha;
    _beta = beta;
}

double SEIRTestFactor::potential( const vector<unsigned int> &state_its) {

    SEIRState out = _states[state_its[0]];

    if( _positive) 
        return out.phase() == SEIRState::I ? 1.0-_alpha : _beta;
    else 
        return out.phase() == SEIRState::I ? _alpha : 1.0-_beta;
}