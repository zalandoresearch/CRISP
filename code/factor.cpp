
#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <any>
#include <cstdint>
#include <functional>

using namespace std;


Factor::Factor( const vector<Node*> &nodes) : 
    _nodes(nodes) 
{
    for( auto n: _nodes) {
        n->_factors.push_back(this);
        n->_messages.push_back( vector<double>(n->_N, 1.0) );
    }
}

vector<double> Factor::message_to( Node* n) {

    size_t i_node = 0; 
    for( ; i_node<_nodes.size() && _nodes[i_node]!=n; i_node++);
    assert( i_node<_nodes.size()); // n needs to be in _nodes


    vector< vector<double>>  incoming_messages;
    vector< vector<std::any>> incoming_states;
    vector< vector<double>::const_iterator>   msg_its;
    vector< vector<std::any>::const_iterator> state_its;

    for( unsigned int i=0; i!=_nodes.size(); i++ ) {
        incoming_messages.push_back(_nodes[i]->message_to(this));
        msg_its.push_back(incoming_messages[i].begin());

        incoming_states.push_back(_nodes[i]->_states);
        state_its.push_back( incoming_states[i].begin());
    }

    vector<double> outgoing_message( n->_N, 0.0);
    vector<double>::iterator out_msg_it = outgoing_message.begin();

    while( true) {

        double p = 1.0;
        for(  int i=0; i<_nodes.size(); i++) {
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
            
            if( msg_its[i] == incoming_messages[i].end()) {
                // if they overrun, wrap them over to the beginning 
                msg_its[i]   = incoming_messages[i].begin();
                state_its[i] = incoming_states[i].begin();  
                if(_nodes[i] == n) 
                    out_msg_it = outgoing_message.begin(); 

                // if the first iterator wrapped around we are done
                if( i==0) goto done;           
            }          
            else break;
        }
    } 
    done:

    return outgoing_message;
}

TableFactor::TableFactor( vector<Node*> &nodes, vector<double> tab) :
    Factor(nodes), _tab(tab)
{
}

double TableFactor::potential( vector<vector<any>::const_iterator> state_its) {
    assert( state_its.size()== _nodes.size());

    int idx = 0;
    for( unsigned int i=0; i<_nodes.size(); i++) {
        if( i>0) idx *= _nodes[i-1]->size();
        idx += std::any_cast<int>(*(state_its[i]));
        //cout << _nodes[i]->size() << " " <<  std::any_cast<int>(*(state_its[i])) << endl;
    }
    //cout << "ïdx: " << idx << endl;
    return _tab[idx];
}


SEIRFactor::SEIRFactor( const vector<double> &qE, const vector<double> &qI, 
                        double p0, double p1, 
                        SEIRNode &in, SEIRNode &out, VirusLoadNode &load) : 
    Factor({&in, &out, &load}), _qE(qE), _qI(qI), _piE(init_pi(qE)), _piI(init_pi(qI)), _p0(p0), _p1(p1) 
{
}


const vector<double> SEIRFactor::init_pi( const vector<double> q) {

    auto res = vector<double>(q.size());
    double scale = 0;
    for( int i=q.size()-1; i>=0; i--) 
        scale += q[i];

    double sum = 0;
    for( int i=q.size()-1; i>=0; i--) {
        sum += q[i]/scale;
        res[i] = q[i]/scale/sum;
    }
    return res;
}

double SEIRFactor::potential( vector<vector<any>::const_iterator> state_its) {

    assert(state_its.size()==3);

    SEIRState in = std::any_cast<SEIRState>(*(state_its[0]));
    SEIRState out = std::any_cast<SEIRState>(*(state_its[1]));
    double load = std::any_cast<double>(*(state_its[2]));

    cout << "("<< in << " " << out << " " << load << ") => ";
    switch( in.phase()) {
        case SEIRState::S: {
            if( in.next(/*change = */ true)  == out) return 1.0 - (1.0-_p0)*pow(1.0-_p1, load);
            if( in.next(/*change = */ false) == out) return (1.0-_p0)*pow(1.0-_p1, load);
        } break;
        case SEIRState::E: {
            if( in.next(/*change = */ true)  == out) return _piE[in.days()];
            if( in.next(/*change = */ false) == out) return 1.0-_piE[in.days()];

        } break;
        case SEIRState::I: {
            if( in.next(/*change = */ true)  == out) return _piI[in.days()];
            if( in.next(/*change = */ false) == out) return 1.0-_piI[in.days()];
        } break;
        case SEIRState::R: {
            return out.phase()==SEIRState::R ? 1.0 : 0.0;
        } break;

    }
    return 0.0;
}


SEIRInitFactor::SEIRInitFactor( SEIRNode &out, bool patient_zero) :
    Factor({&out}), _patient_zero(patient_zero)
{
}

double SEIRInitFactor::potential( vector<vector<any>::const_iterator> state_its) {

    assert(state_its.size()==1);
    SEIRState out = std::any_cast<SEIRState>(*(state_its[0]));

    if( _patient_zero) 
        return out == SEIRState(SEIRState::E,1) ? 1.0 : 0.0;
    else 
        return out == SEIRState(SEIRState::S) ? 1.0 : 0.0;
}