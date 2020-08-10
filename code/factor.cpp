#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <cassert>
#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"

using namespace std;

// Computes the message 'to' to the variable node 'n' using the belief propagation algorithm 
// Note that this implementation enumrates all possible states explicitly
void Factor::message_to_variable(Node* n, MessagePtr to) {
    // Initialize the message to be computed with all zeros
    fill(to->begin(), to->end(), 0.0);

    // Check that the node pointer passed in is actually a neighbor of the factor
    size_t i_node = 0;
    for(; i_node<_nodes.size() && _nodes[i_node] != n; ++i_node) {}
    assert(i_node < _nodes.size()); 


    vector<MessagePtr> incoming_messages;
    vector<Message::const_iterator> msg_its;
    vector<unsigned int> state_its(_nodes.size());

    // Caches all the messages from all incoming variable neighbors
    for(unsigned int i=0; i != _nodes.size(); ++i) {
        incoming_messages.push_back(_nodes[i]->normalized_message_to_factor(this));
        msg_its.push_back(incoming_messages[i]->begin());
    }

    MessagePtr outgoing_message = to;
    Message::iterator out_msg_it = outgoing_message->begin();

    // Loop over all possible state configurations of all connected neighbors (including the variable we sent to)
    while(true) {
        // (1) Compute the summand of factor value times the product of all incoming messages 
        //     (except from the variable we sent to)
        double p = 1.0;
        for(size_t i=0; i<_nodes.size(); ++i) {
            if(i != i_node) {
                p *= *(msg_its[i]);
            }
        }
        double upd = p * potential(state_its);
        *(out_msg_it) += upd;

        // (2) Generate the next state configuration
        for(int i=_nodes.size()-1; i>=0; --i) {
            msg_its[i]++;
            state_its[i]++;
            if(_nodes[i] == n) out_msg_it++;

            if(msg_its[i] == incoming_messages[i]->end()) {
                // if they overrun, wrap them over to the beginning
                msg_its[i]   = incoming_messages[i]->begin();
                state_its[i] = 0;
                if(_nodes[i] == n)
                    out_msg_it = outgoing_message->begin();

                // if the first iterator wrapped around we are done
                if( i==0) goto done;
            }
            else 
                break;
        }
    }
    done:
    return;
}

// Retrieves the value of function from a flattened table
double TableFactor::potential(const vector<unsigned int>& state_its) {
    assert(state_its.size()== _nodes.size());

    int idx = 0;
    for(unsigned int i=0; i<_nodes.size(); ++i) {
        if(i>0) idx *= _nodes[i-1]->size();
        idx += state_its[i];
    }
    return _tab[idx];
}

// Standard constructor of the main CRISP++ infection spread factor
SEIRFactor::SEIRFactor(const Distribution& qE,          // Distribution of patients being n days in E state (see (6) and (7))
                       const Distribution& qI,          // Distribution of patients being n days in I state (see (6) and (8))
                       double p0,                       // Exogeneous propability of infection spread per day (see (5))
                       double p1,                       // Probability of infection spread per contact (see (5))
                       SEIRNode& in,                    // Variable y_{u,t} in (10)
                       SEIRNode& out,                   // Variable y_{u,t+1} in (10)
                       vector<SEIRNode*> contacts) :    // Variables y_{v,t} in (10)
    Factor(init_helper(in,out,contacts)), 
    _piE(init_pi(qE)), 
    _piI(init_pi(qI)), 
    _p0(p0), 
    _p1(p1), 
    _states(in.states())
{
}

// Computes the complete set of neighboring factors for the main SEIR infection spread factor
vector<Node *> SEIRFactor::init_helper(SEIRNode& in, SEIRNode& out, vector<SEIRNode*> contacts) {
    vector<Node*> nodes({ &in, &out});
    nodes.insert(nodes.end(), contacts.begin(), contacts.end());
    return nodes;
}

// Computes the probablity of being n days in a state (E or I) from the logP and lopPTail
const vector<double> SEIRFactor::init_pi(const Distribution& q) {
    auto res = vector<double>(q.getMaxOutcomeValue()+1);
    for(unsigned int i=0; i<res.size(); ++i) {
        res[i] = exp(q.getLogP(i)-q.getLogPTail(i));
    }
    return res;
}

// Computes the forward message (see Subsection 3.2.1 in the CRISP++ paper)
//  Note that A=1 in this implementation due to message normalization
void SEIRFactor::message_forward( MessagePtr to) {
    MessagePtr output_message = to;
    MessagePtr input_message = _nodes[0]->normalized_message_to_factor(this);

    // This code computes (30) in the CRISP++ paper
    double B = (1.0-_p0);
    // Loop over all contacts (_nodes[0] = y_{u,t} and _nodes[1] = y_{u,t+1} so we star at _nodes[2])
    for(size_t i=2; i<_nodes.size(); ++i) {
        auto m = ((SEIRNode*)_nodes[i])->normalized_infection_message_to_factor(this);
        B *= m->at(0) + (1.0-_p1) * m->at(1);
    }

    // Compute the output message by looking at all state transitions from the input states
    auto it_input_message = input_message->cbegin();
    for(auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) {
        auto in = *it_input_states;
        const double& in_value = *it_input_message;

        switch(in.phase()) {
            case SEIRState::S:
                (*output_message)[_states[in.next(true)]] = in_value * (1.0 - B);                   // see (32)
                (*output_message)[_states[in.next(false)]] = in_value * B;                          // see (31)
                break;
            case SEIRState::E:
                (*output_message)[_states[in.next(true)]]  += in_value * _piE[in.days()];           // see (34)
                if(_states.can_continue(in)) {
                    (*output_message)[_states[in.next(false)]] = in_value * (1.0-_piE[in.days()]);  // see (33)
                }
                break;
            case SEIRState::I:
                (*output_message)[_states[in.next(true)]]  += in_value * _piI[in.days()];           // see (36)
                if( _states.can_continue(in)) {
                    (*output_message)[_states[in.next(false)]] = in_value * (1.0-_piI[in.days()]);  // see (35)
                }
                break;
            case SEIRState::R:
                (*output_message)[_states[in]] += in_value;                                         // see (36)
                break;
        }
        it_input_message++;
    }

    return;
}

// Computes the backward message (see Subsection 3.2.2 in the CRISP++ paper)
//  Note that A=1 in this implementation due to message normalization
void SEIRFactor::message_backward(MessagePtr to) {
    MessagePtr input_message = to;
    MessagePtr output_message = _nodes[1]->normalized_message_to_factor(this);

    // This code computes (30) in the CRISP++ paper
    double B = (1.0-_p0);
    // Loop over all contacts (_nodes[0] = y_{u,t} and _nodes[1] = y_{u,t+1} so we star at _nodes[2])
    for(size_t i=2; i<_nodes.size(); ++i) {
        auto m = ((SEIRNode*)_nodes[i])->normalized_infection_message_to_factor(this);
        B *= m->at(0) + (1.0-_p1) * m->at(1);
    }

    // Compute the output message by looking at all state transitions from the input states
    for(auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) {
        auto in = *it_input_states;

        switch(in.phase()) {
            case SEIRState::S:
                (*input_message)[_states[in]] = (*output_message)[_states[in.next(true)]] * (1.0 - B);
                (*input_message)[_states[in]] += (*output_message)[_states[in.next(false)]] * B;                // see (37)
                break;
            case SEIRState::E:
                (*input_message)[_states[in]] = (*output_message)[_states[in.next(true)]] * _piE[in.days()];    // see equation between (37) and (38)
                if( _states.can_continue(in)) {
                    (*input_message)[_states[in]] += (*output_message)[_states[in.next(false)]] * (1.0-_piE[in.days()]);    // see equation between (37) and (38)
                }
                break;
            case SEIRState::I:
                (*input_message)[_states[in]] = (*output_message)[_states[in.next(true)]] * _piI[in.days()];    // see equation between (38) and (39)
                if( _states.can_continue(in)) {
                    (*input_message)[_states[in]] += (*output_message)[_states[in.next(false)]] * (1.0-_piI[in.days()]);    // see equation between (38) and (39)
                }
                break;
            case SEIRState::R:
                (*input_message)[_states[in]] = (*output_message)[_states[in]];                                // see (40)
                break;
        }
    }

    return;
}

// Computes the message to a contact (see Subsection 3.2.3 in the CRISP++ paper)
//  Note that A=1 in this implementation due to message normalization
void SEIRFactor::message_vertical(Node* n, MessagePtr to) {
    auto input_message = _nodes[0]->normalized_message_to_factor(this);
    auto output_message = _nodes[1]->normalized_message_to_factor(this);

    // This code computes (42) in the CRISP++ paper
    double B_k = (1.0-_p0);
    for(size_t i=2; i<_nodes.size(); ++i) {
        if(_nodes[i] != n) {
            auto m = ((SEIRNode*)_nodes[i])->normalized_infection_message_to_factor(this);
            B_k *= m->at(0) + (1.0-_p1) * m->at(1);
        }
    }

    // Since all the messages are the same for the detailed states I_1,...,I_N and S,E_1,...,E_M,R,
    // we only compute these two values
    double p_I = 0.0;
    double p_N = 0.0;

    for(auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) {
        auto in = *it_input_states;
        double p;

        switch(in.phase()) {
            case SEIRState::S:
                p_N += (*input_message)[_states[in]] * (*output_message)[_states[in.next(true)]] * (1.0 - B_k);     // see (44)
                p_N += (*input_message)[_states[in]] * (*output_message)[_states[in.next(false)]] * B_k;            // see (44)
                p_I += (*input_message)[_states[in]] * (*output_message)[_states[in.next(true)]] * (1.0 - B_k * (1.0 - _p1));   // see (43)
                p_I += (*input_message)[_states[in]] * (*output_message)[_states[in.next(false)]] * B_k * (1.0 - _p1);          // see (43)
                break;
            case SEIRState::E:
                p =  (*input_message)[_states[in]] * (*output_message)[_states[in.next(true)]] * _piE[in.days()];           // see (43) and (44)
                p += (*input_message)[_states[in]] * (*output_message)[_states[in.next(false)]] * (1.0 - _piE[in.days()]);  // see (43) and (44)
                p_N += p;
                p_I += p;
                break;
            case SEIRState::I:
                p =  (*input_message)[_states[in]] * (*output_message)[_states[in.next(true)]] * _piI[in.days()];           // see (43) and (44)
                p += (*input_message)[_states[in]] * (*output_message)[_states[in.next(false)]] * (1.0 - _piI[in.days()]);  // see (43) and (44)
                p_N += p;
                p_I += p;
                break;
            case SEIRState::R:
                p = (*input_message)[_states[in]] * (*output_message)[_states[in.next(false)]];     // see (43) and (44)
                p_N += p;
                p_I += p;
                break;
        }
    }

    // Finally, copy the two states into the detailed output message
    MessagePtr outgoing_message = to;
    for(int i=0; i<n->size(); ++i) {
        (*outgoing_message)[i] = (_states[i].phase()==SEIRState::I ? p_I : p_N);
    }

    return;
}

// Computes the message to the variable node using the belief propagation algorithm in Section 3.2 of the CRISP++ paper
void SEIRFactor::message_to_variable(Node* n, MessagePtr to) {
    fill(to->begin(), to->end(), 0.0);
    if(n==_nodes[0]) return message_backward(to);
    if(n==_nodes[1]) return message_forward(to);
    return message_vertical(n, to);
}

// Standard constructor
SEIRInitFactor::SEIRInitFactor(SEIRNode& out, bool patient_zero, double p0) :
    Factor({&out}), 
    _patient_zero(patient_zero), 
    _p0(p0), 
    _states(out.states()) { }

// Computes the prior message taking into account the person being a patient-zero
void SEIRInitFactor::message_to_variable(Node*, MessagePtr to) {
    fill(to->begin(), to->end(), 0.0);
    if(_patient_zero) {
        (*to)[_states[SEIRState(SEIRState::E,1)]] = 1.0;
    }
    else {
        (*to)[_states[SEIRState(SEIRState::S)]] = 1.0-_p0;
        (*to)[_states[SEIRState(SEIRState::E,1)]] = _p0;
    }
    return;
}

// Standard constructor
SEIRTestFactor::SEIRTestFactor(SEIRNode& out, bool positive, double alpha, double beta) :
    Factor({&out}), 
    _positive(positive),
    _alpha(alpha),
    _beta(beta),
    _states(out.states()) { }

// The value of the factor for all detailed SEIR states
double SEIRTestFactor::potential(const vector<unsigned int>& state_its) {
    SEIRState out = _states[state_its[0]];

    if( _positive) {
        return (out.phase() == SEIRState::I) ? 1.0-_alpha : _beta;
    }
    else {
        return (out.phase() == SEIRState::I) ? _alpha : 1.0-_beta;
    }

    // We should never get here
    return 0.0;
}
