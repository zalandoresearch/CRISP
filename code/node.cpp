#include <cassert>
#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"

// Compenentwise in-place multiplicaton of two messages
void operator*=(Message& a, const Message& b) {
    assert(a.size() == b.size());

    for(unsigned long i=0; i<a.size(); ++i) {
        a[i] *= b[i];
    } 
}

// Componentwise in-place division of two messages
void operator/=(Message& a, const Message& b) {
    assert(a.size() == b.size());

    for(unsigned long i=0; i<a.size(); ++i) {
        if(abs(b[i]) > 1e-10) {
            a[i] /= b[i];
        }
    }
    return;
}

// Normalizes a message and return the pointer to the normalized message
Message& normalize(Message &v) {
    double sum=0.0;
    for(auto x: v) sum += x;
    for(auto& x: v) x /= sum;
    return v;
}

// Normalizes a message and return the new normalized message
Message normalize(const Message &v_) {
    Message v(v_);
    double sum=0.0;
    for(auto x: v) sum += x;
    assert(sum > 1e-10);
    for(auto& x: v) x /= sum;
    return v;
}

// Standard constructor (with a uniform marginal)
Node::Node(int N): _N(N), _messages(), _marginal(N,1.0) {}

// Adds an associates factor to the variable node and creates a new message 
// from that factor to this variable node
void Node::add_factor(Factor* f) {
    _factors.push_back(f);
    _messages.emplace_back(nullptr);
}

// Computes an update of the messages from the i-th neighboring factor
void Node::update_factor(int i) {
    // If the message is not yet allocated the allocate and initialize with a uniform message
    if(!_messages[i]) {
        _messages[i] = MessagePtr(new Message(_N,1.0)); 
    }

    // (1) Remove the current message from the cached marginal. This is done in order to prevent a copy of the message
    //     to remove it before it is updated by the message_to_variable method of the associated i-th factor but it means
    //     that the message_to_variable method of the i-th factor MUST NEVER call the function normalized_message_to_factor because 
    //     until step (3), the member _marginal contains effectively the message to the factor already!
    _marginal /= *(_messages[i]);
    // (2) Delegate the computation of the message from the factor TO this variable node to the factor
    _factors[i]->message_to_variable(this, _messages[i]);
    // (3) Multiply back the new message to have the revised marginal cached in this variable node
    _marginal *= *(_messages[i]);
    // Update the message update counter
    Node::message_update_counter++;

    return;
}

// Computes an update of all messages from the neighboring factors (as well as the marginal)
void Node::update() {
    for(size_t i=0; i<_factors.size(); ++i) {
        update_factor(i);
    }
    return;
}

// Resets the variable node to a uniform distribution (both marginals and all messages)
void Node::reset() {
    fill(_marginal.begin(), _marginal.end(), 1.0);
    for(auto m: _messages) {
        if (m != nullptr) {
            fill(m->begin(), m->end(), 1.0);
        }
    }
    return;
}

// Computes the normalized message (i.e., probability) that the variable node sends to a neighboring factor
MessagePtr Node::normalized_message_to_factor(const Factor *f) const {
    MessagePtr m(new Message(_marginal));
    // TODO: This method performs a linear search - this is O(N) rather than O(log(N)) if we used a hashmap
    for(size_t i=0; i < _factors.size(); ++i) {
        if(_factors[i] == f) {
            if(_messages[i] != nullptr) {
               (*m) /= (*_messages[i]);
            }
        }
    }
    normalize(*m);
    return m;
}

// Initialize the message update counter to zero
int Node::message_update_counter = 0;

// Performs an update of all messages from neighboring factors according to the update policy
void SEIRNode::update(UpdatePolicy upd_policy) {
    for(size_t i=0; i < _factors.size(); ++i) {
        if (auto f = dynamic_cast<SEIRFactor*>(_factors[i])) {
            // UpdatePolicy::forward only updates the message from the factor of the next time step. 
            // NOTE THAT THIS implementation is tightly bound to the variable insertion order in SEIRFactor::init_helper! 
            if(upd_policy == forward && f->get_node(1) != this) continue;
            // UpdatePolicy::forward only updates the message from the factor of the previous time step. 
            // NOTE THAT THIS implementation is tightly bound to the variable insertion order in SEIRFactor::init_helper! 
            if(upd_policy == backward && f->get_node(0) != this) continue;
        }
        update_factor(i);
    }
}

// Computes the probability of being in state I_1,...,I_N or not to be in I_1,...,I_N
MessagePtr SEIRNode::normalized_infection_message_to_factor(const Factor *f) const {
    MessagePtr m = normalized_message_to_factor(f);
    double p_I = 0.0;
    for(SEIRState s=SEIRState(SEIRState::I,1); _states.can_continue(s); s = s.next(/* change = */false)) {
        p_I += (*m)[_states[s]];
    }  
    return MessagePtr(new Message({1.0-p_I, p_I}));
}

// Computes a new message over the 4 states (i.e., S,E,I,R) from a message over the 
// detailed state space (i.e, S,E_1,...,E_M,I_1,...,I_N,R) by summation
Message basic_states(const Message& message, const SEIRStateSpace& states) {
    assert(message.size() == states.size());

    Message out(4);
    for(size_t i=0; i<message.size(); ++i) {
        out[states[i].phase()] += message[i];
    }
    return out;
}
