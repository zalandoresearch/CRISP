#ifndef _NODE_HPP
#define _NODE_HPP

#include <memory>
#include <iostream>
#include <vector>
#include "seir.hpp"

using namespace std;

// A message is defined via the non-negative contributions for each state
typedef vector<double> Message;
// This shared pointer will be heavily used in the factor graph construction
typedef shared_ptr<Message> MessagePtr;

// Computes a new message over the 4 states (i.e., S,E,I,R) from a message over the 
// detailed state space (i.e, S,E_1,...,E_M,I_1,...,I_N,R) by summation
Message basic_states(const Message&, const SEIRStateSpace&);

// Outputs a vector on a stream in Matlab notation
template <class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
    os << "[";
    for(unsigned int i=0; i != v.size(); ++i) {
        os << v[i] << (i==v.size()-1 ? "":", ");
    }
    os << "],";
    return os;
}

// Normalize a message to become a probability distribution
Message& normalize(Message&);
// Normalize a message in place to become a probability distribution
Message normalize(const Message&);

class Factor;

// This class represents a variable node in the factor graph
class Node {
    Node(const Node&) = delete;
protected:
    vector<Factor*>     _factors;           // list of associated factors
    int                 _N;                 // length of message as discrete distribution
    vector<MessagePtr>  _messages;          // cache of messages to neighbouring factors (used for loopy belief propagation)
    Message             _marginal;          // cache of marginals of this variable
 
    // Computes an update of the messages from the i-th neighboring factor
    void update_factor(int);
public:
    Node(int);
    virtual void add_factor(Factor*);

    // Returns the number of states of the associates variable (e.g., the support of the discrete distribution)
    int size() const {return _N;};

    // Computes an update of all messages from the neighboring factors (as well as the marginal)
    void update();
    // Resets the variable node to a uniform distribution (both marginals and all messages)
    virtual void reset();

    // Computes the normalized message (i.e., probability) that the variable node sends to a neighboring factor
    MessagePtr normalized_message_to_factor(const Factor *f = 0) const;

    // Counts the total number of message updates
    static int message_update_counter;
};


// This class represents a variable node over a detailed SEIR state in the factor graph
class SEIRNode: public Node {
    const SEIRStateSpace& _states;
    using Node::update;
public:
    // The type of message passing schedule for this node in the factor graph
    enum UpdatePolicy { forward, backward, full };
    // Standard constructor for variable node over a detailed SEIR state in the factor graph
    SEIRNode(const SEIRStateSpace& all_states) : Node(all_states.size()), _states(all_states) {}
    virtual ~SEIRNode() {};

    // Performs an update of all messages from neighboring factors according to the update policy
    void update(UpdatePolicy);

    // Returns the list of all detailed SEIR states that this variable can take
    const SEIRStateSpace& states() { return _states; };

    // Computes the probability of being in state I_1,...,I_N or not to be in I_1,...,I_N
    MessagePtr normalized_infection_message_to_factor(const Factor *f=0) const;
};

#endif