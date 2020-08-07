#include <vector>
#include <iostream>
#include <map>
#include <cassert>
#include "seir.hpp"
#include "node.hpp"

using namespace std;

// Standard constructor
SEIRState::SEIRState(Phase p, int d) : _p(p), _d(d) {};

// Copy constructor
SEIRState::SEIRState(const SEIRState& other) : _p(other._p), _d(other._d) {};

// Returns the next SEIR state in the state space sequence
//  change == true:     S->E_1, E_m -> I_1, I_n -> R
//  change == false:    S->S, E_m -> E_{m+1}, I_n -> I_{n+1}, R -> R
SEIRState SEIRState::next(bool change) const {
    if(change) {
        switch(_p) {
            case S: return SEIRState(E,1);
            case E: return SEIRState(I,1);
            case I: return SEIRState(R);
            case R: throw std::invalid_argument("Phase R cannot be changed from!");
        }
    }
    else {
        switch(_p) {
            case S: return SEIRState(S);
            case E: return SEIRState(E,_d+1);
            case I: return SEIRState(I,_d+1);
            case R: return SEIRState(R);
        }
    }
    return SEIRState(S);
}

// Compares two detailed SEIR states
bool SEIRState::operator==(const SEIRState& other) const {
    if(_p==S && other._p==S) return true;
    if(_p==R && other._p==R) return true;
    if(_d != other._d) return false;
    if(_p==E && other._p==E) return true;
    if(_p==I && other._p==I) return true;

    return false;
}

// Compares a detailed SEIR state with the major SEIR states, i.e. phases (S, E, I, or R)
bool SEIRState::operator==(SEIRState::Phase other) const {
    return _p == other;
}

// Constructor for a whole state space of detailed SEIR states
SEIRStateSpace::SEIRStateSpace(int dE, int dI):
    dEMax(dE), dIMax(dI)
 {
    push_back(SEIRState::S);
    for(int d=1; d<=dE; push_back(SEIRState(SEIRState::E,d++)));
    for(int d=1; d<=dI; push_back(SEIRState(SEIRState::I,d++)));
    push_back(SEIRState::R);
}

// Returns the index of a detailed SEIR state when arranging them starting from S, to
// E_1,...,E_M,I_1,...,I_N,R (in other words, it's between 0 and dEMax+dIMax+1)
int SEIRStateSpace::operator[](const SEIRState& s) const {
    switch(s.phase()) {
        case SEIRState::S: return 0;
        case SEIRState::E: return s.days();
        case SEIRState::I: return s.days()+dEMax;
        case SEIRState::R: return 1+dEMax+dIMax;
    }
    // we should never get here
    return -1;
}

// Returns wether or not a detailed E or I state has reached the maximum possible days
bool SEIRStateSpace::can_continue(const SEIRState& s) const {
    switch(s.phase()) {
        case SEIRState::S: return true;
        case SEIRState::E: return s.days()<dEMax;
        case SEIRState::I: return s.days()<dIMax;
        case SEIRState::R: return true;
    }
    // we should never get here
    return false;
}

// Outputs a detailed SEIR state on a stream
ostream& operator<<(ostream& os, const SEIRState& s) {
    switch(s.phase()) {
        case SEIRState::S: os << "S"; break;
        case SEIRState::E: os << "E(" << s.days() << ")"; break;
        case SEIRState::I: os << "I(" << s.days() << ")"; break;
        case SEIRState::R: os << "R"; break;
    }
    return os;
}