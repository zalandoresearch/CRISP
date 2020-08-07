#ifndef _SEIR_HPP
#define _SEIR_HPP

#include <vector>
#include <map>

using namespace std;

// This class represents the extended SEIR state where each E and I state is 
// further equipped with a day counter (see Section 2.2 in the CRISP++ paper)
class SEIRState {
public:
    enum Phase { S=0, E=1, I=2, R=3};           // the standard four states of the SEIR model

    SEIRState(Phase p, int d=0);
    SEIRState(const SEIRState& other);

    SEIRState next(bool change=false) const;
    bool operator==(const SEIRState& other) const;
    bool operator==(Phase other) const;

    Phase phase() const { return _p;}
    int days() const { return _d;}

private:
    Phase   _p;                                 // Main phase (S, E, I, or R) of the state
    int     _d;                                 // 1-based number of days being in state E or I
};

ostream& operator<<(ostream& os, const SEIRState& s);

class SEIRStateSpace: public vector<SEIRState> {
    SEIRStateSpace(const SEIRStateSpace&) = delete; //avoid copying this object, one instance should always be sufficient
public:
    const int dEMax;
    const int dIMax;

    SEIRStateSpace(int dE, int dI);
    int operator[](const SEIRState& s) const;
    const SEIRState& operator[](int i) const {return vector<SEIRState>::operator[](i);};

    bool can_continue(const SEIRState& s) const;
};

#endif