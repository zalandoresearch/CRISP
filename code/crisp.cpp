#include "crisp.hpp"


std::ostream &operator<<(std::ostream &os, Contact const &c) {
    return os << "(" << c.getTargetIndividual() << "," << c.getSourceIndividual() << "," << c.getTime() << "," << c.getCount() <<")";
}
