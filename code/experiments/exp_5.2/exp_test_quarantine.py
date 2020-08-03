import numpy as np
from math import exp, log
from enum import Enum
from matplotlib.pyplot import *
from crisp import Distribution, LBPPIS, GibbsPIS
import itertools
import argparse
import csv


class InfectionState(Enum):
    SUSCEPTIBLE = 0
    EXPOSED     = 1
    INFECTIOUS  = 2
    RECOVERED   = 3

class CRISP():
    def __init__(self, S, T, qE, qI, alpha, beta, p0, p1, p_symptomatic):
        self.S = S                              # Total number of individuals
        self.T = T                              # Total number of time steps
        self.qE = qE                            # Distribution of the state of exposed
        self.qI = qI                            # Distribution of the state of infectious
        self.alpha = alpha                      # False-negative rate of the I-test
        self.beta = beta                        # False-positive rate of the I-test
        self.p0 = p0                            # Exogenous probability of infection
        self.p1 = p1                            # Infection probability due to a contact
        self.p_symptomatic = p_symptomatic      # The probability of showing symptoms
        self.t = 0                              # Time step we are currently in
 
        self.pis = None                         # PopulationInfectionStatus object
        self.P_zu = None                        # The actual infection status for everyone at the latest time step
        self.non_infectious_individuals = np.full((self.S,),True)   # Caches all the non-infectious individuals to detect an onset of symptoms

    # get infection statistics from the world mode;
    def stats(self, individuals=None):
        if self.P_zu is None:
            return np.array([0., 0., 0., 0.])
        res = self.P_zu.sum(axis=0) if individuals is None else self.P_zu[individuals].sum(axis=0)
        return res

    # Samples one more step forward
    def advance(self, contacts):
        # Increase time step
        self.t += 1

        if self.pis is None:
            self.pis = GibbsPIS(self.S, 1, contacts, [],
                                Distribution(self.qE), Distribution(self.qI),
                                self.alpha, self.beta,
                                self.p0, self.p1, True)
        else:
            self.pis.advance(contacts, [])

        self.P_zu = self.pis.get_infection_status()
        infectious_individuals = (self.P_zu[:,2]>0.999)
        symptomatic = (np.random.rand(self.S) > self.p_symptomatic)
        individuals_with_symptom_onset = np.where(infectious_individuals & self.non_infectious_individuals & symptomatic)[0]
        self.non_infectious_individuals &= ~infectious_individuals

        return individuals_with_symptom_onset

    # Computes the test outcomes for a given set of test candidates
    def sample_test_outcomes(self, test_candidates):
        if len(test_candidates) == 0:
            return []

        p_positive = np.array([[self.beta, self.beta, 1.0 - self.alpha, self.beta]])
        p_positive = (self.P_zu[test_candidates] * p_positive).sum(axis=1)
        test_outcomes = (np.random.rand(*p_positive.shape) < p_positive).astype(np.int)
        test_times = np.full_like(test_candidates, self.t-1)

        return np.c_[test_candidates, test_times, test_outcomes]

    # Returns the infection trace of all individuals
    def get_world_state(self):
        return self.pis.get_individual_traces()

class Contacts():
    def __init__(self, S, T, qE, qI, p1, R_0 = 2.5):
        # Expected number of days of infectiousness
        qIBar = sum(np.arange(len(qI))*qI)
        self.no_contacts_per_day = R_0/(qIBar * p1)
        print("Average number of contacts per day = {:.3}".format(self.no_contacts_per_day))

        # Precompute all contacts in the constructor
        self.contacts = {}
        l = np.arange(S)
        mask = (l[:,np.newaxis] > l[np.newaxis,:])
        p0 = R_0 / qIBar / p1 / (S - 1)
        p = mask * p0

        for t in range(T):
            c = np.where(np.random.rand(S,S) < p)
            c = np.c_[c[0],c[1],np.full_like(c[0],t),np.ones_like(c[0])]

            self.contacts[t] = np.r_[c,c[:,[1,0,2,3]]]

    # Returns the contact matrix for time step t filtering out the quarantined individuals
    def get_contacts(self, t, quarantined_individuals):
        return [(u,v,t,x) for (u,v,t,x) in self.contacts[t] if not(u in quarantined_individuals) and not(v in quarantined_individuals)] 

class PolicyEvaluator():
    def __init__(self, S, T, qE, qI, alpha, beta, p0, p1, contacts, p_symptomatic = 0.5, policy_start = 30, write_world = False):
        self.crisp = CRISP(S, T, qE, qI, alpha, beta, p0, p1, p_symptomatic)
        self.S = S                                              # Total number of individuals
        self.T = T                                              # Total number of time steps
        self.contacts = contacts                                # A fixed contact generator
        self.policy_start = policy_start                        # The day that the testing & quarantine policy starts
        self.write_world = write_world                          # Stores whether or not the world state is written to disc
        self.quarantine_stats = np.array([[0] * 4] * T)         # Stores the statistics over the infection state of quarantined people each day
        self.test_stats = np.array([[0] * 2] * T)               # Stores the total number of tests and the total number of positive test outcomes

    # Compute the infection score for everyone as a probability distribution over state S,E,I,R at the beginning of the day
    def get_infection_score(self):
        raise NotImplementedError

    # Computes the candidates that should be tested based on the infection score
    def compute_test_candidates(self, infection_score):
        raise NotImplementedError

    # Computes the individuals that should be quarantined based on the infection score
    def compute_quarantine_status(self, infection_score):
        raise NotImplementedError

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        raise NotImplementedError

    # Writes a list as a CSV file
    def write_csv_file(self, l, filename):
        if (self.write_world):
            with open(filename, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(l)

    # Samples one entire infection trace for every individual
    def evaluate(self):
        # Reset all the stats
        self.quarantine_stats = np.array([[0] * 4] * self.T)

        # Reset the testing statistics
        self.test_stats = np.array([[0] * 2] * self.T)

        for t in range(self.T-1):
            print("Time step: t = {}".format(t+1))

            # At the start of the day ...

            if (t+1 > self.policy_start):
                # ... (1) get the infection score for every individual (i.e., an S x 4 matrix of probabilities over P(z_{ut} \in \{S,E,I,R\}))
                infection_score = self.get_infection_score()

                # ... (2) compute both the test candidates and the individuals to be quarantind based on the infection score
                test_candidates = self.compute_test_candidates(infection_score)
                quarantined_individuals = self.compute_quarantine_status(infection_score)
            else:
                test_candidates = []
                quarantined_individuals = []
            
            # ... (3) count one day in quarantine for everyone in quarantine
            self.quarantine_stats[t+1] = self.crisp.stats(quarantined_individuals)

            # ... (4) compute the new contacts based on who is in quarantine
            today_contacts = self.contacts.get_contacts(t, quarantined_individuals)
            self.write_csv_file(today_contacts, "contacts{}.csv".format(t))

            # ... (5) advance the world model by one step based on the new contacts
            individuals_with_symptom_onset = self.crisp.advance(today_contacts)

            # ... (6) compute the test outcomes of the test candidates based on the advanced world model
            today_test_outcomes = self.crisp.sample_test_outcomes(test_candidates)
            self.write_csv_file(today_test_outcomes, "test{}.csv".format(t))

            # if (t+1 >= self.policy_start):
            # ... (7) and finally advance the infection score computation
            self.advance_infection_score_model(today_contacts, today_test_outcomes, individuals_with_symptom_onset)

            # Update the test stats
            self.test_stats[t+1] = [len(today_test_outcomes),sum([1 for (u,t,o) in today_test_outcomes if (o==1)])]


## Implements a strategy of no tests and no quarantining
class NoPolicy(PolicyEvaluator):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # Compute the infection score for everyone as a probability distribution over state S,E,I,R at the beginning of the day
    def get_infection_score(self):
        return []

    # Computes the candidates that should be tested based on the infection score
    def compute_test_candidates(self, infection_score):
        return []

    # Computes the individuals that should be quarantined based on the infection score
    def compute_quarantine_status(self, infection_score):
        return []

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        pass

## Implements a testing strategy where each person is locked down
class LockdownPolicy(PolicyEvaluator):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # Compute the infection score for everyone as a probability distribution over state S,E,I,R at the beginning of the day
    def get_infection_score(self):
        return []

    # Computes the candidates that should be tested based on the infection score
    def compute_test_candidates(self, infection_score):
        return []

    # Computes the individuals that should be quarantined based on the infection score
    def compute_quarantine_status(self, infection_score):
        return np.arange(S)

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        pass

## Implements a testing strategy where each person who has symptom onset is tested and all positively test people are quarantined
class SymptomPolicy(PolicyEvaluator):
    def __init__(self, quarantine_days=14, test_capacity=0.01, **kwargs):
        super().__init__(**kwargs)
        self.daily_test_capcity = int(test_capacity * self.S)   # The capacity of daily testing as a percentage of number of total individuals
        self.quarantine_days = quarantine_days                  # The number of days that people have to stay in quarantine
        
        self.individuals_with_symptom_onset = []                # The people who had a symptom onset
        self.days_in_quarantine = np.array([0]*self.S)          # Stores how many days positively tested people are still in quarantine

        print("Daily test capacity = {} tests".format(self.daily_test_capcity))

    # Compute the infection score for everyone as a probability distribution over state S,E,I,R at the beginning of the day
    def get_infection_score(self):
        return []

    # Computes the candidates that should be tested based on the infection score
    def compute_test_candidates(self, infection_score):
        if len(self.individuals_with_symptom_onset) <= self.daily_test_capcity:
            return self.individuals_with_symptom_onset
        else:
            return np.random.choice(self.individuals_with_symptom_onset, self.daily_test_capcity, replace=False)

    # Computes the individuals that should be quarantined based on the infection score
    def compute_quarantine_status(self, infection_score):
        quarantined_individuals = np.where(self.days_in_quarantine > 0)[0]
        self.days_in_quarantine[quarantined_individuals] -= 1
        return quarantined_individuals

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        # Store everyone with symptom onsets as future test candidate
        self.individuals_with_symptom_onset = individuals_with_symptom_onset

        # Store everyone who tested positive as a quarantine candidate
        for (u,_,o) in new_test_outcomes:
            if (o == 1):
                self.days_in_quarantine[u] = self.quarantine_days

# Implements a contact tracing strategy where each person who has symptom onset is tested and all their contacts are 
# prioritized and tested
class ContactTracingPolicy(PolicyEvaluator):
    def __init__(self, quarantine_days=14, test_capacity=0.01, past_days=7, **kwargs):
        super().__init__(**kwargs)
        self.daily_test_capcity = int(test_capacity * self.S)   # The capacity of daily testing as a percentage of number of total individuals
        self.quarantine_days = quarantine_days                  # The number of days that people have to stay in quarantine
        self.past_days = past_days                              # Stores the number of days in the past

        self.past_contacts = {}                                 # Stores the contacts of the past self.past_days many days
        self.idx = 0                                            # Pointer to the latest contact matrix
        self.test_outcomes = []                                 # Stores the new test outcomes needed for computing the infection score
        self.individuals_with_symptom_onset = []                # The people who had a symptom onset

        self.days_in_quarantine = np.array([0]*self.S)          # Stores how many days an individual is still in quarantine
        self.positive_outcome_contacts = np.array([0] * self.S) # Stores how many infectious contacts individual u had

        print("Daily test capacity = {} tests".format(self.daily_test_capcity))

    # Compute the infection score for everyone; since we use member variables, we do not need the infection score
    def get_infection_score(self):
        for (u,_,o) in self.test_outcomes:
            if (o == 1):
                # If a person tests positive, they will be quarantined ...
                self.days_in_quarantine[u] = self.quarantine_days
                # ... and all the contacts of the past days will also be quarantined, if they are not already
                for t in range(self.past_days):
                    # ... also update the stats that counts the number of positively tested people that someone has 
                    # been in contact with (that will be important for prioritizing testing)
                    for (w, v, _, _) in self.past_contacts[(self.idx - t)%self.past_days]:
                         if (w == u):
                            self.positive_outcome_contacts[v] += 1
                            # ... and then quarantine the person if they were not already quarantined
                            if (self.days_in_quarantine[v] <= 0):
                                self.days_in_quarantine[v] = self.quarantine_days  
            else:
                # If a person test negative, definitely remove them from quarantine ...
                self.days_in_quarantine[u] = 0
                # ... and set the number of positive test outcome contacts to zero (because they definitely did not infect the person)
                self.positive_outcome_contacts[u] = 0
        return []

    # Computes the candidates that should be tested based on the infection score
    def compute_test_candidates(self, infection_score):
        if (len(self.individuals_with_symptom_onset) <= self.daily_test_capcity):
            quarantined_individuals = np.where(self.days_in_quarantine > 0)[0]
            prioritied_quarantined_list = sorted(quarantined_individuals, key=lambda u:-self.positive_outcome_contacts[u])
            return (self.individuals_with_symptom_onset + prioritied_quarantined_list)[:self.daily_test_capcity]
        else:
            return np.random.choice(self.individuals_with_symptom_onset, self.daily_test_capcity, replace=False)

    # Computes the individuals that should be quarantined based on the infection score
    def compute_quarantine_status(self, infection_score):
        quarantined_individuals = np.where(self.days_in_quarantine > 0)[0]
        self.days_in_quarantine[quarantined_individuals] -= 1
        return quarantined_individuals

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        # Store everyone with symptom onsets as future test candidate
        self.individuals_with_symptom_onset = list(individuals_with_symptom_onset)

        # Stores the contact matrix for this time step
        self.idx = (self.idx + 1) % self.past_days
        self.past_contacts[self.idx] = new_contacts

        # Copy the test outcomes for the infection score computation
        self.test_outcomes = new_test_outcomes

# Implements the Gibbs sampling infection score estimation where people get prioritized in testing by P(z_{ut}=I),
# get into quarantine based on P(z_{ut}=E) and out of it by P(z_{ut}=S) or P(z_{ut}=R)
class GibbsScoringPolicy(PolicyEvaluator):
    def __init__(self, test_capacity=0.01, pEI_threshold=0.3, pSR_threshold=0.5, **kwargs):
        super().__init__(**kwargs)
        self.daily_test_capcity = int(test_capacity * self.S)   # The capacity of daily testing as a percentage of number of total individuals
        self.pEI_threshold = pEI_threshold                      # The threshold on P(z_{ut}=E) or P(z_{ut}=I) where people get into quarantine
        self.pSR_threshold = pSR_threshold                      # The threshold on P(z_{ut}=S) or P(z_{ut}=R) where people get out off quarantine

        self.positive_test_outcome = np.array([False] * self.S) # Indicates whether or not an individual already had a positive test outcome
        self.quarantined = np.array([False] * self.S)           # The individuals who are currently quanrantined
        self.individuals_with_symptom_onset = []                # The people who had a symptom onset

        self.scorer = None                                      # The Gibbs sampler used for scoring

        print("Daily test capacity = {} tests".format(self.daily_test_capcity))

    # Compute the infection score for everyone; since we use member variables, we do not need the infection score
    def get_infection_score(self):
        return self.scorer.get_infection_status(N=100)

    # Computes the candidates that should be tested based on the infection score
    def compute_test_candidates(self, infection_score):
        if (len(self.individuals_with_symptom_onset) <= self.daily_test_capcity):
            # Only consider individuals for additional testing who do not have symptoms and who did not test positive before
            candidates = list(set(np.where(self.positive_test_outcome == False)[0]) - 
                              set(self.individuals_with_symptom_onset))
            # Sort all these individuals by P(z_{ut}=I) descending
            prioritied_list = sorted(candidates, key=lambda u:-infection_score[u,2])
            ret = (self.individuals_with_symptom_onset + prioritied_list)[:self.daily_test_capcity]
            print("Generated test candidates: {}".format(ret))
            return ret
        else:
            return np.random.choice(self.individuals_with_symptom_onset, self.daily_test_capcity, replace=False)

    # Computes the individuals that should be quarantined based on the infection score
    def compute_quarantine_status(self, infection_score):
        # Decide on releasing and enforcing the quarantine based on P(z_{ut} \in \{E,I\}) or P(z_{ut} \ in \{S,R\})
        infected_score = infection_score[:,1] + infection_score[:,2]
        release_quarantine = list(np.where((self.quarantined == True) & (infected_score < (1-self.pSR_threshold)))[0])
        enforce_quarantine = list(np.where((self.quarantined == False) & (infected_score > self.pEI_threshold))[0])

        # Change the quarantine status
        self.quarantined[release_quarantine] = False
        self.quarantined[enforce_quarantine] = True

        print("Generated quarantine candidates: {}".format(np.where(self.quarantined == True)[0]))

        return np.where(self.quarantined == True)[0]

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        # Store everyone with symptom onsets as future test candidate; they will definitely be tested
        self.individuals_with_symptom_onset = list(individuals_with_symptom_onset)
        # Get everyone who has tested positive
        positive_tested_individuals = [u for (u,_,o) in new_test_outcomes if (o==1)]
        self.positive_test_outcome[positive_tested_individuals] = True

        print("Found individuals with symptoms onset: {} \nFound individuals with positive test outcome: {}".format(self.individuals_with_symptom_onset,positive_tested_individuals))

        # Advance the Gibbs sampling-based scoring by the sparse array of new contacts and new test outcomes
        if self.scorer is None:
            self.scorer = GibbsPIS(self.S, 1, new_contacts, new_test_outcomes,
                                   Distribution(self.crisp.qE), Distribution(self.crisp.qI),
                                   self.crisp.alpha, self.crisp.beta,
                                   self.crisp.p0*10, self.crisp.p1, False)
        else:
            self.scorer.advance(new_contacts, new_test_outcomes)

class LBPScoringPolicy(GibbsScoringPolicy):

    # Advances the infection score model based on one new contact matrix, test outcomes and people with symptom onset
    def advance_infection_score_model(self, new_contacts, new_test_outcomes, individuals_with_symptom_onset):
        # Store everyone with symptom onsets as future test candidate; they will definitely be tested
        self.individuals_with_symptom_onset = list(individuals_with_symptom_onset)
        # Get everyone who has tested positive
        positive_tested_individuals = [u for (u,_,o) in new_test_outcomes if (o==1)]
        self.positive_test_outcome[positive_tested_individuals] = True

        print("Found individuals with symptoms onset: {} \nFound individuals with positive test outcome: {}".format(self.individuals_with_symptom_onset,positive_tested_individuals))

        # Advance the Gibbs sampling-based scoring by the sparse array of new contacts and new test outcomes
        if self.scorer is None:
            self.scorer = LBPPIS(self.S, 1, new_contacts, new_test_outcomes,
                                 Distribution(self.crisp.qE), Distribution(self.crisp.qI),
                                 self.crisp.alpha, self.crisp.beta,
                                 self.crisp.p0*10, self.crisp.p1, False)
        else:
            self.scorer.advance(new_contacts, new_test_outcomes)



def argument_parser() :
    # Create the parser
    my_parser = argparse.ArgumentParser(description='Simulates testing and quarantining policies for COVID-19')
    my_parser.add_argument('--S',               type=int,  required=False, default=1000,   help="The total number of individuals")
    my_parser.add_argument('--T',               type=int,  required=False, default=150,    help="The total number of time steps")
    my_parser.add_argument('--p0',              type=float,required=False, default=1e-4,   help="The probability of infection without contacts")
    my_parser.add_argument('--p1',              type=float,required=False, default=0.025,  help="The probability of infection of a contact")
    my_parser.add_argument('--alpha',           type=float,required=False, default=0.001,  help="The false negative rate of test I-test")
    my_parser.add_argument('--beta',            type=float,required=False, default=0.01,   help="The false positive rate of the I-test")
    my_parser.add_argument('--R0',              type=float,required=False, default=2.5,    help="The R0 factor of COVID-19")
    my_parser.add_argument('--policy-start',    type=int,  required=False, default=30,     help="The day that the policy starts")
    my_parser.add_argument('--test-capacity',   type=float,required=False, default=0.01,   help="The daily capacity of tests")
    my_parser.add_argument('--qDays',           type=int,  required=False, default=14,     help="The number of days a predictive infected is quarantined")
    my_parser.add_argument('--pDays',           type=int,  required=False, default=5,      help="The number of past days to consider in contact tracing")
    my_parser.add_argument('--pEI-threshold',   type=float,required=False, default=0.3,    help="The threshold on P(z_{ut}=E) or P(z_{ut}=I) above of which a person goes into quarantine")
    my_parser.add_argument('--pSR-threshold',   type=float,required=False, default=0.7,    help="The threshold on P(z_{ut}=S) or P(z_{ut}=R) above of which a person goes out of quarantine")
    my_parser.add_argument('--seed',            type=int,  required=False, default=42,     help="The random seed for contacts generation")
    my_parser.add_argument('--iteration',       type=int,  required=False, default=0,      help="An iteration counter (will be ignored but allows to run for the same random seed for contacts generation)")
    my_parser.add_argument('--interactive',     choices=['on', 'off'], default='on',       help="Run the policy evaluation interactively or with file output")
    my_parser.add_argument('--save-world',      choices=['on', 'off'], default='off',      help="Save the world state to files for each time-step")
    my_parser.add_argument('--policy',          choices=['no','lockdown','symptom','contact','score','lbp'], default='contact', help="The policy to be evaluated")

    return my_parser

if __name__ == "__main__":

    my_parser = argument_parser()
    args = my_parser.parse_args()

    T = args.T
    S = args.S
    alpha = args.alpha
    beta = args.beta
    p0 = args.p0
    p1 = args.p1
    R0 = args.R0
    policy_start = args.policy_start
    test_capacity = args.test_capacity
    write_world = (args.save_world == 'on')

    # Initialize the random seed
    np.random.seed(args.seed)

    # The discrete distributions of the duration of exposure and infectiouness
    qEVec = [0.0000000000, 0.05908981283, 0.1656874653, 0.1819578343, 0.154807057,
             0.1198776096, 0.08938884645, 0.06572939883, 0.04819654533,
             0.03543733758, 0.02620080839, 0.01950646727, 0.01463254844, 
             0.0110616426, 0.008426626119]

    qIVec = [0.000000000000, 0.000000000000, 0.00000000000, 0.000000000000, 0.000000000000,
             0.0001178655952, 0.0006658439543, 0.002319264193, 0.005825713197, 0.01160465163,
             0.01949056696, 0.02877007836, 0.03842711373, 0.04743309657, 0.05496446107,
             0.06050719418, 0.06386313651, 0.065094874, 0.06444537162, 0.06225794729,
             0.0589104177, 0.05476817903, 0.05015542853, 0.0453410888, 0.04053528452,
             0.03589255717, 0.03151878504, 0.02747963753, 0.02380914891, 0.02051758911,
             0.01759822872, 0.01503287457, 0.0127962154, 0.01085910889, 0.009190974483,
             0.007761463001, 0.006541562648, 0.005504277076]

    qEVec = [q/sum(qEVec) for q in qEVec]
    qIVec = [q/sum(qIVec) for q in qIVec]

    contacts = Contacts(S=S, T=T, qE=qEVec, qI=qIVec, p1=p1, R_0 = R0)
    if (args.policy=='no'):
        policy = NoPolicy(S=S, T=T, qE=qEVec, qI=qIVec, alpha=alpha, beta=beta, p0=p0, p1=p1, policy_start=policy_start, contacts=contacts, write_world=write_world)
        fig_title = 'Evaluation of No Policy'
    elif (args.policy=='lockdown'):
        policy = LockdownPolicy(S=S, T=T, qE=qEVec, qI=qIVec, alpha=alpha, beta=beta, p0=p0, p1=p1, policy_start=policy_start, contacts=contacts, write_world=write_world)
        fig_title = 'Evaluation of Lockdown Policy'
    elif (args.policy=='symptom'):
        policy = SymptomPolicy(S=S, T=T, qE=qEVec, qI=qIVec, alpha=alpha, beta=beta, p0=p0, p1=p1, policy_start=policy_start, contacts=contacts, test_capacity=test_capacity, quarantine_days=args.qDays, write_world=write_world)
        fig_title = 'Evaluation of Symptom Policy'
    elif (args.policy=='contact'):
        policy = ContactTracingPolicy(S=S, T=T, qE=qEVec, qI=qIVec, alpha=alpha, beta=beta, p0=p0, p1=p1, policy_start=policy_start, contacts=contacts, test_capacity=test_capacity, quarantine_days=args.qDays, past_days=args.pDays, write_world=write_world)
        fig_title = 'Evaluation of Contact Tracing Policy'
    elif (args.policy=='score'):
        policy = GibbsScoringPolicy(S=S, T=T, qE=qEVec, qI=qIVec, alpha=alpha, beta=beta, p0=p0, p1=p1, policy_start=policy_start, contacts=contacts, test_capacity=test_capacity, pEI_threshold=args.pEI_threshold, pSR_threshold=args.pSR_threshold, write_world=write_world)
        fig_title = 'Evaluation of Gibbs Score Policy'
    elif (args.policy=='lbp'):
        policy = LBPScoringPolicy(S=S, T=T, qE=qEVec, qI=qIVec, alpha=alpha, beta=beta, p0=p0, p1=p1, policy_start=policy_start, contacts=contacts, test_capacity=test_capacity, pEI_threshold=args.pEI_threshold, pSR_threshold=args.pSR_threshold)
        fig_title = 'Evaluation of LBP Score Policy'
    else:
        raise ValueError("unknown policy: {}".format(args.policy))

    policy.evaluate()
    Z = policy.crisp.get_world_state()

    # Plot the infection curves for the whole population
    fig, (ax1, ax2) = subplots(2,1, figsize = [8,10])
    p = (np.cumsum(Z,axis=1)[:,np.newaxis] >= np.arange(T).reshape(1,-1,1)).sum(axis=0); 
    p = np.diff(p,axis=1,prepend=0); 
    p = np.concatenate([p,S-p.sum(1,keepdims=True)],axis=1)
    ax1.plot(p,'.-')
    ax1.set_title('Infection Status over Time')
    ax1.set_xlabel('Days')
    ax1.set_ylabel('Number of Individuals')
    ax1.grid(True)

    # Plot the quarantine curves for the whole population
    bars = np.cumsum(policy.quarantine_stats,axis=1)
    ax2.bar(np.arange(T),policy.quarantine_stats[:,0],color="blue",edgecolor="white",width=1)
    ax2.bar(np.arange(T),policy.quarantine_stats[:,1],bottom=bars[:,0],color="orange",edgecolor="white",width=1)
    ax2.bar(np.arange(T),policy.quarantine_stats[:,2],bottom=bars[:,1],color="green",edgecolor="white",width=1)
    ax2.bar(np.arange(T),policy.quarantine_stats[:,3],bottom=bars[:,2],color="red",edgecolor="white",width=1)
    ax2.set_title('Quarantine Status')
    ax2.set_xlabel('Days')
    ax2.set_ylabel('Number of Individuals')


    fig.suptitle(fig_title, fontsize=16)

    if (args.interactive == 'on'):
        print("Total number of test: {} with {:5.2f}% positive outcomes".format(policy.test_stats.sum(axis=0)[0], policy.test_stats.mean(axis=0)[1]/policy.test_stats.mean(axis=0)[0]*100))
        show()
    else:
        np.savetxt("traces.csv", Z, delimiter=",", fmt="%d")
        np.savetxt("qstats.csv", policy.quarantine_stats, delimiter=",", fmt="%d")
        np.savetxt("tstats.csv", policy.test_stats, delimiter=",", fmt="%d")
        fig.savefig("output.png")