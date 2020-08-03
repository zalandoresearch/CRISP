import sys
import warnings
from typing import Union

from tqdm.auto import trange

sys.path.append('../..')

from crisp import Distribution, GibbsPIS
import argparse

from matplotlib.pyplot import *
from matplotlib import cycler

import numpy as np
import random


def init_contacts(S, T, qIbar=20.0, R0: Union[float, np.array] = 2.5, p1=0.01, decay=0.1,
                  R0_mit=(2.5, 0.5), t_mit=None, H=None, seed=42):

    random.seed(seed)
    np.random.seed(seed+1)

    if type(R0) is float:
        R0 = np.ones(T) * R0
    elif type(R0) is np.ndarray:
        assert len(R0) == T
    else:
        raise ValueError("parameter R0 must be float of np.array, was {}".format(type(R0)))

    # Precompute all contacts in the constructor
    contacts = {}
    l = np.arange(S)
    l0 = l[:,np.newaxis]
    l1 = l[np.newaxis,:]

    # a lower triangular binary mask of potention (unidirectional) contacts
    mask = (l[:,np.newaxis] > l[np.newaxis,:])
    idx = list(zip(*np.where(mask)))

    # if H is set, we have H local cliques
    if H is not None:
        maskb = mask * (l0 - l1 < l0 % H + 1) # the mask for in-clique contatcs
        maska = mask * (~maskb)               # the mask for inter-clique contacts
        pa = R0_mit[1] / qIbar / p1 / (S - H)
        pb = R0_mit[0] / qIbar / p1 / (H - 1)
        if pb > 1.0:
            warnings.warn("Mitigation results in decreased nominal R0, increase H to suppress this warning!")

        idxa = list(zip(*np.where(maska)))
        idxb = list(zip(*np.where(maskb)))

    def sample(idx, p0):
        N = len(idx)
        n = np.random.binomial(N, p0)
        c = np.array(random.sample(idx,n))
        c = np.c_[c, np.full_like(c[:,0], t), np.ones_like(c[:,0])]
        return np.r_[c, c[:, [1, 0, 2, 3]]]

    for t in trange(T):

        if t_mit is None or t<t_mit:
            p0 = R0[t] / qIbar / p1 / (S - 1)
            contacts[t] = sample(idx,p0)
        else:
            contacts[t] = np.r_[sample(idxa,pa),sample(idxb,pb)]

    return contacts


if __name__=="__main__":

    my_parser = argparse.ArgumentParser(description='Simulates testing and quarantining policies for COVID-19')
    my_parser.add_argument('--S', type=int, required=False, default=10000, help="The total number of individuals")
    my_parser.add_argument('--T', type=int, required=False, default=274, help="The total number of time steps")
    my_parser.add_argument('--p0', type=float, required=False, default=0.000001,
                           help="The probability of infection without contacts")
    my_parser.add_argument('--p1', type=float, required=False, default=0.01,
                           help="The probability of infection of a contact")
    my_parser.add_argument('--alpha', type=float, required=False, default=0.001,
                           help="The false negative rate of test I-test")
    my_parser.add_argument('--beta', type=float, required=False, default=0.01,
                           help="The false positive rate of the I-test")
    my_parser.add_argument('--R0', type=float, required=False, default=2.5, help="The R0 factor of COVID-19")
    my_parser.add_argument('--it', type=int, required=False, default=10, help="Numper of iterations to average over")
    my_parser.add_argument('--seed', type=int,  required=False, default=42, help="The random seed for contacts generation")
    args = my_parser.parse_args()

    T = args.T
    S = args.S
    alpha = args.alpha
    beta = args.beta
    p0 = args.p0
    p1 = args.p1
    R0 = args.R0

    It = args.it

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

    qE = Distribution([q/sum(qEVec) for q in qEVec])
    qI = Distribution([q/sum(qIVec) for q in qIVec])


    def make_figure(contacts, t_branch=None, contacts_branch=None):

        P = np.zeros((T,4))
        P_branch = np.zeros((T,4))


        for it in range(It):
            pis_branch = None
            for t in trange(T, desc="iteration {}".format(it)):
                if t==0:
                    pis = GibbsPIS(S, 1, contacts[t], [],
                                                     qE, qI,
                                                     alpha, beta,
                                                     p0, p1, True)
                    if t_branch is not None:
                        pis_branch = GibbsPIS(S, 1, contacts[t] if t<t_branch else contacts_branch[t], [],
                                                        qE, qI,
                                                        alpha, beta,
                                                        p0, p1, True)
                else:
                    pis.advance(contacts[t], [])
                    if t_branch is not None:
                        pis_branch.advance(contacts[t] if t<t_branch else contacts_branch[t], [])

                update = pis.get_infection_status().mean(0)
                P[t] += update
                if t_branch is not None:
                    update = pis_branch.get_infection_status().mean(0)
                P_branch[t] += update

        P /= It
        P_branch /= It

        fig = figure(figsize=(7.5, 4.5))

        ax = fig.gca()
        ax.set_prop_cycle( cycler(color=["orange","red","blue"]))
        for i in range(1, P.shape[1]):
            ax.plot(P[:, i]*S, linestyle='-', linewidth=2)
        if t_branch is not None:
            ax = fig.gca()
            ax.set_prop_cycle(cycler(color=["orange", "red", "blue"]))
            for i in range(1, P.shape[1]):
                ax.plot(np.arange(t_branch,T), P_branch[t_branch:, i] * S,
                        linestyle='--', linewidth=2)

        xlabel('days after patient 0 got infected')
        legend(['E', 'I', 'R'])
        grid(True)


        return fig


    contacts = init_contacts(S=S, T=T, R0=R0, p1=p1, seed=args.seed)
    fig_0 = make_figure(contacts)
    title('no mitigation')

    contacts_mit = init_contacts(S=S, T=T, R0=R0, p1=p1,  R0_mit=(R0-0.5,0.5), t_mit=60, H=20, seed=args.seed)
    fig_4 = make_figure(contacts_mit, t_branch=60, contacts_branch=contacts)
    fig_4.gca().axvline(x=60, color=[0.8, 0.8, 0.8], linestyle='--')
    title('mitigation with localized contact pattern')

    R0_1 = np.array([R0] * 60 + [1.0] * (T - 60))
    contacts_1 = init_contacts(S=S, T=T, R0=R0_1, seed=args.seed )
    fig_1 = make_figure(contacts_1)
    fig_1.gca().set_ylim([None, S*0.06])
    fig_1.gca().axvline(x=60, color=[0.8, 0.8, 0.8], linestyle='--')
    title('mitigation after 60 days')

    R0_2 = np.array([R0] * 60 + [0.5] * (T - 60))
    contacts_2 = init_contacts(S=S, T=T, R0=R0_2, seed=args.seed )
    fig_2 = make_figure(contacts_2)
    fig_2.gca().set_ylim(fig_1.gca().get_ylim())
    fig_2.gca().axvline(x=60, color=[0.8, 0.8, 0.8], linestyle='--')
    fig_2.gca().set_ylim(fig_2.gca().get_ylim())
    title('suppression after 60 days')

    R0_3 = np.array([R0] * 60 + [0.5] * 60 + [R0] * (T - 120))
    contacts_3 = init_contacts(S=S, T=T, R0=R0_3, seed=args.seed )
    fig_3 = make_figure(contacts_2, t_branch=120, contacts_branch=contacts_3)
    fig_3.gca().axvline(x=60, color=[0.8, 0.8, 0.8], linestyle='--')
    fig_3.gca().axvline(x=120, color=[0.8, 0.8, 0.8], linestyle='--')
    fig_3.gca().set_ylim(fig_1.gca().get_ylim())
    title('release after 60 days lockdown')

    fig_0.savefig('experiment51a.png')
    fig_1.savefig('experiment51b.png')
    fig_2.savefig('experiment51c.png')
    fig_3.savefig('experiment51d.png')
    fig_4.savefig('experiment51e.png')

