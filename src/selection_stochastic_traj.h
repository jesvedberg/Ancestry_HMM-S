#ifndef __SELECTION_STOCHASTIC_TRAJECTORY_H
#define __SELECTION_STOCHASTIC_TRAJECTORY_H


// generates vector with allele frequency of selected allele over time/generations
void selection_stochastic_trajectory(vector<double> &trajectory, double s, double m, int generations, int ne, int reps) 
{
    double fixed_freq;
    vector<double> traj_sum(generations+1,0); //// +1??????
    default_random_engine rand_gen;

    double n_sel;
    double n_nonsel;
    double fix_freq;
    int fix_gen = 0;
    double ns_fit;
    double nns_fit;
    double ns_freq;
    double new_ns;

    
    for (int i = 0 ; i < reps ; i ++) {
        n_sel = ne * m;
        n_nonsel = ne * (1 - m);
    
        fixed_freq = 0;
        traj_sum[0] += m;

        for (int g = 1; g <= generations; g++) { //// not sure about start and end
            ns_fit = n_sel * (1 + s);
            nns_fit = n_nonsel;
            ns_freq = ns_fit / (ns_fit + nns_fit);

            if ( ns_freq == 0.0 && ns_freq == 1.0) {
                fix_freq = ns_freq;
                fix_gen = g;
                break;
            }

            binomial_distribution<> repopulate(ne, ns_freq);
            new_ns = repopulate(rand_gen);
            traj_sum[g] += new_ns/ne;

            n_sel = new_ns;
            n_nonsel = ne - new_ns;
        }
        if (fix_gen > 0) {
            for (int f = fix_gen; f <= generations; f++) {
                traj_sum[f] += fix_gen;
            }
        }
    }

    for (int i = 0; i < traj_sum.size(); i++) {
        trajectory.push_back(traj_sum[i]/reps);
    }
}

#endif