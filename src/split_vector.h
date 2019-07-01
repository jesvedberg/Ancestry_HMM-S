#ifndef __SPLIT_VECTOR_H
#define __SPLIT_VECTOR_H

/// splits vector (of recombination rates) into two vectors at the selected site.
/// the back vector is generated in reverse order
void split_vector(int sel_site, vector<double> &whole_vec, vector<double> &back_vec, vector<double> &fwd_vec) 
{
    for (int i = sel_site; i < whole_vec.size(); i++) {
        fwd_vec.push_back(whole_vec[i]) ;
    }

    // Inverts vector from the selected site to the start. Removes the first site, as it has recombination rate 0.5
    for (int i = sel_site; i > 0; i--) {
        back_vec.push_back(whole_vec[i]) ;
    }
}

// test. remove
void split_vector_int(int sel_site, vector<int> &whole_vec, vector<int> &back_vec, vector<int> &fwd_vec) 
{
    for (int i = sel_site; i < whole_vec.size(); i++) {
        fwd_vec.push_back(whole_vec[i]) ;
    }

    // Inverts vector from the selected site to the start. Removes the first site, as it has recombination rate 0.5
    for (int i = sel_site; i > 0; i--) {
        back_vec.push_back(whole_vec[i]) ;
    }
}

#endif
