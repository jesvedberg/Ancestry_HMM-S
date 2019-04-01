#ifndef __FWD_ITER_H
#define __FWD_ITER_H

vector<vector<double>> fwd_iter(vector<double> &recombination_rate, vector<double> &basefreq, int tt, double m, int generations, int n) 
{

    vector<double> freq(basefreq) ;
    vector<double> freq_(basefreq);
    vector<vector<double>> transition_rates ;
    double sum ;
    double h11;
    double h12;
    double h21;
    double h22;
    double r;

    double a1;
    double a1_;

    double h11_;
    double h12_;
    double h21_;
    double h22_;

    double p_coal;

    for ( int site = 0 ; site < recombination_rate.size() ; site ++ ) {
        //for (double site = 0; site < 0.5 ; site += r) {
        r = recombination_rate[site] ;
        vector<double> tr ;
        //freq_[tt+1] = m ;

        h11 = m;
        h12 = 0;
        h21 = 0;
        h22 = 1-m;

        p_coal = 1 ;

        for (int t = tt+1; t < tt + generations; t++) {
            a1 = h11 + h12;
            a1_ = h11 + h21;

            h11_ = h11*(1-r) + a1*r*a1_*p_coal;
            h12_ = h12*(1-r) + a1*r*(1-a1_)*p_coal;
            h21_ = h21*(1-r) + (1-a1)*r*a1_*p_coal;
            h22_ = h22*(1-r) + (1-a1)*r*(1-a1_)*p_coal;

            freq_[t] = ( h11_ + h21_ )/( h11 + h12 + h21 + h22 );

            h11 = h11_ ; 
            h12 = h12_ ; 
            h21 = h21_ ; 
            h22 = h22_ ; 
            
            h11 = h11/freq[t]*freq[t+1]; 
            h12 = h12/freq[t]*freq[t+1]; 
    
            h21 = h21/(1-freq[t])*(1-freq[t+1]); 
            h22 = h22/(1-freq[t])*(1-freq[t+1]); 

            sum = h11 + h12 + h21 + h22; 
            h11 = h11/sum; 
            h12 = h12/sum; 
            h21 = h21/sum; 
            h22 = h22/sum; 

            freq_[t+1] = (h11 + h21) ;
            p_coal *= ( 1 - 1/(2*n) ) ;

            //cerr << freq[t] << "," << freq_[t] << " " ;
        }
        //cerr << endl ;
        /*cout << site << " " << setprecision(15) << h11 << " " << h12 << " " << h21 << " " << h22 << " ";
        cout << h12/((h11+h12)*r) << " ";
        cout << h21/((h21+h22)*r) << " ";
        cout << h12/((h12+h22)*r) << " ";
        cout << h21/((h21+h11)*r) << " " << endl;*/
        //cerr << h12/((h11+h12)*r) << " ";
        tr.push_back(h12/((h11+h12)*r)) ;
        tr.push_back(h21/((h21+h22)*r)) ;
        tr.push_back(h12/((h12+h22)*r)) ;
        tr.push_back(h21/((h21+h11)*r)) ;

        transition_rates.push_back(tr) ;

        //copy(begin(freq_), end(freq_), begin(freq));
        //memcpy(freq, freq_, 300*sizeof(double));
        freq = freq_ ;
    }

    return transition_rates;
}   

#endif
