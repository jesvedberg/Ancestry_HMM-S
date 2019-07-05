#ifndef __CMD_LINE_H
#define __CMD_LINE_H

/// command line information and global parameters
class cmd_line {
public:
    
    /// terms to bound the optimization
    double t_max ;
    double t_min ;
    
    /// to bound proportion search
    double p_max ;
    double p_min ;
    
    /// to create intial simplex points
    double t_length ;
    double p_length ;
    
    /// number of restarts
    int n_restarts ;
    
    /// proportion ancestry for 0-n must sum to 1
    /// these are therefore the final ancestry proportion, not necessary the proportion that fluxed if additional pulses occured closer to the present
    vector<double> ancestry_proportion ;
    
    /// store relevant ancestry information
    vector<pulse> ancestry_pulses ;
    
    /// diploid effective population size ( i.e. 2n )
    double ne ;
    
    /// tolerance for parameter search
    double tolerance ;
    
    /// minimum recombinational distance between markers
    double minimum_distance ;
    
    /// error rates for reads (if read based) or genotypes (if genotype based)
    double error_rate ;
    
    /// bool sample is expressed as genotypes, not read counts
    bool genotype ;
    
    /// ancestral genotype frequencies are fixed
    bool ancestral_fixed ;
    
    /// viterbi output
    /// caution: not recommended for more samples of ploidy > 1
    bool viterbi ;
    
    /// number of digits of precision to include
    int precision ;
    
    /// output actual pulses rather than ancestry states
    bool output_pulses ;
    
    /// error rates specifed
    bool error_rates ; 
    
    /// input file name
    string input_file ;
    
    /// sample file
    string sample_file ;
    
    /// bootstrap
    int n_bootstraps ;
    int block_size ; 

    /// +=+=+=+=+=+=+ selection +=+=+=+=+=+=+
    bool is_limit;
    string limit_chr ;
    int limit_win_start ;
    int limit_win_end ;
    //int sel_site ;
    double sel_min;
    double sel_max;
    int pos_margin;
    int pos_min;
    int pos_max;
    double pos_limit;
    double sel_limit;
    bool is_limitpos;

    bool calc_grid;
    int grid_pstart;
    int grid_pstop;
    int grid_pstep;
    double grid_sstart;
    double grid_sstop;
    double grid_sstep;

    bool test_point;
    int test_pos;
    double test_sel;

    bool run_gss;
    int gs_pstart;
    int gs_pstop;
    int gs_pstep;
    double gs_sstart;
    double gs_sstop;
    int gs_max_iterations;
    double gs_precision;

    // HMM chain window size
    string win_unit;
    double win_morgan;
    double win_percent;
    // int win_bp;
    
    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;

} ;

#endif

