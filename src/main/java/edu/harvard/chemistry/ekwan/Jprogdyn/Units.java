package edu.harvard.chemistry.ekwan.Jprogdyn;

/**
 * This interface holds some unit conversions.
 */
public interface Units {

    /** speed of light, cm/s */
    public static final double C                         = 29979245800.0;      
    
    /** Planck's constant, m^2 * kg / s */
    public static final double H                         = 6.626075E-34;       
    
    /** molar gas constant, J / (mol * K) */
    public static final double R_GAS_J                   = 8.31447;            
    
    /** molar gas constant, kcal / (mol * K) */
    public static final double R_GAS_KCAL                = 0.0019858;          
    
    /** Boltzmann's constant, eV / K */
    public static final double BOLTZMANN_EV              = 8.61733238E-5;      
    
    /** Boltzmann's constant, kcal / K */
    public static final double BOLTZMANN_KCAL            = 3.29983E-27;        
    
    /** converts J to kcal/mol */
    public static final double J_TO_KCAL_PER_MOL         = 6.0221415E23/4184.0;
    
    /** converts kcal/mol to J */
    public static final double KCAL_PER_MOL_TO_J         = 4184.0/6.0221415E23;
    
    /** Avogadro's number */
    public static final double AVOGADROS_NUMBER          = 6.0221415E23;       
    
    /** Angstroms in a bohr */
    public static final double BOHR                      = 0.52917721092;      
    
    /** kcal per mol/hartree */
    public static final double KCAL_PER_HARTREE          = 627.509369;         
    
    /** J per mol/hartree */
    public static final double J_PER_HARTREE             = 2625499.62;         
}
