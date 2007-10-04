#ifndef CRYPTCELLMUTATIONSTATES_HPP_
#define CRYPTCELLMUTATIONSTATES_HPP_

/**
 * Possible types mutation state for TissueCells.
 */
typedef enum CryptCellMutationState_
{
    HEALTHY,				// Wild-type cell
    APC_ONE_HIT,			// APC +/-
    APC_TWO_HIT,			// APC -/-
    BETA_CATENIN_ONE_HIT,	// Beta-catenin with a change at residue 45
    LABELLED,               // To paint a different colour but not actually mutant
} CryptCellMutationState;


#endif /*CRYPTCELLMUTATIONSTATES_HPP_*/
