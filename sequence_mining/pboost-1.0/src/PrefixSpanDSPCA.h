
#ifndef	PREFIXSPANDSPCA_H
#define	PREFIXSPANDSPCA_H

#include <assert.h>

#include <vector>
#include <ext/hash_set>
#include <gmm/gmm.h>

#include <PrefixSpan.h>

namespace SeqMine
{

class
PrefixSpanDSPCA : public PrefixSpan
{
private:
	unsigned int basisNonZero;
	unsigned int basisCount;
	unsigned int sequenceCount;
	unsigned int reportCount;

	typedef __gnu_cxx::hash_set<SequenceT, sequence_hash, sequence_identical>
		forbidden_pattern_list_type;
	forbidden_pattern_list_type forbiddenSequences;

	bool foundPattern;

protected:
	virtual bool TestExtendSequence (const Sequence& seq) const;
	virtual void Report (const Sequence& seq);
	virtual double ExplorationPriority (const Sequence& seq) const;
	virtual double ExplorationPriorityBound (const Sequence& seq) const;

public:
	/* The principal covariance submatrix of the n'th principal component.
	 */
	std::vector<gmm::dense_matrix<double> > Cv;
	unsigned int components;

	/* Cpatterns[n] contains all patterns of the n'th principal component.
	 * The order in each Cpatterns[n] is the same as in Cv[n].
	 */
	typedef std::vector<SequenceT> component_patterns_type;
	std::vector<component_patterns_type> Cpatterns;
	typedef std::vector<double> component_eval_type;
	std::vector<component_eval_type> Ceval;

	SequenceT bestSeq;
	double bestSeqLambdaMax;
	gmm::dense_matrix<double> bestC;	// new covariance matrix

	/* Direct Sparse PCA Subsequence Mining Algorithm
	 *
	 * basisNonZero: The number of most non-zero coefficients in each basis
	 *    vector.
	 * basisCount: The number of basis vectors to find.
	 * options: The general PrefixSpan options, defining the feature space.
	 */
	PrefixSpanDSPCA (unsigned int basisNonZero, unsigned int basisCount,
		PrefixSpanOptions options)
		: PrefixSpan (options),
			basisNonZero (basisNonZero), basisCount (basisCount),
			foundPattern (false), components (0)
	{
		assert (basisNonZero > 0);
		assert (basisCount > 0);
	}

	/* Finds the sparse basis vectors capturing the most variance in a given
	 * set of sequences.
	 */
	void PCA (std::vector<const Sequence*>& S);
};

};

#endif

