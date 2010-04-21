
#ifndef	PREFIXSPANLPBOOST15_H
#define	PREFIXSPANLPBOOST15_H

#include <vector>

#include <PrefixSpan.h>
#include <PrefixSpanLPBoost.h>
#include <BestPatternList.h>

namespace SeqMine
{

class
PrefixSpanLPBoost15 : PrefixSpan
{
private:
	mutable unsigned int extendTestCount;
	unsigned int maximumIterations;
	bool memoryConserving;

	/* Compute the 1.5-class 1/0 subseteq stump gain given sequence 'seq'.
	 * That is, compute
	 *
	 *   gain = \sum_{i=1}^{N} d_i y_i (2I(seq \subseteq x_i) - 1)
	 */
	double gain (const Sequence& seq) const;

	/* Compute the maximum gain possible for any decision stumps based on any
	 * supersequence of the given sequence 'seq'.
	 *
	 * This is different for the 1.5-class case.
	 */
	double gainbound (const Sequence& seq) const;

protected:
	/* Boosting parameters and variables
	 */
	std::vector<double>& Ylabels;
	std::vector<double>& W;

	double minimumRequiredGain;

	/* Minimum support constraint in the positive subset of the training set.
	 */
	unsigned int pos_minimum_support;

	virtual bool TestExtendSequence (const Sequence& seq) const;
	virtual void FilterFrequentItems (
		std::map<unsigned int, std::set<unsigned int> >& occ) const;
	virtual void Report (const Sequence& seq);
	virtual double ExplorationPriority (const Sequence& seq) const;
	virtual double ExplorationPriorityBound (const Sequence& seq) const;

public:
	PrefixSpanLPBoost::SearchStrategy searchStrategy;
	BestPatternList bestPatterns;

	/* Ylabels: (1,n) vector of sample labels in { -1, 1 }.
	 * W: (1,n) vector of real positive weights, one for each sample.
	 * minimumRequiredGain: lower necessary bound (theta) on new pattern's
	 *    gain.
	 * keepCount: the number patterns to extract, >= 1.
	 * maximumIterations: If non-zero, denotes a limit on the iteration count.
	 * minsup: minimum required support of the found patterns or zero if no
	 *    minimum support is required.
	 * pos_minimum_support: minimum support in the positive training set.
	 * minlen: minimum required sequence length or zero if no minimum length
	 *    is required.
	 */
	PrefixSpanLPBoost15 (PrefixSpanOptions options,
		std::vector<double>& Ylabels, std::vector<double>& W,
		double minimumRequiredGain = 0.0, unsigned int keepCount = 1,
		unsigned int maximumIterations = 0,
		unsigned int pos_minimum_support = 0,
		PrefixSpanLPBoost::SearchStrategy searchStrategy =
			PrefixSpanLPBoost::DepthFirstSearch)
		: PrefixSpan (options),
			maximumIterations (maximumIterations),
			Ylabels (Ylabels), W (W), minimumRequiredGain (minimumRequiredGain),
			pos_minimum_support (pos_minimum_support),
			searchStrategy (searchStrategy),
			bestPatterns (keepCount)
	{
	}

	virtual void Mine (std::vector<const Sequence*>& S);
};

};

#endif

