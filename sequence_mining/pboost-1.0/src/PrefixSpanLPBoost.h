
#ifndef	PREFIXSPANLPBOOST_H
#define	PREFIXSPANLPBOOST_H

#include <vector>

#include <PrefixSpan.h>
#include <BestPatternList.h>

namespace SeqMine
{

class
PrefixSpanLPBoost : PrefixSpan
{
private:
	mutable unsigned int extendTestCount;
	unsigned int maximumIterations;
	bool memoryConserving;

	/* Compute the 2-class +1/-1 subseteq stump gain given sequence 'seq' and
	 * the positive stump reply of 'ybase'.  That is, compute
	 *
	 *   gain = \sum_{i=1}^{N} d_i y_i ybase (2I(seq \subseteq x_i) - 1)
	 */
	double gain (const Sequence& seq, double ybase) const;

	/* Compute the 2-class one-sided stump gain given sequence 'seq' and
	 * the positive stump reply of 'ybase'.  That is, compute
	 *
	 *   gain = \sum_{i=1}^{N} d_i y_i ybase I(seq \subseteq x_i)
	 */
	double gainDirectional (const Sequence& seq, double ybase) const;

	/* Compute the maximum gain possible for any decision stumps based on any
	 * supersequence of the given sequence 'seq'.
	 */
	double gainbound (const Sequence& seq) const;

	/* The gainbound for the one-sided stumps.
	 */
	double gainboundDirectional (const Sequence& seq) const;

protected:
	/* Boosting parameters and variables
	 */
	std::vector<double>& Ylabels;
	std::vector<double>& W;

	double minimumRequiredGain;

	/* The constant term in the gainbound computation.  We precompute it once.
	 */
	double boostWeightSum;

	virtual bool TestExtendSequence (const Sequence& seq) const;
	virtual void Report (const Sequence& seq);
	virtual double ExplorationPriority (const Sequence& seq) const;
	virtual double ExplorationPriorityBound (const Sequence& seq) const;

public:
	enum SearchStrategy {
		AStar = 1,
		DepthFirstSearch = 2,
		IDAStar = 3,
	};

	/* Use one-sided directional decision stumps?
	 *
	 * If true, the stumps have the form:
	 *     h(x ; t, \omega) = \omega I(t \subseteq x)
	 */
	bool directional;

	/* Search traversal strategy/type.
	 *
	 * DepthFirstSearch is a naive one-step-lookahead strategy.
	 * AStar is the most efficient in runtime, but bears a giant memory
	 *    requirement.
	 * IDAStar is a tradeoff, using few memory but a better traversal than
	 *    DFS.
	 */
	SearchStrategy searchStrategy;

	/* The list of mined patterns.  After calling ::Mine, this list contains
	 * the ordered set of maximum gain patterns, best gains first.
	 */
	BestPatternList bestPatterns;

	/* Ylabels: (1,n) vector of sample labels in { -1, 1 }.
	 * W: (1,n) vector of real positive weights, one for each sample.
	 * minimumRequiredGain: lower necessary bound (theta) on new pattern's
	 *    gain.
	 * keepCount: the number patterns to extract, >= 1.
	 * maximumIterations: If non-zero, denotes a limit on the iteration count.
	 * minsup: minimum required support of the found patterns or zero if no
	 *    minimum support is required.
	 * minlen: minimum required sequence length or zero if no minimum length
	 *    is required.
	 * memoryConserving: If true, only few memory is used but the search space
	 *    is generated on demand, which is slow.  If false, a global search
	 *    space tree is used.
	 */
	PrefixSpanLPBoost (PrefixSpanOptions options,
		std::vector<double>& Ylabels, std::vector<double>& W,
		double minimumRequiredGain = 0.0, unsigned int keepCount = 1,
		unsigned int maximumIterations = 0,
		SearchStrategy searchStrategy = DepthFirstSearch)
		: PrefixSpan (options),
			maximumIterations (maximumIterations),
			Ylabels (Ylabels), W (W), minimumRequiredGain (minimumRequiredGain),
			directional (false),
			searchStrategy (searchStrategy),
			bestPatterns (keepCount)
	{
		boostWeightSum = 0.0;
		for (unsigned int n = 0 ; n < Ylabels.size() ; ++n)
			boostWeightSum += W[n] * Ylabels[n];
	}

	virtual void Mine (std::vector<const Sequence*>& S);
};

};

#endif

