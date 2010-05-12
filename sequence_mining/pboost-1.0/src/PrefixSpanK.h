
#ifndef	PREFIXSPANK_H
#define	PREFIXSPANK_H

#include <map>
#include <functional>
#include <utility>
#include <algorithm>

#include <PrefixSpan.h>
#include <BestPatternList.h>

namespace SeqMine
{
using namespace ::boost::multi_index;

/* Top-K frequent subsequence mining.
 */
class
PrefixSpanK : public PrefixSpan
{
private:
	unsigned int K;

protected:
	virtual bool TestExtendSequence (const Sequence& seq) const;
	virtual void Report (const Sequence& seq);
	virtual double ExplorationPriority (const Sequence& seq) const;
	virtual double ExplorationPriorityBound (const Sequence& seq) const;

public:
	/* The list of mined patterns.  After calling ::Mine, this list contains
	 * the ordered set of maximum frequent patterns, most frequent ones first.
	 */
	typedef std::multimap<unsigned int, SequenceT, std::greater<unsigned int> >
		frequent_map_type;
	frequent_map_type mostFrequent;

	/* Top-K frequent PrefixSpan algorithm
	 *
	 * K: The number of most frequent subsequences to find.
	 * options: The general PrefixSpan options, defining the feature space.
	 */
	PrefixSpanK (unsigned int K, PrefixSpanOptions options)
		: PrefixSpan (options), K (K)
	{
	}

	/* Finds the P top frequent subsequences in the set of sequences 'S'.
	 * P >= K, such that the P most frequent subsequences are mined.
	 *
	 * A given K will usually not be exactly at a frequency boundary, hence
	 * the next boundary after the top K sequences is chosen.
	 *
	 * Note that P could be less than K if no more frequent subsequences
	 * exist, for example with tight minlen/maxlen constraints.
	 */
	virtual void Mine (std::vector<const Sequence*>& S);

	/* After mining: return the minimum support threshold found.
	 */
	unsigned int MinimumSupportThreshold (void) const;
};

};

#endif


