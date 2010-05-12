
#ifndef	PREFIXSPAN_H
#define	PREFIXSPAN_H

#include <map>
#include <functional>
#include <vector>
#include <list>
#include <set>

#include <Sequence.h>
#include <ExplorationGuide2.h>

namespace SeqMine
{

class
PrefixSpanOptions
{
public:
	unsigned int minsup;	/* Minimum absolute support */
	unsigned int minlen;	/* Minimum sequence length (counted in items) */
	unsigned int maxlen;	/* Maximum sequence length (counted in items) */
	int maxgap;				/* Maximum allowed gap or -1 for infinite gaps */

	PrefixSpanOptions (void)
		: minsup (0), minlen (0), maxlen (0), maxgap (-1)
	{
	}

	PrefixSpanOptions (unsigned int minsup, unsigned int minlen = 0,
		unsigned int maxlen = 0, int maxgap = -1) :
		minsup (minsup), minlen (minlen), maxlen (maxlen),
		maxgap (maxgap)
	{
	}
};


class
PrefixSpan
{
private:
	class ExplorationPriorityElement {
	public:
		SequenceT seq;
		unsigned int item;
		bool itemJoined;

		ExplorationPriorityElement (SequenceT& seq, unsigned int item,
			bool itemJoined)
			: seq (seq), item (item), itemJoined (itemJoined)
		{
		}
	};

	typedef std::set<unsigned int, std::less<unsigned int> >::iterator set_iterator;
	typedef std::vector<const Sequence*>::iterator seq_iterator;
	typedef std::map<unsigned int, unsigned int>::iterator freq_iterator;

	/* Mine the frequent items across the set of sequences 'S'.  Report the
	 * frequent items in 'frequent'.
	 *
	 * lastItem: The last item of the last element of the current sequence.
	 * S: The set of sequences to mine frequent items in.
	 * maximumGapLength: If -1 arbitrary gap lengths are allowed, if >0, it
	 *    specifies the maximum allowed gap length for a subsequence to match.
	 *
	 * Output:
	 * 'occ' and 'occ_': The maps [item_idx] -> set of sequence indices
	 */
	void MineFrequentItems (unsigned int lastItem,
		const std::vector<unsigned int>* lastPrefixElement,
		std::map<unsigned int, std::set<unsigned int> >& occ,
		std::map<unsigned int, std::set<unsigned int> >& occ_,
		std::vector<const Sequence*>& S,
		int maximumGapLength = -1) const;

	/* Merge multiple mapped index sets into a joint output set.
	 *
	 * alpha: The current prefix sequence.
	 * S: The sequences, we are interested in the .Occurence members.
	 * setMap: The indices of active sets.
	 * setOut: The joint output set.
	 */
	void indexSetJoin (const SequenceT& alpha,
		const std::vector<const Sequence*>& S,
		const std::set<unsigned int>& setMap, std::set<unsigned int>& setOut) const;

	std::vector<const Sequence*> prefixProject (unsigned int item, bool itemJoined,
		std::vector<const Sequence*> S,
		const std::vector<unsigned int>* lastPrefixElement,
		bool maximumGapMining = false) const;

	std::vector<const Sequence*> prefixProjectGap (unsigned int item,
		std::vector<const Sequence*> S, const Sequence& subseq,
		const std::vector<unsigned int>* lastPrefixElement, unsigned int ahead) const;

	/* Clean the non frequent items from the occurence map 'occ'.
	 *
	 * Note that this is the only point where the minimum support condition is
	 * implemented.  As the condition is tied into the algorithm there is no
	 * way of allowing TestExtendSequence to decide about frequency.
	 */
	void cleanNonFrequent
		(std::map<unsigned int, std::set<unsigned int> >& occ) const;

protected:
	typedef std::map<unsigned int, std::set<unsigned int> >::iterator occ_iterator;

	/* Iterative deepening internal variable: in order to avoid adding the
	 * same (good) element twice to the overall best pattern list, in
	 * successive mining iterations, we do not report patterns which we
	 * reported before.
	 *
	 * To actually do this is easy: we only report patterns at the last level.
	 */
	bool usingIterativeDeepening;
	unsigned int iterativeDeepeningLevel;
	bool seenLeafLevelExtensions;

	/* Standard PrefixSpan options
	 */
	PrefixSpanOptions options;

	ExplorationGuide eguide;

	unsigned int max_item_index;

	std::vector<const Sequence*>* Sall;

	/* Debug method to verify that PrefixSpan's internal 'Occurence' list is
	 * in a consistent state.  This method is very expensive to call, use it
	 * only for debugging runs.
	 *
	 * In case an inconsistency is found, the details are given and we throw
	 * an assertion.
	 */
	bool verifyOccurenceCorrectness (const Sequence& alpha) const;

	/* Predicate: shall the given sequence be mined further (extended)?
	 * This can be used to plug in a gain-bound based decision for weighted
	 * subsequence mining.
	 */
	virtual bool TestExtendSequence (const Sequence& seq) const;

	/* Filter function: can remove (but not add) items from occ.
	 */
	virtual void FilterFrequentItems (
		std::map<unsigned int, std::set<unsigned int> >& occ) const;

	/* Report a frequent subsequence 'seq'.  The frequency map is available
	 * through seq.Occurence.
	 */
	virtual void Report (const Sequence& seq);

	/* The priority for ordered tree exploration.  A higher priority will lead
	 * to have the sequence visited earlier.  For frequent subsequence mining
	 * this is of no use and an overhead only, but for weighted subsequence
	 * mining this can lead to substantial search space cuts.
	 */
	virtual double ExplorationPriority (const Sequence& seq) const;
	virtual double ExplorationPriorityBound (const Sequence& seq) const;

	/* Iterative prioritized mining: execute one iteration.
	 * The central data structure used is eguide.
	 */
	void MineIteration (void);

	/* Main recursive mining function.
	 *
	 * alpha: The current prefix sequence.
	 * S: The projection set S|alpha.
	 * lengthTail: Length extensions still allowed.
	 */
	void Mine (SequenceT& alpha, std::vector<const Sequence*>& S,
		unsigned int lengthTail);

	bool doVerify;
	mutable unsigned long long verifyCount;

public:
	std::list<std::pair<unsigned int, SequenceT> > subsequences;

	PrefixSpan (PrefixSpanOptions options)
		: usingIterativeDeepening (false),
			options (options), doVerify (false), verifyCount (0)
	{
	}

	virtual ~PrefixSpan (void)
	{
	}

	/* Sequence mining, starting with the empty sequence.
	 * Finds all frequent subsequences in the set of sequences 'S'.
	 * Also sets max_item_index properly.
	 */
	virtual void Mine (std::vector<const Sequence*>& S);

	/* Iterative prioritized mining.
	 */
	virtual void Mine2 (std::vector<const Sequence*>& S);

	/* Iterative-deepening prioritized mining.  (Jason says: "sounds like
	 * id-Astar to me").
	 */
	virtual void Mine3 (std::vector<const Sequence*>& S);

	void SetVerification (bool doVerify = true)
	{
		this->doVerify = doVerify;
	}
};

};

#endif

