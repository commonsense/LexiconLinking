
#include <limits>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>
#include <boost/lambda/algorithm.hpp>

#include <PrefixSpanK.h>


namespace SeqMine
{

using namespace boost::lambda;

bool
PrefixSpanK::TestExtendSequence (const Sequence& seq) const
{
	/* This is subtle: normally we would only generate frequent subsequences,
	 * but it can happen that we are in a recursion and min_support just got
	 * increased, so we still see non-frequent subsequences
	 */
	return (seq.Occurence.size () >= options.minsup);
}

double
PrefixSpanK::ExplorationPriority (const Sequence& seq) const
{
	/* For plain frequent mining we don't care about priorities, but here it
	 * is beneficial to harvest high-frequency subbranches first.
	 */
	return (seq.Occurence.size ());
}

double
PrefixSpanK::ExplorationPriorityBound (const Sequence& seq) const
{
	/* The bound is just the frequency here.
	 */
	return (seq.Occurence.size () + 0.5);
}


void
PrefixSpanK::Report (const Sequence& seq)
{
	if (doVerify)
		verifyOccurenceCorrectness (seq);

	unsigned int freq = seq.Occurence.size ();

	if (freq < options.minsup)
		return;

	/* Check minimum length requirement.
	 */
	if (seq.length () < options.minlen)
		return;

	/* Add frequent sequence to most-frequent list
	 */
	mostFrequent.insert (std::pair<unsigned int, SequenceT> (freq, seq));

	if (mostFrequent.size () <= K)
		return;

	frequent_map_type::iterator newMFend = --mostFrequent.end ();
	for (unsigned int n = (mostFrequent.size () - K) ; n > 0 ; --n) {
		// TODO: this is inefficient, and we should use something like
		// boost::multi_index, however it does not support absolute random
		// access indices.
		//std::cerr << "   n = " << n << ", newMFend->first = " << newMFend->first << std::endl;
		--newMFend;
	}

	if (newMFend->first == options.minsup)
		return;

	/* Set new minimum support and find new end of sequence.
	 */
	options.minsup = newMFend->first;
	std::cerr << "   minsup increased to " << options.minsup;
	while (newMFend != mostFrequent.end () && newMFend->first == options.minsup)
		++newMFend;

	unsigned int old_count = mostFrequent.size ();
	mostFrequent.erase (newMFend, mostFrequent.end ());

	std::cerr << ", pruned " << (old_count-mostFrequent.size()) << " sequences, "
		<< mostFrequent.size () << " sequences remaining." << std::endl;
}


void
PrefixSpanK::Mine (std::vector<const Sequence*>& S)
{
	Mine3 (S);
}

unsigned int
PrefixSpanK::MinimumSupportThreshold (void) const
{
	return (options.minsup);
}

}

