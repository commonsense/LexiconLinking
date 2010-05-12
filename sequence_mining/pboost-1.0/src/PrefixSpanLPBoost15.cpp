
#include <vector>

#include <PrefixSpanLPBoost15.h>


namespace SeqMine
{

double
PrefixSpanLPBoost15::ExplorationPriority (const Sequence& seq) const
{
	return (gain (seq));
}

double
PrefixSpanLPBoost15::ExplorationPriorityBound (const Sequence& seq) const
{
	return (gainbound (seq));
}

bool
PrefixSpanLPBoost15::TestExtendSequence (const Sequence& seq) const
{
	extendTestCount += 1;

	if (extendTestCount % 10000 == 0) {
		std::cout << "[" << extendTestCount << "] (" << seq << ")" << std::endl;
		std::cout << "    gainbound " << gainbound (seq) << ", needed "
			<< bestPatterns.LeastBestGain ();

		if (usingIterativeDeepening) {
			std::cout << ", ID level " << iterativeDeepeningLevel;
			std::cout << " (LLext: " << seenLeafLevelExtensions << ")";
		}

		if (eguide.Size () > 0)
			std::cout << ", eguide.Size () = " << eguide.Size ();

		std::cout << std::endl;
	}

	if (maximumIterations != 0 && extendTestCount >= maximumIterations)
		return (false);

	if (options.maxlen > 0 && seq.length () >= options.maxlen)
		return (false);

	/* "Do we still have not enough patterns?" or
	 * "Is the upcoming possible improvement larger than the worst of the best
	 *  patterns so far?"
	 */
	return (bestPatterns.IsCapacityLeft () ||
		gainbound (seq) > bestPatterns.LeastBestGain ());
}


void
PrefixSpanLPBoost15::FilterFrequentItems (
	std::map<unsigned int, std::set<unsigned int> >& occ) const
{
	if (pos_minimum_support == 0)
		return;

	for (occ_iterator oi = occ.begin () ; oi != occ.end () ; ) {
		/* Check set size for minimum support
		 */
		unsigned int pos_support = 0;
		for (std::set<unsigned int>::const_iterator iv = oi->second.begin () ;
			iv != oi->second.end () ; ++iv)
		{
			if (Ylabels[*iv] >= 0.0)
				pos_support += 1;
		}

		if (pos_support < pos_minimum_support) {
			occ.erase (oi++);
		} else {
			++oi;
		}
	}
}


double
PrefixSpanLPBoost15::gainbound (const Sequence& seq) const
{
	double gainb = 0.0;

	for (std::set<unsigned int>::const_iterator si = seq.Occurence.begin () ;
		si != seq.Occurence.end () ; ++si)
	{
		if (Ylabels[*si] > 0.0)
			gainb += W[*si];
	}

#ifdef	DEBUG
	std::cout << "::gainbound (" << seq << ") = " << gainb << std::endl;
#endif

	return (gainb);
}

void
PrefixSpanLPBoost15::Report (const Sequence& seq)
{
	/* Compute gain for h(. ; seq)
	 */
	double gainmax = gain (seq);

	if (bestPatterns.Insert (gainmax, seq, 1.0)) {
		std::cout << "   added sequence, gain " << gainmax << ":" << std::endl;
		std::cout << seq << std::endl;
	}
}


double
PrefixSpanLPBoost15::gain (const Sequence& seq) const
{
	double sum = 0.0;

	for (std::set<unsigned int>::const_iterator si = seq.Occurence.begin () ;
		si != seq.Occurence.end () ; ++si)
	{
		sum += W[*si]*Ylabels[*si];
	}

	return (sum);
}


void
PrefixSpanLPBoost15::Mine (std::vector<const Sequence*>& S)
{
	extendTestCount = 0;

	switch (searchStrategy) {
	case PrefixSpanLPBoost::DepthFirstSearch:
		PrefixSpan::Mine (S);
		break;
	case PrefixSpanLPBoost::AStar:
		PrefixSpan::Mine2 (S);
		break;
	case PrefixSpanLPBoost::IDAStar:
		PrefixSpan::Mine3 (S);
		break;
	default:
		std::cerr << "Unknown search strategy type." << std::endl;
		exit (EXIT_FAILURE);
		break;
	}
}

}


