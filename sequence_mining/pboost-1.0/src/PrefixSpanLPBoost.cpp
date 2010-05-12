/* PrefixSpan based stump finding
 * Copyright (C) 2007  Sebastian Nowozin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <vector>
#include <iostream>

#include <stdlib.h>

#include <PrefixSpanLPBoost.h>


namespace SeqMine
{

double
PrefixSpanLPBoost::ExplorationPriority (const Sequence& seq) const
{
	double curgain_pos;
	double curgain_neg;

	if (directional) {
		curgain_pos = gainDirectional (seq, 1.0);
		curgain_neg = gainDirectional (seq, -1.0);
	} else {
		curgain_pos = gain (seq, 1.0);
		curgain_neg = gain (seq, -1.0);
	}

	double gainmax = curgain_pos;
	if (curgain_neg > curgain_pos)
		gainmax = curgain_neg;

#if 0
	std::cout << "PrefixSpanLPBoost::ExplorationPriority (" << seq << ") = "
		<< gainmax << std::endl;
#endif

	return (gainmax);
}

double
PrefixSpanLPBoost::ExplorationPriorityBound (const Sequence& seq) const
{
	return (gainbound (seq));
}

bool
PrefixSpanLPBoost::TestExtendSequence (const Sequence& seq) const
{
	extendTestCount += 1;

	if (extendTestCount % 10000 == 0) {
		std::cout << "[" << extendTestCount << "] (" << seq << ")" << std::endl;
		std::cout << "    gainbound "
			<< (directional ? gainboundDirectional (seq) : gainbound (seq))
			<< ", needed " << bestPatterns.LeastBestGain ();

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
		(directional ? gainboundDirectional (seq) : gainbound (seq))
			> bestPatterns.LeastBestGain ());
}

double
PrefixSpanLPBoost::gainbound (const Sequence& seq) const
{
	/* <t^,y^> = argmax_{t \in F,y \in {+/- 1}} d_i y_i h_{<t,y>}, where
	 * F = \unionset_{i=1}^{L} { t | t \subseteq x_i }.
	 */
	double gain_pos = 0.0;
	double gain_neg = 0.0;

#ifdef	DEBUG
	std::cout << "   following gainbound stemming from:" << std::endl;
	std::cout << "   ";
#endif
	for (std::set<unsigned int>::const_iterator si = seq.Occurence.begin () ;
		si != seq.Occurence.end () ; ++si)
	{
#ifdef	DEBUG
		std::cout << *si << ", ";
#endif

		if (Ylabels[*si] <= 0.0)
			gain_neg += W[*si];
		else
			gain_pos += W[*si];
	}
#ifdef DEBUG
	std::cout << std::endl;
#endif

	gain_neg = 2.0*gain_neg + boostWeightSum;
	gain_pos = 2.0*gain_pos - boostWeightSum;
	double gainmax = (gain_neg >= gain_pos) ? gain_neg : gain_pos;

#ifdef	DEBUG
	std::cout << "::gainbound (" << seq << ") = " << gainmax << std::endl;
#endif

	return (gainmax);
}


double
PrefixSpanLPBoost::gainboundDirectional (const Sequence& seq) const
{
	/* gain(t',\omega) <= max{ \sum_{\{n|y_n=1\}} \lambda_n I(t \subseteq x_n),
	 *                         \sum_{\{n|y_n=-1\}} \lamdba_n I(t \subseteq x_n) }
	 */
	double gain_pos = 0.0;
	double gain_neg = 0.0;

	/* Loop only over I(t \subseteq x_n) samples
	 */
	for (std::set<unsigned int>::const_iterator si = seq.Occurence.begin () ;
		si != seq.Occurence.end () ; ++si)
	{
		/* Add \lambda_n to the respective gain bound
		 */
		if (Ylabels[*si] <= 0.0)
			gain_neg += W[*si];
		else
			gain_pos += W[*si];
	}

	double gainmax = (gain_neg >= gain_pos) ? gain_neg : gain_pos;

	return (gainmax);
}


void
PrefixSpanLPBoost::Report (const Sequence& seq)
{
#ifdef DEBUG
	/* Expensive verification of correctness
	 */
	verifyOccurenceCorrectness (seq);
#endif

	/* 1. Compute gain for h(. ; 1, seq) and h(. ; -1, seq)
	 */
	double curgain_pos = gain (seq, 1.0);
	double curgain_neg = gain (seq, -1.0);

	double ybase = 1.0;
	double gainmax = curgain_pos;
	if (curgain_neg > curgain_pos) {
		ybase = -1.0;
		gainmax = curgain_neg;
	}

	if (bestPatterns.Insert (gainmax, seq, ybase)) {
		std::cout << "   added sequence, gain " << gainmax << ":" << std::endl;
		std::cout << seq << std::endl;
	}
}


double
PrefixSpanLPBoost::gain (const Sequence& seq, double ybase) const
{
	//std::cout << "PrefixSpanLPBoost (" << seq << ", " << ybase << ")" << std::endl;

	double sum = 0.0;

	/* Compute actual gain (eq. (10) in [Dimiriz2002]) and Problem 1 in graph
	 * boost paper.
	 *
	 * For L samples x_i, calculate
	 *   gain = \sum_{i=1}^{L} d_i y_i h(x_i)
	 * where h is implicitly defined by this current subgraph pattern.
	 */
	unsigned int last_sample_index = 0;
	for (std::set<unsigned int>::const_iterator si = seq.Occurence.begin () ;
		si != seq.Occurence.end () ; ++si)
	{
		/* All samples in between the last one and this one
		 */
		for (unsigned int n = last_sample_index ; n < *si ; ++n)
			sum -= W[n]*Ylabels[n]*ybase;

		/* The current sample
		 */
		sum += W[*si]*Ylabels[*si]*ybase;
		last_sample_index = *si + 1;
	}

	/* Remaining samples at the end
	 */
	for (unsigned int n = last_sample_index ; n < W.size () ; ++n)
		sum -= W[n]*Ylabels[n]*ybase;

#ifdef	DEBUG
	std::cout << "::gain (" << seq << ") is " << sum << std::endl;
#endif
	return (sum);
}


double
PrefixSpanLPBoost::gainDirectional (const Sequence& seq, double ybase) const
{
	/* Compute actual gain (eq. (10) in [Dimiriz2002]) and Problem 1 in graph
	 * boost paper.
	 *
	 * For L samples x_i, calculate
	 *   gain = \sum_{i=1}^{L} d_i y_i h(x_i)
	 *        = \sum_{i=1}^{L} d_i y_i I(t \subseteq x_i)
	 *
	 * where h is implicitly defined by this current subgraph pattern.
	 */
	double sum = 0.0;
	for (std::set<unsigned int>::const_iterator si = seq.Occurence.begin () ;
		si != seq.Occurence.end () ; ++si)
	{
		sum += W[*si]*Ylabels[*si]*ybase;
	}

	return (sum);
}


void
PrefixSpanLPBoost::Mine (std::vector<const Sequence*>& S)
{
	extendTestCount = 0;

	switch (searchStrategy) {
	case DepthFirstSearch:
		PrefixSpan::Mine (S);
		break;
	case AStar:
		PrefixSpan::Mine2 (S);
		break;
	case IDAStar:
		PrefixSpan::Mine3 (S);
		break;
	default:
		std::cerr << "Unknown search strategy type." << std::endl;
		exit (EXIT_FAILURE);
		break;
	}
}

}


