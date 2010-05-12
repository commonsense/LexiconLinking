/* pspan PrefixSpan implementation
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

/* PrefixSpan.cpp - PrefixSpan implementation for frequent and weighted
 *    substructure mining.
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 21st May 2007
 */

#include <iostream>
#include <algorithm>
#include <map>
#include <iterator>

#include <PrefixSpan.h>
#include <ExplorationGuide.h>

namespace SeqMine
{

/* Convert 'setMap' indices which point to S, into absolute sequence indices.
 */
void
PrefixSpan::indexSetJoin (const SequenceT& alpha,
	const std::vector<const Sequence*>& S,
	const std::set<unsigned int>& setMap, std::set<unsigned int>& setOut) const
{
	std::set<unsigned int> setTemp;

	for (std::set<unsigned int>::const_iterator smi = setMap.begin () ;
		smi != setMap.end () ; ++smi)
	{
		setTemp.insert (S[*smi]->Occurence.begin (), S[*smi]->Occurence.end ());
	}

	setOut.clear ();

	std::set_intersection (alpha.Occurence.begin (), alpha.Occurence.end (),
		setTemp.begin (), setTemp.end (),
		std::inserter (setOut, setOut.begin ()));
}
	
void
PrefixSpan::Mine (std::vector<const Sequence*>& S)
{
	this->Sall = &S;

	/* Set max_item_index correctly
	 */
	max_item_index = 0;

	for (std::vector<const Sequence*>::iterator iv = S.begin () ;
		iv != S.end () ; ++iv)
	{
		for (unsigned int elem = 0 ; elem < (*iv)->lengthElements () ; ++elem) {
			for (unsigned int iidx = 0 ; iidx < (*iv)->lengthSet (elem) ; ++iidx) {
				unsigned int item = (*iv)->get (elem, iidx);
				if (item > max_item_index)
					max_item_index = item;
			}
		}
	}

	max_item_index += 1;
#ifdef DEBUG
	std::cout << "Maximum item index: " << max_item_index << std::endl;
#endif

	/* Start the mining
	 */
	SequenceT alphaStart;
	for (unsigned int n = 0 ; n < S.size () ; ++n)
		alphaStart.Occurence.insert (S[n]->Occurence.begin (),
			S[n]->Occurence.end ());

	Mine (alphaStart, S, (options.maxlen == 0) ?
		std::numeric_limits<unsigned int>::max () : (options.maxlen-1));

	this->Sall = NULL;
}


void
PrefixSpan::Mine2 (std::vector<const Sequence*>& S)
{
	this->Sall = &S;

	/* Set max_item_index correctly
	 */
	max_item_index = 0;

	for (std::vector<const Sequence*>::iterator iv = S.begin () ;
		iv != S.end () ; ++iv)
	{
		for (unsigned int elem = 0 ; elem < (*iv)->lengthElements () ; ++elem) {
			for (unsigned int iidx = 0 ; iidx < (*iv)->lengthSet (elem) ; ++iidx) {
				unsigned int item = (*iv)->get (elem, iidx);
				if (item > max_item_index)
					max_item_index = item;
			}
		}
	}

	max_item_index += 1;
#ifdef DEBUG
	std::cout << "Maximum item index: " << max_item_index << std::endl;
#endif

	/* Start the mining
	 */
	SequenceT alphaStart;
	for (unsigned int n = 0 ; n < S.size () ; ++n)
		alphaStart.Occurence.insert (S[n]->Occurence.begin (),
			S[n]->Occurence.end ());

	eguide = ExplorationGuide ();
	//eguide.Add (0.0, 0.0, alphaStart, S, false);
	ExplorationElement eelem (0.0, 1.0, alphaStart, S, false);
	eguide.Add (eelem);

	while (eguide.IsEmpty () == false)
		MineIteration ();

	this->Sall = NULL;
}


/* Iterative-deepening variant of ::Mine2.  Thus, we do not keep a large
 * search tree in memory, but still achieve a runtime complexity similar to
 * the fast full-search-tree variant.
 */
void
PrefixSpan::Mine3 (std::vector<const Sequence*>& S)
{
	this->Sall = &S;

	/* Set max_item_index correctly
	 */
	max_item_index = 0;

	for (std::vector<const Sequence*>::iterator iv = S.begin () ;
		iv != S.end () ; ++iv)
	{
		for (unsigned int elem = 0 ; elem < (*iv)->lengthElements () ; ++elem) {
			for (unsigned int iidx = 0 ; iidx < (*iv)->lengthSet (elem) ; ++iidx) {
				unsigned int item = (*iv)->get (elem, iidx);
				if (item > max_item_index)
					max_item_index = item;
			}
		}
	}

	max_item_index += 1;

	/* The maximum number of deepening steps
	 */
	unsigned int maximum_depth = (options.maxlen == 0) ?
		std::numeric_limits<unsigned int>::max () : (options.maxlen - 1);

	/* Perform iterative deepening pattern search.
	 */
	usingIterativeDeepening = true;
	for (unsigned int id_depth = 0 ; id_depth < maximum_depth ; ++id_depth) {
		iterativeDeepeningLevel = id_depth+1;

		/* One mining iteration with a maximum depth:
		 *
		 * 1. Find the best pattern up to the given depth/length.
		 * 2. Use the best gain observed in previous iterations to prune on
		 *    future iterations.
		 *
		 * Note that (2.) happens automagically in the PrefixSpanLPBoost
		 * variants as the optimal pattern list is not cleared.
		 */
		SequenceT alphaStart;
		for (unsigned int n = 0 ; n < S.size () ; ++n)
			alphaStart.Occurence.insert (S[n]->Occurence.begin (),
				S[n]->Occurence.end ());

		seenLeafLevelExtensions = false;
		Mine (alphaStart, S, id_depth);

		/* Optimality check: observed possible extensions at the leaf-level?
		 *
		 *   yes -> we have to continue
		 *   no -> we can guarantee global optimality
		 */
		if (seenLeafLevelExtensions == false)
			break;
	}

	this->Sall = NULL;
}


void
PrefixSpan::MineIteration (void)
{
	if (eguide.IsEmpty ()) {
		std::cout << "eguide is empty." << std::endl;
		return;
	}

	ExplorationElement el = eguide.PopTop ();

#ifdef DEBUG
	std::cout << "MineIteration: gain " << el.gain << ", gainbound " << el.gainbound
		<< ", el.alpha = " << el.alpha << std::endl;
#endif

	/* 1. Find all frequent items across the sequences in S by a linear scan
	 * [item_index] holds the set of all active sequences in S
	 */
	std::map<unsigned int, std::set<unsigned int> > occ;
	std::map<unsigned int, std::set<unsigned int> > occ_;

	MineFrequentItems (el.alpha.lastItem (), el.alpha.lastElement (),
		occ, occ_, el.Sproj);

	FilterFrequentItems (occ);
	FilterFrequentItems (occ_);

	/* 2. Insert extensions into the global priority queue.
	 */
	for (occ_iterator fi = occ.begin () ; fi != occ.end () ; ++fi) {
		/* First handle all simple non joined extensions
		 */
		SequenceT alphaP = SequenceT (el.alpha);

		std::vector<unsigned int> elem;
		elem.push_back (fi->first);
		alphaP.appendElement (elem);

		indexSetJoin (el.alpha, el.Sproj, fi->second, alphaP.Occurence);

		Report (alphaP);
		if (TestExtendSequence (alphaP) == false)
			continue;

		/* Check gain/gainbound
		 * TODO: gainbound is computed twice here (TestExtendSequence also
		 * computes it)
		 */
		double gain = ExplorationPriority (alphaP);
		double gainbound = ExplorationPriorityBound (alphaP);

		// FIXME: here, gainbound = gain with a good gain will never reach
		// this point
		std::vector<const Sequence*> Sproj =
			prefixProject (fi->first, false, el.Sproj, alphaP.lastElement ());

		ExplorationElement eelem (gain, gainbound, alphaP, Sproj, true);
		eguide.Add (eelem);
	}

	for (occ_iterator fi = occ_.begin () ; fi != occ_.end () ; ++fi) {
		/* First handle all simple non joined extensions
		 */
		SequenceT alphaP = SequenceT (el.alpha);

		std::vector<unsigned int> elem;
		alphaP.appendItem (fi->first);

		indexSetJoin (el.alpha, el.Sproj, fi->second, alphaP.Occurence);

		Report (alphaP);
		if (TestExtendSequence (alphaP) == false)
			continue;

		/* Check gain/gainbound
		 * TODO: gainbound is computed twice here (TestExtendSequence also
		 * computes it)
		 */
		double gain = ExplorationPriority (alphaP);
		double gainbound = ExplorationPriorityBound (alphaP);

		// FIXME: here, gainbound = gain with a good gain will never reach
		// this point
		std::vector<const Sequence*> Sproj =
			prefixProject (fi->first, true, el.Sproj, alphaP.lastElement ());

		ExplorationElement eelem (gain, gainbound, alphaP, Sproj, true);
		eguide.Add (eelem);
	}

	if (el.cleanupS)
		el.Cleanup ();
}


/* At the time ::Mine is called, the sequence 'alpha' has been seen at least
 * min_support times.  alpha.Occurence.size() is the count of 'alpha' in the
 * projected set 'S'.
 *
 * This is the central PrefixSpan mining method, used for iterative-deepening.
 */
void
PrefixSpan::Mine (SequenceT& alpha, std::vector<const Sequence*>& S,
	unsigned int lengthTail)
{
	/* 1. Find all frequent items across the sequences in S by a linear scan
	 * [item_index] holds the set of all active sequences in S
	 */
	std::map<unsigned int, std::set<unsigned int> > occ;
	std::map<unsigned int, std::set<unsigned int> > occ_;

	MineFrequentItems (alpha.lastItem (), alpha.lastElement (),
		occ, occ_, S, alpha.lengthElements () == 0 ? -1 : options.maxgap);

	/* 2. Additional filtering step, unused for frequent subsequence mining.
	 *    Here additional exploration criterions, for example minimum support
	 * constraints in a sample subset can be implemented.
	 */
	FilterFrequentItems (occ);
	FilterFrequentItems (occ_);

	/* Ordered exploration.  This is not at all useful for complete frequent
	 * subsequence mining.  For weighted subsequence mining however, a greedy
	 * strategy of maximizing the gain can avoid exploration of the search
	 * space significantly.
	 */
	std::multimap<double, ExplorationPriorityElement, std::greater<double> >
		orderedExploration;

	for (occ_iterator fi = occ.begin () ; fi != occ.end () ; ++fi) {
		/* First handle all simple non joined extensions
		 */
		SequenceT alphaP = SequenceT (alpha);

		std::vector<unsigned int> elem;
		elem.push_back (fi->first);
		alphaP.appendElement (elem);

		indexSetJoin (alpha, S, fi->second, alphaP.Occurence);
		if (doVerify)
			verifyOccurenceCorrectness (alphaP);

		orderedExploration.insert (
			std::pair<double, ExplorationPriorityElement> (
				ExplorationPriority (alphaP),
				ExplorationPriorityElement (alphaP, fi->first, false)));
	}

	for (occ_iterator fi = occ_.begin () ; fi != occ_.end () ; ++fi) {
		/* We now deal with the case of frequent items of kind '_b'.
		 */
		if (alpha.checkAppendPossible (fi->first) == false)
			continue;

		SequenceT alphaP = SequenceT (alpha);
		alphaP.appendItem (fi->first);

		indexSetJoin (alpha, S, fi->second, alphaP.Occurence);
		if (doVerify)
			verifyOccurenceCorrectness (alphaP);

		orderedExploration.insert (
			std::pair<double, ExplorationPriorityElement> (
				ExplorationPriority (alphaP),
				ExplorationPriorityElement (alphaP, fi->first, true)));
	}

	/* Prioritized exploration of subsequences
	 */
	for (std::multimap<double, ExplorationPriorityElement,
			std::greater<double> >::iterator iv = orderedExploration.begin () ;
		iv != orderedExploration.end () ; ++iv)
	{
		/* This check is trivially false whenever min_support is fixed, but
		 * for the top-K scheme, this can be important.
		 */
		if (iv->second.seq.Occurence.size () < options.minsup)
			continue;

		/* If we use iterative deepening, only report leaf level patterns.
		 * Otherwise, always report.
		 */
		if (usingIterativeDeepening == false || lengthTail == 0)
			Report (iv->second.seq);

		bool testExtendSequenceResult;
		if (usingIterativeDeepening) {
			/* Iterative deepening: check for optimality
			 */
			testExtendSequenceResult = TestExtendSequence (iv->second.seq);

			if (lengthTail == 0 && testExtendSequenceResult)
				seenLeafLevelExtensions = true;
		} else {
			/* DFS mining
			 */
			if (lengthTail == 0)
				continue;

			testExtendSequenceResult = TestExtendSequence (iv->second.seq);
		}
		if (lengthTail == 0 || testExtendSequenceResult == false)
			continue;

		/* It is clear we should recurse into the projections, so explicitly
		 * create them and mine them.
		 */
		std::vector<const Sequence*> Sproj;
		if (options.maxgap < 0 || iv->second.itemJoined) {
			Sproj = prefixProject (iv->second.item,
				iv->second.itemJoined, S, iv->second.seq.lastElement (),
				options.maxgap < 0 ? false : true);
		} else {
			/* For standard b extensions and maximum gap subsequence mining,
			 * project on all occurences of the item in the next maxgap
			 * elements of the sequence set.
			 *
			 * An exception is the first step, where we are required to
			 * project on _all_ occurences in the sequence
			 */
			Sproj = prefixProjectGap (iv->second.item, S, iv->second.seq,
				iv->second.seq.lastElement (),
				alpha.lengthElements () == 0 ?
					std::numeric_limits<int>::max () : options.maxgap);
		}

		Mine (iv->second.seq, Sproj, lengthTail-1);

		for (seq_iterator si = Sproj.begin () ; si != Sproj.end () ; ++si)
			delete (*si);
	}
}


bool
PrefixSpan::verifyOccurenceCorrectness (const Sequence& alpha) const
{
	++verifyCount;
	if (verifyCount % 1000 == 1) {
		std::cout << "[" << verifyCount
			<< "] PrefixSpan::verifyOccurenceCorrectness" << std::endl;
	}

	/* Expensive verification for debugging:
	 * We explicitly check if our test function agrees on the occurence function.
	 */
	for (unsigned int n = 0 ; n < (*Sall).size () ; ++n) {
		bool occContains = (alpha.Occurence.find (n) != alpha.Occurence.end ());
		bool seqContains;
		if (options.maxgap < 0)
			seqContains = ((*Sall)[n])->contains (alpha);
		else
			seqContains = ((*Sall)[n])->containsGap (alpha, options.maxgap);

		/* They agree: consistent
		 */
		if (occContains == seqContains)
			continue;

		/* FATAL: our test function disagrees with PrefixSpan's internal
		 * state.
		 */
		std::cout << "FATAL: Inconsistent subseteq relation (maxgap="
			<< options.maxgap << ")" << std::endl;
		std::cout << "sub sequence: " << alpha << std::endl;
		std::cout << "sup sequence (#" << n << "): " << *(*Sall)[n] << std::endl;
		std::cout << "occContains: " << (occContains ? "TRUE" : "FALSE") << std::endl;
		std::cout << "seqContains: " << (seqContains ? "TRUE" : "FALSE") << std::endl;
		std::cout << std::endl;
		std::cout << "alpha.Occurence:" << std::endl;

		for (std::set<unsigned int>::iterator iv = alpha.Occurence.begin () ;
			iv != alpha.Occurence.end () ; ++iv)
		{
			std::cout << *iv << " ";
		}
		std::cout << std::endl;

		std::cout << "seq.contains occurences:" << std::endl;
		for (unsigned int n = 0 ; n < (*Sall).size () ; ++n) {
			bool seqContains;
			if (options.maxgap < 0)
				seqContains = ((*Sall)[n])->contains (alpha);
			else
				seqContains = ((*Sall)[n])->containsGap (alpha, options.maxgap);

			if (seqContains)
				std::cout << n << " ";
		}
		std::cout << std::endl;

		assert (0);
		return (false);
	}
	return (true);
}


/* Project a set of sequences by projecting each sequence individually, if
 * possible.
 */

/*
bool
PrefixSpan::verifyOccurenceCorrectnessBackup (const Sequence& alpha) const
{
	++verifyCount;
	if (verifyCount % 1000 == 1) {
		std::cout << "[" << verifyCount
			<< "] PrefixSpan::verifyOccurenceCorrectness" << std::endl;
	}

	for (unsigned int n = 0 ; n < (*Sall).size () ; ++n) {
		bool occContains = (alpha.Occurence.find (n) != alpha.Occurence.end ());
		bool seqContains;
		if (options.maxgap < 0)
			seqContains = ((*Sall)[n])->contains (alpha);
		else
			seqContains = ((*Sall)[n])->containsGap (alpha, options.maxgap);

		if (occContains == seqContains)
			continue;

		std::cout << "FATAL: Inconsistent subseteq relation (maxgap="
			<< options.maxgap << ")" << std::endl;
		std::cout << "sub sequence: " << alpha << std::endl;
		std::cout << "sup sequence (#" << n << "): " << *(*Sall)[n] << std::endl;
		std::cout << "occContains: " << (occContains ? "TRUE" : "FALSE") << std::endl;
		std::cout << "seqContains: " << (seqContains ? "TRUE" : "FALSE") << std::endl;
		std::cout << std::endl;
		std::cout << "alpha.Occurence:" << std::endl;

		for (std::set<unsigned int>::iterator iv = alpha.Occurence.begin () ;
			iv != alpha.Occurence.end () ; ++iv)
		{
			std::cout << *iv << " ";
		}
		std::cout << std::endl;

		std::cout << "seq.contains occurences:" << std::endl;
		for (unsigned int n = 0 ; n < (*Sall).size () ; ++n) {
			bool seqContains;
			if (options.maxgap < 0)
				seqContains = ((*Sall)[n])->contains (alpha);
			else
				seqContains = ((*Sall)[n])->containsGap (alpha, options.maxgap);

			if (seqContains)
				std::cout << n << " ";
		}
		std::cout << std::endl;

		assert (0);
		return (false);
	}
	return (true);
}
*/

/* Project a set of sequences by projecting each sequence individually, if
 * possible.
 */
std::vector<const Sequence*>
PrefixSpan::prefixProject (unsigned int item, bool itemJoined,
	std::vector<const Sequence*> S,
	const std::vector<unsigned int>* lastPrefixElement,
	bool maximumGapMining) const
{
	std::vector<const Sequence*> res = std::vector<const Sequence*> ();

	for (std::vector<const Sequence*>::iterator si = S.begin () ;
		si != S.end () ; ++si)
	{
		/* At most a single projection
		 */
		const Sequence* proj = (*si)->prefixProject (item, itemJoined,
			lastPrefixElement, maximumGapMining);

		if (proj == NULL)
			continue;

		if (proj->length () == 0) {
			delete proj;
			continue;
		}

		res.push_back (proj);
	}

	return (res);
}

std::vector<const Sequence*>
PrefixSpan::prefixProjectGap (unsigned int item,
	std::vector<const Sequence*> S, const Sequence& subseq,
	const std::vector<unsigned int>* lastPrefixElement, unsigned int ahead) const
{
	std::vector<const Sequence*> res;

	for (std::vector<const Sequence*>::iterator si = S.begin () ;
		si != S.end () ; ++si)
	{
		/* Imagine a sequence like [a][b][c][d][e][d][f][g] and a maximum gap
		 * length of 2 (so gaps of length 0, 1 or 2 are allowed).  If we have
		 * the current prefix of [a][b][c] and we project on [d], we should
		 * obtain:
		 *     [a][b][c][d].[e][d][f][g]  _and_
		 *     [a][b][c][d][e][d].[f][g],
		 * because otherwise we might miss the item [g] in a follow-up
		 * projection.  If we only have the first sequence, we would not be
		 * able to project on [g], whereas with the second one we can project.
		 *
		 * Hence: the projections are _all_ the sequences which can be
		 *    projected within the maximum gap constraint.
		 */
		std::vector<const Sequence*> proj =
			(*si)->prefixProjectGap (item, lastPrefixElement, ahead);

		if (proj.size () == 0)
			continue;

#if 0
		// FIXME:
		if ((*si)->Occurence.size () == 1 && *((*si)->Occurence.begin ()) == 8
			&& (item == 671 || item == 752))
		{
			// DEBUG
			std::cout << "Projecting sequence #8, current form:" << std::endl;
			std::cout << "  sub: " << subseq << std::endl;
			std::cout << "  item: " << item << ", ahead " << ahead << std::endl;

			std::cout << "  lastPrefixElement: ";
			for (std::vector<unsigned int>::const_iterator civ = lastPrefixElement->begin () ;
				civ != lastPrefixElement->end () ; ++civ)
				std::cout << *civ << ", ";
			std::cout << std::endl;

			std::cout << "  #8: " << *(*si) << std::endl;
			std::cout << std::endl;

			std::cout << "Resulting in " << proj.size () << " projections:" << std::endl;
			for (std::vector<const Sequence*>::iterator pi = proj.begin () ;
				pi != proj.end () ; ++pi)
			{
				std::cout << "   (len " << (*pi)->length () << "): "
					<< *(*pi) << std::endl;
			}
		}
#endif

		for (std::vector<const Sequence*>::iterator pi = proj.begin () ;
			pi != proj.end () ; )
		{
			if ((*pi)->length () == 0) {
				delete (*pi);
				pi = proj.erase (pi);
			} else
				++pi;
				
						
		}

        /*std::cout << "Resulting in " << proj.size () << " projections:" << std::endl;
		for (std::vector<const Sequence*>::iterator pi = proj.begin () ;
			pi != proj.end () ; ++pi)
		{
			std::cout << "   (len " << (*pi)->length () << "): "
				<< *(*pi) << std::endl;
		}*/
		res.insert (res.end (), proj.begin (), proj.end ());
	}

	return (res);
}


bool
PrefixSpan::TestExtendSequence (const Sequence& seq) const
{
	/* This is just frequent subsequence mining, hence we want to extend all
	 * sequences.
	 */
	return (true);
}


void
PrefixSpan::FilterFrequentItems (
	std::map<unsigned int, std::set<unsigned int> >& occ) const
{
	return;
}


void
PrefixSpan::Report (const Sequence& seq)
{
	/* Expensive verification of correctness
	 */
	if (doVerify)
		verifyOccurenceCorrectness (seq);

	/* Check minimum length requirement.
	 */
	if (seq.length () < options.minlen)
		return;

	unsigned int freq = seq.Occurence.size ();
	subsequences.push_back (std::pair<unsigned int, SequenceT> (freq, seq));
}


double
PrefixSpan::ExplorationPriority (const Sequence& seq) const
{
	/* Frequent mining: we don't care about priorities
	 */
	return (0.0);
}

double
PrefixSpan::ExplorationPriorityBound (const Sequence& seq) const
{
	/* Frequent mining: trivial bound, examine everything
	 */
	return (1.0);
}

/* This function takes one pass through the entire set S.
 *
 * Note that the output indices in 'occ' and 'occ_' are indices to S, not
 * absolute sequence indices.
 */
void
PrefixSpan::MineFrequentItems (unsigned int lastItem,
	const std::vector<unsigned int>* lastPrefixElement,
	std::map<unsigned int, std::set<unsigned int> >& occ,
	std::map<unsigned int, std::set<unsigned int> >& occ_,
	std::vector<const Sequence*>& S,
	int maximumGapLength) const
{
	occ.clear ();
	occ_.clear ();

	/* All sequences in 'S' are projected off the same prefix.  Hence, our
	 * goal here is to find frequent items in this shared-prefix sequences.
	 *
	 * There are two kinds of items,
	 *
	 *    i) "vanilla items", which are occuring _after_ the element of the
	 *       sequence projection point,
	 *   ii) "_b items", which are occuring in the same element as the
	 *       sequence's projection point.
	 *
	 * Frequent vanilla items produce a new element in the new prefix
	 * sequence, such as if "c" is a new frequent vanilla item and the current
	 * prefix is [a][b], the resulting prefix is [a][b][c].  "_c items" join
	 * into the last element of the current prefix, such that the resulting
	 * prefix would be [a][bc].
	 *
	 * Here we also handle maximum-gap constraints (maxgap G) by
	 *    i) projecting in the first iteration after _each_ occurence, not
	 *       just the first
	 *   ii) in all later iterations only check the next G elements, and
	 *       project on all occurences there
	 */
	unsigned int sidx = 0;
	for (std::vector<const Sequence*>::iterator iv = S.begin () ;
		iv != S.end () ; ++iv, ++sidx)
	{
		/* Special handling for (_b) elements: in case this sequence has its
		 * split point in the interior of its first element
		 * (firstElementSplitted).  Then, all remaining items in the first
		 * element are in the _b occurence set.
		 */
		if ((*iv)->firstElementSplitted ()) {
			for (unsigned int iidx = 0 ; iidx < (*iv)->lengthSet (0) ; ++iidx) {
				if ((*iv)->get (0, iidx) > lastItem)
					occ_[(*iv)->get (0, iidx)].insert (sidx);
			}
		}

		/* Normal vanilla items, the easiest, as they open a new element.
		 * Hence we only need to be careful for the very first element in the
		 * sequence.
		 *
		 * We also handle _b items here.  A linear scan is performed for each
		 * element and if the current prefix's last element is a subset of the
		 * current element, then all remaining items are added as _b items.
		 */
		unsigned int eidx = (*iv)->firstElementSplitted () ? 1 : 0;

		/* Maximum look ahead (infinite if arbitrary gap length are allowed).
		 * If we are dealing with a split element, then there is no gap, even
		 * if we start with eidx=1, so we compensate for that here.
		 */
		unsigned int ahead = (maximumGapLength < 0) ?
			std::numeric_limits<int>::max () : (maximumGapLength + eidx);

		/* Stopping condition: maximum look-ahead or sequence length exceeded.
		 */
		for ( ; eidx <= ahead && eidx < (*iv)->lengthElements () ; ++eidx) {
			bool matchedPrefixElement = false;
			unsigned int p_iidx = 0;

			for (unsigned int iidx = 0 ; iidx < (*iv)->lengthSet (eidx) ; ++iidx) {
				unsigned int item = (*iv)->get (eidx, iidx);

				occ[item].insert (sidx);

				/* This is again a little subtle: for the normal subsequence
				 * mining with arbitrary long gaps, we have to check
				 * subsequence elements for _b matches.  For maximum gap
				 * mining, these have already been handled by multiple
				 * projections.
				 */
				if (maximumGapLength < 0) {
					/* If we already found the last prefix element within this
					 * element, then we copy all remaining items into the _b set.
					 */
					if (matchedPrefixElement)
						occ_[item].insert (sidx);

					/* Continue linear match in case it has not been matched
					 * yet.
					 */
					if (matchedPrefixElement == false &&
						lastPrefixElement != NULL &&
						item == (*lastPrefixElement)[p_iidx])
					{
						p_iidx += 1;

						if (p_iidx == lastPrefixElement->size ())
							matchedPrefixElement = true;
					}
				}
			}
		}
	}

	/* Discard non frequent occurence maps.
	 */
	cleanNonFrequent (occ);
	cleanNonFrequent (occ_);
}



void
PrefixSpan::cleanNonFrequent
	(std::map<unsigned int, std::set<unsigned int> >& occ) const
{
	// TODO: replace with
	// occ.erase (remove_if (occ.begin (), occ.end (), removePredicate), occ.end ());

	for (occ_iterator oi = occ.begin () ; oi != occ.end () ; ) {
		/* Check set size for minimum support
		 */
		if (oi->second.size () < options.minsup) {
			occ.erase (oi++);
		} else {
			++oi;
		}
	}
}

}


