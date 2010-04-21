
#ifndef	SEQUENCE_H
#define	SEQUENCE_H

#include <vector>
#include <set>
#include <functional>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include <string>

#include <assert.h>
#include <unistd.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace SeqMine {

class
Sequence
{
public:
	virtual ~Sequence (void) { }

	/* Public members
	 */

	/* The set of sequence indices where this particular sequence appears in.
	 * Semantics are as follows: i) the empty sequence appears in every
	 * sequence, ii) at start, the sequence appears only in itself.
	 */
	std::set<unsigned int> Occurence;

	/* lengthElements: The number of elements.
	 * lengthSet: The number of elements in the given set.
	 */
	virtual unsigned int lengthElements () const = 0;
	virtual unsigned int lengthSet (unsigned int eidx) const = 0;

	virtual unsigned int get (unsigned int eidx, unsigned int iidx) const = 0;

#if 0
	/* Sort by, i) total length, ii) lexicographically
	 */
	virtual bool operator< (const Sequence& seq) const
	{
		unsigned int l1 = length ();
		unsigned int l2 = seq.length ();

		if (l1 != l2)
			return (l1 < l2);

		/* Case: equal length sequences
		 */
		unsigned int eidx2 = 0;
		unsigned int iidx2 = 0;
		for (unsigned int eidx1 = 0 ; eidx1 < lengthElements () ; ++eidx1, ++eidx2)
		{
			if (lengthSet (eidx1) != seq.lengthSet (eidx2))
				return (lengthSet (eidx1) < seq.lengthSet (eidx2));

			iidx2 = 0;
			for (unsigned int iidx1 = 0 ; iidx1 < lengthSet (eidx1) ;
				++iidx1, ++iidx2)
			{
				if (get (eidx1, iidx1) != seq.get (eidx2, iidx2))
					return (get (eidx1, iidx1) < seq.get (eidx2, iidx2));
			}
		}

		return (true);
	}
#endif

	virtual void writeSequenceFile (const char* filename) const
	{
		std::ofstream pout (filename);

		for (unsigned int eidx = 0 ; eidx < lengthElements () ; ++eidx) {
			for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
				if (iidx > 0)
					pout << " ";

				pout << get (eidx, iidx);
			}
			pout << std::endl;
		}
		pout.close ();
	}

	/* length: The total length of the sequence, linearized.
	 */
	virtual unsigned int length () const
	{
		unsigned int el = lengthElements ();
		unsigned int count = 0;

		for (unsigned int eidx = 0 ; eidx < el ; ++eidx)
			count += lengthSet (eidx);

		return (count);
	}

	virtual unsigned int lastItem (void) const
	{
		assert (lengthElements () > 0);
		assert (lengthSet (lengthElements () - 1) > 0);

		return (get (lengthElements () - 1,
			lengthSet (lengthElements () - 1) - 1));
	}

	virtual unsigned int indexAfter (unsigned int eidx, unsigned int item) const
	{
		for (unsigned int n = 0 ; n < lengthSet (eidx) ; ++n)
			if (get (eidx, n) == item)
				return (n + 1);

		return (lengthSet (eidx));
	}

	virtual bool contains (unsigned int eidx, unsigned int item,
		unsigned int start_iidx) const
	{
		for (unsigned int n = start_iidx ; n < lengthSet (eidx) ; ++n) {
			if (get (eidx, n) == item)
				return (true);
		}

		return (false);
	}

	virtual bool contains (unsigned int eidx, const Sequence& pat,
		unsigned int peidx) const
	{
		for (unsigned int n = 0 ; n < pat.lengthSet (peidx) ; ++n) {
			if (contains (eidx, pat.get (peidx, n), 0) == false)
				return (false);
		}

		return (true);
	}

	/* Recursive function checking subsequence presence with maximum gap
	 * constraint
	 */
	virtual bool containsGap_ (unsigned int beidx, unsigned bss_eidx,
		const Sequence& subseq, unsigned int maxgap) const
	{
		/* Recursion termination: the complete subsequence has been matched
		 */
		if (bss_eidx >= subseq.lengthElements () && bss_eidx >= subseq.lengthElements ())
			return (true);

		/* Two cases
		 *    i) beidx = 0, then we recurse into all element matches
		 *    ii) beidx > 0, we recurse only into the next maxgap+1 element
		 *        matches
		 *
		 * In any case, this method terminates early in case a match has been
		 * found.  TODO: modify this so a density of matches is returned over
		 * the elements.
		 */
		unsigned int ss_eidx = bss_eidx;
		if (beidx == 0) {
			for (unsigned int eidx = 0 ; eidx < lengthElements () ; ++eidx) {
				unsigned int ss_iidx = 0;

				for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
					if (get (eidx, iidx) != subseq.get (ss_eidx, ss_iidx))
						continue;

					/* One single item is matched, advance
					 */
					ss_iidx += 1;
					if (ss_iidx == subseq.lengthSet (ss_eidx))
						goto recurse_seq;
				}
				continue;

			recurse_seq:
				/* Recurse into each first match
				 */
				bool sres = containsGap_ (eidx+1, 1, subseq, maxgap);
				if (sres)
					return (true);
			}
		} else {
			for (unsigned int eidx = beidx ;
				(eidx-beidx) <= maxgap && eidx < lengthElements () ; ++eidx)
			{
				unsigned int ss_iidx = 0;

				for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
					if (get (eidx, iidx) != subseq.get (ss_eidx, ss_iidx))
						continue;

					/* One single item is matched, advance
					 */
					ss_iidx += 1;
					if (ss_iidx == subseq.lengthSet (ss_eidx))
						goto recurse_seq2;
				}
				continue;

			recurse_seq2:
				/* Recurse into each maxgap-constrained match
				 */
				bool sres = containsGap_ (eidx+1, ss_eidx+1, subseq, maxgap);
				if (sres)
					return (true);
			}
		}

		return (false);
	}

	/* 'contains' subsequence predicate with maximum gap length constraints
	 */
	virtual bool containsGap (const Sequence& subseq, unsigned int maxgap) const
	{
		return (containsGap_ (0, 0, subseq, maxgap));
	}

	/* 'contains' subsequence predicate, matches arbitrary long gaps.
	 *
	 * If 'pmatch' is non-NULL, it will contain the index-to-index mapping
	 */
	virtual bool contains (const Sequence& subseq,
		std::vector<unsigned int>* pmatch = NULL) const
	{
		if (pmatch != NULL)
			pmatch->clear ();

		/* Check if 'subseq' appears in this instance.
		 */
		unsigned int ss_eidx = 0;

		for (unsigned int eidx = 0 ; eidx < lengthElements () ; ++eidx) {
			unsigned int ss_iidx = 0;

			/* Compare two single sets
			 */
			for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
				if (get (eidx, iidx) != subseq.get (ss_eidx, ss_iidx))
					continue;

				/* One single item is matched, advance
				 */
				ss_iidx += 1;
				if (ss_iidx == subseq.lengthSet (ss_eidx))
					break;
			}

			if (ss_iidx == subseq.lengthSet (ss_eidx)) {
				/* Record index matches
				 */
				if (pmatch != NULL)
					pmatch->push_back (eidx);

				/* This elment was matched, advance in sequence
				 */
				ss_eidx += 1;

				if (ss_eidx == subseq.lengthElements ())
					return (true);	// matched
			}
		}

		if (pmatch != NULL)
			pmatch->clear ();

		return (false);
	}

	/* Test for identity of two sequences.
	 */
	virtual bool identical (const Sequence& s2) const
	{
		if (lengthElements () != s2.lengthElements () ||
			length () != s2.length ())
			return (false);

		for (unsigned int n = 0 ; n < lengthElements () ; ++n) {
			if (lengthSet (n) != s2.lengthSet (n))
				return (false);

			for (unsigned int k = 0 ; k < lengthSet (n) ; ++k)
				if (get (n, k) != s2.get (n, k))
					return (false);
		}

		return (true);
	}

	virtual bool identicalSet (unsigned int eidx1, const Sequence& s2,
		unsigned int eidx2) const
	{
		if (lengthSet (eidx1) != s2.lengthSet (eidx2))
			return (false);

		for (unsigned int k = 0 ; k < lengthSet (eidx1) ; ++k)
			if (get (eidx1, k) != s2.get (eidx2, k))
				return (false);

		return (true);
	}

	virtual bool firstElementSplitted (void) const
	{
		return (false);
	}

	virtual const Sequence* prefixProject (unsigned int item, bool itemJoined,
		const std::vector<unsigned int>* lastPrefixElement,
		bool maximumGapMining = false) const;
	virtual std::vector<const Sequence*> prefixProjectGap (
		unsigned int item, const std::vector<unsigned int>* lastPrefixElement,
		unsigned int ahead) const;
};


/* Real sequence, storing the actual elements.
 * The sequence is mutable through appendElement and deleteElement.
 */
class
SequenceT : public Sequence
{
private:
	typedef std::vector<unsigned int>::iterator vec_iterator;
	typedef std::vector<unsigned int>::const_iterator cvec_iterator;

	std::vector<std::vector<unsigned int> > seq;

	bool _firstElementSplitted;

	void copySequence (const Sequence& source)
	{
		seq.clear ();
		for (unsigned int eidx = 0 ; eidx < source.lengthElements () ; ++eidx) {
			std::vector<unsigned int> curel (source.lengthSet (eidx));

			for (unsigned int iidx = 0 ; iidx < source.lengthSet (eidx) ; ++iidx)
				curel[iidx] = source.get (eidx, iidx);

			seq.push_back (curel);
		}
	}

	/* Serialization
	 */
    friend class boost::serialization::access;

    template<class Archive>
    void serialize (Archive & ar, const unsigned int version)
    {
		ar & seq;
		ar & _firstElementSplitted;
    }

public:
	SequenceT (void)
		: _firstElementSplitted (false)
	{
	}

	SequenceT (unsigned int seqIndex)
		: _firstElementSplitted (false)
	{
		Occurence.insert (seqIndex);
	}

	SequenceT (const Sequence& source, unsigned int seqIndex)
		: _firstElementSplitted (false)
	{
		copySequence (source);

		Occurence.clear ();
		Occurence.insert (seqIndex);
	}

	/* Produce a vanilla T sequence from any kind of sequence.  (That is,
	 * project it onto a canonical representation.)
	 */
	SequenceT (const Sequence& source)
	{
		Occurence = source.Occurence;

		copySequence (source);
		_firstElementSplitted = source.firstElementSplitted ();
	}

	SequenceT (const Sequence& source, unsigned int eidx, unsigned int iidx)
	{	// TODO: optimize this, this is very expensive!
		seq.clear ();
		_firstElementSplitted = (iidx > 0);

		for ( ; eidx < source.lengthElements () ; ++eidx) {
			std::vector<unsigned int> curel (source.lengthSet (eidx));

			for ( ; iidx < source.lengthSet (eidx) ; ++iidx)
				curel[iidx] = source.get (eidx, iidx);

			iidx = 0;
			seq.push_back (curel);
		}

		Occurence = source.Occurence;
	}

	virtual ~SequenceT (void) { }

	unsigned int lengthElements () const;
	unsigned int lengthSet (unsigned int eidx) const;
	unsigned int lastItem (void) const;
	unsigned int get (unsigned int eidx, unsigned int iidx) const;

	/* Override, as we can do better with a set.
	 */
	bool contains (unsigned int eidx, unsigned int item,
		unsigned int start_iidx) const;
	bool contains (unsigned int eidx,
		std::vector<unsigned int> elem,
		unsigned int start_iidx) const;

	/* Append a new element (item set) to the existing sequence.
	 *
	 * Note: All items in 'elem' must be in lexicographic order and no
	 *       duplications are allowed.
	 */
	void appendElement (std::vector<unsigned int>& elem);

	/* Delete the last element (item set) from the sequence.
	 */
	void deleteElement (void);

	/* Append a single item to the last element in this sequence.
	 *
	 * Note: we assume only items are inserted which pass the
	 *       'checkAppendPossible', that is, those which come in lexicographic
	 *       order after the last current item.
	 */
	void appendItem (unsigned int item);

	/* Delete the last item in the last element of this sequence.
	 */
	void deleteItem (void);

	/* Return a pointer to the last element vector.  (Used for frequent mining
	 * and prefix projection).
	 */
	const std::vector<unsigned int>* lastElement (void) const;

	/* Check whether alpha' would be a valid sequence, where alpha' is the
	 * current sequence with 'item' appended to the last element.
	 */
	bool checkAppendPossible (unsigned int item) const;

	virtual bool firstElementSplitted (void) const
	{
		return (_firstElementSplitted);
	}

	/* Read a sequence from a text file, in the set-row format.
	 */
	static SequenceT* readFromFile (const char* filename);
};


/* Pseudoprojection sequence, storing only a reference to the original
 * sequence and an offset.
 */
// TODO: declare this class const
class
SequenceP : public Sequence
{
protected:
	const Sequence& seqbase;	// reference sequence
	unsigned int ebase;	// base element index
	unsigned int ibase;	// base item index in first element

public:
	/* Create a new pseudoprojected sequence from the real sequence baseseq at
	 * position (eidx,iidx).
	 *
	 * Use underlying SequenceT type so as to only use one indirection
	 * step.
	 */
	SequenceP (const SequenceP& baseseq, unsigned int eidx, unsigned int iidx)
		: seqbase (baseseq.seqbase)
	{
		/* Adjust the indices depending on movement size
		 */
		ebase = baseseq.ebase;
		if (eidx == 0) {
			/* Moving within the same element.
			 */
			ibase = baseseq.ibase + iidx;
		} else {
			/* Moving one element or more.
			 */
			ebase += eidx;
			ibase = iidx;
		}

		/* Note: the Occurence set is copied over to SequenceP, which may not
		 * be what you want in general.
		 */
		Occurence = baseseq.Occurence;
	}

	SequenceP (const Sequence& baseseq, unsigned int eidx, unsigned int iidx)
		: seqbase (baseseq), ebase (eidx), ibase (iidx)
	{
		/* Note: the Occurence set is copied over to SequenceP, which may not
		 * be what you want in general.
		 */
		Occurence = baseseq.Occurence;
	}

	unsigned int lengthElements () const
	{
		return (seqbase.lengthElements () - ebase);
	}

	unsigned int lengthSet (unsigned int eidx) const
	{
		assert (eidx < lengthElements ());

		/* The first set has to be dealt with differently: we have to consider
		 * the starting index.
		 */
		if (eidx == 0)
			return (seqbase.lengthSet (ebase) - ibase);

		return (seqbase.lengthSet (ebase + eidx));
	}

	unsigned int get (unsigned int eidx, unsigned int iidx) const
	{
		/*std::cout << "SequenceP::get (eidx: "
			<< eidx << ", iidx: " << iidx << ")" << std::endl;
		*/

		if (eidx == 0)
			return (seqbase.get (ebase, iidx + ibase));

		return (seqbase.get (ebase + eidx, iidx));
	}

	virtual bool firstElementSplitted (void) const
	{
		return (ibase > 0);
	}
};

};

std::ostream& operator<< (std::ostream& os, const SeqMine::Sequence& seq);
std::ostream& operator<< (std::ostream& os,
	const std::vector<const SeqMine::Sequence*>& seqs);

std::istream& operator>> (std::istream& is, SeqMine::SequenceT& seq);

struct sequence_hash
{
	size_t operator() (const SeqMine::Sequence& seq) const
	{
		unsigned long long code = 2147483647;

		for (unsigned int eidx = 0 ; eidx < seq.lengthElements () ; ++eidx) {
			unsigned int leidx = seq.lengthSet (eidx);

			for (unsigned int iidx = 0 ; iidx < leidx ; ++iidx) {
				unsigned int citem = seq.get (eidx, iidx);

				code *= citem;
				code += ~citem;
				code %= 2147383647;
			}
		}

		return ((size_t) code);
	}
};

struct sequence_identical
{
	bool operator() (const SeqMine::Sequence& s1, const SeqMine::Sequence& s2) const
	{
		return (s1.identical (s2));
	}
};

#endif


