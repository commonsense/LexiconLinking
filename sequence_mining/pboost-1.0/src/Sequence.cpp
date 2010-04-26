
#include <algorithm>
#include <iostream>
#include <cstring>
#include <iterator>
#include <strstream>	// TODO someday: replace with sstream (stringstring)
#include <list>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>

#include <stdlib.h>
//#include <Lexicon.cpp>
#include <Sequence.h>

namespace SeqMine
{

/* Do a prefix projection, if possible.
 *
 * item: The new prefix item T,
 * itemJoined: If true, the prefix item is not standing freely but is the last
 *    item of the prefix element (...T),
 * lastPrefixElement: contains the entire _new_ prefix, _including_ item.
 *
 * Return NULL if no projection was possible/found.
 */
const Sequence*
Sequence::prefixProject (unsigned int item, bool itemJoined,
	const std::vector<unsigned int>* lastPrefixElement,
	bool maximumGapMining) const
{
	assert (lastPrefixElement != NULL);

#ifdef DEBUG
	std::cout << "[" << *this << "].prefixProject (" << item << ", "
		<< itemJoined << ")" << std::endl;
#endif

	/* Case A: a _b projection.
	 */
	if (itemJoined) {
		/* The prefix is already contained in this element, hence we only need
		 * to match 'item'.
		 */
		if (firstElementSplitted ()) {
			for (unsigned int iidx = 0 ; iidx < lengthSet (0) ; ++iidx) {
				unsigned int cur_item = get (0, iidx);

				if (cur_item < item)
					continue;	// -> a match is still possible
				else if (cur_item > item)
					break;	// -> no match is possible

				/* We got a correct matching, do the projection.
				 *
				 * There are, in principle, two prefix projections points '|':
				 *   i) (...T)|rest
				 *  ii) (...T|cdef)rest
				 */
				if ((iidx+1) == lengthSet (0)) {	// case i)
					if (lengthElements () > 1)		// is there a rest?
						return (new SequenceT (*this, 1, 0));
					else
						return (NULL);	// end of sequence
				} else
					return (new SequenceT (*this, 0, iidx+1));
			}
		}

		/* For maximum-gap mining, we can safely ignore the remaining part, as
		 * the projections were carried out explicitly
		 */
		if (maximumGapMining)
			return (NULL);

		/* All remaining elements: first find last prefix element in it, then
		 * item.
		 */
		for (unsigned int eidx = firstElementSplitted () ? 1 : 0 ;
			eidx < lengthElements () ; ++eidx)
		{
			unsigned int pp_iidx = 0;

			for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
				unsigned int cur_item = get (eidx, iidx);

				if (cur_item < (*lastPrefixElement)[pp_iidx])
					continue;
				else if (cur_item > (*lastPrefixElement)[pp_iidx])
					break;

				/* Matched one more element of the prefix's last element.
				 */
				pp_iidx += 1;
				if (pp_iidx < lastPrefixElement->size ())
					continue;

				/* Final match, do projection
				 *
				 * There are, in principle, two prefix projections points '|':
				 *   i) (...T)|rest
				 *  ii) (...T|cdef)rest
				 */
				if ((iidx+1) == lengthSet (eidx)) {	// case i)
					if (lengthElements () > (eidx+1))		// is there a rest?
						return (new SequenceT (*this, eidx+1, 0));
					else
						return (NULL);	// end of sequence
				} else	// case ii)
					return (new SequenceT (*this, eidx, iidx+1));
			}
		}
		return (NULL);
	}

	/* Case B: a vanilla b projection.
	 */
	unsigned int eidx = firstElementSplitted () ? 1 : 0;
	for ( ; eidx < lengthElements () ; ++eidx) {
		for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
			unsigned int citem = get (eidx, iidx);

			if (citem < item)
				continue;	// keep searching

			if (citem > item)
				break;

			/* Found sequence, return simple projection
			 */
			iidx += 1;
			if (iidx == lengthSet (eidx)) {
				eidx += 1;
				iidx = 0;
			}
			if (eidx >= lengthElements ())
				return (NULL);	// end of sequence reached

			/* A simple projection
			 */
			return (new SequenceT (*this, eidx, iidx));
		}
	}

	return (NULL);
}





std::vector<const Sequence*>
Sequence::prefixProjectGap (unsigned int item,
	const std::vector<unsigned int>* lastPrefixElement, unsigned int ahead) const
{
	assert (lastPrefixElement != NULL);

	/* Here we only need to handle vanilla b projections, however, there is a
	 * twist: we need to project on all sequences which match in a given
	 * look-ahead (ahead).
	 */
	unsigned int eidx = firstElementSplitted () ? 1 : 0;
	ahead += eidx;

	std::vector<const Sequence*> res;

	for ( ; eidx <= ahead && eidx < lengthElements () ; ++eidx) {
		unsigned int leidx = lengthSet (eidx);

		for (unsigned int iidx = 0 ; iidx < leidx ; ++iidx) {
			unsigned int citem = get (eidx, iidx);

			if (citem < item)
				continue;	// keep searching

			if (citem > item)
				break;

			/* Found sequence, carry out simple projection
			 */
			{
				unsigned int p_iidx = iidx + 1;
				unsigned int p_eidx = eidx;

				if (p_iidx == lengthSet (p_eidx)) {
					p_eidx += 1;
					p_iidx = 0;
				}
				if (p_eidx >= lengthElements ())
					goto bail;

				res.push_back (new SequenceT (*this, p_eidx, p_iidx));
			}
		}
	}

bail:
	return (res);
}






std::vector<const Sequence*>
Sequence::prefixFindGap (unsigned int item,
	const std::vector<unsigned int>* lastPrefixElement, unsigned int ahead) const
{
	assert (lastPrefixElement != NULL);

	/* Here we only need to handle vanilla b projections, however, there is a
	 * twist: we need to project on all sequences which match in a given
	 * look-ahead (ahead).
	 */
	unsigned int eidx = firstElementSplitted () ? 1 : 0;
	ahead += eidx;

	std::vector<const Sequence*> res;

	for ( ; eidx <= ahead && eidx < lengthElements () ; ++eidx) {
		unsigned int leidx = lengthSet (eidx);

		for (unsigned int iidx = 0 ; iidx < leidx ; ++iidx) {
			unsigned int citem = get (eidx, iidx);

			if (citem < item)
				continue;	// keep searching

			if (citem > item)
				break;

			/* Found sequence, carry out simple projection
			 */
			{
				unsigned int p_iidx = iidx + 1;
				unsigned int p_eidx = eidx;

				if (p_iidx == lengthSet (p_eidx)) {
					p_eidx += 1;
					p_iidx = 0;
				}
				if (p_eidx >= lengthElements ())
					goto bail;

				res.push_back (new SequenceT (*this, p_eidx, p_iidx));
			}
		}
	}

bail:
	return (res);
}

/*** SequenceT class
 ***/
SequenceT*
SequenceT::readFromFile (const char* filename,  Lexicon & L)
{
	SeqMine::SequenceT* seq = new SeqMine::SequenceT ();

	std::ifstream in (filename);
	if (in.fail ()) {
		delete seq;

		return (NULL);
	}

	/* Read in the sequence: each line a bag, each bag a set of ordered
	 * indices.
	 */
	while (in.eof () == false) {
		//std::cout << "    new elem" << std::endl;
		std::string line;
		std::getline (in, line);

		/* Skip over empty elements: they are not important to our sequence
		 * mining.
		 */
		if (line.size () == 0)
			continue;

		std::istringstream is (line);

		std::set<unsigned int> elemS;
		while (is.eof () == false) {
			unsigned int item;

			is >> item;
			elemS.insert (L.getId(item));
		}

		/* Sort the element (set)
		 */
		std::vector<unsigned int> elem (elemS.size ());
		std::copy (elemS.begin (), elemS.end (), elem.begin ());

		seq->appendElement (elem);
	}
	in.close ();

	return (seq);
}


/* Print with Lexicon */
void SequenceT::decodeAndPrint(Lexicon & lex)
{
	for (unsigned int eidx = 0 ; eidx < lengthElements () ; ++eidx) {
		if (eidx == 0)
			std::cout << "[ ";
		else
			std::cout << " [ ";

		for (unsigned int iidx = 0 ; iidx < lengthSet (eidx) ; ++iidx) {
			if (eidx == 0 && iidx == 0 && firstElementSplitted ())
				std::cout << "_";

			std::cout << lex.getString(get (eidx, iidx)) << " ";
		}
		std::cout << "]";
	}

    return;
}




unsigned int
SequenceT::lengthElements () const
{
	return (seq.size ());
}

unsigned int
SequenceT::lengthSet (unsigned int eidx) const
{
	assert (eidx < lengthElements ());

	return (seq[eidx].size ());
}

unsigned int
SequenceT::lastItem (void) const
{
	// TODO: document the meaning of ::max() as return value.
	if (seq.size () == 0)
		return (std::numeric_limits<unsigned int>::max ());

	assert (seq[seq.size () - 1].size () > 0);
	return (*(--seq[seq.size () - 1].end ()));
}


const std::vector<unsigned int>*
SequenceT::lastElement (void) const
{
	if (seq.size () == 0)
		return (NULL);

	return (&seq[seq.size()-1]);
}


unsigned int
SequenceT::get (unsigned int eidx, unsigned iidx) const
{
	return (seq[eidx][iidx]);
}

bool
SequenceT::contains (unsigned int eidx, unsigned int item,
	unsigned int start_iidx) const
{
	/* Compute minimum required position
	 */
	cvec_iterator sb = seq[eidx].begin ();
	advance (sb, start_iidx);
	cvec_iterator fi = find (sb, seq[eidx].end (), item);

	return (fi != seq[eidx].end ());
}

void
SequenceT::appendElement (std::vector<unsigned int>& elem)
{
	seq.push_back (elem);
}

void
SequenceT::deleteElement (void)
{
	seq.pop_back ();
}

void
SequenceT::appendItem (unsigned int item)
{
	assert (seq.size () > 0);
	seq[seq.size () - 1].push_back (item);
}

void
SequenceT::deleteItem (void)
{
	seq[seq.size () - 1].pop_back ();
}

bool
SequenceT::checkAppendPossible (unsigned int item) const
{
	if (seq.size () == 0)
		return (true);

	if (seq[seq.size () - 1].size () == 0)
		return (true);

	//std::cout << "lastItem: " << lastItem () << ", item: " << item << std::endl;
	return (lastItem () < item);
}


}






std::ostream&
operator<< (std::ostream& os, const SeqMine::Sequence& seq)
{
	for (unsigned int eidx = 0 ; eidx < seq.lengthElements () ; ++eidx) {
		if (eidx == 0)
			os << "[ ";
		else
			os << " [ ";

		for (unsigned int iidx = 0 ; iidx < seq.lengthSet (eidx) ; ++iidx) {
			if (eidx == 0 && iidx == 0 && seq.firstElementSplitted ())
				os << "_";

			os << seq.get (eidx, iidx) << " ";
		}
		os << "]";
	}

	return (os);
}

std::ostream&
operator<< (std::ostream& os, const std::vector<const SeqMine::Sequence*>& seqs)
{
	os << "{ ";
	for (std::vector<const SeqMine::Sequence*>::const_iterator si = seqs.begin () ;
		si != seqs.end () ; ++si)
	{
		os << "[" << *(*si) << "], ";
	}
	os << " }";

	return (os);
}


template <class T, class Iterator>
static void tokenize (const char *str, Iterator iterator)
{
	std::istrstream is (str, std::strlen (str));
	std::copy (std::istream_iterator <T> (is), std::istream_iterator <T> (), iterator);
}

std::istream&
operator>> (std::istream& is, SeqMine::SequenceT& seq)
{
	char line[32768];
	std::vector<std::string> result;


	if (is.getline (line, sizeof(line)-1) == false)
		return (is);

	line[sizeof (line) - 1] = '\0';

	result.clear ();
	tokenize<std::string> (line, std::back_inserter (result));

	if (result.empty ())
		return (is);

#ifdef DEBUG
	for (std::vector<std::string>::iterator ri = result.begin () ;
		ri != result.end () ; ++ri)
	{
		std::cout << "res: " << *ri << std::endl;
	}
#endif

	//std::vector<unsigned int> pset;
	std::list<unsigned int> pset;
	for (unsigned int i = 0 ; i < result.size () ; ++i) {
		if (result[i] == "[")
			continue;

		if (result[i] == "]") {
			pset.unique ();
			std::vector<unsigned int> pvec =
				std::vector<unsigned int> (pset.size ());
			std::copy (pset.begin (), pset.end (), pvec.begin ());
			seq.appendElement (pvec);
			pset.clear ();
			continue;
		}

		if (result[i] == ";") {
			if (pset.empty () == false) {
				pset.unique ();

				std::vector<unsigned int> pvec =
					std::vector<unsigned int> (pset.size ());
				std::copy (pset.begin (), pset.end (), pvec.begin ());
				seq.appendElement (pvec);
			}

			break;
		}

		unsigned int item = atoi (result[i].c_str ());
		pset.push_back (item);
	}

	return (is);
}



