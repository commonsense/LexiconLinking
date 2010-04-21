
#ifndef	EXPLORATIONGUIDE_H
#define	EXPLORATIONGUIDE_H

#include <algorithm>
#include <vector>
#include <map>
#include <limits>
#include <functional>
#include <math.h>

#include <Sequence.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/ordered_index.hpp>

using namespace ::boost::multi_index;

namespace SeqMine
{

struct ExplorationElement
{
	double gain;
	double gainbound;

	SequenceT alpha;
	mutable std::vector<const Sequence*> Sproj;
	bool cleanupS;

	ExplorationElement (double gain, double gainbound, SequenceT& alpha,
		std::vector<const Sequence*>& Sproj, bool cleanupS)
		: gain (gain), gainbound (gainbound), alpha (alpha),
			Sproj (Sproj), cleanupS (cleanupS)
	{
	}

	bool operator< (const ExplorationElement& ee) const
	{
		/* Ordering relation: if gain is almost equal, prefer smaller
		 * sequences, if gain is different, choose higher gain first.
		 */
		if (fabs (gain - ee.gain) < 1e-5)
			return (alpha.length () < ee.alpha.length ());

		return (gain > ee.gain);
	}

	bool Cleanup (void) const
	{
		if (cleanupS == false)
			return (true);

		for (std::vector<const Sequence*>::iterator si = Sproj.begin () ;
			si != Sproj.end () ; ++si)
		{
			delete (*si);
		}

		return (true);
	}
};

typedef multi_index_container<
	ExplorationElement,
	indexed_by<
		/* First ordering relation: by gain (tie break by size), in order to
		 * be able to quickly pop the best gain sequence
		 */
		ordered_non_unique<identity<ExplorationElement> >,

		/* Second ordering relation: by gainbound, increasing
		 */
		ordered_non_unique<member<ExplorationElement, double,
			&ExplorationElement::gainbound> >
	>
> ExplorationList;


class
ExplorationGuide
{
private:
	ExplorationList elist;

	void discardByNewGain (double gain)
	{
		/* Obtain the iterator that points to the first element with a higher
		 * gainbound than the current gain
		 */
		ExplorationList::nth_index_iterator<1>::type firstgood
			= elist.get<1>().upper_bound (gain);

		/* Remove all elements with a gainbound below the current gain
		 */
		for (ExplorationList::nth_index_iterator<1>::type iv =
			elist.get<1>().begin () ; iv != firstgood ; ++iv)
		{
//			elist.modify (iv, Cleanup);
			iv->Cleanup ();
		}

		elist.get<1>().erase (elist.get<1>().begin (), firstgood);
	}

public:
	bool IsEmpty (void) const
	{
		return (elist.empty ());
	}

	double NextGain (void) const
	{
		if (IsEmpty ())
			return (std::numeric_limits<double>::min ());

		return (elist.begin ()->gain);
	}

	void Add (ExplorationElement& el)
	{
		elist.insert (el);
		discardByNewGain (el.gain);
	}

	size_t Size (void) const
	{
		return (elist.size ());
	}

	ExplorationElement PopTop (void)
	{
		ExplorationElement res = *elist.get<0>().begin ();
		elist.get<0>().erase (elist.get<0>().begin ());

		return (res);
	}
};

};

#endif

