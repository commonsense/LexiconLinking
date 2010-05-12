
#ifndef	EXPLORATIONGUIDE_H
#define	EXPLORATIONGUIDE_H

// TODO: someday replace with
// http://www.boost.org/libs/multi_index/doc/index.html

#include <algorithm>
#include <vector>
#include <map>
#include <limits>
#include <functional>

#include <Sequence.h>

namespace SeqMine
{

class ExplorationElement
{
public:
	SequenceT alpha;
	std::vector<const Sequence*> Sproj;
	bool cleanupS;

	ExplorationElement (void)
	{
	}

	ExplorationElement (SequenceT& alpha,
		std::vector<const Sequence*>& Sproj, bool cleanupS)
		: alpha (alpha), Sproj (Sproj), cleanupS (cleanupS)
	{
	}

	void Cleanup (void)
	{
		for (std::vector<const Sequence*>::iterator si = Sproj.begin () ;
			si != Sproj.end () ; ++si)
		{
			delete (*si);
		}
	}
};

class
ExplorationGuide
{
private:
	unsigned long long uid;

	/* Two redundant heaps in order to realize the following two operations:
	 *
	 * 1. Deletion of an element with a given id.
	 * 2. Deletion of any element with a gainbound below or equal to a given
	 *    gain.
	 */
	std::multimap<double, unsigned long long, std::less<double> > gainboundMap1;
	std::map<unsigned long long, double, std::less<double> > gainboundMap2;

	std::multimap<double, unsigned long long, std::greater<double> > gainMap1;
	std::map<unsigned long long, double, std::greater<double> > gainMap2;

	std::map<unsigned long long, ExplorationElement> explorationStore;

	unsigned int informNewGain (double gain)
	{
		/* Nothing there -> nothing to remove ;-)
		 */
		if (IsEmpty ())
			return (0);

		/* Now, the remaining patterns can be safely removed if their
		 * gainbound is less than or equal to the current gain.
		 */
		unsigned int removeCount = 0;
		for (std::multimap<double, unsigned long long, std::less<double> >::iterator
			iv = gainboundMap1.begin () ; iv != gainboundMap1.end () ; ++iv)
		{
			/* gainbound exceeds the gain -> nothing to remove left.
			 */
			if (iv->first > gain)
				break;

			/* The current sequence can be dropped.  Do so safely.
			 */
			removeElement (iv->second, true);
			removeCount += 1;
		}

		return (removeCount);
	}

	ExplorationElement removeElement (unsigned long long id, bool cleanup)
	{
		double gainbound = gainboundMap2[id];
		gainboundMap2.erase (id);

		for (std::multimap<double, unsigned long long, std::less<double> >::iterator
			iv = gainboundMap1.find (gainbound) ; iv != gainboundMap1.end () ; ++iv)
		{
			if (iv->first > gainbound)
				break;

			if (iv->second == id) {
				gainboundMap1.erase (iv);
				break;
			}
		}

		double gain = gainMap2[id];
		gainMap2.erase (id);

		for (std::multimap<double, unsigned long long, std::greater<double> >::iterator
			iv = gainMap1.find (gain) ; iv != gainMap1.end () ; ++iv)
		{
			if (iv->first > gain)
				break;

			if (iv->second == id) {
				gainMap1.erase (iv);
				break;
			}
		}

		ExplorationElement res = explorationStore[id];
		if (cleanup)
			res.Cleanup ();

		explorationStore.erase (id);

		return (res);
	}

public:
	bool IsEmpty (void) const
	{
		return (gainMap1.empty ());
	}

	double NextGain (void) const
	{
		if (IsEmpty ())
			return (std::numeric_limits<double>::min ());

		return (gainMap1.begin ()->first);
	}

	void Add (double gain, double gainbound, SequenceT& alpha,
		std::vector<const Sequence*>& Sproj, bool cleanupS)
	{
		unsigned long long newId = uid++;

		gainboundMap1.insert (std::pair<double, unsigned long long>
			(gainbound, newId));
		gainboundMap2.insert (std::pair<unsigned long long, double>
			(newId, gainbound));

		gainMap1.insert (std::pair<double, unsigned long long> (gain, newId));
		gainMap2.insert (std::pair<unsigned long long, double> (newId, gain));
		explorationStore[newId] = ExplorationElement (alpha, Sproj, cleanupS);
	}

	ExplorationElement PopTop (void)
	{
		unsigned long long id = gainMap1.begin ()->second;

		ExplorationElement res = removeElement (id, false);

		return (res);
	}

	ExplorationGuide (void)
	{
		uid = 0;
	}
};

};

#endif

