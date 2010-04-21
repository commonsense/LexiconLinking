
#ifndef	BESTPATTERNLIST_H
#define	BESTPATTERNLIST_H

#include <algorithm>
#include <map>
#include <limits>
#include <functional>

#include <Sequence.h>

namespace SeqMine
{

/* Ordered (distance, (Sequence, ybase)) tupples
 */
class
BestPatternList :
	public std::multimap<double, std::pair<SequenceT, double>, std::greater<double> >
{
private:
	unsigned int max_size;

public:
	BestPatternList (void)
	{
		max_size = std::numeric_limits<unsigned int>::max ();
	}

	BestPatternList (unsigned int maximum_size)
		: max_size (maximum_size)
	{
	}

	double LeastBestGain (void) const
	{
		if (empty ())
			return (std::numeric_limits<double>::min ());

		return ((--end ())->first);
	}

	bool IsCapacityLeft (void) const
	{
		return (size () < max_size);
	}

	/* Return true if the pattern was actually inserted
	 */
	bool Insert (double gain, const Sequence& seq, double ybase)
	{
		if (LeastBestGain () >= gain)
			return (false);

		std::multimap<double, std::pair<SequenceT, double>,
			std::greater<double> >::insert
				(std::pair<double, std::pair<SequenceT, double> > (gain,
					std::pair<SequenceT, double> (SequenceT (seq), ybase)));

		if (size () > max_size)
			erase (--end ());

		return (true);
	}
};

};

#endif

