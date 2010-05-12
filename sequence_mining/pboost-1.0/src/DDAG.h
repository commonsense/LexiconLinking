
#ifndef	DDAG_H
#define	DDAG_H

#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

namespace Classifier
{

class
DDAG
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize (Archive & ar, const unsigned int version)
    {
		ar & K;
		ar & table;
    }

	void constructTable (void);

public:
	/* Number of classes
	 */
	unsigned int K;

	/* The DDAG table, encoding the i-vs-j classes.  Two identical class
	 * labels denote a final multiclass decision for the respective class.
	 *
	 * So it usually looks like:
	 * <1:4>, <1:3>, <1:2>, <1:1>, <2:4>, <2:3>, <2:2>, <3:4>, <3:3>, <4:4>
	 *
	 * Evaluation starts at the left hand side and is navigated by the
	 * relativeOffset method.
	 */
	std::vector<std::pair<int, int> > table;

	/* Build the DDAG
	 */
	DDAG (unsigned int K)
		: K (K)
	{
		constructTable ();
	}

	DDAG (void)
	{
	}

	/* Relative offset in DDAG table.  Our current position is the 'i'-vs-'j'
	 * classifier and a decision 'd' (+1 or -1) has been made.  Then the
	 * relative offset tells us how to move in the DDAG array.
	 */
	int relativeOffset (int i, int j, int d) const;
};

};

BOOST_CLASS_TRACKING(Classifier::DDAG, boost::serialization::track_never)

#endif

