
#ifndef	CLASSIFIER_H
#define	CLASSIFIER_H

/* Classifier result class.
 *
 * This stores a boosting-type sequence classifier, single and multiclass
 * versions.
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 30th March 2007
 */

#include <iostream>
#include <fstream>
#include <vector>

// libboost-serialization in text format, input/output
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <Sequence.h>
#include <DDAG.h>

namespace
Classifier
{

class TwoClassClassifier
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize (Archive & ar, const unsigned int version)
    {
		ar & S;
		ar & alpha;
		ar & classthresh;
		ar & cltype;
		ar & directional;
		ar & maxgap;
    }

public:
	/* Stump sequences
	 */
	std::vector<SeqMine::SequenceT*> S;

	/* Stump weights
	 */
	std::vector<double> alpha;

	/* Class decision threshold
	 */
	double classthresh;

	/* Classifier type
	 */
	unsigned int cltype;

	/* Directional stump flag
	 */
	bool directional;

	/* Maximum gap constraint (-1 for no constraint on the gap length)
	 */
	int maxgap;

	TwoClassClassifier (void)
	{
	}

	TwoClassClassifier (const std::vector<SeqMine::SequenceT*>& S,
		const std::vector<double>& alpha, double classthresh,
		unsigned int cltype, bool directional, int maxgap)
		: S (S), alpha (alpha), classthresh (classthresh), cltype (cltype),
			directional (directional), maxgap (maxgap)
	{
	}
};

class MulticlassClassifier
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize (Archive & ar, const unsigned int version)
    {
		ar & ddag;
		ar & classifiers;
    }

public:
	/* The multiclass decision dag
	 */
	DDAG ddag;

	/* The corresponding 1-vs-1 classifiers
	 */
	std::vector<Classifier::TwoClassClassifier*> classifiers;

	MulticlassClassifier (void)
	{
	}

	MulticlassClassifier (DDAG& ddag,
		std::vector<TwoClassClassifier*>& classifiers)
		: ddag (ddag), classifiers (classifiers)
	{
	}
};

};

BOOST_CLASS_TRACKING(Classifier::MulticlassClassifier, boost::serialization::track_never)
BOOST_CLASS_TRACKING(Classifier::TwoClassClassifier, boost::serialization::track_never)

#endif

