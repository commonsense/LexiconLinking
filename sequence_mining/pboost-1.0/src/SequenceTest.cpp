
#include <iostream>

#include "Sequence.h"

using namespace SeqMine;

int
main (int argc, char* argv[])
{
	SequenceT st = SequenceT ();
	std::vector<unsigned int> iset = std::vector<unsigned int> ();

	iset.push_back (1);
	iset.push_back (2);
	iset.push_back (4);
	iset.push_back (6);

	std::vector<unsigned int> iset2 = std::vector<unsigned int> ();
	iset2.push_back (3);
	iset2.push_back (4);

	st.appendElement (iset);
	st.appendElement (iset2);

	std::cout << "st" << std::endl;
	std::cout << "length: " << st.length () << std::endl;
	std::cout << "length elements: " << st.lengthElements () << std::endl;
	for (unsigned int n = 0 ; n < st.lengthElements () ; ++n)
		std::cout << "   elem " << n << ": " << st.lengthSet (n) << std::endl;

	std::cout << std::endl;

	/* Test projection
	 */
	std::vector<unsigned int> alpha = std::vector<unsigned int> ();
	alpha.push_back (4);

	SequenceT stalpha = SequenceT ();
	stalpha.appendElement (alpha);
	//SequenceTSearch stalphaSearch = SequenceTSearch (stalpha);

	const Sequence* proj = stalphaSearch.prefixProject (st, false);

	std::cout << "proj" << std::endl;
	std::cout << "length: " << proj->length () << std::endl;
	std::cout << "length elements: " << proj->lengthElements () << std::endl;
	for (unsigned int n = 0 ; n < proj->lengthElements () ; ++n) {
		std::cout << "   elem " << n << ": " << proj->lengthSet (n) << std::endl;
		for (unsigned int j = 0 ; j < proj->lengthSet (n) ; ++j)
			std::cout << "     item " << j << ": "
				<< proj->get (n, j) << std::endl;
	}

	std::cout << std::endl;
}


