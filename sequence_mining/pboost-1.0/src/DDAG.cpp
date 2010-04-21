
#include <DDAG.h>

namespace Classifier
{

int
DDAG::relativeOffset (int i, int j, int d) const
{
	if (d == 1)
		return (1);

	return (K-i);
}

void
DDAG::constructTable (void)
{
	for (int i = 0 ; i < (int) K ; ++i)
		for (int j = K-1 ; j >= i ; --j)
			table.push_back (std::pair<int, int> (i, j));
}

}

