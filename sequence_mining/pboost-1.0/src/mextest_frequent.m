
% Test matlab wrapper for weighted subsequence mining.

S1 = {uint32([1 2 3]), uint32([2 4]), uint32([5]), uint32([5]), uint32([1 2 6])};
S2 = {uint32([4]), uint32([2 4]), uint32([7]), uint32([1 3 5]), uint32([6]), uint32([3])};
S3 = {uint32([2 3]), uint32([1 3 5]), uint32([5]), uint32([7])};
S4 = {uint32([1 6]), uint32([1 3 5]), uint32([1 2 3]), uint32([4]), uint32([7])};

S = {S1, S2, S3, S4};

[subseqs, occurence] = mexprefixspan (S, 2, [0 0]);

