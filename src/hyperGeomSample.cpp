#include <random>
#include <iostream>
#include <array>

unsigned int* hyperGeomSample( std::vector<unsigned int> binCounts,  unsigned int totalDrawSize ) {
    unsigned int totalPopnSize = 0;

    unsigned int numBins = binCounts.size();
    unsigned int* cumBinCounts = new unsigned int[numBins];
    unsigned int* sample = new unsigned int[numBins];


    for (unsigned int idx = 0; idx < numBins; idx++ ) {
	totalPopnSize += binCounts[idx];
	cumBinCounts[idx] = totalPopnSize;
	sample[idx] = 0;
    }

    std::random_device rd;
    unsigned int randVal = 0;
    unsigned int mid,lb,ub,midPrev,popnSampledFrom = 0;

    for (int drawNumber = totalDrawSize; drawNumber > 0; drawNumber-- ) {	
	std::uniform_int_distribution<> unifIntDist(1, totalPopnSize);

	randVal = unifIntDist(rd);
	
	popnSampledFrom = numBins - 1;

	if ( randVal <= cumBinCounts[0] ) 
	    popnSampledFrom = 0;
	else if ( randVal > cumBinCounts[numBins-2] & numBins > 2 )
	    popnSampledFrom = numBins - 1;
	else {
	    lb = 0;
	    ub = numBins-1;

	    mid = 0;
	    midPrev = numBins;
	    while (lb < ub ) {
		mid = (lb + ub)/2;
		if ( midPrev == mid ) {
     		    popnSampledFrom = mid + 1;
		    break;
		}
		midPrev = mid;
		if ( cumBinCounts[mid] < randVal )
		    lb = mid;
		else if ( cumBinCounts[mid] > randVal )
		    ub = mid;
		else {
		    popnSampledFrom = mid;
		    break;
		}
	    }
	}

	++sample[popnSampledFrom];
	--totalPopnSize;

	for ( int j = popnSampledFrom; j < numBins; j++ )
	    cumBinCounts[j]--;

    }

    for (int j = 0; j < numBins; j++ )
	std::cout << sample[j] << std::endl;

    delete[] cumBinCounts;
    return sample;
}

int main() {
    unsigned int totalDrawSize = 20;
    std::vector<unsigned int> binCounts {20,1,3,5,1,2};

    hyperGeomSample( binCounts, totalDrawSize );
}

