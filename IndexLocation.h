/*
 * IndexLocation.h
 *
 *  Created on: Oct 13, 2013
 *      Author: ivan
 */

#ifndef INDEXLOCATION_H_
#define INDEXLOCATION_H_

#include <ostream>



class IndexLocation
{
public:
	/** ID of the reference sequence that has been hit. */
	unsigned long long int sequenceId;
	/** Start of the seed hit on the reference sequence (location given in the number of bases from the beginning of the reference sequence). */
	unsigned long long int start;
	/** If equal to zero, the sequence is reversed, otherwise it is on the forward strand (don't ask me why I chose this way of denoting it, instead of the other way around. I don't know myself.). */
	unsigned char isReverse;

	IndexLocation()
	{
		sequenceId = start = 0;
		isReverse = 0;
	}

	IndexLocation(unsigned long long int newSequenceId, unsigned long long int newStart, unsigned char newIsReverse)
	{
		sequenceId = newSequenceId;
		start = newStart;
		isReverse = newIsReverse;
	}

	void clear()
	{
		sequenceId = start = 0;
		isReverse = 0;
	}

	void verbose(std::ostream &outStream)
	{
		outStream << "(sequenceId, start) = (" << sequenceId << ", " << start << ", ";
		if (isReverse)
			outStream << "forward";
		else
			outStream << "reverse";
		outStream << ")";
	}
};

#endif /* INDEXLOCATION_H_ */
