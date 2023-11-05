#include "bitset.h"
#include <cstring>
#include <algorithm>
#include <iostream>

Bitset::Bitset(size_t bitsetSize)
{
    blockSizeBytes = sizeof(u_int64_t);
    blockSizeBites = blockSizeBytes * 8;
    size = bitsetSize;
    blockCount = (size / (blockSizeBites) + 1) + 2;
    blocks = new u_int64_t[blockCount];

    std::memset(blocks, 0, sizeof(u_int64_t) * (blockCount));
}

Bitset::Bitset(const std::string &bitset, size_t bitsetSize, bool reverseString)
{
    blockSizeBytes = sizeof(u_int64_t);
    blockSizeBites = blockSizeBytes * 8;
    size = bitsetSize;
    blockCount = (size / (blockSizeBites) + 1) + 2;
    blocks = new u_int64_t[blockCount];

    std::memset(blocks, 0, blockSizeBytes * blockCount);

    std::string currentBitset = bitset;

    u_int64_t iteratorMask = 1;

    if (reverseString)
    {
        std::reverse(currentBitset.begin(), currentBitset.end());
    }

    size_t globalBitIterator = 0;
    int currentBlockIterator = 1;
    int currentIterator = 0;
    u_int64_t * currentBlock = blocks + currentBlockIterator;

    while (globalBitIterator < size)
    {
        if (currentBitset[globalBitIterator] == '1')
        {
            *currentBlock |= iteratorMask;
        }
        iteratorMask <<= 1;

        ++currentIterator;
        ++globalBitIterator;
        if (currentIterator >= blockSizeBites)
        {
            iteratorMask = 1;
            currentIterator = 0;
            ++currentBlockIterator;
            currentBlock = blocks + currentBlockIterator;
        }
    }
}

std::string Bitset::printDebug() const
{
    return printDebug(blocks, blockSizeBites);
}

std::string Bitset::printDebug(u_int64_t *localBlocks, int localBlockSizeBites) const
{
    std::string response = "";
    u_int64_t iteratorMask = 1;

    for (size_t currentBlockIterator = 0; currentBlockIterator < blockCount; ++currentBlockIterator)
    {
        u_int64_t currentElement = localBlocks[currentBlockIterator];
        for (size_t j = 0; j < localBlockSizeBites; ++j)
        {
            response += ('0' + (currentElement & iteratorMask));
            currentElement >>= 1;
        }
    }

    std::cout << response;
    std::cout << " | block offset : " << insideBlockOffset;
    std::cout << std::endl;
    return response;
}

std::string Bitset::print()
{
    std::string response = "";
    u_int64_t iteratorMask = 1;
    int currentBlockIterator = 1;
    int currentIterator = insideBlockOffset;
    size_t globalBitIterator = 0;
    if (insideBlockOffset < 0)
    {
        currentBlockIterator = 0;
        currentIterator = blockSizeBites + insideBlockOffset;
    }

    u_int64_t currentBlock = blocks[currentBlockIterator];

    if (insideBlockOffset > 0)
    {
        currentBlock >>= insideBlockOffset;
    }

    while (globalBitIterator < size)
    {
        response += ('0' + (currentBlock & iteratorMask));
        currentBlock >>= 1;

        ++currentIterator;
        ++globalBitIterator;
        if (currentIterator >= blockSizeBites)
        {
            currentIterator = 0;
            ++currentBlockIterator;
            currentBlock = blocks[currentBlockIterator];
        }
    }
    std::cout << response << std::endl;
    return response;
}

Bitset &Bitset::operator>>=(int shiftNumber)
{
    int offsetBlockShift = (insideBlockOffset - shiftNumber) / blockSizeBites; // Take care of offsetBlockShift int overflow
    insideBlockOffset = (insideBlockOffset - shiftNumber) % blockSizeBites;

    std::memmove(blocks, blocks + offsetBlockShift, (blockCount + offsetBlockShift) * blockSizeBytes);

    std::memset(blocks, 0, (-offsetBlockShift) * blockSizeBytes);
    /*
    for (int currentBlock = blockCount - 1; currentBlock + offsetBlockShift >= 0; currentBlock--)
    {
        // std::cout << currentBlock << "<-" << currentBlock + offsetBlockShift << std::endl;
        blocks[currentBlock] = blocks[currentBlock + offsetBlockShift];
    }

    for (size_t currentBlock = 0; currentBlock < -offsetBlockShift; currentBlock++)
    {
        blocks[currentBlock] = 0;
    }
    */
    return *this;
}

Bitset &Bitset::operator<<=(int shiftNumber)
{
    int offsetBlockShift = (insideBlockOffset + shiftNumber) / blockSizeBites; // Take care of offsetBlockShift int overflow
    insideBlockOffset = (insideBlockOffset + shiftNumber) % blockSizeBites;

    std::memmove(blocks, blocks + offsetBlockShift, (blockCount - offsetBlockShift) * blockSizeBytes);

    std::memset(blocks - (blockCount - offsetBlockShift), 0, (offsetBlockShift)*blockSizeBytes);
    /*
    for (int currentBlock = blockCount - 1; currentBlock + offsetBlockShift >= 0; currentBlock--)
    {
        // std::cout << currentBlock << "<-" << currentBlock + offsetBlockShift << std::endl;
        blocks[currentBlock] = blocks[currentBlock + offsetBlockShift];
    }

    for (size_t currentBlock = 0; currentBlock < -offsetBlockShift; currentBlock++)
    {
        blocks[currentBlock] = 0;
    }
    */
    return *this;
}

bool Bitset::compare(int fromIndex, const Bitset &compareTo, int toIndex, int length) const
{
    int lengthBlock = length / blockSizeBites + 1;
    u_int64_t *copyCompareFromBlocks = new u_int64_t[lengthBlock];
    u_int64_t *copyCompareToBlocks = new u_int64_t[lengthBlock];

    alignBlocks(copyCompareFromBlocks, fromIndex, copyCompareToBlocks, compareTo, toIndex, length, lengthBlock);

    doReset(copyCompareFromBlocks, blockSizeBites, length, lengthBlock * blockSizeBites);
    doReset(copyCompareToBlocks, compareTo.getBlockSizeBites(), length, lengthBlock * compareTo.getBlockSizeBites());

    // std::cout << "Reseted blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);

    int response = memcmp(copyCompareToBlocks, copyCompareFromBlocks, lengthBlock * blockSizeBytes);

    delete copyCompareFromBlocks;
    delete copyCompareToBlocks;

    return (response == 0);
}

int Bitset::compareDistance(int fromIndex, const Bitset &compareTo, int toIndex, int length) const
{
    int lengthBlock = length / blockSizeBites + 1;
    u_int64_t *copyCompareFromBlocks = new u_int64_t[lengthBlock];
    u_int64_t *copyCompareToBlocks = new u_int64_t[lengthBlock];

    alignBlocks(copyCompareFromBlocks, fromIndex, copyCompareToBlocks, compareTo, toIndex, length, lengthBlock);

    
    // std::cout << "Aligned blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);
    

    int currentBlockIterator = 0;
    int globalBitIterator = 0;
    int currentIterator = 0;

    u_int64_t currentBlockFrom = copyCompareFromBlocks[currentBlockIterator];
    u_int64_t currentBlockTo = copyCompareToBlocks[currentBlockIterator];

    u_int64_t iteratorMask = 1;
    int count = 0;

    while (globalBitIterator < length)
    {
        if ((currentBlockFrom & iteratorMask) != (currentBlockTo & iteratorMask))
        {
            ++count;
        }
        iteratorMask <<= 1;
        ++globalBitIterator;
        ++currentIterator;
        if (currentIterator >= blockSizeBites)
        {
            iteratorMask = 1;
            currentIterator = 0;
            ++currentBlockIterator;
            currentBlockFrom = copyCompareFromBlocks[currentBlockIterator];
            currentBlockTo = copyCompareToBlocks[currentBlockIterator];
        }
    }

    delete copyCompareFromBlocks;
    delete copyCompareToBlocks;

    return count;
}

/*
 * fromIndex in virtual space
 */
void Bitset::alignBlocks(u_int64_t *copyCompareFromBlocks, int fromIndex, u_int64_t *copyCompareToBlocks, const Bitset &compareTo, int toIndex, int length, int lengthBlock) const
{
    u_int64_t *const compareFromBlocks = blocks;
    u_int64_t *const compareToBlocks = compareTo.getBlocks();

    /*
    std::cout << "Enter bitsets: " << std::endl;
    this->printDebug(compareFromBlocks, blockSizeBites);
    compareTo.printDebug(compareToBlocks, compareTo.getBlockSizeBites());
    */

    toIndex += compareTo.getInsideBlockOffset();
    int toBlockIndex = toIndex / compareTo.getBlockSizeBites() + 1;
    int toBitIndex = toIndex % compareTo.getBlockSizeBites();
    int toLengthBlock = lengthBlock;

    if (toBitIndex < 0)
    {
        toBlockIndex -= 1;
        toBitIndex = compareTo.getBlockSizeBites() + toBitIndex;
        ++toLengthBlock;
    }

    fromIndex += insideBlockOffset;
    int fromBlockIndex = fromIndex / blockSizeBites + 1;
    int fromBitIndex = fromIndex % blockSizeBites;
    int fromLengthBlock = lengthBlock;

    if (fromBitIndex < 0)
    {
        fromBlockIndex -= 1;
        fromBitIndex = blockSizeBites + fromBitIndex;
        ++fromLengthBlock;
    }

    // std::cout << "fromBlockIndex:" << fromBlockIndex << std::endl;
    // std::cout << "fromBitIndex:" << fromBitIndex << std::endl;
    // std::cout << "FromInsideBlockOffset:" << insideBlockOffset << std::endl;
    // std::cout << "fromLengthBlock:" << fromLengthBlock << std::endl;
    // std::cout << "toBlockIndex:" << toBlockIndex << std::endl;
    // std::cout << "toBitIndex:" << toBitIndex << std::endl;
    // std::cout << "lengthBlock:" << lengthBlock << std::endl;

    memcpy(copyCompareFromBlocks, compareFromBlocks + fromBlockIndex, fromLengthBlock * blockSizeBytes);

    memcpy(copyCompareToBlocks, compareToBlocks + toBlockIndex, toLengthBlock * compareTo.getBlockSizeBytes());

    // std::cout << "Copied bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // this->printDebug(copyCompareToBlocks, blockSizeBites);

    // Shift few elements (boost/dynamic_bitset)
    if (fromBitIndex > 0)
    {
        const int leftBitShift = blockSizeBites - fromBitIndex;

        for (int i = 0; i < fromLengthBlock - 1; ++i)
        {
            copyCompareFromBlocks[i] = (copyCompareFromBlocks[i] >> fromBitIndex) | (copyCompareFromBlocks[i + 1] << leftBitShift);
        }

        copyCompareFromBlocks[fromLengthBlock - 1] = copyCompareFromBlocks[fromLengthBlock - 1] >> fromBitIndex;
    }

    if (toBitIndex > 0)
    {
        const int leftBitShift = compareTo.getBlockSizeBites() - toBitIndex;

        for (int i = 0; i < toLengthBlock - 1; ++i)
        {
            copyCompareToBlocks[i] = (copyCompareToBlocks[i] >> toBitIndex) | (copyCompareToBlocks[i + 1] << leftBitShift);
        }

        copyCompareToBlocks[toLengthBlock - 1] = copyCompareToBlocks[toLengthBlock - 1] >> toBitIndex;
    }

    // std::cout << "Copied blocks from bitsets: " << std::endl;
    // this->printDebug(copyCompareFromBlocks, blockSizeBites);
    // compareTo.printDebug(copyCompareToBlocks, blockSizeBites);
}

// fromBitIndex inclusive, ToBitIndex exclusive
// It doesnt take into account right and left offset inside block
void Bitset::doReset(u_int64_t *blocks, int blockBiteSize, int fromBiteIndex, int toBiteIndex)
{
    u_int64_t constMask = 1;

    int distance = toBiteIndex - fromBiteIndex;
    int completeBlocksSize = distance / blockBiteSize;

    int currentBlockBiteIndex = fromBiteIndex % blockBiteSize;
    int currentBlockIndex = fromBiteIndex / blockBiteSize;

    u_int64_t *currentBlock = blocks + currentBlockIndex;

    // std::cout << "distance:" << distance << std::endl;
    // std::cout << "completeBlocksSize:" << completeBlocksSize << std::endl;
    // std::cout << "currentBlockIndex:" << currentBlockIndex << std::endl;
    // std::cout << "currentBlockBiteIndex:" << currentBlockBiteIndex << std::endl;

    constMask <<= currentBlockBiteIndex;

    while (fromBiteIndex < toBiteIndex)
    {

        (*currentBlock) &= (~constMask);

        constMask <<= 1;
        ++fromBiteIndex;
        ++currentBlockBiteIndex;
        if (currentBlockBiteIndex >= blockBiteSize)
        {
            ++currentBlock;
            for (int completeBlocksIndex = 0; completeBlocksIndex < completeBlocksSize; ++completeBlocksIndex)
            {
                (*currentBlock) = 0;
                ++currentBlock;
                fromBiteIndex += blockBiteSize;
            }
            currentBlockBiteIndex = 0;
            constMask = 1;
        }
    }
}

bool Bitset::operator[](int targetNumber)
{
    targetNumber += insideBlockOffset;
    int currentBlockIndex = targetNumber / blockSizeBites + 1;
    int currentBlockBiteIndex = targetNumber % blockSizeBites;

    u_int64_t mask = 1;
    bool value = (blocks[currentBlockIndex] >> currentBlockBiteIndex) & mask;
    return value;
}

void Bitset::set(int targetNumber, bool value)
{
    targetNumber += insideBlockOffset;
    int currentBlockIndex = targetNumber / blockSizeBites + 1;
    int currentBlockBiteIndex = targetNumber % blockSizeBites;
    u_int64_t removeMask = 1;
    removeMask <<= currentBlockBiteIndex;

    if (value == false)
    {
        blocks[currentBlockIndex] &= (~removeMask);
    }
    else
    {
        blocks[currentBlockIndex] |= (removeMask);
    }
}

void Bitset::flip(int targetNumber)
{
    targetNumber += insideBlockOffset;
    int currentBlockIndex = targetNumber / blockSizeBites + 1;
    int currentBlockBiteIndex = targetNumber % blockSizeBites;
    u_int64_t flipMask = 1;
    flipMask <<= currentBlockBiteIndex;

    blocks[currentBlockIndex] ^= flipMask;
}

int Bitset::getInsideBlockOffset() const
{
    return insideBlockOffset;
}

size_t Bitset::bitsetSize() const
{
    return size;
}

u_int64_t *Bitset::getBlocks() const
{
    return blocks;
}

int Bitset::getBlockSizeBites() const
{
    return blockSizeBites;
}

int Bitset::getBlockSizeBytes() const
{
    return blockSizeBytes;
}

int Bitset::getBlockCount() const
{
    return blockCount;
}

Bitset::~Bitset()
{
    // delete[] blocks;
}